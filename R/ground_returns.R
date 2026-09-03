## SusNRunPointCloud_Workflow
## Standalone R script assembled from the code chunks in
## SusNRun_PointCloudWorkflow_blank_082626.qmd
##
## Point Cloud Pre-Processing Workflow
##
## To keep any single output from growing unmanageably large on big lidar collections, the full
## set of input tiles is split into batches (a configurable number of tiles each). Steps 2-4
## below run once per batch -- each batch gets its own merged LAZ and (optionally) its own
## GeoPackage -- and the per-batch count rasters from step 4 are mosaicked into one whole-site
## raster before the final reprojection in step 5. Every batch is always filtered down to ground
## returns (Classification == 2) in step 2, whether or not the optional GeoPackage in step 3 is
## created -- the GeoPackage toggle only controls whether those already-ground-only points also
## get written out as a vector layer, not whether the filtering itself happens.
##
## 1. Split the full set of LAZ tiles into batches of a configurable size
## 2. Filter ground returns (Classification == 2) and merge LAZ tiles, one merged file per batch
##    -- always happens, regardless of the GeoPackage option below
## 3. (Optional) Write each batch's ground-only points to its own GeoPackage as a point feature
##    class
## 4. Compute a point-count raster per batch (target pixel size = 3 m), then mosaic all batches
##    together into one whole-site count raster and derive density from it -- the count raster
##    itself is only an internal step here and is never written out as a final output
## 5. Reproject the combined density raster to EPSG:2856 (the workflow's only raster output)
##
## Requires: lasR, terra, sf, future, future.apply
##
## For terra, sf, future, and future.apply: install.packages(c("terra", "sf", "future", "future.apply"))
##
## lasR is not on CRAN (by design, not due to archival)
## Install with: install.packages("lasR", repos = "https://r-lidar.r-universe.dev")
##
## All input parameters below (input_dir, output_dir, n_workers, batch_size_tiles, etc.) are
## read from a plain-text parameter file rather than hardcoded in this script -- see the
## "Parameter file" section just below for its format, and
## ground_returns_params_template.txt for a ready-to-copy example. Pass that file's
## path as a command-line argument when running via Rscript:
##   Rscript ground_returns.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

## ---- Inputs -----------------------------------------------------------------------------

# Marks the start of the entire workflow (first thing that runs), so total run time can be
# reported at the end alongside the start/end timestamps -- see the "Workflow timing" section
# at the bottom of this file.
workflow_start_time <- Sys.time()

library(lasR)          # point cloud engine: filters/merges LAZ tiles and computes hulls
library(terra)         # raster creation, math (density calc), and reprojection
library(sf)            # vector I/O (GeoPackage) and CRS handling for the point layer
library(future)        # backend for running batches in parallel worker processes (see plan() below)
library(future.apply)  # future_lapply() -- parallel dispatch across batches, used in "Process batches"

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input (input_dir, output_dir, n_workers, batch_size_tiles, gpkg_layer,
## export_ground_points_gpkg, target_pixel_size_m, target_crs, fill_empty_cells_with_zero,
## cleanup_intermediate_files, stream_chunk_size) is read from a plain-text parameter file
## using "keyword: value" lines, one per line -- e.g.:
##
##   input_dir: c:\work\data\site1\cascades_north_wali_2023\laz
##   output_dir: c:\work\data\site1\cascades_north_wali_2023\grnd_den\
##   n_workers: 6
##   batch_size_tiles: 50
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types,
## and their defaults. Blank lines and lines starting with "#" are ignored, and a trailing
## "# comment" after a value is stripped. A value may itself contain a colon (e.g. a Windows
## drive letter, "c:\...", or "EPSG:2856") -- only the FIRST colon on a line splits the keyword
## from its value, so that's safe. input_dir and output_dir have no default and must be present
## in the file; every other keyword falls back to the default shown below if omitted.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below (and other scripts' references to a sibling script or template
# file) is found by this script's location on disk rather than by the caller's working
# directory -- e.g. `Rscript R/ground_returns.R` from the repo root and `Rscript ground_returns.R`
# from inside R/ both find R/ground_returns_params_template.txt. Falls back to "." (assume cwd)
# when there's no "--file=" argument to read, e.g. when this script is source()'d from RStudio
# instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript ground_returns.R path/to/params.txt` -- the line below is only used as a
# fallback when no such argument is given (e.g. sourcing this script from RStudio), so set it
# there in that case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "ground_returns_params_template.txt")

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) config_path <- cli_args[1]  # Rscript ... params.txt takes precedence over the hardcoded fallback above

if (!file.exists(config_path)) {
  stop("Parameter file not found: ", config_path, ". Set config_path in this script, or pass ",
       "the parameter file's path as a command-line argument (Rscript ... params.txt).")
}

# Reads a "keyword: value" parameter file into a named list of raw (still character) values.
# Order doesn't matter -- every line is self-labeled by its keyword -- so this just builds a
# lookup table; build_params() below is what applies types/defaults/required-ness on top of it.
read_param_file <- function(path) {
  raw_lines <- readLines(path, warn = FALSE)
  raw_lines <- trimws(raw_lines)
  raw_lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]  # drop blank lines and comment-only lines

  keys <- character(0)
  vals <- character(0)
  for (line in raw_lines) {
    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the FIRST colon splits keyword from value -- keeps values like "c:\..." or "EPSG:2856" intact
    if (colon_pos < 0) {
      stop("Malformed line in parameter file (expected \"keyword: value\"): ", line)
    }
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(substr(line, colon_pos + 1, nchar(line)))
    val <- trimws(sub("\\s+#.*$", "", val))  # strip a trailing " # comment", if any, from the value
    keys <- c(keys, key)
    vals <- c(vals, val)
  }

  if (anyDuplicated(keys)) {
    stop("Duplicate keyword(s) in parameter file ", path, ": ",
         paste(unique(keys[duplicated(keys)]), collapse = ", "))
  }

  setNames(as.list(vals), keys)
}

# The full set of recognized parameter-file keywords: expected type, and default value used
# when a keyword is omitted from the file (required = TRUE means there is no default and the
# file must supply it). This table is the single place that documents what each keyword means;
# comments here mirror the ones that used to sit next to each hardcoded assignment.
param_specs <- list(
  # Folder containing the LAZ/LAS tiles for this site/date. No default -- must be set in the
  # parameter file.
  input_dir = list(type = "character", required = TRUE),

  # Output folder for all intermediate and final files. Created (with any missing parent
  # folders) if it doesn't already exist. No default -- must be set in the parameter file.
  output_dir = list(type = "character", required = TRUE),

  # Maximum number of batches processed at once, in separate worker processes (see "Process
  # batches" below, where this is capped to n_batches -- no point spinning up idle workers if
  # there are fewer batches than n_workers). multisession spawns background Rscript processes
  # rather than forking -- required on Windows, and also works cross-platform -- so each worker
  # starts with none of the packages from this session loaded; process_one_batch() reloads what
  # it needs itself.
  # NOTE: lasR's exec() calls already multithread internally (no explicit ncores is set below,
  # so each call may use all available cores on its own). Running n_workers batches at once on
  # top of that can oversubscribe the CPU -- if that's a problem, pass an explicit, lower ncores
  # into the exec() calls inside filter_merge_batch/detect_native_crs/rasterize_batch_count
  # (e.g. total cores / n_workers each), or lower n_workers.
  n_workers = list(type = "integer", default = 6),

  # Number of LAS/LAZ tiles processed together in one batch. Splitting the full tile set into
  # batches this size keeps each batch's merged LAZ, GeoPackage, and native raster small and
  # bounded, no matter how many tiles input_dir holds overall. Lower this if a single batch's
  # GeoPackage is still too large, or if a batch's point count is still bumping into lasR's
  # 32-bit point-index ceiling (see stream_chunk_size below); raise it if batches come out small
  # enough that merging more tiles together per batch would still be manageable.
  batch_size_tiles = list(type = "integer", default = 50),

  # Layer name used inside every batch's GeoPackage.
  gpkg_layer = list(type = "character", default = "ground_points"),

  # If TRUE, each batch's ground-only points are also streamed out to their own GeoPackage (see
  # the "Export" section below). This is the slowest step in the per-batch loop and isn't needed
  # if all you want is the density raster, so it's off by default -- set TRUE to get a vector
  # point layer per batch alongside the raster.
  export_ground_points_gpkg = list(type = "logical", default = FALSE),

  # Target pixel size, in meters, for the (internal, per-batch) count raster and the density
  # raster derived from it.
  target_pixel_size_m = list(type = "numeric", default = 3),

  # Target CRS the final density raster is reprojected to.
  target_crs = list(type = "character", default = "EPSG:2856"),  # NAD83(HARN) / Washington South (Meters)

  # If TRUE, cells with no ground points get density/count = 0 instead of NA.
  fill_empty_cells_with_zero = list(type = "logical", default = FALSE),

  # If TRUE, every per-batch intermediate file (merged LAZ, per-batch GeoPackage, per-batch
  # native-CRS count raster) is deleted once the workflow finishes, keeping only the final
  # EPSG:2856 density raster (density_tif_2856_path). Set FALSE to keep the per-batch files
  # around -- e.g. to inspect them or resume a run without redoing batches (see
  # filter_merge_batch below). See the "Clean up intermediate files" section near the end of
  # this file.
  cleanup_intermediate_files = list(type = "logical", default = TRUE),

  # Spatial chunk size (native CRS units) used whenever we stream a batch's merged LAZ file.
  # lasR's internal point index is a 32-bit integer (max 2,147,483,647 points), so reading a
  # very large file in one unchunked pass can overflow it. Chunking keeps each pass well under
  # that ceiling and bounds memory use. This value doesn't need to be precise for correctness --
  # shrink it if a chunk is still dense enough to overflow (or shrink batch_size_tiles instead,
  # if it's a whole batch, not just a chunk within it, that's overflowing).
  stream_chunk_size = list(type = "integer", default = 100000)
)

# Applies param_specs on top of the raw (character) values read_param_file() returned: converts
# each recognized keyword to its declared type, falls back to its default when the keyword is
# absent from the file, and stops with a clear error if a required keyword (no default) is
# missing or a value can't be parsed as its declared type. Keywords in the file that aren't in
# param_specs are ignored with a warning, rather than silently accepted, so a typo'd keyword in
# the parameter file doesn't just vanish unnoticed.
build_params <- function(raw_params, specs, path) {
  unknown <- setdiff(names(raw_params), names(specs))
  if (length(unknown) > 0) {
    warning("Ignoring unrecognized parameter(s) in ", path, ": ", paste(unknown, collapse = ", "))
  }

  out <- list()
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    if (nm %in% names(raw_params)) {
      raw_val   <- raw_params[[nm]]
      converted <- switch(spec$type,
        character = as.character(raw_val),
        numeric   = suppressWarnings(as.numeric(raw_val)),
        integer   = suppressWarnings(as.integer(raw_val)),
        logical   = suppressWarnings(as.logical(raw_val)),
        stop("Internal error: unknown parameter type '", spec$type, "' for '", nm, "'")
      )
      if (spec$type != "character" && is.na(converted)) {
        stop("Could not parse parameter '", nm, "' (value \"", raw_val, "\") as ", spec$type,
             " -- check ", path, ".")
      }
      out[[nm]] <- converted
    } else if (!is.null(spec$default)) {
      out[[nm]] <- spec$default
    } else if (isTRUE(spec$required)) {
      stop("Required parameter '", nm, "' is missing from parameter file: ", path)
    }
  }
  out
}

raw_params <- read_param_file(config_path)
params     <- build_params(raw_params, param_specs, config_path)
list2env(params, envir = globalenv())  # makes input_dir, output_dir, n_workers, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)  # create output_dir (and any missing parent folders) if it doesn't exist yet

# Builds a batch-specific output path, e.g. batch_path("ground_merged.laz", 3) ->
# ".../b003_ground_merged.laz". Namespacing every intermediate by batch number keeps batches
# from overwriting each other's outputs, and lets an interrupted run resume without redoing
# batches that already finished (see filter_merge_batch below).
batch_path <- function(name, batch_id) {
  file.path(output_dir, sprintf("b%03d_%s", batch_id, name))
}

# Final, whole-site raster output -- assembled by mosaicking each batch's native count raster
# together once all batches finish, deriving density from that combined total, then
# reprojecting the density raster to EPSG:2856 (see the "Combine batches" and reprojection
# sections below). The native-CRS count/density mosaics are only kept in memory as
# intermediates on the way to this EPSG:2856 raster, not written to disk, and the count raster
# itself is never written out in any CRS -- this path is the only whole-site raster file the
# workflow produces.
# There is no single whole-site GeoPackage by design -- that's the file batching exists to keep
# from growing unmanageably large -- so, when export_ground_points_gpkg is TRUE, ground points
# stay split across the per-batch GeoPackages built in the "Export" section.
density_tif_2856_path <- file.path(output_dir, "ground_density_epsg2856.tif")

## ---- Discover and batch LAZ tiles --------------------------------------------------------
## Lists every input tile once, then splits that list into batch_size_tiles-sized groups.
## Everything from here on operates on one batch at a time.

las_files <- list.files(input_dir, pattern = "\\.(laz|las)$",
                        full.names = TRUE, ignore.case = TRUE)  # gather every LAS/LAZ tile in input_dir (case-insensitive extension match), with full paths

if (length(las_files) == 0) {
  stop("No LAS/LAZ files found in input_dir. Check the path.")  # fail fast with a clear message if the folder is empty or misconfigured
}

n_tiles   <- length(las_files)
n_batches <- ceiling(n_tiles / batch_size_tiles)
tile_batches <- split(las_files, ceiling(seq_len(n_tiles) / batch_size_tiles))  # list of character vectors, one per batch, in file-list order

message("Found ", n_tiles, " tiles; split into ", n_batches,
        " batch(es) of up to ", batch_size_tiles, " tiles each.")

## ---- Filter and merge LAZ -----------------------------------------------------------------
## Filters ground returns by tile, then merges one batch's tiles into a single per-batch file.
## This filtering always runs -- it's independent of export_ground_points_gpkg below, so every
## batch's merged LAZ is ground-only whether or not a GeoPackage is ever created for it. Defined
## here as a function; called once per batch in the loop below.

# Filters ground returns (Classification == 2) from one batch's tiles and merges them into a
# single per-batch LAZ file. Skips the work if that batch's merged file already exists, so a
# failed/interrupted run can be restarted without redoing batches that already completed.
filter_merge_batch <- function(las_files_batch, merged_path) {
  if (file.exists(merged_path)) {
    message("Merged file already exists, skipping filter/merge: ", merged_path)
    return(invisible(merged_path))
  }

  read_stage  <- reader(filter = keep_ground())  # keep_ground() filters at read time, not a separate stage
  write_stage <- write_las(merged_path)               # no "*" -> single merged file for this batch

  filter_merge_pipeline <- read_stage + write_stage   # chain the two stages into one lasR pipeline: read+filter, then write
  exec(filter_merge_pipeline, on = las_files_batch)    # run the pipeline over this batch's tiles only

  message("Ground-only merged point cloud written to: ", merged_path)
  invisible(merged_path)
}

## ---- Coordinate system ---------------------------------------------------------------------
## Detects the native CRS and units of a batch's merged LAZ file, needed for that batch's
## GeoPackage and count raster. Defined here as a function; called once per batch below.

# Detects the native CRS/units of a batch's merged LAZ file and derives the raster resolution
# (in native-CRS units) corresponding to target_pixel_size_m. Returned as a list so each
# batch's result can be checked for consistency with the others afterward -- all batches come
# from the same site/collection, so they should all report the same CRS, units, and resolution.
detect_native_crs <- function(merged_path) {
  probe_pipeline <- reader() + hulls()  # hulls() computes a hull polygon from the actual point coordinates -- it does read every point, but that's still far cheaper than a full rasterize() pass
  probe_hull <- exec(probe_pipeline, on = merged_path,
                     with = list(chunk = stream_chunk_size))  # stream in spatial chunks rather than reading the whole merged file at once -- see stream_chunk_size above

  native_crs   <- st_crs(probe_hull)         # CRS embedded in the LAZ file, carried through onto the hull
  native_units <- native_crs$units_gdal      # linear units of that CRS (e.g. "metre", "foot", "us survey foot")
  if (is.null(native_units)) {
    stop("Could not determine linear units of the LAS file's CRS for ", merged_path, ". ",
         "Check that the point cloud has a valid projected CRS assigned.")  # a missing or geographic (lat/long) CRS has no linear units, so bail out here
  }

  is_us_survey_feet <- grepl("us survey foot", native_units, ignore.case = TRUE)  # check the more specific "us survey foot" unit first...
  is_feet <- !is_us_survey_feet && grepl("foot|feet|ft", native_units, ignore.case = TRUE)  # ...then treat any other "foot"-like unit as an international foot

  if (is_us_survey_feet) {
    native_res <- target_pixel_size_m * (3937 / 1200)   # US survey foot: 1 metre = 3937/1200 US survey feet
  } else if (is_feet) {
    native_res <- target_pixel_size_m / 0.3048            # international foot: 1 ft = 0.3048 m exactly, so metres -> feet is division by 0.3048
  } else {
    native_res <- target_pixel_size_m                     # already metres, no conversion needed
  }

  list(crs = native_crs, units = native_units, res = native_res)
}

## ---- Export ground points to GeoPackage ----------------------------------------------------
## Optional (controlled by export_ground_points_gpkg in Setup). Streams one batch's merged
## ground-only points (Classification == 2) into that batch's own GeoPackage. Defined here as a
## function; called once per batch below, only when export_ground_points_gpkg is TRUE.

# NOTE: a batch's merged point cloud can still hold more points than lasR's internal point
# index can address in a single unchunked read (a 32-bit integer, capped at 2^31 - 1 =~ 2.1
# billion points). On top of that, a file this large can't be pulled into R as one data.frame/sf object
# anyway (it would need tens to hundreds of GB of RAM, and a GeoPackage layer with billions of
# point features isn't something GIS software can realistically open). So instead of reading
# everything at once, stream each batch's file in spatial chunks (via stream_chunk_size)
# and write each chunk straight to that batch's GeoPackage as it's read -- only ever holding
# one chunk in memory. Splitting the whole tile set into batches (batch_size_tiles) is the
# other half of this: it bounds how large any single GeoPackage can grow in the first place.
export_to_gpkg <- function(merged_path, gpkg_path, gpkg_layer, native_crs) {
  first_chunk  <- TRUE   # becomes FALSE after the first chunk is written, so later chunks append instead of overwriting
  total_points <- 0      # running total across all chunks (a plain double, so it can't overflow like a 32-bit integer could)
  class_tally  <- integer(0)  # running classification-code tally across all chunks, combined at the end

  write_chunk_to_gpkg <- function(data) {
    # `data` is one spatial chunk's worth of points (X, Y, Z, Classification) -- never the whole
    # point cloud. lasR calls this once per chunk as it streams through merged_path.
    if (nrow(data) == 0) return(nrow(data))  # lasR can call back with an empty chunk near tile edges; nothing to do

    tab <- table(data$Classification)   # tally classification codes seen in this chunk
    class_tally <<- c(class_tally, tab) # fold into the running tally (collapsed/checked below)

    # Sanity check: keep_ground() should mean everything here is class 2, but confirm rather
    # than assume -- check per chunk so a bad chunk fails fast instead of after a long run.
    if (!all(names(tab) == "2")) {
      stop("Merged ground file contains classes other than 2 -- check ",
           "that keep_ground() behaved as expected. File: ", merged_path)
    }

    chunk_sf_native <- st_as_sf(data, coords = c("X", "Y"), crs = native_crs)  # this chunk's points as native-CRS geometries; Z stays as an attribute column
    chunk_sf_2856   <- st_transform(chunk_sf_native, target_crs)              # reproject just this chunk to the target CRS

    st_write(chunk_sf_2856, gpkg_path, layer = gpkg_layer, driver = "GPKG",
             delete_layer = first_chunk, append = !first_chunk, quiet = TRUE)
    # the first chunk (re)creates this batch's layer; every later chunk appends onto it
    first_chunk  <<- FALSE
    total_points <<- total_points + nrow(data)

    nrow(data)  # return just a count, not the data, so exec()'s combined result stays tiny
  }

  cb       <- callback(write_chunk_to_gpkg, expose = "xyzc", no_las_update = TRUE)  # expose X, Y, Z, Classification per chunk; no_las_update = TRUE means we're only reading, not modifying the point cloud
  read_pipeline <- reader() + cb  # chain reader + callback into one pipeline

  exec(read_pipeline, on = merged_path,
       with = list(chunk = stream_chunk_size))  # stream through this batch's file in spatial chunks, writing to its GeoPackage as we go

  if (total_points == 0) {
    stop("No ground points found after filtering in ", merged_path,
         ". Is the point cloud already ground-classified?")
  }

  message("Read ", total_points, " ground points from ", merged_path, ".")
  message("Classification breakdown:")
  print(tapply(class_tally, names(class_tally), sum))  # collapse the per-chunk tallies into one breakdown for this batch

  message("Ground return points written to: ", gpkg_path,
          " (layer: ", gpkg_layer, ")")

  total_points
}

## ---- Create Point Density Raster -----------------------------------------------------------
## Ground points only. Rasterizes one batch's ground-only points using lasR's own streaming
## rasterize() stage (not terra::rasterize()), which tallies points per cell while reading and
## never needs the point cloud loaded into R -- necessary at this file's scale (see note in the
## Export chunk above). Only the raw point COUNT is computed per batch; density for the whole
## site is derived once, after all batches' count rasters are combined (see "Combine batches"
## below) -- computing density per batch first and then mosaicking density values would
## double-count or misweight cells that straddle a batch boundary, whereas summing raw counts
## first and dividing once is exact. The count raster is only ever an internal step toward that
## density value -- the workflow's only raster output is the density raster (see "Project
## Density Raster to EPSG:2856" below). Defined here as a function; called once per batch below.

# lasR::rasterize(res, "count") is a lasR pipeline stage: it counts points per cell directly
# inside the C++ streaming engine and writes the result straight to `ofile`, so no point data
# ever has to pass through R. Cells the point cloud doesn't cover come out as NA. Namespaced
# explicitly as lasR::rasterize() because terra (loaded after lasR in Setup) also exports a
# `rasterize` generic that masks lasR's -- calling the bare name here would instead dispatch to
# terra::rasterize() and fail with "unable to find an inherited method ... x = 'numeric', y =
# 'character'", since terra's rasterize() expects a SpatVector/SpatRaster, not a lasR pipeline
# stage's arguments.
rasterize_batch_count <- function(merged_path, native_res, count_tif_path) {
  count_stage    <- lasR::rasterize(native_res, "count", ofile = count_tif_path)
  count_pipeline <- reader() + count_stage
  exec(count_pipeline, on = merged_path,
       with = list(chunk = stream_chunk_size))  # stream in the same spatial chunks used for the GeoPackage export, for the same reason

  count_rast <- rast(count_tif_path)  # read the count raster lasR just wrote back in
  names(count_rast) <- "ground_point_count"
  message("Batch native-CRS ground count raster written to: ", count_tif_path)
  count_rast
}

## ---- Run the workflow across all batches ---------------------------------------------------
## Calls the four functions defined above once per batch, in order, using batch-specific paths.
## Batches run in parallel, up to n_workers at a time (set in Setup), via
## future.apply::future_lapply().

# Runs the full per-batch pipeline (filter/merge -> detect CRS -> optional GeoPackage export ->
# rasterize) for one batch. Called once per batch below via future_lapply(), so each call may
# run in its own worker process (see plan(multisession, ...) in Setup) rather than in this
# session. Returns a plain list of paths/CRS info rather than the count raster itself --
# terra::SpatRaster objects hold an external C++ pointer that doesn't survive being passed back
# from a multisession worker, so the raster is re-read from disk in the main session instead
# (see below), now that every worker has finished writing its .tif.
process_one_batch <- function(b) {
  library(lasR)   # each multisession worker is a separate R process, with none of this
  library(terra)  # session's packages loaded -- so the packages these functions need have to
  library(sf)     # be (re)loaded here, inside the worker.

  message("=== Batch ", b, " of ", n_batches, " (", length(tile_batches[[b]]), " tiles) ===")

  merged_path <- batch_path("ground_merged.laz", b)
  gpkg_path   <- batch_path("ground_points.gpkg", b)
  count_path  <- batch_path("ground_count_native_crs.tif", b)

  filter_merge_batch(tile_batches[[b]], merged_path)

  crs_info <- detect_native_crs(merged_path)

  if (export_ground_points_gpkg) {
    export_to_gpkg(merged_path, gpkg_path, gpkg_layer, crs_info$crs)
  }

  rasterize_batch_count(merged_path, crs_info$res, count_path)  # writes count_path; the SpatRaster this returns is discarded here -- see note above

  list(merged_path = merged_path, gpkg_path = gpkg_path, count_path = count_path,
       crs_info = crs_info)
}

# Cap workers at n_batches -- e.g. n_workers = 4 with only 2 batches would otherwise spin up 2
# worker processes that never get any work.
plan(multisession, workers = min(n_workers, n_batches))

# future_lapply() dispatches one process_one_batch() call per batch across the worker processes
# set up just above. future.seed = FALSE because nothing here draws random numbers -- it just
# silences future's "no parallel-safe RNG seed" advisory.
batch_results <- future_lapply(seq_along(tile_batches), process_one_batch, future.seed = FALSE)

# Batches are done with parallel workers -- everything from here runs sequentially in this
# session, so release the worker processes.
plan(sequential)

# Per-batch intermediate file paths, tracked here so they can be deleted afterward if
# cleanup_intermediate_files is TRUE (see the "Clean up intermediate files" section near the
# end of this file).
batch_merged_paths <- vapply(batch_results, function(r) r$merged_path, character(1))
batch_gpkg_paths   <- vapply(batch_results, function(r) r$gpkg_path, character(1))
batch_count_paths  <- vapply(batch_results, function(r) r$count_path, character(1))
batch_native_crs   <- lapply(batch_results, function(r) r$crs_info)  # checked for consistency below

# Re-read each batch's native-CRS count raster from disk now that every worker has finished --
# see the note on process_one_batch() above for why the SpatRaster itself isn't passed back
# directly from a worker.
batch_count_rasts <- lapply(batch_count_paths, function(p) {
  r <- rast(p)
  names(r) <- "ground_point_count"
  r
})

# All batches come from the same site/collection, so they should share one native CRS, units,
# and resolution -- confirm that rather than assume it, since a mismatch (e.g. mixed-CRS tiles
# slipping into input_dir) would otherwise silently corrupt the mosaic step below.
native_units_all <- vapply(batch_native_crs, function(x) x$units, character(1))
native_res_all   <- vapply(batch_native_crs, function(x) x$res, numeric(1))
if (length(unique(native_units_all)) > 1 || max(native_res_all) - min(native_res_all) > 1e-6) {
  stop("Batches do not share a consistent native CRS/resolution -- check for mixed-CRS tiles ",
       "in input_dir.")
}
native_crs <- batch_native_crs[[1]]$crs
native_res <- batch_native_crs[[1]]$res
message("Native-CRS resolution corresponding to ", target_pixel_size_m,
        " m pixels: ", round(native_res, 6), " (", native_units_all[1], ")")

## ---- Combine batches -----------------------------------------------------------------------
## Mosaics every batch's native-CRS count raster into one whole-site count raster, then derives
## density from that combined total. Both are kept in memory only -- the count raster is purely
## an intermediate on the way to the density value and is never written to disk in any CRS, and
## the density raster is only written to disk once, after reprojection, in the next step.

# Mosaic every batch's native count raster into one whole-site count raster. fun = "sum"
# (rather than the default "first") means that if any batches' tiles happen to overlap in
# space, their counts are added together instead of one silently overwriting the other; for
# the ordinary case of non-overlapping tiles this behaves just like a plain mosaic. This
# assumes every batch's raster shares the same origin/grid alignment, which holds here because
# lasR's rasterize() snaps cell edges to multiples of native_res regardless of a batch's
# extent, and native_res was confirmed identical across batches above.
count_rast <- if (length(batch_count_rasts) == 1) {
  batch_count_rasts[[1]]
} else {
  do.call(mosaic, c(batch_count_rasts, fun = "first"))
}
names(count_rast) <- "ground_point_count"

if (fill_empty_cells_with_zero) {
  count_rast <- subst(count_rast, NA, 0)  # replace "no points" cells (NA) with 0, per fill_empty_cells_with_zero
}
# Native-CRS count_rast is kept in memory only -- it's an intermediate on the way to
# density_rast below (see "Project Density Raster to EPSG:2856"), never a final output itself,
# so it's not written to disk.

n_na <- global(count_rast, fun = "isNA")[1, 1]  # number of cells with no ground points (NA), for reporting
n_total <- ncell(count_rast)                    # total number of cells in the raster
message(sprintf("%d of %d cells (%.1f%%) have no ground points.",
                 n_na, n_total, 100 * n_na / n_total))  # report the fraction of empty cells so coverage gaps are visible

# Normalized density (ground points per square meter), computed once from the combined counts
cell_area    <- target_pixel_size_m^2   # true physical cell area, in m^2 instead of native_res^2, which would be in native (e.g. square feet) units and give density in the wrong units
density_rast <- count_rast / cell_area  # convert combined counts to density (points per m^2) using the true physical cell area
names(density_rast) <- "ground_point_density"

# Likewise, native-CRS density_rast is kept in memory only -- see the note on count_rast above.
# count_rast itself is discarded after this point -- it was only ever a step toward density_rast
# and is never written to disk or reprojected.

## ---- Project Density Raster to EPSG:2856 -----------------------------------------------------
## Bilinear resampling is appropriate here since density is a continuous field. This is the
## workflow's only raster output -- no count raster is written, in any CRS.

density_rast_2856 <- project(density_rast, target_crs, method = "bilinear",
                             res = target_pixel_size_m)  # reproject the combined density raster; bilinear interpolation is appropriate for a continuous field
writeRaster(density_rast_2856, density_tif_2856_path, overwrite = TRUE)
message("Reprojected ground density raster written to: ", density_tif_2856_path)

## ---- Quick check to see if it worked ---------------------------------------------------------

print(density_rast_2856)  # inspect CRS, extent, and resolution of the reprojected density raster

# GeoPackage export is optional (export_ground_points_gpkg in Setup); only look for/sample them
# if that step actually ran.
if (export_ground_points_gpkg) {
  # There is no single whole-site GeoPackage (that's the whole reason tiles were batched --
  # see batch_size_tiles above), so list the per-batch GeoPackages and sample from the first one.
  # The full point layer within a batch is still too large to hold in R at once (that's why it
  # was streamed straight to the GeoPackage instead), so pull back only a handful of rows with a
  # SQL LIMIT rather than reading the layer in and then subsetting it.
  batch_gpkgs <- vapply(seq_len(n_batches), function(b) batch_path("ground_points.gpkg", b),
                        character(1))
  message("Batch GeoPackages written:")
  print(batch_gpkgs)

  ground_sample <- st_read(batch_gpkgs[1], layer = gpkg_layer,
                           query = paste0("SELECT * FROM \"", gpkg_layer, "\" LIMIT 5"),
                           quiet = TRUE)
  print(ground_sample)  # peek at a handful of reprojected ground points from the first batch
} else {
  message("export_ground_points_gpkg is FALSE -- no per-batch GeoPackages were written.")
}

## ---- Clean up intermediate files -------------------------------------------------------------
## Optional (controlled by cleanup_intermediate_files in Setup). Deletes every per-batch
## intermediate file -- merged LAZ, per-batch GeoPackage (if export_ground_points_gpkg was
## TRUE), per-batch native-CRS count raster, and any .lax spatial-index sidecar files lasR left
## behind in output_dir -- now that they've been folded into the final output. Only the final
## EPSG:2856 density raster (density_tif_2856_path) is kept; the native-CRS count/density
## mosaics were never written to disk in the first place (see "Combine batches" above), and no
## count raster is written out in any CRS.

if (cleanup_intermediate_files) {
  # lasR builds a .lax spatial-index sidecar (same basename, .lax extension) next to any
  # LAZ/LAS file it reads or writes -- including each batch's merged LAZ -- to speed up spatial
  # queries. Rather than assuming exactly one .lax per merged file, just glob output_dir for
  # any ".lax" so nothing gets left behind regardless of naming.
  lax_paths <- list.files(output_dir, pattern = "\\.lax$", full.names = TRUE, ignore.case = TRUE)

  intermediate_paths <- c(batch_merged_paths, batch_count_paths, lax_paths,
                          if (export_ground_points_gpkg) batch_gpkg_paths else character(0))
  intermediate_paths <- intermediate_paths[file.exists(intermediate_paths)]  # only try to delete files that actually exist

  if (length(intermediate_paths) == 0) {
    message("No intermediate files found to clean up.")
  } else {
    removed <- file.remove(intermediate_paths)
    message(sum(removed), " of ", length(intermediate_paths), " intermediate file(s) deleted.")
    if (!all(removed)) {
      warning("Could not delete: ", paste(intermediate_paths[!removed], collapse = ", "))
    }
  }
} else {
  message("cleanup_intermediate_files is FALSE -- per-batch intermediate files were left in ",
          output_dir)
}

## ---- Workflow timing -----------------------------------------------------------------------
## Reports the start time, end time, and total elapsed time for the entire workflow, using the
## workflow_start_time timestamp captured at the very top of the Setup chunk.

workflow_end_time <- Sys.time()
workflow_run_time <- workflow_end_time - workflow_start_time

message("Workflow start time: ", format(workflow_start_time, "%Y-%m-%d %H:%M:%S"))
message("Workflow end time:   ", format(workflow_end_time, "%Y-%m-%d %H:%M:%S"))
message("Total run time:      ", format(unclass(workflow_run_time), digits = 4),
        " ", units(workflow_run_time))
