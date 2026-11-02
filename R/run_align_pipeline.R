## run_align_pipeline.R
##
## Driver script: gets one reference and one to-be-aligned DTM and DSM (both DSMs required) -- each
## either mosaicked from a folder of GeoTIFF tiles, processed from one raster file, or read from
## an already-mosaicked raster -- then counts ground returns per reference DTM cell for two sets
## of LAZ/LAS tiles (a "reference" set and the set to be aligned to it -- the same per-set work
## run_ground_returns_batches.R does, here fixed to exactly these two named sets instead of an
## arbitrary list), and runs the Fortran program Align via TWutils (the same call align_dtms.R
## makes) on them, using the two count rasters as refGrnd/alignGrnd.
##
## Everything can be limited to one study-area polygon ([study_area] clip_polygon): only tiles
## intersecting it are used, and every raster built from tiles is masked to it. With
## clip_polygon NOFILE (the default), every tile in each tile folder is used, unmasked.
##
## Stage 1 (mosaic_dtm.R, up to x4): each of [ref_dtm], [align_dtm], [ref_dsm] and [align_dsm]
## that sets tile_dir or input_dtm/input_dsm is run as its own
## `Rscript mosaic_dtm.R <params.txt>` subprocess, one at a time -- see
## run_ground_returns_batches.R for why each run gets its own subprocess (isolated future::plan
## and global-environment state) rather than being source()'d into this session. The tiles
## (intersecting clip_polygon, if set) -- or the single input_dtm/input_dsm raster -- are
## mosaicked, converted to meters, reprojected to target_crs if they aren't already in it
## (keeping their own pixel size), masked to the polygon, and written to <output_dtm>.flt /
## <output_dsm>.flt. The reference DTM is built first, because Stage 2 counts on its grid; the
## align DTM and the DSMs are built in a background R session while Stage 2 runs (logging to
## <output>_mosaic_dtm.log), and Stage 3 waits for them.
## Alternatively, a section can name an already-mosaicked raster with existing_dtm/existing_dsm;
## that raster is used as-is (no mosaic_dtm.R run, no clipping/reprojection/unit conversion), so
## it must already be an ESRI GridFloat .flt in meters and in target_crs -- checked before
## anything runs. Both DSM sections are required, like the DTM sections: Align uses canopy
## height for its canopy classes, stable ground and the level of detection's open-ground rules.
##
## Stage 2 (ground_returns.R x2): each of [ref_set] and [align_set] below is run the same way,
## with template_raster set to the reference DTM, so the output grid -- extent, cell corners,
## and cellsize -- is fixed before any points are rasterized: each ground point is counted in
## the reference DTM cell containing it, giving ground_count_ref_grid.flt (points per reference
## DTM cell) inside each set's own output_dir.
##
## Stage 3 (Align): once both count rasters and all DTMs/DSMs exist, Align is run with
## refDTM/alignDTM/refDSM/alignDSM set to Stage 1's rasters and refGrnd/alignGrnd set to Stage
## 2's two count rasters. align_supervisor.R runs it (as its own Rscript subprocess): one Align
## run that, with tileNx/tileNy: auto, chooses the tile count and prior by cross-validation
## (TILES SELECT = AUTO), stops each loop once the aligned DEM has converged, fits the density
## LoD, writes the gradient raster LShunter needs, and moves the outputs to their final paths --
## see that script's header for its keywords and defaults.
##
## Stage 4 (strip_noise.R x2, only with a [strip_noise] section): each survey's noise measured
## where its own flight passes overlap (see strip_noise.R), on a sample of its LAZ tiles, classed by
## canopy height (that survey's DSM - DTM) and Align's gradient raster, with the reference DTM as
## the template grid and each set's pass_tally.csv (from Stage 2) when it exists. Outputs
## <label>_noise_*.csv in the section's output_dir (default: strip_noise\ in the run folder); a
## survey is skipped once its <label>_noise_classes.csv exists. It doesn't feed Align, so it runs
## after it -- even when Align was already done. Flags: --no-strip-noise runs Stages 1-3 only,
## --strip-noise-only Stage 4 only (run_landslide_pipelines.R uses both, to run Stage 4 beside
## LShunter).
##
## Resumability: every file this pipeline produces is checked for existence before it's
## (re)created, and skipped if it's already there -- same convention used throughout this repo
## (see ground_returns.R's per-batch merged-file skip and align_dtms.R's outOutlier skip). In
## particular, since outOutlier is the very last thing the whole pipeline produces, if its output
## already exists the entire pipeline -- both ground_returns.R runs included -- is skipped
## outright; see "Skip everything" below. Existence is the only check (except that a count raster
## older than the reference DTM stops the pipeline, since it may be on a stale grid): after
## changing clip_polygon (or any other parameter), delete the outputs it affects so they're
## rebuilt.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment),
## terra, future (for the background mosaics), plus whatever ground_returns.R and mosaic_dtm.R
## require (lasR, terra, sf, lidR, future, future.apply) in the R installation Rscript runs them
## with.
##
## ---- Pipeline config file --------------------------------------------------------------------
##
## All inputs for every stage come from one plain-text config file (see
## run_align_pipeline_params_template.txt for a ready-to-copy example), read from
## pipeline_config_path below or a command-line argument:
##   Rscript run_align_pipeline.R path/to/your_pipeline_config.txt
##
## The file uses "keyword: value" lines (same rules as ground_returns.R's and align_dtms.R's own
## parameter files -- order within a section doesn't matter, blank/"#"-only lines are ignored, a
## trailing "# comment" is stripped, only the first colon on a line splits keyword from value)
## grouped under these section headers:
##
##   [study_area] -- optional; clip_polygon (polygon file, any CRS; default NOFILE = no clipping)
##                   and target_crs (default EPSG:2856). Passed to every ground_returns.R and
##                   mosaic_dtm.R run, so they are set only here; the same keywords in any other
##                   section are ignored (with a warning). Also water_polygons (with water_layer
##                   for a multi-layer source such as a geodatabase, and water_where, an SQL
##                   condition that picks features): lakes, rasterized onto the reference DTM's
##                   grid as water_mask_ref_grid.flt in the run folder and given to Align as its
##                   WATER MASK. Not rivers: they move between surveys, so their polygons match
##                   neither survey's channel.
##   [defaults]   -- optional; ground_returns.R parameters shared by both [ref_set] and
##                   [align_set] below, unless a set repeats the same keyword to override it.
##                   extra_counts defaults to water=9 here: the surveys' water returns are
##                   counted and given to Align as WATER POINTS (extra_counts: NONE turns it off).
##   [ref_set]    -- ground_returns.R parameters (input_dir, output_dir, ...) for the reference
##                   LAZ tiles (pre-change).
##   [align_set]  -- ground_returns.R parameters for the LAZ tiles to be aligned (post-change).
##   [ref_dtm]    -- mosaic_dtm.R parameters (tile_dir, output_dtm, z_units, ...) for the
##                   reference DTM tiles (pre-change); or input_dtm (one raster file) with
##                   output_dtm, to clip/convert/reproject a single DTM the same way; or
##                   existing_dtm alone to use an already-mosaicked DTM as-is.
##   [align_dtm]  -- the same, for the DTM to be aligned (post-change).
##   [ref_dsm]    -- the same as [ref_dtm] but for the reference DSM, with the keywords
##                   output_dsm/input_dsm/existing_dsm. Required, like the DTMs.
##   [align_dsm]  -- the same, for the DSM to be aligned. Required.
##   [strip_noise] -- optional; strip_noise.R's keywords (noise_method, tile_fraction, n_workers,
##                   ...), passed on as they are, plus run (default TRUE), output_dir, ref_label and
##                   align_label (output prefixes, default ref/align). input_dir, template_raster,
##                   dtm, dsm, gradient_raster, pass_tally and clip_polygon are set here.
##   [align]      -- align_supervisor.R's keywords (iterations, convergence, k, tileNx, ...,
##                   outOutlier, outExceed, outLod, outCenDif, scratch_dir, executable_dir,
##                   program_name; see param_specs there), passed on as they are.
##                   refDTM/alignDTM, refDSM/alignDSM, and refGrnd/alignGrnd are NOT set here --
##                   they are always taken from the sections above; any found in [align] are
##                   ignored (with a warning).
##
## input_dir/output_dir (per LAZ set), tile_dir or input_dtm with output_dtm, or existing_dtm
## (per DTM), and outDTM/outbins/outDif/outOutlier/scratch_dir/executable_dir (in [align]) have
## no default and must be set. Everything else falls back to ground_returns.R's, mosaic_dtm.R's, or
## align_supervisor.R's own built-in default -- see param_specs in each of those scripts.

library(TWutils)

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so
# sibling scripts and template files are found by their location on disk rather than by the
# caller's working directory -- e.g. `Rscript R/run_align_pipeline.R ...` from the repo root and
# `Rscript run_align_pipeline.R ...` from inside R/ both correctly find R/ground_returns.R and
# R/run_align_pipeline_params_template.txt. Falls back to "." (assume cwd) when there's no
# "--file=" argument to read, e.g. when this script is source()'d from RStudio instead of run
# via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

pipeline_config_path <- file.path(get_script_dir(), "run_align_pipeline_params_template.txt")  # fallback; overridden by a command-line argument below

cli_args <- commandArgs(trailingOnly = TRUE)
no_strip_noise   <- "--no-strip-noise" %in% cli_args    # Stages 1-3 only
strip_noise_only <- "--strip-noise-only" %in% cli_args  # Stage 4 only
cli_args <- setdiff(cli_args, c("--no-strip-noise", "--strip-noise-only"))
if (length(cli_args) >= 1) pipeline_config_path <- cli_args[1]  # Rscript ... pipeline_config.txt takes precedence over the hardcoded fallback above

if (!file.exists(pipeline_config_path)) {
  stop("Pipeline config file not found: ", pipeline_config_path, ". Set pipeline_config_path ",
       "in this script, or pass its path as a command-line argument ",
       "(Rscript ... pipeline_config.txt).")
}

ground_returns_script <- file.path(get_script_dir(), "ground_returns.R")  # sibling script, same folder as this one
mosaic_dtm_script     <- file.path(get_script_dir(), "mosaic_dtm.R")      # likewise
supervisor_script     <- file.path(get_script_dir(), "align_supervisor.R")  # likewise
strip_noise_script    <- file.path(get_script_dir(), "strip_noise.R")       # likewise
param_dir <- file.path(tempdir(), "run_align_pipeline_params")  # generated ground_returns.R/mosaic_dtm.R parameter files are written here
dir.create(param_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Parse the pipeline config file -----------------------------------------------------------

# Reads a [section]-headed config file into a named list of sections, each itself a named list
# of "keyword: value" pairs (raw character strings -- build_params()/ground_returns.R's own
# build_params() do the actual type conversion). Unlike run_ground_returns_batches.R's [set]
# sections (which repeat, one per LAZ set), every section name here is expected exactly once.
read_sectioned_config <- function(path) {
  raw_lines <- readLines(path, warn = FALSE)
  raw_lines <- trimws(raw_lines)
  raw_lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]  # drop blank lines and comment-only lines

  sections <- list()
  active <- NULL              # NULL until the first section header
  current_keys <- character(0)  # keys seen so far in the section currently being filled, to catch duplicates within it

  for (line in raw_lines) {
    section_match <- regmatches(line, regexec("^\\[(.+)\\]$", line))[[1]]
    if (length(section_match) == 2) {
      section_name <- tolower(trimws(section_match[2]))
      if (section_name %in% names(sections)) {
        stop("Duplicate [", section_name, "] section in ", path)
      }
      sections[[section_name]] <- list()
      active <- section_name
      current_keys <- character(0)
      next
    }

    if (is.null(active)) {
      stop("Parameter line found before any [section] header in ", path, ": ", line)
    }

    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the FIRST colon splits keyword from value -- keeps values like "c:\..." or "EPSG:2856" intact
    if (colon_pos < 0) {
      stop("Malformed line in ", path, " (expected \"keyword: value\" or a [section] header): ",
           line)
    }
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(substr(line, colon_pos + 1, nchar(line)))
    val <- trimws(sub("\\s+#.*$", "", val))  # strip a trailing " # comment", if any, from the value

    if (key %in% current_keys) {
      stop("Duplicate keyword '", key, "' within [", active, "] of ", path)
    }
    current_keys <- c(current_keys, key)
    sections[[active]][[key]] <- val
  }

  sections
}

config <- read_sectioned_config(pipeline_config_path)

ground_returns_defaults <- config[["defaults"]]
if (is.null(ground_returns_defaults)) ground_returns_defaults <- list()
# Water handling is on by default: both surveys' lidar water returns (LAS class 9) are counted with
# the ground returns and left out of Align as WATER POINTS. extra_counts: NONE in [defaults] (or a
# set) turns it off; a survey without water points simply marks nothing.
if (is.null(ground_returns_defaults$extra_counts)) ground_returns_defaults$extra_counts <- "water=9"

study_area <- config[["study_area"]]
if (is.null(study_area)) study_area <- list()
ref_set    <- config[["ref_set"]]
align_set  <- config[["align_set"]]
ref_dtm    <- config[["ref_dtm"]]
align_dtm  <- config[["align_dtm"]]
ref_dsm    <- config[["ref_dsm"]]    # required, like the DTMs: Align needs both surveys' DSMs
align_dsm  <- config[["align_dsm"]]
align_section <- config[["align"]]
strip_noise_section <- config[["strip_noise"]]  # optional (Stage 4)

required_sections <- c("ref_set", "align_set", "ref_dtm", "align_dtm", "ref_dsm", "align_dsm", "align")
missing_sections  <- setdiff(required_sections, names(config))
if (length(missing_sections) > 0) {
  stop("Pipeline config file is missing section(s) ",
       paste0("[", missing_sections, "]", collapse = ", "), ": ", pipeline_config_path)
}

## ---- Study area ------------------------------------------------------------------------------
## clip_polygon and target_crs come only from [study_area] and are passed on to every
## ground_returns.R/mosaic_dtm.R run, so every product of the pipeline covers the same area in
## the same CRS. clip_polygon NOFILE (the default) uses every tile in each tile folder, unmasked.

if (is.null(study_area$clip_polygon) || !nzchar(study_area$clip_polygon)) {
  study_area$clip_polygon <- "NOFILE"
}
if (toupper(study_area$clip_polygon) == "NOFILE") {
  message("No clip_polygon -- every LAZ/DTM/DSM tile in each tile folder is used, unmasked.")
} else if (!file.exists(study_area$clip_polygon)) {
  stop("clip_polygon not found: ", study_area$clip_polygon)
}
if (is.null(study_area$target_crs)) study_area$target_crs <- "EPSG:2856"
# Optional water polygons (any vector format sf reads, any CRS; water_layer names the layer of a
# multi-layer source such as a geodatabase; water_where, an SQL condition, picks features, e.g.
# "FCode = 39000"), rasterized onto the reference DTM's grid for Align's WATER MASK (see "Water
# mask" below). Lakes only: rivers move between surveys, so river polygons match neither survey's
# channel -- Align finds rivers from the surveys' own water returns (WATER POINTS) instead.
water_polygons <- if (is.null(study_area$water_polygons)) "NOFILE" else study_area$water_polygons
water_layer    <- study_area$water_layer
water_where    <- study_area$water_where
if (toupper(water_polygons) != "NOFILE" && !file.exists(water_polygons)) {
  stop("water_polygons not found: ", water_polygons)
}
unknown_study_keys <- setdiff(names(study_area), c("clip_polygon", "target_crs", "water_polygons",
                                                   "water_layer", "water_where"))
if (length(unknown_study_keys) > 0) {
  warning("Ignoring unrecognized parameter(s) in [study_area] of ", pipeline_config_path, ": ",
          paste(unknown_study_keys, collapse = ", "))
}
study_area_params <- study_area[c("clip_polygon", "target_crs")]

# Drops clip_polygon/target_crs from one section (with a warning) so [study_area]'s values are
# the only ones ever used.
drop_study_area_keys <- function(section, section_name) {
  overridden <- intersect(names(section), names(study_area_params))
  if (length(overridden) > 0) {
    warning(paste(overridden, collapse = "/"), " in [", section_name, "] ignored -- set in ",
            "[study_area] instead.")
    section[overridden] <- NULL
  }
  section
}
ground_returns_defaults <- drop_study_area_keys(ground_returns_defaults, "defaults")
ref_set   <- drop_study_area_keys(ref_set,   "ref_set")
align_set <- drop_study_area_keys(align_set, "align_set")
ref_dtm   <- drop_study_area_keys(ref_dtm,   "ref_dtm")
align_dtm <- drop_study_area_keys(align_dtm, "align_dtm")
ref_dsm   <- drop_study_area_keys(ref_dsm,   "ref_dsm")
align_dsm <- drop_study_area_keys(align_dsm, "align_dsm")

# template_raster is always the reference DTM (set in Stage 2 below), never a user parameter.
drop_template_key <- function(section, section_name) {
  if (!is.null(section$template_raster)) {
    warning("template_raster in [", section_name, "] ignored -- ground points are always ",
            "counted on the reference DTM's grid.")
    section$template_raster <- NULL
  }
  section
}
ground_returns_defaults <- drop_template_key(ground_returns_defaults, "defaults")
ref_set   <- drop_template_key(ref_set,   "ref_set")
align_set <- drop_template_key(align_set, "align_set")

ground_returns_defaults <- modifyList(ground_returns_defaults, study_area_params)

# <output>.flt -- the file mosaic_dtm.R writes (it accepts its output path with or without ".flt").
flt_path <- function(output) paste0(sub("\\.flt$", "", output, ignore.case = TRUE), ".flt")

# Checks an already-mosaicked raster given as existing_dtm/existing_dsm and returns its .flt
# path. It's handed to Align as-is -- not clipped, unit-converted, or reprojected -- so it must
# already be an ESRI GridFloat (.flt with an ESRI-keyword .hdr, as mosaic_dtm.R writes) in meters
# and in target_crs. kind is "dtm" or "dsm".
check_existing_raster <- function(existing, section_name, kind) {
  key <- paste0("existing_", kind)
  alt <- paste0("tile_dir/output_", kind)
  flt <- flt_path(existing)
  if (!file.exists(flt)) {
    stop(key, " in [", section_name, "] not found: ", flt)
  }
  hdr <- sub("\\.flt$", ".hdr", flt, ignore.case = TRUE)
  if (!file.exists(hdr) ||
      !any(grepl("^\\s*ncols\\s", readLines(hdr, warn = FALSE), ignore.case = TRUE))) {
    stop(key, " in [", section_name, "] has no ESRI GridFloat header (ncols, xllcorner, ",
         "cellsize, ...) at ", hdr, " -- Align can't read it. Build it with ", alt, " instead.")
  }
  raster_crs <- terra::crs(terra::rast(flt))
  if (raster_crs == "") {
    warning(key, " in [", section_name, "] has no CRS (.prj) -- assuming it is already in ",
            study_area_params$target_crs, ": ", flt)
  } else if (!terra::same.crs(raster_crs, study_area_params$target_crs)) {
    stop(key, " in [", section_name, "] is not in target_crs ", study_area_params$target_crs,
         ": ", flt, ". Reproject it, or build it with ", alt, " instead.")
  }
  flt
}

# Each DTM/DSM section gets its raster from exactly one of:
#   existing_dtm/existing_dsm -- an already-mosaicked raster, handed to Align as-is;
#   input_dtm/input_dsm       -- one raster file, run through mosaic_dtm.R (as its input_file) to
#                                be clipped, converted to meters, and reprojected into
#                                output_dtm/output_dsm;
#   tile_dir                  -- a folder of tiles, mosaicked by mosaic_dtm.R into
#                                output_dtm/output_dsm.
# Returns the .flt path Align will read as path, whether Stage 1 must run mosaic_dtm.R for it as
# from_tiles, and if so that run's mosaic_dtm.R parameters (output_dsm is passed on as
# mosaic_dtm.R's output_dtm, input_dtm/input_dsm as its input_file). A source set to NOFILE
# counts as not set.
resolve_raster_section <- function(section, section_name, kind) {
  existing_key <- paste0("existing_", kind)
  input_key    <- paste0("input_", kind)
  output_key   <- paste0("output_", kind)

  source_keys <- c(existing_key, input_key, "tile_dir")
  sources <- Filter(function(k) !is.null(section[[k]]) && nzchar(section[[k]]) &&
                      toupper(section[[k]]) != "NOFILE",
                    source_keys)
  section[setdiff(intersect(source_keys, names(section)), sources)] <- NULL  # drop NOFILE sources

  if (length(sources) > 1) {
    stop("[", section_name, "] sets ", paste(sources, collapse = " and "), " -- set only one of ",
         existing_key, ", ", input_key, ", or tile_dir: ", pipeline_config_path)
  }

  if (identical(sources, existing_key)) {
    ignored <- intersect(names(section), c(output_key, "z_units", "file_pattern",
                                           "mosaic_fun", "clip_buffer_cells"))
    if (length(ignored) > 0) {
      warning(paste(ignored, collapse = "/"), " in [", section_name, "] ignored -- ",
              existing_key, " is set, so nothing is mosaicked for it.")
    }
    return(list(from_tiles = FALSE,
                path = check_existing_raster(section[[existing_key]], section_name, kind)))
  }

  if (length(sources) == 0 || is.null(section[[output_key]])) {
    stop("[", section_name, "] must set either ", existing_key, ", or ", input_key, " or ",
         "tile_dir together with ", output_key, ": ", pipeline_config_path)
  }

  params <- section
  if (identical(sources, input_key)) {
    if (!file.exists(section[[input_key]])) {
      stop(input_key, " in [", section_name, "] not found: ", section[[input_key]])
    }
    # Stage 1 skips mosaic_dtm.R whenever output_<kind> exists, so an output naming the input
    # itself would hand the unprocessed input to Align.
    if (tolower(normalizePath(section[[input_key]], winslash = "/", mustWork = FALSE)) ==
        tolower(normalizePath(flt_path(section[[output_key]]), winslash = "/", mustWork = FALSE))) {
      stop(output_key, " in [", section_name, "] must differ from ", input_key, " (",
           section[[input_key]], ") -- it names the processed copy written from it.")
    }
    ignored <- intersect(names(section), c("file_pattern", "mosaic_fun"))
    if (length(ignored) > 0) {
      warning(paste(ignored, collapse = "/"), " in [", section_name, "] ignored -- ", input_key,
              " is a single raster, so there are no tiles to match or combine.")
      params[ignored] <- NULL
    }
    params[[input_key]] <- NULL
    params$input_file <- section[[input_key]]
  }
  params[[output_key]] <- NULL
  params$output_dtm <- section[[output_key]]
  list(from_tiles = TRUE, path = flt_path(section[[output_key]]),
       params = modifyList(params, study_area_params))
}

ref_dtm   <- resolve_raster_section(ref_dtm,   "ref_dtm",   "dtm")
align_dtm <- resolve_raster_section(align_dtm, "align_dtm", "dtm")
ref_dsm   <- resolve_raster_section(ref_dsm,   "ref_dsm",   "dsm")
align_dsm <- resolve_raster_section(align_dsm, "align_dsm", "dsm")

## ---- Resolve Align's parameters up front -----------------------------------------------------
## Done before Stages 1-2 run (not after, as in align_dtms.R) so outOutlier is known early
## enough to decide, below, whether the whole pipeline can be skipped -- refGrnd/alignGrnd and
## refDTM/alignDTM/refDSM/alignDSM only need ref_set/align_set's output_dir and the DTM/DSM
## sections' output_dtm/output_dsm (or existing_dtm/existing_dsm), which are already known from
## the config without actually having run anything yet.

derived_align_keys <- intersect(names(align_section),
                                c("refGrnd", "alignGrnd", "refDTM", "alignDTM", "refDSM",
                                  "alignDSM", "refWater", "alignWater", "waterMask"))
if (length(derived_align_keys) > 0) {
  warning(paste(derived_align_keys, collapse = "/"), " in [align] ignored -- refGrnd/alignGrnd ",
          "(and refWater/alignWater) are always Stage 2's count rasters (on the reference DTM ",
          "grid), waterMask is [study_area]'s water_polygons on that grid, and ",
          "refDTM/alignDTM/refDSM/alignDSM come from ",
          "[ref_dtm]/[align_dtm]/[ref_dsm]/[align_dsm].")
  align_section[derived_align_keys] <- NULL
}

# [align] goes to align_supervisor.R as it is: that script holds the keywords and their
# defaults and checks them. Only the keywords this script itself needs are checked here.
for (key in c("outDTM", "outbins", "outDif", "outOutlier", "scratch_dir", "executable_dir")) {
  if (is.null(align_section[[key]]) || !nzchar(align_section[[key]])) {
    stop("Required parameter '", key, "' is missing from [align] in ", pipeline_config_path)
  }
}
align_params <- align_section

# Expected paths only -- none of these files needs to exist yet at this point; whether they do
# is checked in Stages 1-2 below (which skip (re)building each output once it already exists).
# With template_raster set, ground_returns.R always writes its count raster to
# ground_count_ref_grid.flt in the set's output_dir.
align_params$refGrnd   <- file.path(ref_set$output_dir,   "ground_count_ref_grid.flt")
align_params$alignGrnd <- file.path(align_set$output_dir, "ground_count_ref_grid.flt")
align_params$refDTM    <- ref_dtm$path
align_params$alignDTM  <- align_dtm$path
align_params$refDSM    <- ref_dsm$path
align_params$alignDSM  <- align_dsm$path

# Water-return counts: a set whose ground_returns.R run counts "water" (extra_counts, e.g.
# water=9) writes water_count_ref_grid.flt beside its ground counts, and Align's WATER POINTS
# reads it.
counts_water <- function(set_params) {
  ec <- modifyList(ground_returns_defaults, set_params)$extra_counts
  !is.null(ec) && any(grepl("^\\s*water\\s*=", strsplit(ec, ",", fixed = TRUE)[[1]]))
}
water_count_path <- function(set_params) file.path(set_params$output_dir, "water_count_ref_grid.flt")
if (counts_water(ref_set))   align_params$refWater   <- water_count_path(ref_set)
if (counts_water(align_set)) align_params$alignWater <- water_count_path(align_set)
# [study_area]'s water polygons, rasterized onto the reference DTM's grid in the run folder
water_mask_flt <- if (toupper(water_polygons) == "NOFILE") NULL else
  file.path(dirname(align_params$outOutlier), "water_mask_ref_grid.flt")
if (!is.null(water_mask_flt)) align_params$waterMask <- water_mask_flt

refGrnd    <- align_params$refGrnd
alignGrnd  <- align_params$alignGrnd
outDTM     <- align_params$outDTM
outOutlier <- sub("\\.flt$", "", align_params$outOutlier, ignore.case = TRUE)

## ---- Skip everything if the final output already exists --------------------------------------
## outOutlier is the last file the whole pipeline produces (every ground_returns.R/mosaic_dtm.R
## run feeds into it via Align). If it's already there, there is nothing left for any stage to
## do, so skip Stages 1-3 entirely rather than just checking right before the Align call, as
## align_dtms.R does -- that would still go through every (potentially slow) earlier run for no
## reason once the pipeline has already finished.

outOutlier_path <- paste0(outOutlier, ".flt")

# Writes a named list out as a "keyword: value" parameter file and returns that file's path.
write_param_file <- function(params, path) {
  lines <- vapply(names(params), function(nm) paste0(nm, ": ", params[[nm]]), character(1))
  writeLines(lines, path)
  path
}

# Runs one sibling script (ground_returns.R, mosaic_dtm.R, ...) as its own Rscript subprocess on
# an already-written parameter file and returns the path to the output it's expected to produce.
# Skips the Rscript call entirely if that output already exists. label is just for messages.
# With log_file set, the subprocess's console output goes to that file instead of this
# session's -- needed for the background mosaics (see Stage 1), whose output would otherwise
# be lost.
run_script_step <- function(script, label, param_path, expected_output, log_file = NULL) {
  if (all(file.exists(expected_output))) {
    message("\n==== ", label, ": output already exists, skipping ", basename(script), ": ",
            paste(expected_output, collapse = ", "), " ====")
    return(expected_output[1])
  }

  message("\n==== ", label, " (", basename(script), ") ====",
          if (!is.null(log_file)) paste0("\n  log: ", log_file))

  start_time <- Sys.time()
  status <- if (is.null(log_file)) {
    system2("Rscript", args = shQuote(c(script, param_path)))
  } else {
    system2("Rscript", args = shQuote(c(script, param_path)),
            stdout = log_file, stderr = log_file)
  }
  elapsed <- Sys.time() - start_time

  if (status != 0) {
    stop(basename(script), " failed for ", label, " (parameters: ", param_path,
         ") -- exit status ", status, ".",
         if (!is.null(log_file) && file.exists(log_file)) {
           paste0(" Last lines of its log (", log_file, "):\n",
                  paste(tail(readLines(log_file, warn = FALSE), 10), collapse = "\n"))
         })
  }

  message(label, " finished in ", format(unclass(elapsed), digits = 4), " ", units(elapsed))

  if (!all(file.exists(expected_output))) {
    stop(basename(script), " reported success for ", label, " but expected output is missing: ",
         paste(expected_output[!file.exists(expected_output)], collapse = ", "))
  }
  expected_output[1]
}

if (strip_noise_only) {
  message("--strip-noise-only: Stages 1-3 skipped.")
} else if (file.exists(outOutlier_path)) {
  message("Final outlier output already exists -- pipeline already complete, nothing to do: ",
          outOutlier_path)
} else {

  ## ---- Stage 1: mosaic_dtm.R for each DTM/DSM built from tiles or one input file --------------
  ## One Rscript subprocess per raster (see run_ground_returns_batches.R for why each run gets
  ## its own subprocess); mosaic_dtm.R handles DSM tiles exactly like DTM tiles, and an
  ## input_dtm/input_dsm raster as a single tile. Each raster is skipped (not rebuilt) if its
  ## .flt already exists, and one given as existing_dtm/existing_dsm -- or a DSM that is NOFILE
  ## -- is never mosaicked at all.
  ##
  ## Only the reference DTM is built here, before anything else: Stage 2 counts ground points on
  ## its grid. The align DTM and both DSMs aren't needed until Align (Stage 3), so they're built
  ## in a background R session (one after another) while Stage 2 counts in the foreground, and
  ## Stage 3 waits for them. Each background mosaic writes its console output to
  ## <output>_mosaic_dtm.log next to its output, since a background session's subprocess output
  ## isn't relayed back here; their "====" messages appear only once they've all finished.

  run_mosaic <- function(label, raster, log_file = NULL) {
    if (identical(raster$path, "NOFILE")) {
      message("\n==== ", label, ": none (NOFILE) ====")
      return(raster$path)
    }
    if (!raster$from_tiles) {
      message("\n==== ", label, ": using existing raster, skipping mosaic_dtm.R: ", raster$path,
              " ====")
      return(raster$path)
    }
    param_path <- write_param_file(raster$params,
                                   file.path(param_dir, paste0("params_", label, ".txt")))
    run_script_step(mosaic_dtm_script, label, param_path, raster$path, log_file)
  }

  ref_dtm_flt <- run_mosaic("ref_dtm", ref_dtm)

  later_rasters <- list(align_dtm = align_dtm, ref_dsm = ref_dsm, align_dsm = align_dsm)
  needs_mosaic  <- vapply(later_rasters, function(r) isTRUE(r$from_tiles) && !file.exists(r$path),
                          logical(1))
  mosaic_log <- function(raster) paste0(sub("\\.flt$", "", raster$path, ignore.case = TRUE),
                                        "_mosaic_dtm.log")

  if (any(needs_mosaic)) {
    message("\n==== Starting ", paste(names(later_rasters)[needs_mosaic], collapse = ", "),
            " in a background R session (mosaic_dtm.R), to run during Stage 2 ====")
    for (label in names(later_rasters)[needs_mosaic]) {
      message("  ", label, " log: ", mosaic_log(later_rasters[[label]]))
    }
    future::plan(future::multisession, workers = 1)
    background_mosaics <- future::future({
      lapply(names(later_rasters), function(label) {
        raster <- later_rasters[[label]]
        run_mosaic(label, raster, log_file = if (isTRUE(raster$from_tiles)) mosaic_log(raster))
      })
    }, seed = TRUE)
  } else {
    background_mosaics <- NULL
    later_flts <- lapply(names(later_rasters),
                         function(label) run_mosaic(label, later_rasters[[label]]))
  }

  message("\nStage 1: reference DTM ready: ", ref_dtm_flt)

  ## ---- Stage 2: ground_returns.R, once each for ref_set and align_set ------------------------
  ## Sequential, one Rscript subprocess at a time -- parallelism happens *within* each run,
  ## across that run's own batches, via that run's own n_workers. Each run gets
  ## template_raster = the reference DTM, so its output is the ground-point count per reference
  ## DTM cell (ground_count_ref_grid.flt), counted directly on that grid: same extent, cell
  ## corners, and cellsize. Each set's count raster is skipped (not regenerated) if it already
  ## exists -- unless it's older than the reference DTM, in which case the DTM was rebuilt after
  ## the points were counted on it and the counts may be on a stale grid, so this stops rather
  ## than hand them to Align.

  # Merges one set's own parameters over ground_returns_defaults (set-specific values win), adds
  # the reference DTM as template_raster, writes the result out as a "keyword: value"
  # ground_returns.R parameter file, and returns that file's path. Same helper as
  # write_param_file() in run_ground_returns_batches.R.
  write_ground_returns_param_file <- function(set_params, defaults, path) {
    params <- modifyList(defaults, set_params)
    if (is.null(params$input_dir) || is.null(params$output_dir)) {
      stop("Both [ref_set] and [align_set] must set input_dir and output_dir.")
    }
    params$template_raster <- ref_dtm_flt
    write_param_file(params, path)
  }

  # Runs ground_returns.R for one named set and returns the path to the count raster it produces.
  # With water counts asked for, the set counts as done only once the water raster exists too.
  run_ground_returns_for_set <- function(label, set_params, defaults, expected_count_flt) {
    expected <- c(expected_count_flt, if (counts_water(set_params)) water_count_path(set_params))
    for (f in expected[file.exists(expected)]) {
      if (file.mtime(f) < file.mtime(ref_dtm_flt)) {
        stop(label, " count raster (", f, ") is older than the reference DTM (", ref_dtm_flt,
             "), so it may have been counted on a different grid. Delete it and rerun this ",
             "pipeline to recount.")
      }
    }
    param_path <- write_ground_returns_param_file(set_params, defaults,
                                                  file.path(param_dir, paste0("params_", label, ".txt")))
    run_script_step(ground_returns_script, label, param_path, expected)
  }

  ## ---- Water mask: [study_area]'s water polygons on the reference DTM's grid -----------------
  ## Rasterized with terra (cells whose centre falls in a polygon), masked to the DTM's data, and
  ## written as an ESRI GridFloat (1 water, 0 not) for Align's WATER MASK. Rebuilt when older than
  ## the reference DTM or the polygon file.
  if (!is.null(water_mask_flt)) {
    # The source (file, layer, condition) is recorded beside the mask, so changing any of them
    # rebuilds it, as a newer DTM or polygon file does
    source_txt <- sub("\\.flt$", "_source.txt", water_mask_flt, ignore.case = TRUE)
    water_source <- paste(c(normalizePath(water_polygons, winslash = "/"), water_layer %||% "",
                            water_where %||% ""), collapse = " | ")
    stale <- !file.exists(water_mask_flt) || !file.exists(source_txt) ||
      !identical(readLines(source_txt, warn = FALSE), water_source) ||
      file.mtime(water_mask_flt) < max(file.mtime(c(ref_dtm_flt, water_polygons)))
    if (stale) {
      message("\n==== Water mask: rasterizing ", water_polygons,
              if (!is.null(water_layer)) paste0(" (layer ", water_layer, ")"),
              if (!is.null(water_where)) paste0(" where ", water_where), " ====")
      template <- terra::rast(ref_dtm_flt)
      if (terra::crs(template) == "") terra::crs(template) <- study_area_params$target_crs
      layer <- if (is.null(water_layer)) sf::st_layers(water_polygons)$name[1] else water_layer
      polys <- if (is.null(water_where)) sf::st_read(water_polygons, layer, quiet = TRUE) else
        sf::st_read(water_polygons, query = paste0("SELECT * FROM \"", layer, "\" WHERE ", water_where),
                    quiet = TRUE)
      if (nrow(polys) == 0) stop("No water polygons in ", water_polygons, " (layer ", layer, ")",
                                 if (!is.null(water_where)) paste0(" where ", water_where), ".")
      polys <- sf::st_zm(sf::st_make_valid(sf::st_transform(polys, terra::crs(template))))
      wm <- terra::rasterize(terra::vect(polys), template, field = 1, background = 0)
      wm <- terra::mask(wm, template)
      n_water <- terra::global(wm == 1, "sum", na.rm = TRUE)[1, 1]
      # Written under a "_partial" name and renamed into place, the .flt last, as ground_returns.R
      # and mosaic_dtm.R do, so an interrupted write never leaves a file that looks finished
      dir.create(dirname(water_mask_flt), showWarnings = FALSE, recursive = TRUE)
      out_base     <- sub("\\.flt$", "", water_mask_flt, ignore.case = TRUE)
      partial_base <- paste0(out_base, "_partial")
      raster_exts  <- c(".flt", ".hdr", ".prj", ".stx", ".flt.aux.xml")
      unlink(paste0(c(out_base, partial_base), rep(raster_exts, each = 2)))
      wm <- terra::writeRaster(wm, paste0(partial_base, ".flt"), filetype = "EHdr",
                               datatype = "FLT4S", NAflag = -9999, gdal = "BYTEORDER=LSB",
                               overwrite = TRUE)
      e <- terra::ext(wm)
      writeLines(sprintf("%-14s %s",  # ESRI GridFloat header for the Fortran grid reader (see mosaic_dtm.R)
                         c("ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "NODATA_value",
                           "byteorder"),
                         c(terra::ncol(wm), terra::nrow(wm), sprintf("%.6f", terra::xmin(e)),
                           sprintf("%.6f", terra::ymin(e)), sprintf("%.12f", terra::res(wm)[1]),
                           "-9999", "LSBFIRST")),
                 paste0(partial_base, ".hdr"))
      rm(wm); invisible(gc())
      for (x in c(setdiff(raster_exts, ".flt"), ".flt")) {
        if (file.exists(paste0(partial_base, x)) && !file.rename(paste0(partial_base, x), paste0(out_base, x))) {
          stop("Could not rename ", partial_base, x, " to ", out_base, x)
        }
      }
      writeLines(water_source, source_txt)
      message("Water mask written to ", water_mask_flt, ": ",
              round(n_water * prod(terra::res(template)) / 1e6, 3), " km2 of water")
    } else {
      message("\n==== Water mask: already up to date, ", water_mask_flt, " ====")
    }
  }

  # If counting fails, wait for the background mosaics to finish before stopping: stopping now
  # would end the background session and could leave its current mosaic unfinished.
  tryCatch({
    run_ground_returns_for_set("ref_set",   ref_set,   ground_returns_defaults, refGrnd)
    run_ground_returns_for_set("align_set", align_set, ground_returns_defaults, alignGrnd)
  }, error = function(e) {
    if (!is.null(background_mosaics)) {
      message("\nStage 2 failed -- waiting for the background mosaics to finish before stopping.")
      try(future::value(background_mosaics), silent = TRUE)
      future::plan(future::sequential)
    }
    stop(e)
  })

  message("\nStage 2 complete. Ground-point count rasters on the reference DTM grid:")
  message("  ref_set:   ", refGrnd)
  message("  align_set: ", alignGrnd)

  if (!is.null(background_mosaics)) {
    message("\n==== Waiting for the background mosaics (",
            paste(names(later_rasters)[needs_mosaic], collapse = ", "), ") ====")
    later_flts <- future::value(background_mosaics)  # re-raises a background mosaic's error, if any
    future::plan(future::sequential)
  }
  align_dtm_flt <- later_flts[[1]]
  ref_dsm_flt   <- later_flts[[2]]
  align_dsm_flt <- later_flts[[3]]

  message("\nAll DTMs/DSMs ready:")
  message("  ref_dtm:   ", ref_dtm_flt)
  message("  align_dtm: ", align_dtm_flt)
  message("  ref_dsm:   ", ref_dsm_flt)
  message("  align_dsm: ", align_dsm_flt)

  ## ---- Stage 3: Align --------------------------------------------------------------------------
  ## Reaching here means outOutlier_path did not exist when this run started, so Align has not
  ## already produced this pipeline's final output -- run it.

  # align_supervisor.R reads the same keywords as [align], plus the six raster paths.
  supervisor_params <- write_param_file(align_params,
                                        file.path(param_dir, "params_align_supervisor.txt"))
  run_script_step(supervisor_script, "Align", supervisor_params, outOutlier_path)
  message("Align finished. Aligned DTM written to: ", outDTM)
}

## ---- Stage 4: strip_noise.R, once per survey ---------------------------------------------------
## Each survey's own flight-pass overlaps, classed by its own canopy height and Align's gradient.
if (!is.null(strip_noise_section) && !no_strip_noise &&
    !identical(toupper(strip_noise_section$run %||% "TRUE"), "FALSE")) {
  noise_dir <- strip_noise_section$output_dir %||% file.path(dirname(align_params$outOutlier), "strip_noise")
  noise_labels <- c(ref = strip_noise_section$ref_label %||% "ref",
                    align = strip_noise_section$align_label %||% "align")
  set_keys <- c("input_dir", "output_dir", "label", "template_raster", "dtm", "dsm", "gradient_raster",
                "pass_tally", "clip_polygon")
  ignored <- intersect(names(strip_noise_section), set_keys[-2])
  if (length(ignored)) warning(paste(ignored, collapse = "/"), " in [strip_noise] ignored -- set from the other sections.")
  passed <- strip_noise_section[setdiff(names(strip_noise_section), c(set_keys, "run", "ref_label", "align_label"))]
  gradient <- flt_path(TWutils::align_gradient_raster(ref_dtm$path, as.numeric(align_params$radius %||% 15)))
  surveys <- list(ref = list(set = ref_set, dtm = ref_dtm$path, dsm = ref_dsm$path),
                  align = list(set = align_set, dtm = align_dtm$path, dsm = align_dsm$path))
  for (sv in names(surveys)) {
    x <- surveys[[sv]]
    tally <- file.path(x$set$output_dir, "pass_tally.csv")
    params <- c(list(input_dir = x$set$input_dir, output_dir = noise_dir, label = noise_labels[[sv]],
                     template_raster = ref_dtm$path, dtm = x$dtm, dsm = x$dsm,
                     gradient_raster = if (file.exists(gradient)) gradient else "NOFILE",
                     pass_tally = if (file.exists(tally)) tally else "NOFILE",
                     clip_polygon = study_area_params$clip_polygon), passed)
    run_script_step(strip_noise_script, paste0("strip_noise_", sv),
                    write_param_file(params, file.path(param_dir, paste0("params_strip_noise_", sv, ".txt"))),
                    file.path(noise_dir, paste0(noise_labels[[sv]], "_noise_classes.csv")))
  }
}

message("\nPipeline complete.")
