## run_align_pipeline.R
##
## Driver script: builds ground-point-density rasters for two sets of LAZ/LAS tiles (a
## "reference" set and the set to be aligned to it -- the same per-set work
## run_ground_returns_batches.R does, here fixed to exactly these two named sets instead of an
## arbitrary list), gets one reference and one to-be-aligned DTM (and, optionally, DSM) -- each
## either mosaicked from a folder of GeoTIFF tiles or read from an already-mosaicked raster --
## then runs the Fortran program Align via TWutils (the same call align_dtms.R makes) on them,
## using the two density rasters as refGrnd/alignGrnd.
##
## Everything can be limited to one study-area polygon ([study_area] clip_polygon): only tiles
## intersecting it are used, and every raster built from tiles is masked to it. With
## clip_polygon NOFILE (the default), every tile in each tile folder is used, unmasked.
##
## Stage 1 (ground_returns.R x2): each of [ref_set] and [align_set] below is run as its own
## `Rscript ground_returns.R <params.txt>` subprocess, one at a time -- see
## run_ground_returns_batches.R for why each run gets its own subprocess (isolated future::plan
## and global-environment state) rather than being source()'d into this session.
##
## Stage 2 (mosaic_dtm.R, up to x4): each of [ref_dtm], [align_dtm], and the optional [ref_dsm]
## and [align_dsm] that sets tile_dir is run the same way, as its own `Rscript mosaic_dtm.R
## <params.txt>` subprocess: the tiles (intersecting clip_polygon, if set) are mosaicked,
## converted to meters, reprojected to target_crs if they aren't already in it (keeping their own
## pixel size), masked to the polygon, and written to <output_dtm>.flt / <output_dsm>.flt.
## Alternatively, a section can name an already-mosaicked raster with existing_dtm/existing_dsm;
## that raster is used as-is (no mosaic_dtm.R run, no clipping/reprojection/unit conversion), so
## it must already be an ESRI GridFloat .flt in meters and in target_crs -- checked before
## anything runs. A DSM section that's left out (or sets existing_dsm: NOFILE) means no DSM.
##
## Stage 3 (Align): once both density rasters and all DTMs/DSMs exist, TWutils::align() is called
## directly in this session (as in align_dtms.R) with refDTM/alignDTM/refDSM/alignDSM set to
## Stage 2's rasters and refGrnd/alignGrnd set to Stage 1's two rasters --
## ground_density_epsg2856.tif inside each set's own output_dir, the fixed filename
## ground_returns.R always writes its final raster to.
##
## Resumability: every file this pipeline produces is checked for existence before it's
## (re)created, and skipped if it's already there -- same convention used throughout this repo
## (see ground_returns.R's per-batch merged-file skip and align_dtms.R's outOutlier skip). In
## particular, since outOutlier is the very last thing the whole pipeline produces, if its output
## already exists the entire pipeline -- both ground_returns.R runs included -- is skipped
## outright; see "Skip everything" below. Existence is the only check: after changing
## clip_polygon (or any other parameter), delete the outputs it affects so they're rebuilt.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment),
## terra, plus whatever ground_returns.R and mosaic_dtm.R require (lasR, terra, sf, lidR, future,
## future.apply) in the R installation Rscript runs them with.
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
##                   section are ignored (with a warning).
##   [defaults]   -- optional; ground_returns.R parameters shared by both [ref_set] and
##                   [align_set] below, unless a set repeats the same keyword to override it.
##   [ref_set]    -- ground_returns.R parameters (input_dir, output_dir, ...) for the reference
##                   LAZ tiles (pre-change).
##   [align_set]  -- ground_returns.R parameters for the LAZ tiles to be aligned (post-change).
##   [ref_dtm]    -- mosaic_dtm.R parameters (tile_dir, output_dtm, z_units, ...) for the
##                   reference DTM tiles (pre-change), or existing_dtm alone to use an
##                   already-mosaicked DTM instead.
##   [align_dtm]  -- the same, for the DTM to be aligned (post-change).
##   [ref_dsm]    -- optional; the same as [ref_dtm] but for a reference DSM, with the keywords
##                   output_dsm/existing_dsm. Leave the section out, or set existing_dsm: NOFILE,
##                   for no DSM.
##   [align_dsm]  -- optional; the same, for the DSM to be aligned.
##   [align]      -- align_dtms.R / TWutils::align() parameters (iterations, k, ..., outOutlier,
##                   scratch_dir, executable_dir, program_name). refDTM/alignDTM,
##                   refDSM/alignDSM, and refGrnd/alignGrnd are NOT set here -- they are always
##                   taken from the sections above; any found in [align] are ignored (with a
##                   warning).
##
## input_dir/output_dir (per LAZ set), tile_dir/output_dtm or existing_dtm (per DTM), and
## outDTM/outbins/outDif/outOutlier/scratch_dir/executable_dir (in [align]) have no default and
## must be set. Everything else falls back to ground_returns.R's, mosaic_dtm.R's, or
## align_dtms.R's own built-in default -- see param_specs in each of those scripts.

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
if (length(cli_args) >= 1) pipeline_config_path <- cli_args[1]  # Rscript ... pipeline_config.txt takes precedence over the hardcoded fallback above

if (!file.exists(pipeline_config_path)) {
  stop("Pipeline config file not found: ", pipeline_config_path, ". Set pipeline_config_path ",
       "in this script, or pass its path as a command-line argument ",
       "(Rscript ... pipeline_config.txt).")
}

ground_returns_script <- file.path(get_script_dir(), "ground_returns.R")  # sibling script, same folder as this one
mosaic_dtm_script     <- file.path(get_script_dir(), "mosaic_dtm.R")      # likewise
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

study_area <- config[["study_area"]]
if (is.null(study_area)) study_area <- list()
ref_set    <- config[["ref_set"]]
align_set  <- config[["align_set"]]
ref_dtm    <- config[["ref_dtm"]]
align_dtm  <- config[["align_dtm"]]
ref_dsm    <- config[["ref_dsm"]]    # optional -- no reference DSM when absent
align_dsm  <- config[["align_dsm"]]  # likewise
align_section <- config[["align"]]

required_sections <- c("ref_set", "align_set", "ref_dtm", "align_dtm", "align")
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
unknown_study_keys <- setdiff(names(study_area), c("clip_polygon", "target_crs"))
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

# Each DTM/DSM section either names an already-mosaicked raster (existing_dtm/existing_dsm) or
# the tiles to build one from (tile_dir + output_dtm/output_dsm). Returns the .flt path Align
# will read as path, whether Stage 2 must run mosaic_dtm.R for it as from_tiles, and if so that
# run's mosaic_dtm.R parameters (output_dsm is passed on as mosaic_dtm.R's output_dtm). With
# optional = TRUE (the DSM sections), a missing/empty section or existing_dsm: NOFILE resolves to
# path NOFILE -- no DSM.
resolve_raster_section <- function(section, section_name, kind, optional = FALSE) {
  existing_key <- paste0("existing_", kind)
  output_key   <- paste0("output_", kind)
  no_raster    <- list(from_tiles = FALSE, path = "NOFILE")

  if (optional && length(section) == 0) return(no_raster)

  if (!is.null(section[[existing_key]])) {
    ignored <- intersect(names(section), c("tile_dir", output_key, "z_units", "file_pattern",
                                           "mosaic_fun", "clip_buffer_cells"))
    if (length(ignored) > 0) {
      warning(paste(ignored, collapse = "/"), " in [", section_name, "] ignored -- ",
              existing_key, " is set, so nothing is mosaicked for it.")
    }
    if (optional && toupper(section[[existing_key]]) == "NOFILE") return(no_raster)
    return(list(from_tiles = FALSE,
                path = check_existing_raster(section[[existing_key]], section_name, kind)))
  }

  if (is.null(section$tile_dir) || is.null(section[[output_key]])) {
    stop("[", section_name, "] must set either ", existing_key, ", or tile_dir and ", output_key,
         if (optional) paste0(" (or ", existing_key, ": NOFILE, or leave the section out)"),
         ": ", pipeline_config_path)
  }
  params <- section
  params[[output_key]] <- NULL
  params$output_dtm <- section[[output_key]]
  list(from_tiles = TRUE, path = flt_path(section[[output_key]]),
       params = modifyList(params, study_area_params))
}

ref_dtm   <- resolve_raster_section(ref_dtm,   "ref_dtm",   "dtm")
align_dtm <- resolve_raster_section(align_dtm, "align_dtm", "dtm")
ref_dsm   <- resolve_raster_section(ref_dsm,   "ref_dsm",   "dsm", optional = TRUE)
align_dsm <- resolve_raster_section(align_dsm, "align_dsm", "dsm", optional = TRUE)

## ---- Resolve Align's parameters up front -----------------------------------------------------
## Done before Stages 1-2 run (not after, as in align_dtms.R) so outOutlier is known early
## enough to decide, below, whether the whole pipeline can be skipped -- refGrnd/alignGrnd and
## refDTM/alignDTM/refDSM/alignDSM only need ref_set/align_set's output_dir and the DTM/DSM
## sections' output_dtm/output_dsm (or existing_dtm/existing_dsm), which are already known from
## the config without actually having run anything yet.

derived_align_keys <- intersect(names(align_section),
                                c("refGrnd", "alignGrnd", "refDTM", "alignDTM", "refDSM",
                                  "alignDSM"))
if (length(derived_align_keys) > 0) {
  warning(paste(derived_align_keys, collapse = "/"), " in [align] ignored -- refGrnd/alignGrnd ",
          "are always Stage 1's density rasters, and refDTM/alignDTM/refDSM/alignDSM come from ",
          "[ref_dtm]/[align_dtm]/[ref_dsm]/[align_dsm].")
  align_section[derived_align_keys] <- NULL
}

# Same table as align_dtms.R's param_specs, minus refGrnd/alignGrnd/refDTM/alignDTM/refDSM/
# alignDSM (those are set directly from Stage 1's and Stage 2's outputs below, never read from
# [align]).
align_param_specs <- list(

  iterations = list(type = "integer", default = 5),
  k          = list(type = "numeric", default = 0),
  dampener   = list(type = "numeric", default = 1),

  outDTM = list(type = "character", required = TRUE),

  tileNx  = list(type = "integer", default = 0),
  tileNy  = list(type = "integer", default = 0),
  overlap = list(type = "numeric", default = 0.5),

  radius = list(type = "numeric", default = 15),

  nslope   = list(type = "integer", default = 7),
  maxSlope = list(type = "numeric", default = 1),
  nAzimuth = list(type = "integer", default = 8),

  outbins    = list(type = "character", required = TRUE),
  outDif     = list(type = "character", required = TRUE),
  outOutlier = list(type = "character", required = TRUE),

  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),

  program_name = list(type = "character", default = "align")
)

# Same as align_dtms.R's build_params(): converts each recognized keyword to its declared type,
# falls back to its default when omitted, and stops on a missing required keyword or an
# unparsable value. Keywords not in align_param_specs are ignored with a warning.
build_params <- function(raw_params, specs, path) {
  unknown <- setdiff(names(raw_params), names(specs))
  if (length(unknown) > 0) {
    warning("Ignoring unrecognized parameter(s) in [align] of ", path, ": ",
            paste(unknown, collapse = ", "))
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
             " -- check [align] in ", path, ".")
      }
      out[[nm]] <- converted
    } else if (!is.null(spec$default)) {
      out[[nm]] <- spec$default
    } else if (isTRUE(spec$required)) {
      stop("Required parameter '", nm, "' is missing from [align] in ", path)
    }
  }
  out
}

align_params <- build_params(align_section, align_param_specs, pipeline_config_path)

# Expected paths only -- none of these files needs to exist yet at this point; whether they do
# is checked in Stages 1-2 below (which skip (re)running ground_returns.R/mosaic_dtm.R once
# their own output already exists).
align_params$refGrnd   <- file.path(ref_set$output_dir,   "ground_density_epsg2856.tif")
align_params$alignGrnd <- file.path(align_set$output_dir, "ground_density_epsg2856.tif")
align_params$refDTM    <- ref_dtm$path
align_params$alignDTM  <- align_dtm$path
align_params$refDSM    <- ref_dsm$path    # NOFILE when there's no [ref_dsm]
align_params$alignDSM  <- align_dsm$path  # NOFILE when there's no [align_dsm]

invisible(list2env(align_params, envir = globalenv()))  # makes refDTM, alignDTM, refDSM, alignDSM, refGrnd, alignGrnd, iterations, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

## ---- Skip everything if the final output already exists --------------------------------------
## outOutlier is the last file the whole pipeline produces (every ground_returns.R/mosaic_dtm.R
## run feeds into it via Align). If it's already there, there is nothing left for any stage to
## do, so skip Stages 1-3 entirely rather than just checking right before the Align call, as
## align_dtms.R does -- that would still go through every (potentially slow) earlier run for no
## reason once the pipeline has already finished.

outOutlier_path <- paste0(outOutlier, ".flt")

if (file.exists(outOutlier_path)) {
  message("Final outlier output already exists -- pipeline already complete, nothing to do: ",
          outOutlier_path)
} else {

  ## ---- Stage 1: ground_returns.R, once each for ref_set and align_set ------------------------
  ## Sequential, one Rscript subprocess at a time (see run_ground_returns_batches.R for why each
  ## set gets its own subprocess) -- parallelism happens *within* each run, across that run's own
  ## batches, via that run's own n_workers. Each set's own density raster is skipped (not
  ## regenerated) if it already exists, same as ground_returns.R's own per-batch resumability.

  # Merges one set's own parameters over ground_returns_defaults (set-specific values win),
  # writes the result out as a "keyword: value" ground_returns.R parameter file, and returns that
  # file's path. Same helper as write_param_file() in run_ground_returns_batches.R.
  write_ground_returns_param_file <- function(set_params, defaults, path) {
    params <- modifyList(defaults, set_params)
    if (is.null(params$input_dir) || is.null(params$output_dir)) {
      stop("Both [ref_set] and [align_set] must set input_dir and output_dir.")
    }
    write_param_file(params, path)
  }

  # Writes a named list out as a "keyword: value" parameter file and returns that file's path.
  write_param_file <- function(params, path) {
    lines <- vapply(names(params), function(nm) paste0(nm, ": ", params[[nm]]), character(1))
    writeLines(lines, path)
    path
  }

  # Runs one sibling script (ground_returns.R or mosaic_dtm.R) as its own Rscript subprocess on an
  # already-written parameter file and returns the path to the output it's expected to produce.
  # Skips the Rscript call entirely if that output already exists. label is just for messages.
  run_script_step <- function(script, label, param_path, expected_output) {
    if (file.exists(expected_output)) {
      message("\n==== ", label, ": output already exists, skipping ", basename(script), ": ",
              expected_output, " ====")
      return(expected_output)
    }

    message("\n==== ", label, " (", basename(script), ") ====")

    start_time <- Sys.time()
    status <- system2("Rscript", args = shQuote(c(script, param_path)))
    elapsed <- Sys.time() - start_time

    if (status != 0) {
      stop(basename(script), " failed for ", label, " (parameters: ", param_path,
           ") -- exit status ", status, ".")
    }

    message(label, " finished in ", format(unclass(elapsed), digits = 4), " ", units(elapsed))

    if (!file.exists(expected_output)) {
      stop(basename(script), " reported success for ", label, " but expected output is missing: ",
           expected_output)
    }
    expected_output
  }

  # Runs ground_returns.R for one named set and returns the path to the density raster it
  # produces.
  run_ground_returns_for_set <- function(label, set_params, defaults, expected_density_tif) {
    param_path <- write_ground_returns_param_file(set_params, defaults,
                                                  file.path(param_dir, paste0("params_", label, ".txt")))
    run_script_step(ground_returns_script, label, param_path, expected_density_tif)
  }

  ref_density_tif   <- run_ground_returns_for_set("ref_set",   ref_set,   ground_returns_defaults,
                                                   refGrnd)
  align_density_tif <- run_ground_returns_for_set("align_set", align_set, ground_returns_defaults,
                                                   alignGrnd)

  message("\nStage 1 complete. Ground-return density rasters:")
  message("  ref_set:   ", ref_density_tif)
  message("  align_set: ", align_density_tif)

  ## ---- Stage 2: mosaic_dtm.R for each DTM/DSM built from tiles ---------------------------------
  ## Same one-subprocess-per-run pattern as Stage 1 (mosaic_dtm.R handles DSM tiles exactly like
  ## DTM tiles). Each raster is skipped (not rebuilt) if its .flt already exists, and one given
  ## as existing_dtm/existing_dsm -- or a DSM that is NOFILE -- is never mosaicked at all.

  run_mosaic <- function(label, raster) {
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
    run_script_step(mosaic_dtm_script, label, param_path, raster$path)
  }

  ref_dtm_flt   <- run_mosaic("ref_dtm",   ref_dtm)
  align_dtm_flt <- run_mosaic("align_dtm", align_dtm)
  ref_dsm_flt   <- run_mosaic("ref_dsm",   ref_dsm)
  align_dsm_flt <- run_mosaic("align_dsm", align_dsm)

  message("\nStage 2 complete. DTMs/DSMs:")
  message("  ref_dtm:   ", ref_dtm_flt)
  message("  align_dtm: ", align_dtm_flt)
  message("  ref_dsm:   ", ref_dsm_flt)
  message("  align_dsm: ", align_dsm_flt)

  ## ---- Stage 3: Align --------------------------------------------------------------------------
  ## Reaching here means outOutlier_path did not exist when this run started, so Align has not
  ## already produced this pipeline's final output -- run it.

  message("\n==== Align ====")

  returnCode <- TWutils::align(refDTM,
                               alignDTM,
                               refGrnd,
                               alignGrnd,
                               refDSM,
                               alignDSM,
                               iterations,
                               k,
                               dampener,
                               outDTM,
                               tileNx,
                               tileNy,
                               overlap,
                               radius,
                               nslope,
                               maxSlope,
                               nAzimuth,
                               outbins,
                               outDif,
                               outOutlier,
                               scratch_dir,
                               executable_dir,
                               program_name)

  if (returnCode != 0) {
    stop("Error in Align")
  }

  message("Align finished. Aligned DTM written to: ", outDTM)
}

message("\nPipeline complete.")
