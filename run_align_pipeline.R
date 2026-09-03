## run_align_pipeline.R
##
## Driver script: builds ground-point-density rasters for two sets of LAZ tiles (a "reference"
## set and the set to be aligned to it -- the same per-set work run_ground_returns_batches.R
## does, here fixed to exactly these two named sets instead of an arbitrary list), then runs the
## Fortran program Align via TWutils (the same call align_dtms.R makes) using those two
## just-produced density rasters as refGrnd/alignGrnd.
##
## Stage 1 (ground_returns.R x2): each of [ref_set] and [align_set] below is run as its own
## `Rscript ground_returns.R <params.txt>` subprocess, one at a time -- see
## run_ground_returns_batches.R for why each run gets its own subprocess (isolated future::plan
## and global-environment state) rather than being source()'d into this session.
##
## Stage 2 (Align): once both density rasters exist, TWutils::align() is called directly in this
## session (as in align_dtms.R) with refGrnd/alignGrnd set to those two rasters'
## paths -- ground_density_epsg2856.tif inside each set's own output_dir, the fixed filename
## ground_returns.R always writes its final raster to.
##
## Resumability: every file this pipeline produces is checked for existence before it's
## (re)created, and skipped if it's already there -- same convention used throughout this repo
## (see ground_returns.R's per-batch merged-file skip and align_dtms.R's outOutlier skip). In
## particular, since outOutlier is the very last thing the whole pipeline produces, if its output
## already exists the entire pipeline -- both ground_returns.R runs included -- is skipped
## outright; see "Skip everything" below.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## ---- Pipeline config file --------------------------------------------------------------------
##
## All inputs for both stages come from one plain-text config file (see
## run_align_pipeline_params_template.txt for a ready-to-copy example), read from
## pipeline_config_path below or a command-line argument:
##   Rscript run_align_pipeline.R path/to/your_pipeline_config.txt
##
## The file uses "keyword: value" lines (same rules as ground_returns.R's and align_dtms.R's own
## parameter files -- order within a section doesn't matter, blank/"#"-only lines are ignored, a
## trailing "# comment" is stripped, only the first colon on a line splits keyword from value)
## grouped under four section headers:
##
##   [defaults]   -- optional; ground_returns.R parameters shared by both [ref_set] and
##                   [align_set] below, unless a set repeats the same keyword to override it.
##   [ref_set]    -- ground_returns.R parameters (input_dir, output_dir, ...) for the reference
##                   LAZ tiles (pre-change).
##   [align_set]  -- ground_returns.R parameters for the LAZ tiles to be aligned (post-change).
##   [align]      -- align_dtms.R / TWutils::align() parameters (refDTM, alignDTM, iterations,
##                   k, ..., outOutlier, scratch_dir, executable_dir, program_name). refGrnd and
##                   alignGrnd are NOT set here -- they are always taken from the density rasters
##                   Stage 1 produces for [ref_set] and [align_set]; any refGrnd/alignGrnd found
##                   in [align] is ignored (with a warning).
##
## input_dir/output_dir (per set) and refDTM/alignDTM/outDTM/outbins/outDif/outOutlier/
## scratch_dir/executable_dir (in [align]) have no default and must be set. Everything else
## falls back to ground_returns.R's or align_dtms.R's own built-in default -- see param_specs in
## each of those scripts for the full list.

library(TWutils)

pipeline_config_path <- "run_align_pipeline_params_template.txt"  # fallback; overridden by a command-line argument below

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) pipeline_config_path <- cli_args[1]  # Rscript ... pipeline_config.txt takes precedence over the hardcoded fallback above

if (!file.exists(pipeline_config_path)) {
  stop("Pipeline config file not found: ", pipeline_config_path, ". Set pipeline_config_path ",
       "in this script, or pass its path as a command-line argument ",
       "(Rscript ... pipeline_config.txt).")
}

ground_returns_script <- "ground_returns.R"  # path to ground_returns.R -- edit if this driver isn't run from the same folder
param_dir <- file.path(tempdir(), "run_align_pipeline_params")  # generated ground_returns.R parameter files are written here
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

ref_set   <- config[["ref_set"]]
align_set <- config[["align_set"]]
align_section <- config[["align"]]

if (is.null(ref_set) || is.null(align_set) || is.null(align_section)) {
  stop("Pipeline config file must contain [ref_set], [align_set], and [align] sections: ",
       pipeline_config_path)
}

## ---- Resolve Align's parameters up front -----------------------------------------------------
## Done before Stage 1 runs (not after, as in align_dtms.R) so outOutlier is known early enough
## to decide, below, whether the whole pipeline can be skipped -- refGrnd/alignGrnd only need
## ref_set$output_dir/align_set$output_dir, which are already known from the config without
## actually having run ground_returns.R yet.

if (!is.null(align_section$refGrnd) || !is.null(align_section$alignGrnd)) {
  warning("refGrnd/alignGrnd in [align] are ignored -- they are always taken from the ",
          "ref_set/align_set density rasters Stage 1 produces.")
  align_section$refGrnd   <- NULL
  align_section$alignGrnd <- NULL
}

# Same table as align_dtms.R's param_specs, minus refGrnd/alignGrnd (those two are set directly
# from Stage 1's output below, never read from the config file).
align_param_specs <- list(
  refDTM   = list(type = "character", required = TRUE),
  alignDTM = list(type = "character", required = TRUE),

  refDSM   = list(type = "character", default = "NOFILE"),
  alignDSM = list(type = "character", default = "NOFILE"),

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

# Expected paths only -- neither file needs to exist yet at this point; whether they do is
# checked in Stage 1 below (which skips (re)running ground_returns.R once its own output already
# exists).
align_params$refGrnd   <- file.path(ref_set$output_dir,   "ground_density_epsg2856.tif")
align_params$alignGrnd <- file.path(align_set$output_dir, "ground_density_epsg2856.tif")

list2env(align_params, envir = globalenv())  # makes refDTM, alignDTM, refGrnd, alignGrnd, iterations, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

## ---- Skip everything if the final output already exists --------------------------------------
## outOutlier is the last file the whole pipeline produces (both ground_returns.R runs feed into
## it via Align). If it's already there, there is nothing left for either stage to do, so skip
## Stage 1 (both ground_returns.R runs) and Stage 2 (Align) entirely rather than just checking
## right before the Align call, as align_dtms.R does -- that would still repeat both
## (potentially slow) ground_returns.R runs for no reason once the pipeline has already finished.

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
    lines <- vapply(names(params), function(nm) paste0(nm, ": ", params[[nm]]), character(1))
    writeLines(lines, path)
    path
  }

  # Runs ground_returns.R for one named set (label is just for messages/file-naming) and returns
  # the path to the density raster it produces. Skips the Rscript call entirely if that raster
  # already exists.
  run_ground_returns_for_set <- function(label, set_params, defaults, expected_density_tif) {
    if (file.exists(expected_density_tif)) {
      message("\n==== ", label, ": density raster already exists, skipping ground_returns.R: ",
              expected_density_tif, " ====")
      return(expected_density_tif)
    }

    message("\n==== ", label, ": ", set_params$input_dir, " ====")

    param_path <- file.path(param_dir, paste0("params_", label, ".txt"))
    write_ground_returns_param_file(set_params, defaults, param_path)

    start_time <- Sys.time()
    status <- system2("Rscript", args = shQuote(c(ground_returns_script, param_path)))
    elapsed <- Sys.time() - start_time

    if (status != 0) {
      stop("ground_returns.R failed for ", label, " (input_dir: ", set_params$input_dir,
           ") -- exit status ", status, ".")
    }

    message(label, " finished in ", format(unclass(elapsed), digits = 4), " ", units(elapsed))

    if (!file.exists(expected_density_tif)) {
      stop("ground_returns.R reported success for ", label, " but expected output is missing: ",
           expected_density_tif)
    }
    expected_density_tif
  }

  ref_density_tif   <- run_ground_returns_for_set("ref_set",   ref_set,   ground_returns_defaults,
                                                   refGrnd)
  align_density_tif <- run_ground_returns_for_set("align_set", align_set, ground_returns_defaults,
                                                   alignGrnd)

  message("\nStage 1 complete. Ground-return density rasters:")
  message("  ref_set:   ", ref_density_tif)
  message("  align_set: ", align_density_tif)

  ## ---- Stage 2: Align --------------------------------------------------------------------------
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
