## align_dtms.R
## Standalone R script assembled from the first code chunk in Coregistration.qmd (@sec-coregister).
##
## Runs the Fortran program Align (part of the Netstream suite) via the local TWutils package to
## co-register one DTM ("align") to another ("ref"), estimating a horizontal/vertical shift
## (corrected for hillslope gradient and canopy-height change) and writing the aligned DTM plus
## the difference/outlier diagnostics used to derive a minimum level of detection (LoD).
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## All input parameters below (refDTM, alignDTM, iterations, k, ..., program_name) are read from
## a plain-text parameter file rather than hardcoded in this script -- see the "Parameter file"
## section just below for its format, and align_dtms_params_template.txt for a ready-to-copy
## example. Pass that file's path as a command-line argument when running via Rscript:
##   Rscript align_dtms.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input (refDTM, alignDTM, refGrnd, alignGrnd, refDSM, alignDSM, iterations, k,
## dampener, outDTM, tileNx, tileNy, overlap, radius, nslope, maxSlope, nAzimuth, outbins,
## outDif, outOutlier, scratch_dir, executable_dir, program_name) is read from a plain-text
## parameter file using "keyword: value" lines, one per line -- e.g.:
##
##   refDTM: c:\work\data\site1\elev_2017.flt
##   alignDTM: c:\work\data\site1\elev_2023.flt
##   iterations: 4
##   k: 1.5
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types,
## and their defaults. Blank lines and lines starting with "#" are ignored, and a trailing
## "# comment" after a value is stripped. A value may itself contain a colon (e.g. a Windows
## drive letter, "c:\...") -- only the FIRST colon on a line splits the keyword from its value,
## so that's safe. refDTM, alignDTM, refGrnd, alignGrnd, outDTM, outbins, outDif, outOutlier,
## scratch_dir, and executable_dir have no default and must be present in the file; every other
## keyword falls back to the default shown below if omitted.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory -- e.g. `Rscript R/align_dtms.R` from the repo root and
# `Rscript align_dtms.R` from inside R/ both find R/align_dtms_params_template.txt. Falls back
# to "." (assume cwd) when there's no "--file=" argument to read, e.g. when this script is
# source()'d from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript align_dtms.R path/to/params.txt` -- the line below is only used as a
# fallback when no such argument is given (e.g. sourcing this script from RStudio), so set it
# there in that case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "align_dtms_params_template.txt")

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
    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the FIRST colon splits keyword from value -- keeps values like "c:\..." intact
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
# file must supply it). This table is the single place that documents what each keyword means,
# mirroring TWutils::align()'s own arguments/defaults.
param_specs <- list(
  # Reference DTM (pre-change), and the DTM to be shifted/aligned to it (post-change). No
  # defaults -- must be set in the parameter file.
  refDTM   = list(type = "character", required = TRUE),
  alignDTM = list(type = "character", required = TRUE),

  # Ground-point density rasters corresponding to refDTM/alignDTM (e.g. as produced by
  # ground_returns.R), used to weight/restrict the alignment. No defaults -- must be set in the
  # parameter file.
  refGrnd   = list(type = "character", required = TRUE),
  alignGrnd = list(type = "character", required = TRUE),

  # Optional reference/align digital surface models. "NOFILE" (align()'s own default) means not
  # used.
  refDSM   = list(type = "character", default = "NOFILE"),
  alignDSM = list(type = "character", default = "NOFILE"),

  # Number of times to repeat the alignment calculation.
  iterations = list(type = "integer", default = 5),

  # Tukey-fence multiplier used to identify and exclude outliers (e.g. landslide sites) from the
  # alignment calculation; k <= 0 filters more points as outliers.
  k = list(type = "numeric", default = 0),

  # Proportion of each iteration's calculated shift to actually apply (1.0 = apply it in full).
  dampener = list(type = "numeric", default = 1),

  # Output path (no extension) for the aligned DTM. No default -- must be set in the parameter
  # file.
  outDTM = list(type = "character", required = TRUE),

  # Number of tiles in the x (east-west) and y (north-south) directions used to split the
  # analysis area, and the fractional overlap between adjacent tiles.
  tileNx  = list(type = "integer", default = 0),
  tileNy  = list(type = "integer", default = 0),
  overlap = list(type = "numeric", default = 0.5),

  # Radius (m) used when calculating local slope and aspect at each grid point.
  radius = list(type = "numeric", default = 15),

  # Number of slope bins, the maximum slope gradient included in those bins, and the number of
  # aspect (azimuth) bins, used when tabulating elevation differences by slope/aspect for outlier
  # detection.
  nslope    = list(type = "integer", default = 7),
  maxSlope  = list(type = "numeric", default = 1),
  nAzimuth  = list(type = "integer", default = 8),

  # Output path (no extension) for the per-slope/aspect-bin csv summary, the difference raster,
  # and the outlier-flagged raster. No defaults -- must be set in the parameter file (outOutlier
  # is also used below to decide whether Align has already been run -- see "Run Align").
  outbins    = list(type = "character", required = TRUE),
  outDif     = list(type = "character", required = TRUE),
  outOutlier = list(type = "character", required = TRUE),

  # Scratch folder for Align's intermediate files, and the folder containing the Align
  # executable. No defaults -- must be set in the parameter file.
  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),

  # Name of the Align executable (without extension) inside executable_dir.
  program_name = list(type = "character", default = "align")
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
list2env(params, envir = globalenv())  # makes refDTM, alignDTM, iterations, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

## ---- Run Align -------------------------------------------------------------------------------
## Skips the work if outOutlier's output file already exists, so a failed/interrupted run can be
## restarted without redoing an alignment that already completed (same resumability convention
## used elsewhere in this repo, e.g. ground_returns.R's per-batch merged files).

if (!file.exists(paste0(outOutlier, ".flt"))) {
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
} else {
  message("Outlier output already exists, skipping Align: ", paste0(outOutlier, ".flt"))
}
