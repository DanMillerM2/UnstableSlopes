## LSstats.R
##
## Runs the Fortran program LSstats via the local TWutils package's LSstats() wrapper.
## LSstats summarizes a set of rasters (e.g. elevation derivatives, DoD k-values) over the
## candidate landslide polygons delineated by LShunter and over the landform classes delineated
## by RIL. For each listed raster i it writes two csv files:
##   - <out_csv>_<i>.csv: one row per landslide polygon -- dominant landform, area,
##     min/5/25/50/75/95 percent quantiles/max of the raster's values, and the shortest distance
##     to a floodplain edge (landform class 1);
##   - <out_csv>_landforms_<i>.csv: the same quantiles for each landform class.
## Distance to the floodplain edge is either read from an existing raster (distance_raster) or
## computed by LSstats out to radius meters (and optionally written to out_distance_raster).
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## This script supports both of TWutils::LSstats()'s modes:
##   - Mode 1: a user-supplied, already-written LSstats input file (input_file parameter below)
##     is handed straight to LSstats.exe, bypassing LSstats_input() entirely. Every other
##     parameter below (polygon_raster, landform_raster, rasters, out_csv, ...) is ignored in
##     this mode -- the input file already carries everything the program needs -- except
##     executable_dir, which is still used as a fallback if the input file has no working
##     EXECUTABLE DIR: line of its own (see TWutils::resolve_executable_dir()). The output csv
##     base name lives in the file's own OUTPUT CSV keyword, so this script's "skip if already
##     done" check (below) does not apply in mode 1 -- that mode always runs.
##   - Mode 2 (the default): this script builds the input file itself, via
##     TWutils::LSstats_input() (called by TWutils::LSstats()), from the parameters below.
##
## All input parameters are read from a plain-text parameter file rather than hardcoded in this
## script -- see the "Parameter file" section just below for its format, and
## LSstats_params_template.txt for a ready-to-copy mode-2 example (or
## LSstats_params_mode1_example.txt for a mode-1 one, pointing at an existing input file). Pass
## that file's path as a command-line argument when running via Rscript:
##   Rscript LSstats.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input is read from a plain-text parameter file using "keyword: value" lines,
## one per line -- e.g.:
##
##   polygon_raster: c:\work\data\site1\patch_2006_2023
##   rasters: c:\work\data\site1\grad_2006, c:\work\data\site1\outlier_2006_2023
##   scratch_dir: c:\work\scratch
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types, and
## their defaults. Blank lines and lines starting with "#" are ignored, and a trailing
## "# comment" after a value is stripped. A value may itself contain a colon (e.g. a Windows
## drive letter, "c:\..."); only the FIRST colon on a line splits the keyword from its value, so
## that's safe. rasters is a comma-separated list of raster paths, numbered 1, 2, ... in the
## order given (that number is the <i> in the output csv names). polygon_raster,
## landform_raster, rasters, out_csv, scratch_dir, and executable_dir have no default and must
## be present in the file (in mode 2); every other keyword falls back to the default shown in
## param_specs if omitted.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory. Falls back to "." (assume cwd) when there's no "--file=" argument
# to read, e.g. when this script is source()'d from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript LSstats.R path/to/params.txt` -- the line below is only used as a fallback
# when no such argument is given (e.g. sourcing this script from RStudio), so set it there in
# that case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "LSstats_params_template.txt")

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

# Returns TRUE for a raw text value that means "not supplied" -- same sentinel words
# TWutils::is_missing_path() recognizes (case-insensitive "nofile"/"none"/"na", blank).
is_nofile_text <- function(raw_val) {
  tolower(trimws(raw_val)) %in% c("nofile", "none", "na", "")
}

# The full set of recognized parameter-file keywords: expected type, and default value used
# when a keyword is omitted from the file (required = TRUE means there is no default and the
# file must supply it). Mirrors TWutils::LSstats()/LSstats_input()'s arguments.
param_specs <- list(
  # Optional: path to an already-written LSstats input file (selects TWutils::LSstats()'s
  # mode 1). When set (anything other than NOFILE), this script skips LSstats_input() entirely
  # and hands the file straight to LSstats.exe -- every mode-2-only keyword below is ignored.
  # executable_dir (below) still applies, as a fallback used only if this file has no working
  # "EXECUTABLE DIR:" line of its own. Leave as NOFILE (the default) to build an input file from
  # the parameters below instead (mode 2).
  input_file = list(type = "character", default = "NOFILE"),

  # Landslide polygon (patch) raster from LShunter, and landform raster from RIL (classes 1-10,
  # 1 = floodplain). No defaults -- must be set (mode 2 only; ignored in mode 1).
  polygon_raster  = list(type = "character", required = TRUE),
  landform_raster = list(type = "character", required = TRUE),

  # Comma-separated list of rasters to summarize, numbered 1, 2, ... in the order given. No
  # default (mode 2 only).
  rasters = list(type = "character", required = TRUE),

  # Base name for the output csv files; LSstats appends "_<i>.csv" and "_landforms_<i>.csv"
  # itself (any ".csv" extension given here is dropped). No default (mode 2 only -- in mode 1
  # it comes from the input file's own OUTPUT CSV keyword).
  out_csv = list(type = "character", required = TRUE),

  # Scratch directory (must already exist; LSstats's input file is written here, mode 2 only)
  # and the folder containing LSstats.exe (both modes -- see input_file above).
  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),

  # Distance to the floodplain edge: either an existing distance raster, or, when
  # distance_raster is NOFILE, a search radius (m) out to which LSstats computes it itself --
  # radius is then required and must be > 0. out_distance_raster optionally saves the computed
  # distance for reuse as distance_raster in later runs; both radius and out_distance_raster
  # are ignored when distance_raster is supplied.
  distance_raster     = list(type = "character", default = "NOFILE"),
  radius              = list(type = "numeric", default = NA_real_, allow_na = TRUE),
  out_distance_raster = list(type = "character", default = "NOFILE")
)

# Applies param_specs on top of the raw (character) values read_param_file() returned: converts
# each recognized keyword to its declared type, falls back to its default when the keyword is
# absent from the file, and stops with a clear error if a required keyword (no default) is
# missing or a value can't be parsed as its declared type. Keywords in the file that aren't in
# param_specs are ignored with a warning, so a typo'd keyword doesn't just vanish unnoticed. A
# numeric spec with allow_na = TRUE may be given as NA/NONE/NOFILE (or blank) to mean "not
# applicable".
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
      if (isTRUE(spec$allow_na) && is_nofile_text(raw_val)) {  # optional numeric explicitly given as NA/NONE/NOFILE/blank
        out[[nm]] <- NA_real_
        next
      }
      converted <- switch(spec$type,
        character = as.character(raw_val),
        numeric   = suppressWarnings(as.numeric(raw_val)),
        integer   = suppressWarnings(as.integer(raw_val)),
        logical   = suppressWarnings(as.logical(raw_val)),
        stop("Internal error: unknown parameter type '", spec$type, "' for '", nm, "'")
      )
      if (spec$type %in% c("numeric", "integer", "logical") && anyNA(converted)) {
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

# Mode 1 (existing input file) vs. mode 2 (build one from the parameters below) -- decide this
# before build_params() runs so mode 2's required keywords don't force the parameter file to
# carry values that TWutils::LSstats() would just ignore anyway when input_file is supplied.
# executable_dir stays optional in mode 1 too -- it's only used there as a fallback.
raw_input_file <- if (!is.null(raw_params[["input_file"]])) raw_params[["input_file"]] else "NOFILE"
use_existing_input_file <- !is_nofile_text(raw_input_file)
if (use_existing_input_file) {
  for (nm in c("polygon_raster", "landform_raster", "rasters", "out_csv", "scratch_dir",
               "executable_dir")) {
    param_specs[[nm]]$required <- FALSE
    param_specs[[nm]]$default  <- "NOFILE"
  }
}

params <- build_params(raw_params, param_specs, config_path)
list2env(params, envir = globalenv())  # makes polygon_raster, rasters, out_csv, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

## ---- Run LSstats -----------------------------------------------------------------------------

if (use_existing_input_file) {

  # Mode 1: hand the existing input file straight to LSstats.exe. Every mode-2-only parameter
  # above is ignored by TWutils::LSstats() in this mode -- the file already carries everything
  # the program needs. executable_dir is passed through as a fallback only (NOFILE if left
  # unset), used only if the input file has no working "EXECUTABLE DIR:" line of its own.
  message("Using existing LSstats input file: ", input_file)
  TWutils::LSstats(input_file = input_file,
                   executable_dir = executable_dir)

  message("LSstats finished using existing input file: ", input_file)

} else {

  # Split the comma-separated raster list into a character vector.
  raster_files <- trimws(strsplit(rasters, ",", fixed = TRUE)[[1]])
  raster_files <- raster_files[nzchar(raster_files)]
  if (length(raster_files) == 0) {
    stop("Parameter 'rasters' lists no raster files -- check ", config_path, ".")
  }

  # TWutils::LSstats() only checks that radius is a number when distance_raster is omitted, and
  # its own default (-9999) would pass that check and only fail inside LSstats.exe ("RADIUS not
  # specified") -- so require a positive radius here instead.
  if (is_nofile_text(distance_raster) && (is.na(radius) || radius <= 0)) {
    stop("distance_raster is NOFILE, so radius (> 0, in meters) is required for LSstats to ",
         "compute distance to the floodplain edge -- check ", config_path, ".")
  }

  # Skips the work if the last raster's landform csv already exists, so a failed/interrupted
  # pipeline can be restarted without redoing an LSstats run that already completed (same
  # resumability convention used elsewhere in this repo, e.g. LShunter.R's out_patch skip).
  # LSstats writes <out_csv>_<i>.csv then <out_csv>_landforms_<i>.csv for each raster in turn,
  # so the last raster's landform csv is the final file written.
  out_csv_base <- sub("\\.csv$", "", out_csv, ignore.case = TRUE)
  last_csv <- paste0(out_csv_base, "_landforms_", length(raster_files), ".csv")

  if (!file.exists(last_csv)) {

    TWutils::LSstats(polygon_raster = polygon_raster,
                     landform_raster = landform_raster,
                     rasters = raster_files,
                     out_csv = out_csv,
                     scratch_dir = scratch_dir,
                     distance_raster = distance_raster,
                     radius = if (is.na(radius)) -9999 else radius,
                     out_distance_raster = out_distance_raster,
                     executable_dir = executable_dir)

    message("LSstats finished. Output csv files written to: ", out_csv_base, "_<i>.csv and ",
            out_csv_base, "_landforms_<i>.csv")
  } else {
    message("Output csv already exists, skipping LSstats: ", last_csv)
  }
}
