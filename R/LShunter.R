## LShunter.R
##
## Runs the Fortran program LShunter via the local TWutils package's LShunter() wrapper.
## LShunter finds candidate landslide patches on a DEM of difference (DoD), using a raster created
## by program Align: its exceedance probability raster (exceedance; recommended, thresholds in
## S = -log10(p)) or its outlier (k-value) raster (outlier). It works in two rounds -- a strict set of
## thresholds finds patch cores, then a looser set grows them -- so most thresholds come in
## pairs (threshold1/threshold2, min1/min2, max_accum1/max_accum2). Patches are also screened by
## gradient, flow accumulation and (optionally) proximity to roads, and the result is written as
## a patch raster.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## This script supports both of TWutils::LShunter()'s modes:
##   - Mode 1: a user-supplied, already-written LShunter input file (input_file parameter below)
##     is handed straight to LShunter.exe, bypassing LShunterInput() entirely. Every other
##     parameter below (outlier, accum, out_patch, scratch_dir, thresholds, ...) is ignored in
##     this mode -- the input file already carries everything the program needs -- except
##     executable_dir, which is still used as a fallback if the input file has no working
##     EXECUTABLE DIR: line of its own (see TWutils::resolve_executable_dir()). The output
##     raster's path lives in the file's own OUTPUT PATCH RASTER keyword, so this script's
##     "skip if already done" check (below) does not apply in mode 1 -- that mode always runs.
##   - Mode 2 (the default): this script builds the input file itself, via
##     TWutils::LShunterInput() (called by TWutils::LShunter()), from the outlier/accum/
##     out_patch/... parameters below.
##
## Gradient (mode 2) must be supplied as a raster: the one Align uses (Align's OUTPUT GRADIENT
## RASTER: YES writes <reference DEM>_grad<RADIUS>m, as sin(slope) -- gradient_units SIN, the
## default). LShunter no longer computes gradient itself, so it shares Align's slope and length
## scale; min1/min2 stay tan(slope).
##
## All input parameters are read from a plain-text parameter file rather than hardcoded in this
## script -- see the "Parameter file" section just below for its format, and
## LShunter_params_template.txt for a ready-to-copy mode-2 example (or
## LShunter_params_mode1_example.txt for a mode-1 one, pointing at an existing input file). Pass
## that file's path as a command-line argument when running via Rscript:
##   Rscript LShunter.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input is read from a plain-text parameter file using "keyword: value" lines,
## one per line -- e.g.:
##
##   outlier: c:\work\data\site1\outlier_2023
##   threshold1: -5.0
##   scratch_dir: c:\work\scratch
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types, and
## their defaults. Blank lines and lines starting with "#" are ignored, and a trailing
## "# comment" after a value is stripped. A value may itself contain a colon (e.g. a Windows
## drive letter, "c:\..."); only the FIRST colon on a line splits the keyword from its value, so
## that's safe. gradient, accum, out_patch, scratch_dir, and executable_dir have no default and
## must be present in the file (in mode 2), as must exactly one of exceedance and outlier; every
## other keyword falls back to the default shown in param_specs if omitted.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory. Falls back to "." (assume cwd) when there's no "--file=" argument
# to read, e.g. when this script is source()'d from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript LShunter.R path/to/params.txt` -- the line below is only used as a fallback
# when no such argument is given (e.g. sourcing this script from RStudio), so set it there in
# that case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "LShunter_params_template.txt")

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
# file must supply it). Mirrors TWutils::LShunter()/LShunterInput()'s arguments; the threshold
# defaults are the values used in MappingLandslides.qmd.
param_specs <- list(
  # Optional: path to an already-written LShunter input file (selects TWutils::LShunter()'s
  # mode 1). When set (anything other than NOFILE), this script skips LShunterInput() entirely
  # and hands the file straight to LShunter.exe -- every mode-2-only keyword below is ignored.
  # executable_dir (below) still applies, as a fallback used only if this file has no working
  # "EXECUTABLE DIR:" line of its own. Leave as NOFILE (the default) to build an input file from
  # the parameters below instead (mode 2).
  input_file = list(type = "character", default = "NOFILE"),

  # Raster created by Align: exactly one of exceedance (OUTPUT EXCEEDANCE PROBABILITY RASTER,
  # recommended) and outlier (OUTPUT OUTLIER RASTER). Flow accumulation raster (.flt) created by
  # bldgrds. accum has no default -- must be set (mode 2 only; ignored in mode 1).
  exceedance = list(type = "character", default = "NOFILE"),
  outlier    = list(type = "character", default = "NOFILE"),
  accum      = list(type = "character", required = TRUE),

  # Output patch raster (.flt), given without extension. No default (mode 2 only -- in mode 1
  # the output path comes from the input file's own OUTPUT PATCH RASTER keyword).
  out_patch = list(type = "character", required = TRUE),

  # Scratch directory (must already exist; LShunter's input file is written here, mode 2 only)
  # and the folder containing LShunter.exe (both modes -- see input_file above).
  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),

  # Maximum value for a patch, 1st (strict) and 2nd (looser) rounds: outlier (k) values with
  # outlier, or S = -log10(p) values with exceedance (e.g. -4 and -2, i.e. p <= 0.01% and 1%).
  # The defaults suit outlier; set both for exceedance.
  threshold1 = list(type = "numeric", default = -5.0),
  threshold2 = list(type = "numeric", default = -1.5),

  # Minimum gradient for a patch, 1st and 2nd rounds, as tan(slope).
  min1 = list(type = "numeric", default = 0.2),
  min2 = list(type = "numeric", default = 0.1),

  # Maximum flow accumulation for a patch, 1st and 2nd rounds.
  max_accum1 = list(type = "numeric", default = 1000),
  max_accum2 = list(type = "numeric", default = 1000),

  # Minimum patch size in square meters.
  min_size = list(type = "numeric", default = 10),

  # Gradient raster (.flt) that Align uses (Align's OUTPUT GRADIENT RASTER: YES writes
  # <reference DEM>_grad<RADIUS>m). Required (mode 2). gradient_units: SIN (default, as Align
  # writes it) or TAN; LShunter converts to tan(slope), the units of min1/min2.
  gradient       = list(type = "character", required = TRUE),
  gradient_units = list(type = "character", default = "SIN"),

  # Optional road-layer polyline shapefile (.shp) used to exclude road-related earth movement,
  # and the buffer (m) around it. road_buffer is required when road_shapefile is supplied.
  road_shapefile = list(type = "character", default = "NOFILE"),
  road_buffer    = list(type = "numeric", default = NA_real_, allow_na = TRUE),

  # --- Optional (LShunter 1.1): patch table, patch significance and controls. All from Align's
  # density-LoD outputs; leave at NOFILE/NA to leave each out. ---
  # Align's centred difference raster (OUTPUT CENTERED DIFFERENCE RASTER): patch depths and
  # volumes; needed by out_table and min_patch_significance.
  cen_dif = list(type = "character", default = "NOFILE"),
  # Align's IQR raster (<OUTPUT LOD RASTER>_iqr) and noise correlogram
  # (<BINS OUTPUT>_densityLodCorr.csv), together: each patch's volume uncertainty. noise_cover:
  # ALL (LShunter's default), OPEN or VEGETATED.
  iqr         = list(type = "character", default = "NOFILE"),
  noise_corr  = list(type = "character", default = "NOFILE"),
  noise_cover = list(type = "character", default = "NOFILE"),
  # RIL landform raster and distance-to-channel raster (m), reported per patch in the table.
  landform     = list(type = "character", default = "NOFILE"),
  dist_channel = list(type = "character", default = "NOFILE"),
  # Output patch table (.csv): area, centroid, depth and volume, and with iqr/noise_corr the
  # volume SD and significance.
  out_table = list(type = "character", default = "NOFILE"),
  # Drop patches whose -volume/SD(volume) is below this z value (needs cen_dif, iqr, noise_corr).
  min_patch_significance = list(type = "numeric", default = NA_real_, allow_na = TRUE),
  # NULL SIMULATION: realizations of simulated correlated noise run through the same rules (needs
  # noise_corr); null_tile is the FFT tile in points (a power of 2, at least 512).
  null_realizations = list(type = "integer", default = NA_integer_, allow_na = TRUE),
  null_tile         = list(type = "integer", default = NA_integer_, allow_na = TRUE),
  null_seed         = list(type = "integer", default = NA_integer_, allow_na = TRUE),
  # SEED DENSITY: seed patches only where both surveys' local ground-return density is at least
  # seed_den_min returns/m2 -- the two count rasters Align uses (ground_count_ref_grid.flt of each
  # LAZ set); seed_den_window is the window half-width in points.
  seed_den_ref    = list(type = "character", default = "NOFILE"),
  seed_den_align  = list(type = "character", default = "NOFILE"),
  seed_den_min    = list(type = "numeric", default = NA_real_, allow_na = TRUE),
  seed_den_window = list(type = "integer", default = NA_integer_, allow_na = TRUE),
  # SIGN FLIP CONTROL: also run the same rules on the negated input -- the gain patches, a control
  # for false alarms.
  sign_flip = list(type = "logical", default = FALSE),

  # If TRUE, allow overwriting an existing LShunter input file in scratch_dir.
  overwrite = list(type = "logical", default = TRUE)
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
# before build_params() runs so mode 2's required keywords (outlier, accum, out_patch,
# scratch_dir, executable_dir) don't force the parameter file to carry values that
# TWutils::LShunter() would just ignore anyway when input_file is supplied. executable_dir stays
# optional in mode 1 too -- it's only used there as a fallback.
raw_input_file <- if (!is.null(raw_params[["input_file"]])) raw_params[["input_file"]] else "NOFILE"
use_existing_input_file <- !is_nofile_text(raw_input_file)
if (use_existing_input_file) {
  for (nm in c("gradient", "accum", "out_patch", "scratch_dir", "executable_dir")) {
    param_specs[[nm]]$required <- FALSE
    param_specs[[nm]]$default  <- "NOFILE"
  }
}

params <- build_params(raw_params, param_specs, config_path)
list2env(params, envir = globalenv())  # makes outlier, accum, out_patch, threshold1, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

## ---- Run LShunter ----------------------------------------------------------------------------

if (use_existing_input_file) {

  # Mode 1: hand the existing input file straight to LShunter.exe. Every mode-2-only parameter
  # above is ignored by TWutils::LShunter() in this mode -- the file already carries everything
  # the program needs. executable_dir is passed through as a fallback only (NOFILE if left
  # unset), used only if the input file has no working "EXECUTABLE DIR:" line of its own.
  message("Using existing LShunter input file: ", input_file)
  TWutils::LShunter(input_file = input_file,
                    Executable_dir = executable_dir)

  message("LShunter finished using existing input file: ", input_file)

} else {

  # Skips the work if out_patch's output raster already exists, so a failed/interrupted
  # pipeline can be restarted without redoing an LShunter run that already completed (same
  # resumability convention used elsewhere in this repo, e.g. RIL.R's out_RIL skip).
  out_patch_flt <- if (grepl("\\.flt$", out_patch, ignore.case = TRUE)) out_patch else paste0(out_patch, ".flt")
  # With a patch table asked for, the run is complete only once that is written too.
  out_table_csv <- if (is_nofile_text(out_table)) NULL else
    paste0(sub("\\.csv$", "", out_table, ignore.case = TRUE), ".csv")

  if (!file.exists(out_patch_flt) || (!is.null(out_table_csv) && !file.exists(out_table_csv))) {

    # TWutils::LShunter()/LShunterInput() check that exactly one of Outlier and Exceedance is
    # given, that Gradient is, that road_buffer comes with Roads, and that the optional patch
    # table and control inputs come with what they need; LShunterInput() wants NULL (not a
    # sentinel) for an optional number or choice that isn't given, so translate NA/NOFILE.
    na_null <- function(x) if (is.na(x)) NULL else x
    TWutils::LShunter(Outlier = outlier,
                      Exceedance = exceedance,
                      Gradient = gradient,
                      Accum = accum,
                      OutPatch = out_patch,
                      ScratchDir = scratch_dir,
                      threshold1 = threshold1,
                      threshold2 = threshold2,
                      GradientUnits = gradient_units,
                      min1 = min1,
                      min2 = min2,
                      maxAccum1 = max_accum1,
                      maxAccum2 = max_accum2,
                      MinSize = min_size,
                      Roads = road_shapefile,
                      road_buffer = na_null(road_buffer),
                      CenDif = cen_dif,
                      IQR = iqr,
                      NoiseCorr = noise_corr,
                      NoiseCover = if (is_nofile_text(noise_cover)) NULL else noise_cover,
                      Landform = landform,
                      DistChannel = dist_channel,
                      OutTable = out_table,
                      MinPatchSignificance = na_null(min_patch_significance),
                      NullRealizations = na_null(null_realizations),
                      NullTile = na_null(null_tile),
                      NullSeed = na_null(null_seed),
                      SeedDenRef = seed_den_ref,
                      SeedDenAlign = seed_den_align,
                      SeedDenMin = na_null(seed_den_min),
                      SeedDenWindow = na_null(seed_den_window),
                      SignFlip = sign_flip,
                      overwrite = overwrite,
                      Executable_dir = executable_dir)

    message("LShunter finished. Output patch raster written to: ", out_patch)
  } else {
    message("Output raster already exists, skipping LShunter: ", out_patch_flt)
  }
}
