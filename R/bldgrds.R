## bldgrds.R
##
## Runs the Fortran program bldgrds via the local TWutils package's bldgrds() wrapper. bldgrds
## computes flow direction and D-infinity contributing area for a DEM, then traces the channel
## network downstream and writes it out as a node-list database (consumed by downstream programs
## such as netrace), optionally excavating the DEM along mapped road-crossing/culvert lines,
## masking out mapped water bodies before tracing channels, and writing a node point shapefile.
## Channel initiation is by area-slope threshold, separately calibrated for low- and high-gradient
## terrain, refined against a local relief raster and plan curvature.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## IMPORTANT: TWutils::bldgrds()/bldgrds_input() are verified against a real working reference
## input file (from the Sprague River project), NOT against GridUtilities/bldGrds2.f90 as checked
## out on this machine -- that checkout has been confirmed to NOT be the source bldgrds.exe was
## actually built from (see the "bldgrds" entries in TWutils's own CLAUDE.md). This script's
## parameter set matches bldgrds_input()'s current signature, which is the one to trust; if you
## hit a Fortran-side error, don't try to reconcile it against bldGrds2.f90's SELECT CASE.
##
## Unlike RIL.R/align_dtms.R's Fortran programs, bldgrds has no single output raster/keyword to
## read back afterwards, and -- as of this rewrite -- bldgrds_input() no longer takes path/dem_id
## arguments at all (dropped along with the rest of the older, untrustworthy bldGrds2.f90-derived
## keyword set), so this script can no longer derive or predict bldgrds' output file names. There
## is consequently no "skip if already done" resumability check here (see "Run bldgrds" below) --
## use_existing_files (a bldgrds_input() argument, exposed below) is the current API's own
## equivalent, telling bldgrds to reuse existing intermediate files rather than recalculating them.
##
## All input parameters below (dem, aspect_length, ..., overwrite) are read from a plain-text
## parameter file rather than hardcoded in this script -- see the "Parameter file" section just
## below for its format, and bldgrds_params_template.txt for a ready-to-copy example. Pass that
## file's path as a command-line argument when running via Rscript:
##   Rscript bldgrds.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input is read from a plain-text parameter file using "keyword: value" lines,
## one per line -- e.g.:
##
##   dem: c:\work\data\site1\elev_2023aligned.flt
##   aspect_length: 15
##   plan_length: 15
##   gradient_length_scale: 7.5
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types, and
## their defaults (the same defaults bldgrds_input() itself uses). Blank lines and lines starting
## with "#" are ignored, and a trailing "# comment" after a value is stripped. A value may itself
## contain a colon (e.g. a Windows drive letter, "c:\..."); only the FIRST colon on a line splits
## the keyword from its value, so that's safe. See param_specs below for exactly which keywords
## have no default and must be present in the file -- bldgrds_input() has no physically meaningful
## default for most of its channel-initiation-criteria arguments, so that list is long.
##
## Several optional keywords (excavate_line_buffer, water_mask_min_patch_size,
## water_mask_min_gradient, water_mask_buffer_radius, node_splits) have no meaningful numeric
## default -- leave them as NOFILE (or omit the line) to leave that keyword out of bldgrds' input
## file entirely.
##
## attribute_list (the per-node attributes written to the node-list database and, if requested,
## node_shapefile) is NOT set from this parameter file directly -- it's an R object
## (attribute_spec()/equation_term() calls), not expressible as flat text, same as RIL.R's
## attribute_list. Instead this script always calls
## TWutils::bldgrds_default_attributes(precip_raster) (see the precip_raster parameter below,
## same convention RIL.R uses for its own precip_raster); to use a genuinely custom attribute_list,
## edit this script directly.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory -- e.g. `Rscript R/bldgrds.R` from the repo root and
# `Rscript bldgrds.R` from inside R/ both find R/bldgrds_params_template.txt. Falls back to "."
# (assume cwd) when there's no "--file=" argument to read, e.g. when this script is source()'d
# from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript bldgrds.R path/to/params.txt` -- the line below is only used as a fallback
# when no such argument is given (e.g. sourcing this script from RStudio), so set it there in
# that case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "R/bldgrds_params_template.txt")

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
# TWutils::is_missing_path() recognizes (case-insensitive "nofile"/"none"/"na", blank), used
# here for the "optional_numeric"/"optional_integer" parameter types below.
is_nofile_text <- function(raw_val) {
  tolower(trimws(raw_val)) %in% c("nofile", "none", "na", "")
}

# The full set of recognized parameter-file keywords: expected type, and default value used
# when a keyword is omitted from the file (required = TRUE means there is no default and the
# file must supply it). This table is the single place that documents what each keyword means,
# mirroring TWutils::bldgrds()'s/bldgrds_input()'s own arguments/defaults.
param_specs <- list(
  # Input DEM. No default -- must be set in the parameter file.
  dem = list(type = "character", required = TRUE),

  # Length (m) over which aspect is smoothed (USE SMOOTHED ASPECT: LENGTH SCALE) and plan
  # curvature is measured. No defaults -- must be set in the parameter file.
  aspect_length = list(type = "numeric", required = TRUE),
  plan_length   = list(type = "numeric", required = TRUE),

  # Length (m) over which gradient is measured. No default -- must be set in the parameter file.
  gradient_length_scale = list(type = "numeric", required = TRUE),

  # Weights on plan curvature and aspect in the D8 flow-direction calculation (D8 COEFFICIENTS),
  # and the length scales (m) for the aspect/plan curvature used in that same calculation
  # (D8 LENGTH SCALES). No defaults -- must be set in the parameter file.
  d8_plan_coefficient   = list(type = "numeric", required = TRUE),
  d8_aspect_coefficient = list(type = "numeric", required = TRUE),
  d8_aspect_length      = list(type = "numeric", required = TRUE),
  d8_plan_length        = list(type = "numeric", required = TRUE),

  # Inner and outer buffer distances (m) around a candidate initiation point (INITIATION BUFFER),
  # and the contributing area (sq m) above which channel initiation is forced regardless of
  # other thresholds (INITIATION BUFFER: AREA OVERRIDE). No defaults -- must be set in the
  # parameter file.
  initiation_buffer_inner  = list(type = "numeric", required = TRUE),
  initiation_buffer_outer  = list(type = "numeric", required = TRUE),
  initiation_area_override = list(type = "numeric", required = TRUE),

  # High and low local-relief thresholds for channel initiation (LOCAL RELIEF THRESHOLD). No
  # defaults -- must be set in the parameter file.
  local_relief_threshold_high = list(type = "numeric", required = TRUE),
  local_relief_threshold_low  = list(type = "numeric", required = TRUE),

  # Contributing-area threshold for channel initiation in low- and high-gradient terrain. No
  # defaults -- must be set in the parameter file.
  area_slope_threshold_low_gradient  = list(type = "numeric", required = TRUE),
  area_slope_threshold_high_gradient = list(type = "numeric", required = TRUE),

  # Plan curvature threshold for channel initiation in low- and high-gradient terrain. No
  # defaults -- must be set in the parameter file.
  plan_curvature_threshold_low_gradient  = list(type = "numeric", required = TRUE),
  plan_curvature_threshold_high_gradient = list(type = "numeric", required = TRUE),

  # Minimum flow length (m) for a channel-initiation threshold to apply. No default -- must be
  # set in the parameter file.
  minimum_threshold_flow_length = list(type = "numeric", required = TRUE),

  # Minimum channel length (m).
  minimum_channel_length = list(type = "numeric", default = 0),

  # Optional: precomputed local relief raster (e.g. from DEV()/LocalRelief), reused instead of
  # bldgrds calculating its own.
  local_relief_raster = list(type = "character", default = "NOFILE"),

  # If TRUE, reuse existing intermediate files (USE EXISTING FILES) rather than recalculating
  # them -- the current API's equivalent of a resumability shortcut (see the note near the top
  # of this script on why there's no file-existence skip check here any more).
  use_existing_files = list(type = "logical", default = FALSE),

  # If TRUE, run in calibration mode.
  calibrate = list(type = "logical", default = FALSE),

  # Optional: polyline shapefile (e.g. road crossings) to excavate through the DEM, clearing
  # culvert blockages, and a buffer (m) around it. excavate_line_buffer is only meaningful with
  # excavate_line; leave as NOFILE/omit otherwise.
  excavate_line        = list(type = "character", default = "NOFILE"),
  excavate_line_buffer = list(type = "optional_numeric", default = NA_real_),

  # Optional: water-body polygon/raster mask, and its associated options (all only meaningful
  # with water_mask; leave as NOFILE/FALSE/omit otherwise).
  water_mask                      = list(type = "character", default = "NOFILE"),
  water_mask_min_patch_size       = list(type = "optional_numeric", default = NA_real_),
  water_mask_set_to_min_elevation = list(type = "logical", default = FALSE),
  water_mask_incise_to_center     = list(type = "logical", default = FALSE),
  water_mask_min_gradient         = list(type = "optional_numeric", default = NA_real_),
  water_mask_buffer_radius        = list(type = "optional_numeric", default = NA_real_),
  water_mask_preclude_initiation  = list(type = "logical", default = FALSE),

  # Optional output node point shapefile -- supplying this turns on OUTPUT NODE POINT SHAPEFILE.
  node_shapefile = list(type = "character", default = "NOFILE"),

  # Optional: number of pieces to split each channel segment into for node_shapefile. Only
  # meaningful together with node_shapefile; leave as NOFILE/omit otherwise.
  node_splits = list(type = "optional_integer", default = NA_integer_),

  # Optional mean-annual-precipitation raster, passed to TWutils::bldgrds_default_attributes()
  # to build attribute_list (see the note on attribute_list above). NOFILE means the
  # precip-dependent attributes (mean annual precip/flow, channel width/depth) are omitted --
  # same convention RIL.R uses for its own precip_raster.
  precip_raster = list(type = "character", default = "NOFILE"),

  # Scratch directory (bldgrds' input file is written here) and the folder containing the
  # bldgrds executable. No defaults -- must be set in the parameter file.
  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),

  # If TRUE, allow overwriting an existing bldgrds input file in scratch_dir.
  overwrite = list(type = "logical", default = TRUE)
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
        character        = as.character(raw_val),
        numeric          = suppressWarnings(as.numeric(raw_val)),
        integer          = suppressWarnings(as.integer(raw_val)),
        logical          = suppressWarnings(as.logical(raw_val)),
        # "NOFILE" (or blank/none/na) means "leave this keyword out entirely" -- same sentinel
        # convention as the character NOFILE defaults above -- anything else must parse cleanly.
        optional_numeric = if (is_nofile_text(raw_val)) NA_real_ else suppressWarnings(as.numeric(raw_val)),
        optional_integer = if (is_nofile_text(raw_val)) NA_integer_ else suppressWarnings(as.integer(raw_val)),
        stop("Internal error: unknown parameter type '", spec$type, "' for '", nm, "'")
      )
      if (spec$type %in% c("numeric", "integer", "logical") && anyNA(converted)) {
        stop("Could not parse parameter '", nm, "' (value \"", raw_val, "\") as ", spec$type,
             " -- check ", path, ".")
      }
      if (spec$type %in% c("optional_numeric", "optional_integer") &&
          anyNA(converted) && !is_nofile_text(raw_val)) {
        stop("Could not parse parameter '", nm, "' (value \"", raw_val, "\") as ", spec$type,
             " -- check ", path, ".")
      }
      out[[nm]] <- converted
    } else if (!is.null(spec$default) || spec$type %in% c("optional_numeric", "optional_integer")) {
      out[[nm]] <- spec$default
    } else if (isTRUE(spec$required)) {
      stop("Required parameter '", nm, "' is missing from parameter file: ", path)
    }
  }
  out
}

raw_params <- read_param_file(config_path)
params     <- build_params(raw_params, param_specs, config_path)
list2env(params, envir = globalenv())  # makes dem, aspect_length, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

## ---- Run bldgrds -------------------------------------------------------------------------------
## No "skip if already done" check here -- see the note near the top of this script: bldgrds_input()
## no longer takes path/dem_id arguments, so this script has no way to predict bldgrds' output
## file names and can't tell whether a previous run already finished. use_existing_files (above)
## is the current API's own mechanism for avoiding redundant recalculation.

TWutils::bldgrds(dem = dem,
                 scratch_dir = scratch_dir,
                 aspect_length = aspect_length,
                 plan_length = plan_length,
                 gradient_length_scale = gradient_length_scale,
                 d8_plan_coefficient = d8_plan_coefficient,
                 d8_aspect_coefficient = d8_aspect_coefficient,
                 d8_aspect_length = d8_aspect_length,
                 d8_plan_length = d8_plan_length,
                 initiation_buffer_inner = initiation_buffer_inner,
                 initiation_buffer_outer = initiation_buffer_outer,
                 initiation_area_override = initiation_area_override,
                 local_relief_threshold_high = local_relief_threshold_high,
                 local_relief_threshold_low = local_relief_threshold_low,
                 area_slope_threshold_low_gradient = area_slope_threshold_low_gradient,
                 area_slope_threshold_high_gradient = area_slope_threshold_high_gradient,
                 plan_curvature_threshold_low_gradient = plan_curvature_threshold_low_gradient,
                 plan_curvature_threshold_high_gradient = plan_curvature_threshold_high_gradient,
                 minimum_threshold_flow_length = minimum_threshold_flow_length,
                 minimum_channel_length = minimum_channel_length,
                 local_relief_raster = local_relief_raster,
                 use_existing_files = use_existing_files,
                 calibrate = calibrate,
                 excavate_line = excavate_line,
                 excavate_line_buffer = if (is.na(excavate_line_buffer)) NULL else excavate_line_buffer,
                 water_mask = water_mask,
                 water_mask_min_patch_size = if (is.na(water_mask_min_patch_size)) NULL else water_mask_min_patch_size,
                 water_mask_set_to_min_elevation = water_mask_set_to_min_elevation,
                 water_mask_incise_to_center = water_mask_incise_to_center,
                 water_mask_min_gradient = if (is.na(water_mask_min_gradient)) NULL else water_mask_min_gradient,
                 water_mask_buffer_radius = if (is.na(water_mask_buffer_radius)) NULL else water_mask_buffer_radius,
                 water_mask_preclude_initiation = water_mask_preclude_initiation,
                 node_shapefile = node_shapefile,
                 node_splits = if (is.na(node_splits)) NULL else node_splits,
                 attribute_list = TWutils::bldgrds_default_attributes(precip_raster),
                 overwrite = overwrite,
                 executable_dir = executable_dir)

# TWutils::bldgrds() -> run_program() stops with an error on a nonzero exit code, so reaching
# this line means the run succeeded.
message("bldgrds finished. See scratch_dir (", scratch_dir, ") for its input file and log, and ",
        "the node-list database bldgrds writes alongside the DEM.")
