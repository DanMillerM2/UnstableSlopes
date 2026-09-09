## bldgrds_enforce.R
##
## Runs the Fortran program bldgrds via the local TWutils package's bldgrds_enforce() wrapper --
## the third of bldgrds' three usage scenarios in this repo (alongside bldgrds_nochannels() and
## bldgrds()). Rather than initiating new channels from area-slope/plan-curvature/local-relief
## thresholds (what TWutils::bldgrds() does), this scenario excavates an existing, previously-
## mapped channel-network raster/shapefile (CHANNEL MASK) into the DEM and traces the node-list
## database from that network alone -- NO NEW CHANNELS is always written, precluding any channel
## initiation outside the provided stream layer. None of the channel-initiation-criteria
## parameters TWutils::bldgrds() needs (aspect/plan-curvature/D8 length scales, initiation
## buffer, area-slope/plan-curvature/local-relief thresholds, minimum flow length, ...) apply
## here, since no new initiation happens -- this script's parameter set is correspondingly much
## smaller than bldgrds.R's.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## bldgrds_enforce_input() reproduces the keyword grammar of a working reference run (Skykomish
## project): DEM FILE, SCRATCH, CHANNEL MASK (w/ FILE, DIG, RADIUS, DIRECTIONAL, INIT ALL),
## NO NEW CHANNELS, OUTPUT NODE POINT SHAPEFILE (w/ SPLITS), DRAINAGE WING RASTER, HAND RASTER
## (w/ FLOW THRESHOLD, NORMALIZE), TWI RASTER (w/ GRADIENT LENGTH SCALE), and an ATTRIBUTE LIST
## block closed with END LIST -- NOT END ATTRIBUTE LIST like bldgrds_input() (regular bldgrds.R)
## writes. That's a real discrepancy between the two reference files bldgrds_input() and
## bldgrds_enforce_input() are each built from, not a mistake in one of them -- see the "bldgrds"
## entries in TWutils's own CLAUDE.md; the two haven't been reconciled.
##
## bldgrds has no single output raster/keyword to read back afterwards, and -- as with regular
## bldgrds.R -- this script has no way to derive or predict bldgrds' output file names, so there
## is no "skip if already done" resumability check here.
##
## All input parameters below (dem, channel_mask, channel_mask_dig, ..., overwrite) are read from
## a plain-text parameter file rather than hardcoded in this script -- see the "Parameter file"
## section just below for its format, and bldgrds_enforce_params_template.txt for a ready-to-copy
## example. Pass that file's path as a command-line argument when running via Rscript:
##   Rscript bldgrds_enforce.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input is read from a plain-text parameter file using "keyword: value" lines,
## one per line -- e.g.:
##
##   dem: c:\work\data\skykomish\elev_3m.flt
##   channel_mask: c:\work\data\skykomish\channel_network
##   channel_mask_dig: 15.0
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types, and
## their defaults. Blank lines and lines starting with "#" are ignored, and a trailing
## "# comment" after a value is stripped. A value may itself contain a colon (e.g. a Windows
## drive letter, "c:\..."); only the FIRST colon on a line splits the keyword from its value, so
## that's safe. See param_specs below for exactly which keywords have no default and must be
## present in the file.
##
## Several optional keywords (hand_flow_threshold, twi_gradient_length_scale, node_splits) have
## no meaningful numeric default -- leave them as NOFILE (or omit the line) to leave that keyword
## out of bldgrds' input file entirely.
##
## attribute_list (the per-node attributes written to the node-list database and, if requested,
## node_shapefile) is NOT set from this parameter file directly -- it's an R object
## (attribute_spec()/equation_term() calls), not expressible as flat text. Instead, same
## convention as RIL.R:
##   - attribute_list_file (below), when set, names a separate text file holding a standalone
##     "ATTRIBUTE LIST:" / "END LIST:" block -- the same grammar bldgrds_enforce_input() itself
##     writes -- read with TWutils::read_attribute_list_file(). See
##     bldgrds_enforce_attributes_example.txt for a ready-to-copy example (transcribed from the
##     Skykomish reference run this script's keyword set is built from), including a
##     mean-annual-precipitation-dependent chain given via that block's own MEAN ANNUAL PRECIP:
##     FILE = ... entry -- there is no separate precip_raster keyword anywhere in this parameter
##     file.
##   - When attribute_list_file is NOFILE, this script falls back to
##     TWutils::bldgrds_default_attributes() instead -- the bare elevation/area attributes only,
##     with no precipitation-dependent chain; prefer attribute_list_file when you have a real
##     project-specific attribute list to match.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory -- e.g. `Rscript R/bldgrds_enforce.R` from the repo root and
# `Rscript bldgrds_enforce.R` from inside R/ both find R/bldgrds_enforce_params_template.txt.
# Falls back to "." (assume cwd) when there's no "--file=" argument to read, e.g. when this
# script is source()'d from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript bldgrds_enforce.R path/to/params.txt` -- the line below is only used as a
# fallback when no such argument is given (e.g. sourcing this script from RStudio), so set it
# there in that case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "bldgrds_enforce_params_template.txt")

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
# mirroring TWutils::bldgrds_enforce()'s/bldgrds_enforce_input()'s own arguments/defaults.
param_specs <- list(
  # Input DEM. No default -- must be set in the parameter file.
  dem = list(type = "character", required = TRUE),

  # Existing channel-network raster/shapefile to enforce (CHANNEL MASK: FILE). No default --
  # must be set in the parameter file.
  channel_mask = list(type = "character", required = TRUE),

  # Depth (DEM elevation units) to excavate/burn channel_mask into the DEM (CHANNEL MASK: DIG).
  # No default -- must be set in the parameter file.
  channel_mask_dig = list(type = "numeric", required = TRUE),

  # Radius used when excavating channel_mask into the DEM (CHANNEL MASK: RADIUS).
  channel_mask_radius = list(type = "numeric", default = 0),

  # If TRUE, treat channel_mask as directional (CHANNEL MASK: DIRECTIONAL).
  channel_mask_directional = list(type = "logical", default = FALSE),

  # If TRUE, seed channel initiation at every channel_mask cell, not just its ends (CHANNEL
  # MASK: INIT ALL).
  channel_mask_init_all = list(type = "logical", default = FALSE),

  # Optional output node point shapefile -- supplying this turns on OUTPUT NODE POINT SHAPEFILE
  # (and requires the attribute_list this script builds; see the note above). Leave as NOFILE
  # if not wanted.
  node_shapefile = list(type = "character", default = "NOFILE"),

  # Optional: number of pieces to split each channel segment into for node_shapefile. Only
  # meaningful together with node_shapefile; leave as NOFILE/omit otherwise.
  node_splits = list(type = "optional_integer", default = NA_integer_),

  # Optional: output drainage wing raster. Leave as NOFILE if not wanted.
  drainage_wing_raster = list(type = "character", default = "NOFILE"),

  # Optional: output HAND (height above nearest drainage) raster, with its flow-accumulation
  # threshold (as a proportion) and whether to normalize it. The latter two are only meaningful
  # together with hand_raster; leave as NOFILE/omit otherwise.
  hand_raster          = list(type = "character", default = "NOFILE"),
  hand_flow_threshold  = list(type = "optional_numeric", default = NA_real_),
  hand_normalize       = list(type = "logical", default = FALSE),

  # Optional: output topographic wetness index raster, with the length scale (m) over which
  # gradient is measured for it. The latter is only meaningful together with twi_raster; leave
  # as NOFILE/omit otherwise.
  twi_raster                = list(type = "character", default = "NOFILE"),
  twi_gradient_length_scale = list(type = "optional_numeric", default = NA_real_),

  # Text file holding a standalone "ATTRIBUTE LIST:" / "END LIST:" block (see the note on
  # attribute_list above, and bldgrds_enforce_attributes_example.txt). When NOFILE, this script
  # falls back to TWutils::bldgrds_default_attributes() (bare elevation/area attributes, no
  # precipitation-dependent chain).
  attribute_list_file = list(type = "character", default = "NOFILE"),

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
list2env(params, envir = globalenv())  # makes dem, channel_mask, channel_mask_dig, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

attribute_list <- if (toupper(attribute_list_file) != "NOFILE") {
  message("Reading attribute list from: ", attribute_list_file)
  TWutils::read_attribute_list_file(attribute_list_file)
} else {
  TWutils::bldgrds_default_attributes()
}

## ---- Run bldgrds (enforce mode) ---------------------------------------------------------------

TWutils::bldgrds_enforce(dem = dem,
                         scratch_dir = scratch_dir,
                         channel_mask = channel_mask,
                         channel_mask_dig = channel_mask_dig,
                         channel_mask_radius = channel_mask_radius,
                         channel_mask_directional = channel_mask_directional,
                         channel_mask_init_all = channel_mask_init_all,
                         node_shapefile = node_shapefile,
                         node_splits = if (is.na(node_splits)) NULL else node_splits,
                         drainage_wing_raster = drainage_wing_raster,
                         hand_raster = hand_raster,
                         hand_flow_threshold = if (is.na(hand_flow_threshold)) NULL else hand_flow_threshold,
                         hand_normalize = hand_normalize,
                         twi_raster = twi_raster,
                         twi_gradient_length_scale = if (is.na(twi_gradient_length_scale)) NULL else twi_gradient_length_scale,
                         attribute_list = attribute_list,
                         overwrite = overwrite,
                         executable_dir = executable_dir)

# TWutils::bldgrds_enforce() -> run_program() stops with an error on a nonzero exit code, so
# reaching this line means the run succeeded.
message("bldgrds (enforce mode) finished. See scratch_dir (", scratch_dir, ") for its input ",
        "file and log, and the node-list database bldgrds writes alongside the DEM.")
