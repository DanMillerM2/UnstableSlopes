## valleyfloor.R
##
## Runs the Fortran program ValleyFloor via the local TWutils package's valleyfloor() wrapper.
## For each selected channel, ValleyFloor builds a cell-by-cell height/depth-above-channel
## surface across the surrounding valley, then, when requested, measures valley width at each of
## a series of depth-above-channel thresholds and/or (with method = 4) fits a TIN-based flood-
## inundation surface. ValleyFloor's real output is a binary per-channel data file,
## valleyfloor_<ID>.dat, written unconditionally *next to the DEM* (not in scratch_dir) -- every
## named output raster below (elevation/depth above channel, flood height/depth, D8) is genuinely
## optional, and a run may request none of them (write_data alone is a valid use, e.g. to
## populate data files for a later read_data pass). TWutils::valleyfloor() itself now returns
## that .dat file's path (via TWutils::valleyfloor_dat_file()), which this script reports at the
## end; the optional output rasters, if any were requested, get their own paths reported too.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## valleyfloor_input() reproduces the keyword grammar of a working reference run (Cherry
## project); see TWutils's own CLAUDE.md and valleyfloor_input()'s roxygen docs for the several
## legacy keywords/mode switches that a working input_valleyfloor.txt file still carries but
## ValleyFloor.f90's readInputFile() silently drops as no-ops (not written here, since they'd do
## nothing) -- notably, height-above is triggered by requesting an output raster
## (out_elev/out_depth/out_elev_bil/out_depth_bil), not by a "calculate height above" switch.
##
## Although the .dat file's path is now predictable in advance (TWutils::valleyfloor_dat_file()),
## this script does not add a RIL.R/align_dtms.R-style "skip if already done" check: ValleyFloor
## has its own overwrite_data/"OVERWRITE EXISTING DATA FILES" keyword controlling whether it's OK
## to overwrite that file, defaulting to TRUE, so a normal repeat run is expected to refresh it
## rather than be skipped. Note also that DATA ID (and so data_id below) is only ever honored in
## ALL CHANNELS mode (all_channels = TRUE) -- with all_channels = FALSE, ValleyFloor always falls
## back to the DEM-derived ID regardless of data_id, and TWutils::valleyfloor()'s own return
## value accounts for that.
##
## All input parameters below (dem, scratch_dir, ..., overwrite) are read from a plain-text
## parameter file rather than hardcoded in this script -- see the "Parameter file" section just
## below for its format, and valleyfloor_params_template.txt for a ready-to-copy example. Pass
## that file's path as a command-line argument when running via Rscript:
##   Rscript valleyfloor.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input is read from a plain-text parameter file using "keyword: value" lines,
## one per line -- e.g.:
##
##   dem: c:\work\data\cherry\elev_cherry.flt
##   scratch_dir: c:\work\scratch
##   method: 4
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types, and
## their defaults (the same defaults valleyfloor_input() itself uses -- the Cherry project
## reference run). Blank lines and lines starting with "#" are ignored, and a trailing
## "# comment" after a value is stripped. A value may itself contain a colon (e.g. a Windows
## drive letter, "c:\..."); only the FIRST colon on a line splits the keyword from its value, so
## that's safe. dem, scratch_dir, and executable_dir have no default and must be present in the
## file; every other keyword falls back to the default shown in param_specs if omitted.
##
## valley_buffer is a *group* of related thresholds (valleyfloor_input()'s named-vector
## argument). It's written on one line as comma-separated "NAME=value" pairs, e.g.:
##   valley_buffer: CHANNEL WIDTHS=150, MIN RADIUS=20, MAX RADIUS=1000
## Names must match valleyfloor_input()'s own sub-names exactly (see the "named_numeric" spec
## below); order within the group doesn't matter, but every name in the group must be present.
##
## channel_list (only used, and required, when all_channels is FALSE) and height_above_steps are
## each written as one line of comma-separated plain values, e.g.:
##   channel_list: 1, 2, 3
##   height_above_steps: 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10
##
## A handful of keywords are meaningful only together with another one and have no sensible
## numeric/logical default of their own (min_chan_area, second_expansion, monotonic,
## minimum_channel_length, fill_dem_holes_max_size, fill_dem_holes_min_elev, upstream_node,
## downstream_node) -- leave them as NOFILE (or omit the line) to leave that keyword out of
## ValleyFloor's input file entirely (i.e. pass NULL to valleyfloor_input()).
##
## attribute_list (the list of per-node attributes ValleyFloor computes -- ValleyFloor.f90 aborts
## with "No attributes specified" if none are given) is read from a *separate* text file named by
## the attribute_list_file parameter -- an arbitrary "ATTRIBUTE LIST:" / "END LIST:" block, the
## same grammar ValleyFloor.exe's own input file uses (see valleyfloor_attributes_example.txt for
## a ready-to-copy example). This script reads it with TWutils::read_attribute_list_file(), so any
## attribute name/argument combination ValleyFloor.exe understands can be requested, in any
## order, without editing this script -- including a mean-annual-precipitation-dependent chain,
## whose precipitation raster is given as that block's own MEAN ANNUAL PRECIP: FILE = ... entry;
## there is no separate precip_raster parameter for that anywhere in this parameter file. When
## attribute_list_file is NOFILE, this script falls back to
## TWutils::valleyfloor_default_attributes() instead -- the bare identifier/area attributes only,
## with no precipitation-dependent chain.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory -- e.g. `Rscript R/valleyfloor.R` from the repo root and
# `Rscript valleyfloor.R` from inside R/ both find R/valleyfloor_params_template.txt. Falls back
# to "." (assume cwd) when there's no "--file=" argument to read, e.g. when this script is
# source()'d from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript valleyfloor.R path/to/params.txt` -- the line below is only used as a
# fallback when no such argument is given (e.g. sourcing this script from RStudio), so set it
# there in that case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "valleyfloor_params_template.txt")

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
# TWutils::is_missing_path() recognizes (case-insensitive "nofile"/"none"/"na", blank), used here
# for the "optional_numeric"/"optional_integer"/"optional_logical" parameter types below.
is_nofile_text <- function(raw_val) {
  tolower(trimws(raw_val)) %in% c("nofile", "none", "na", "")
}

# Parses one grouped ("NAME=value, NAME=value, ...") parameter value into a named numeric
# vector, e.g. "CHANNEL WIDTHS=150, MIN RADIUS=20" -> c(`CHANNEL WIDTHS` = 150, `MIN RADIUS` = 20).
# Used for valleyfloor_input()'s one named-vector argument, valley_buffer.
parse_named_vector <- function(raw_val, nm, path) {
  parts <- trimws(strsplit(raw_val, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]

  entry_names <- character(length(parts))
  entry_vals  <- numeric(length(parts))
  for (i in seq_along(parts)) {
    eq_pos <- regexpr("=", parts[i], fixed = TRUE)
    if (eq_pos < 0) {
      stop("Malformed entry in grouped parameter '", nm, "' (expected \"NAME=value\"): \"",
           parts[i], "\" -- check ", path)
    }
    entry_names[i] <- trimws(substr(parts[i], 1, eq_pos - 1))
    val_str        <- trimws(substr(parts[i], eq_pos + 1, nchar(parts[i])))
    entry_vals[i]  <- suppressWarnings(as.numeric(val_str))
    if (is.na(entry_vals[i])) {
      stop("Could not parse value for '", entry_names[i], "' in grouped parameter '", nm,
           "' (value \"", val_str, "\") as numeric -- check ", path)
    }
  }
  setNames(entry_vals, entry_names)
}

# Parses one plain comma-separated list of numbers/integers into a vector, e.g.
# "0, 0.25, 0.5" -> c(0, 0.25, 0.5). Used for height_above_steps and channel_list.
parse_numeric_vector <- function(raw_val, nm, path) {
  parts <- trimws(strsplit(raw_val, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  vals  <- suppressWarnings(as.numeric(parts))
  if (anyNA(vals)) {
    stop("Could not parse one or more entries of '", nm, "' (\"", raw_val, "\") as numeric -- ",
         "check ", path)
  }
  vals
}

parse_integer_vector <- function(raw_val, nm, path) {
  parts <- trimws(strsplit(raw_val, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  vals  <- suppressWarnings(as.integer(parts))
  if (anyNA(vals)) {
    stop("Could not parse one or more entries of '", nm, "' (\"", raw_val, "\") as integer -- ",
         "check ", path)
  }
  vals
}

# The full set of recognized parameter-file keywords: expected type, and default value used
# when a keyword is omitted from the file (required = TRUE means there is no default and the
# file must supply it). This table is the single place that documents what each keyword means,
# mirroring TWutils::valleyfloor()'s/valleyfloor_input()'s own arguments/defaults (the Cherry
# project reference run).
param_specs <- list(
  # Input DEM. No default -- must be set in the parameter file.
  dem = list(type = "character", required = TRUE),

  # --- Channel selection (ALL CHANNELS / CHANNEL LIST) ---
  # If TRUE (the default), process every channel at least min_chan_width wide. Set FALSE and
  # supply channel_list to process only specific channel numbers instead.
  all_channels = list(type = "logical", default = TRUE),

  # Integer channel numbers to process, e.g. "1, 2, 3". Only used, and required, when
  # all_channels is FALSE.
  channel_list = list(type = "integer_vector", default = integer(0)),

  # Minimum channel width (m) or contributing area for a channel to be processed. Only
  # meaningful with all_channels = TRUE; min_chan_area has no default of its own (NOFILE omits
  # it, using min_chan_width alone).
  min_chan_width = list(type = "numeric", default = 1.0),
  min_chan_area  = list(type = "optional_numeric", default = NA_real_),

  # If TRUE, write/read the per-channel binary data files ValleyFloor manages.
  write_data = list(type = "logical", default = TRUE),
  read_data  = list(type = "logical", default = FALSE),

  # If TRUE, overwrite existing per-channel data files (or, with channel_list, the bare
  # OVERWRITE flag).
  overwrite_data = list(type = "logical", default = TRUE),

  # Optional data-file ID tag. Only meaningful with all_channels = TRUE.
  data_id = list(type = "character", default = "NOFILE"),

  # Height-above algorithm: 1 for the default per-cell weighted-average method, 4 for the
  # TIN-based method (required for out_flood_height/out_flood_depth).
  method = list(type = "integer", default = 4),

  # DEM sampling interval, in cells, for the height-above calculation.
  sampling_interval = list(type = "numeric", default = 1),

  # VALLEY BUFFER: how far from the channel to search for valley cells.
  valley_buffer = list(type = "named_numeric",
                       default = c(`CHANNEL WIDTHS` = 150, `MIN RADIUS` = 20, `MAX RADIUS` = 1000)),

  # Multipliers expanding the search radius around each valley cell. second_expansion has no
  # meaningful default of its own beyond valleyfloor_input()'s built-in 3.0 -- set it to NOFILE to
  # omit it (and rely on expansion_factor alone) instead.
  expansion_factor = list(type = "numeric", default = 1.5),
  second_expansion = list(type = "optional_numeric", default = 3.0),

  # If TRUE, restrict the valley mask to the local watershed.
  mask_by_watershed = list(type = "logical", default = TRUE),

  # Maximum channel-depth difference used to limit the distance-to-channel search.
  max_depth_dif = list(type = "numeric", default = 15),

  # Minimum absolute elevation difference that height-above is forced to extend to.
  min_elev_dif = list(type = "numeric", default = 2.0),

  # Exponent weighting nearby channel cells more heavily in the height-above calculation.
  weighting_exponent = list(type = "numeric", default = 1.0),

  # Smoothing passes applied to the height-above surface.
  smoothing_iterations = list(type = "integer", default = 0),

  # Maximum radius (m) for smoothing the flood-inundation surface. Only meaningful with
  # method = 4.
  inundation_smoothing_max_radius = list(type = "numeric", default = 20),

  # Limits, in channel depths, on how far above/below the channel height-above is mapped.
  max_depths_above = list(type = "numeric", default = 12),
  max_depths_below = list(type = "numeric", default = -12),

  # Maximum allowed difference between distance-to-channel measured by cell and by node.
  max_dif_dist_chan_dist_node = list(type = "numeric", default = 0.25),

  # If not NOFILE, write the MONOTONIC flag forcing each channel profile to increase in
  # elevation upstream (TRUE/FALSE writes YES/NO). NOFILE omits the keyword entirely.
  monotonic = list(type = "optional_logical", default = NA),

  # Optional minimum channel length (m). NOFILE omits the keyword.
  minimum_channel_length = list(type = "optional_numeric", default = NA_real_),

  # If TRUE, fill holes in the DEM before processing. fill_dem_holes_max_size (square meters)
  # and fill_dem_holes_min_elev are only meaningful with fill_dem_holes = TRUE; NOFILE omits
  # whichever of them isn't needed.
  fill_dem_holes           = list(type = "logical", default = FALSE),
  fill_dem_holes_max_size  = list(type = "optional_numeric", default = NA_real_),
  fill_dem_holes_min_elev  = list(type = "optional_numeric", default = NA_real_),

  # Optional precomputed D8 flow-direction raster, reused instead of recalculating it.
  input_d8 = list(type = "character", default = "NOFILE"),

  # --- Optional output rasters. Leave as NOFILE if not wanted. Either out_elev or out_depth
  # (in either .flt or .bil form) triggers the height-above calculation; out_flood_height/
  # out_flood_depth require method = 4. ---
  out_elev          = list(type = "character", default = "NOFILE"),
  out_depth         = list(type = "character", default = "NOFILE"),
  out_elev_bil      = list(type = "character", default = "NOFILE"),
  out_depth_bil     = list(type = "character", default = "NOFILE"),
  out_flood_height  = list(type = "character", default = "NOFILE"),
  out_flood_depth   = list(type = "character", default = "NOFILE"),
  out_d8            = list(type = "character", default = "NOFILE"),

  # If TRUE, measure valley widths at each of height_above_steps.
  measure_valley_widths = list(type = "logical", default = TRUE),

  # Channel-depth steps, in channel depths above the channel, at which to measure valley width,
  # e.g. "0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10".
  height_above_steps = list(type = "numeric_vector",
                            default = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10)),

  # Window length, in channel widths, over which valley width is averaged along the channel.
  valley_width_window = list(type = "numeric", default = 20),

  # Text file holding a standalone "ATTRIBUTE LIST:" / "END LIST:" block (see the note on
  # attribute_list above, and valleyfloor_attributes_example.txt). When NOFILE, this script
  # falls back to TWutils::valleyfloor_default_attributes() (bare identifier/area attributes, no
  # precipitation-dependent chain).
  attribute_list_file = list(type = "character", default = "NOFILE"),

  # Optional: restrict processing to the channel segment between these two node IDs. NOFILE
  # omits either keyword.
  upstream_node   = list(type = "optional_integer", default = NA_integer_),
  downstream_node = list(type = "optional_integer", default = NA_integer_),

  # If TRUE, write the TIME IT flag, timing the height-above calculation.
  time_it = list(type = "logical", default = FALSE),

  # If TRUE, write the DEBUG flag.
  debug = list(type = "logical", default = FALSE),

  # Scratch directory (ValleyFloor's input file is written here) and the folder containing the
  # ValleyFloor executable. No defaults -- must be set in the parameter file.
  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),

  # If TRUE, allow overwriting an existing ValleyFloor input file in scratch_dir.
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

  optional_types <- c("optional_numeric", "optional_integer", "optional_logical")

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
        named_numeric    = parse_named_vector(raw_val, nm, path),
        numeric_vector   = parse_numeric_vector(raw_val, nm, path),
        integer_vector   = if (is_nofile_text(raw_val)) integer(0) else parse_integer_vector(raw_val, nm, path),
        # "NOFILE" (or blank/none/na) means "leave this keyword out entirely" -- same sentinel
        # convention as the character NOFILE defaults above -- anything else must parse cleanly.
        optional_numeric = if (is_nofile_text(raw_val)) NA_real_    else suppressWarnings(as.numeric(raw_val)),
        optional_integer = if (is_nofile_text(raw_val)) NA_integer_ else suppressWarnings(as.integer(raw_val)),
        optional_logical = if (is_nofile_text(raw_val)) NA          else suppressWarnings(as.logical(raw_val)),
        stop("Internal error: unknown parameter type '", spec$type, "' for '", nm, "'")
      )
      if (spec$type %in% c("numeric", "integer", "logical") && anyNA(converted)) {
        stop("Could not parse parameter '", nm, "' (value \"", raw_val, "\") as ", spec$type,
             " -- check ", path, ".")
      }
      if (spec$type %in% optional_types && anyNA(converted) && !is_nofile_text(raw_val)) {
        stop("Could not parse parameter '", nm, "' (value \"", raw_val, "\") as ", spec$type,
             " -- check ", path, ".")
      }
      out[[nm]] <- converted
    } else if (!is.null(spec$default) || spec$type %in% optional_types) {
      out[[nm]] <- spec$default
    } else if (isTRUE(spec$required)) {
      stop("Required parameter '", nm, "' is missing from parameter file: ", path)
    }
  }
  out
}

raw_params <- read_param_file(config_path)
params     <- build_params(raw_params, param_specs, config_path)
list2env(params, envir = globalenv())  # makes dem, scratch_dir, method, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

attribute_list <- if (toupper(attribute_list_file) != "NOFILE") {
  message("Reading attribute list from: ", attribute_list_file)
  TWutils::read_attribute_list_file(attribute_list_file)
} else {
  TWutils::valleyfloor_default_attributes()
}

## ---- Run ValleyFloor ---------------------------------------------------------------------------

dat_file <- TWutils::valleyfloor(dem = dem,
                                 scratch_dir = scratch_dir,
                                 all_channels = all_channels,
                                 channel_list = if (length(channel_list) == 0) NULL else channel_list,
                                 min_chan_width = min_chan_width,
                                 min_chan_area = if (is.na(min_chan_area)) NULL else min_chan_area,
                                 write_data = write_data,
                                 read_data = read_data,
                                 overwrite_data = overwrite_data,
                                 data_id = data_id,
                                 method = method,
                                 sampling_interval = sampling_interval,
                                 valley_buffer = valley_buffer,
                                 expansion_factor = expansion_factor,
                                 second_expansion = if (is.na(second_expansion)) NULL else second_expansion,
                                 mask_by_watershed = mask_by_watershed,
                                 max_depth_dif = max_depth_dif,
                                 min_elev_dif = min_elev_dif,
                                 weighting_exponent = weighting_exponent,
                                 smoothing_iterations = smoothing_iterations,
                                 inundation_smoothing_max_radius = inundation_smoothing_max_radius,
                                 max_depths_above = max_depths_above,
                                 max_depths_below = max_depths_below,
                                 max_dif_dist_chan_dist_node = max_dif_dist_chan_dist_node,
                                 monotonic = if (is.na(monotonic)) NULL else monotonic,
                                 minimum_channel_length = if (is.na(minimum_channel_length)) NULL else minimum_channel_length,
                                 fill_dem_holes = fill_dem_holes,
                                 fill_dem_holes_max_size = if (is.na(fill_dem_holes_max_size)) NULL else fill_dem_holes_max_size,
                                 fill_dem_holes_min_elev = if (is.na(fill_dem_holes_min_elev)) NULL else fill_dem_holes_min_elev,
                                 input_d8 = input_d8,
                                 out_elev = out_elev,
                                 out_depth = out_depth,
                                 out_elev_bil = out_elev_bil,
                                 out_depth_bil = out_depth_bil,
                                 out_flood_height = out_flood_height,
                                 out_flood_depth = out_flood_depth,
                                 out_d8 = out_d8,
                                 measure_valley_widths = measure_valley_widths,
                                 height_above_steps = height_above_steps,
                                 valley_width_window = valley_width_window,
                                 attribute_list = attribute_list,
                                 upstream_node = if (is.na(upstream_node)) NULL else upstream_node,
                                 downstream_node = if (is.na(downstream_node)) NULL else downstream_node,
                                 time_it = time_it,
                                 debug = debug,
                                 overwrite = overwrite,
                                 executable_dir = executable_dir)

# TWutils::valleyfloor() -> run_program() stops with an error on a nonzero exit code, so
# reaching this line means the run succeeded. It returns the full path of the binary
# valleyfloor_<ID>.dat data file ValleyFloor wrote (or, with read_data = TRUE, updated) next to
# the DEM -- ValleyFloor's real output, always written whatever else was requested.
message("ValleyFloor finished. Per-channel data file written: ", dat_file)

# Report whichever of the optional output rasters were actually requested, alongside the path
# each was written to -- unlike the .dat file above, ValleyFloor runs fine with none of these.
requested_rasters <- c(out_elev = out_elev, out_depth = out_depth,
                       out_elev_bil = out_elev_bil, out_depth_bil = out_depth_bil,
                       out_flood_height = out_flood_height, out_flood_depth = out_flood_depth,
                       out_d8 = out_d8)
requested_rasters <- requested_rasters[toupper(requested_rasters) != "NOFILE"]
if (length(requested_rasters) > 0) {
  message("Output raster(s) written:\n  ",
          paste(names(requested_rasters), requested_rasters, sep = " -> ", collapse = "\n  "))
}
