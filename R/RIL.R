## RIL.R
##
## Runs the Fortran program RIL (part of the Netstream suite) via the local TWutils package's
## RIL() wrapper. RIL finds channel initiation points on a DEM, traces the channel network
## downstream, and classifies the terrain around it into valley floor, hollow and inner-gorge
## landforms, writing a raster of those classes (and, optionally, a node point shapefile
## carrying per-node attributes).
##
## Requires: TWutils (local package -- must already be installed/available in the R environment)
##
## All input parameters below (dem, scratch_dir, out_RIL, radius, ..., overwrite) are read from
## a plain-text parameter file rather than hardcoded in this script -- see the "Parameter file"
## section just below for its format, and RIL_params_template.txt for a ready-to-copy example.
## Pass that file's path as a command-line argument when running via Rscript:
##   Rscript RIL.R path/to/your_params.txt
## or, when sourcing/running from R/RStudio, set config_path below before running the script.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
##
## Every workflow input is read from a plain-text parameter file using "keyword: value" lines,
## one per line -- e.g.:
##
##   dem: c:\work\data\site1\elev_2023aligned.flt
##   scratch_dir: c:\work\scratch
##   as2_threshold: 1000
##
## Because each line is self-labeled with its keyword, the order lines appear in the file does
## not matter -- see param_specs below for the full set of recognized keywords, their types,
## and their defaults (the same defaults RIL_input() itself uses -- the Post Mortem reference
## run). Blank lines and lines starting with "#" are ignored, and a trailing "# comment" after a
## value is stripped. A value may itself contain a colon (e.g. a Windows drive letter,
## "c:\..."); only the FIRST colon on a line splits the keyword from its value, so that's safe.
## dem, scratch_dir, out_RIL, and executable_dir have no default and must be present in the
## file; every other keyword falls back to the default shown in param_specs if omitted.
##
## A handful of parameters are themselves *groups* of related thresholds -- RIL_input()'s
## named-vector arguments (e.g. closest_node, hollow_gradient, ...). Those are written on one
## line as comma-separated "NAME=value" pairs, e.g.:
##   hollow_gradient: PRIMARY=0.3, SECONDARY=0.25, DIF1=1.0, DIF2=0.45, MIN POLY GRAD=0.2
## Names must match RIL_input()'s own sub-names exactly (see the "named_numeric" specs below);
## order within the group doesn't matter, but every name in the group must be present.
##
## attribute_list (the list of per-node attributes RIL computes) is read from a *separate* text
## file named by the attribute_list_file parameter -- an arbitrary "ATTRIBUTE LIST:" /
## "END LIST:" block, the same grammar RIL.exe's own input file uses (see
## RIL_attributes_example.txt for a ready-to-copy example). This script reads it with
## TWutils::read_attribute_list_file(), so any attribute name/argument combination RIL.exe
## understands can be requested, in any order, without editing this script -- including a
## precipitation-dependent MEAN ANNUAL PRECIP/FLOW/WIDTH/DEPTH chain, whose mean-annual-
## precipitation raster is given as that MEAN ANNUAL PRECIP entry's own FILE argument (see
## RIL_attributes_example.txt); there is no separate precip_raster parameter for that anywhere
## in this parameter file. When attribute_list_file is NOFILE, this script falls back to
## TWutils::ril_default_attributes() instead -- the bare identifier/area/geometry attributes
## only, with no precipitation-dependent chain.

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory -- e.g. `Rscript R/RIL.R` from the repo root and `Rscript RIL.R`
# from inside R/ both find R/RIL_params_template.txt. Falls back to "." (assume cwd) when
# there's no "--file=" argument to read, e.g. when this script is source()'d from RStudio
# instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

# Path to this run's parameter file. Overridden by a command-line argument when the script is
# run via `Rscript RIL.R path/to/params.txt` -- the line below is only used as a fallback when
# no such argument is given (e.g. sourcing this script from RStudio), so set it there in that
# case.
#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "RIL_params_template.txt")

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

# Parses one grouped ("NAME=value, NAME=value, ...") parameter value into a named numeric
# vector, e.g. "PRIMARY=0.3, SECONDARY=0.25" -> c(PRIMARY = 0.3, SECONDARY = 0.25). Used for
# RIL_input()'s named-vector arguments (closest_node, hollow_gradient, ...).
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

# The full set of recognized parameter-file keywords: expected type, and default value used
# when a keyword is omitted from the file (required = TRUE means there is no default and the
# file must supply it). This table is the single place that documents what each keyword means,
# mirroring TWutils::RIL_input()'s own arguments/defaults (the Post Mortem reference run).
param_specs <- list(
  # Input DEM. No default -- must be set in the parameter file.
  dem = list(type = "character", required = TRUE),

  # Scratch directory (RIL's input file is written here) and the folder containing the RIL
  # executable. No defaults -- must be set in the parameter file.
  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),

  # Output RIL raster (.flt), given without extension. No default -- must be set in the
  # parameter file.
  out_RIL = list(type = "character", required = TRUE),

  # Radius (m) for calculating elevation derivatives.
  radius = list(type = "numeric", default = 7.5),

  # CLOSEST NODE RASTER search parameters: NUM WIDTHS, MAX RADIUS, MIN RADIUS.
  closest_node = list(type = "named_numeric",
                      default = c(`NUM WIDTHS` = 0, `MAX RADIUS` = 250, `MIN RADIUS` = 250)),

  # Channel-initiation thresholds.
  as2_threshold            = list(type = "numeric", default = 1000),
  plan_curvature_threshold = list(type = "numeric", default = 0.005),
  gradient_threshold       = list(type = "numeric", default = 0.5),
  fluvial_area_threshold   = list(type = "numeric", default = 0.07),

  # Valley-floor parameters.
  valley_depth_max  = list(type = "numeric", default = 5.0),
  valley_buffer      = list(type = "numeric", default = 2.5),
  valley_max_hole    = list(type = "numeric", default = 1000),
  valley_min_patch   = list(type = "numeric", default = 20),

  # Hollow parameters.
  hollow_gradient = list(type = "named_numeric",
                        default = c(PRIMARY = 0.3, SECONDARY = 0.25, DIF1 = 1.0, DIF2 = 0.45,
                                    `MIN POLY GRAD` = 0.2)),
  hollow_tangential = list(type = "named_numeric",
                          default = c(PRIMARY = -0.005, SECONDARY = -0.02,
                                      `MIN POLY TAN` = 0.0)),
  hollow_profile = list(type = "named_numeric",
                        default = c(PRIMARY = -0.01, SECONDARY = -0.06,
                                    `MIN POLY PROF` = -0.015)),
  grad_proportion = list(type = "named_numeric", default = c(GRAD = 0.5, PROP = 0.1)),
  tan_proportion  = list(type = "named_numeric", default = c(TAN = 0.0, PROP = 0.5)),
  fill_hollow_embayments = list(type = "logical", default = TRUE),
  hollow_area_threshold   = list(type = "numeric", default = 5),
  hollow_max_hole         = list(type = "numeric", default = 1000),
  hollow_min_patch = list(type = "named_numeric", default = c(INITIAL = 0, FINAL = 300)),

  # Inner-gorge parameters.
  gorge_gradient = list(type = "named_numeric",
                        default = c(PRIMARY = 0.7, SECONDARY = 0.45, DIF1 = 0.95, DIF2 = 0.7)),
  gorge_max_hole  = list(type = "numeric", default = 1000),
  gorge_min_patch = list(type = "numeric", default = 0),

  # Hillslope parameters.
  slope_thresholds = list(type = "named_numeric", default = c(SLOPE = 0.4, STEEP = 0.7)),
  curve_thresholds  = list(type = "named_numeric",
                          default = c(CONVERGENT = 0.01, DIVERGENT = -0.01)),

  # Polygon smoothing and roads.
  edge_smoothing_iterations = list(type = "integer", default = 11),
  road_shapefile             = list(type = "character", default = "NOFILE"),
  road_buffer                = list(type = "numeric", default = 15),

  # Optional precomputed input rasters. "NOFILE" means not used. RIL only skips recomputing
  # elevation derivatives when all six of in_upgrad, in_downgrad, in_grad, in_tangential,
  # in_plan and in_prof are supplied; a partial set is ignored.
  in_closest_node     = list(type = "character", default = "NOFILE"),
  in_dist_to_channel  = list(type = "character", default = "NOFILE"),
  in_drainage_wing    = list(type = "character", default = "NOFILE"),
  in_valley_floor     = list(type = "character", default = "NOFILE"),
  in_upgrad           = list(type = "character", default = "NOFILE"),
  in_downgrad         = list(type = "character", default = "NOFILE"),
  in_grad             = list(type = "character", default = "NOFILE"),
  in_tangential       = list(type = "character", default = "NOFILE"),
  in_plan             = list(type = "character", default = "NOFILE"),
  in_prof             = list(type = "character", default = "NOFILE"),

  # Optional intermediate/output rasters and node shapefile, written out for reuse. "NOFILE"
  # means not wanted.
  out_closest_node    = list(type = "character", default = "NOFILE"),
  out_dist_to_channel = list(type = "character", default = "NOFILE"),
  out_drainage_wing   = list(type = "character", default = "NOFILE"),
  out_valley_floor    = list(type = "character", default = "NOFILE"),
  out_upgrad          = list(type = "character", default = "NOFILE"),
  out_downgrad        = list(type = "character", default = "NOFILE"),
  out_grad            = list(type = "character", default = "NOFILE"),
  out_tangential      = list(type = "character", default = "NOFILE"),
  out_plan            = list(type = "character", default = "NOFILE"),
  out_prof            = list(type = "character", default = "NOFILE"),
  out_nodes           = list(type = "character", default = "NOFILE"),
  out_zero_order      = list(type = "character", default = "NOFILE"),

  # Text file holding an arbitrary "ATTRIBUTE LIST:" / "END LIST:" block (see the note on
  # attribute_list near the top of this script, and RIL_attributes_example.txt). When NOFILE,
  # this script falls back to TWutils::ril_default_attributes() (bare identifier/area/geometry
  # attributes, no precipitation-dependent chain).
  attribute_list_file = list(type = "character", default = "NOFILE"),

  # If FALSE, channel-node drainage wings are built with standard D8 flow paths instead of
  # D8-LTD.
  use_ltd = list(type = "logical", default = TRUE),

  # If TRUE, RIL writes extra diagnostic rasters to hardcoded paths under c:\temp, which must
  # already exist.
  debug = list(type = "logical", default = FALSE),

  # If TRUE, allow overwriting an existing RIL input file in scratch_dir.
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
        character     = as.character(raw_val),
        numeric       = suppressWarnings(as.numeric(raw_val)),
        integer       = suppressWarnings(as.integer(raw_val)),
        logical       = suppressWarnings(as.logical(raw_val)),
        named_numeric = parse_named_vector(raw_val, nm, path),
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
params     <- build_params(raw_params, param_specs, config_path)
list2env(params, envir = globalenv())  # makes dem, scratch_dir, out_RIL, radius, ... ordinary top-level variables, exactly as if they'd been assigned by hand below

message("Loaded parameters from: ", config_path)

attribute_list <- if (toupper(attribute_list_file) != "NOFILE") {
  message("Reading attribute list from: ", attribute_list_file)
  TWutils::read_attribute_list_file(attribute_list_file)
} else {
  TWutils::ril_default_attributes()
}

## ---- Run RIL ---------------------------------------------------------------------------------
## Skips the work if out_RIL's output raster already exists, so a failed/interrupted run can be
## restarted without redoing a RIL run that already completed (same resumability convention used
## elsewhere in this repo, e.g. align_dtms.R's outOutlier skip).

out_RIL_flt <- if (grepl("\\.flt$", out_RIL, ignore.case = TRUE)) out_RIL else paste0(out_RIL, ".flt")

if (!file.exists(out_RIL_flt)) {
  ril_raster <- TWutils::RIL(dem = dem,
                             scratch_dir = scratch_dir,
                             out_RIL = out_RIL,
                             radius = radius,
                             closest_node = closest_node,
                             as2_threshold = as2_threshold,
                             plan_curvature_threshold = plan_curvature_threshold,
                             gradient_threshold = gradient_threshold,
                             fluvial_area_threshold = fluvial_area_threshold,
                             valley_depth_max = valley_depth_max,
                             valley_buffer = valley_buffer,
                             valley_max_hole = valley_max_hole,
                             valley_min_patch = valley_min_patch,
                             hollow_gradient = hollow_gradient,
                             hollow_tangential = hollow_tangential,
                             hollow_profile = hollow_profile,
                             grad_proportion = grad_proportion,
                             tan_proportion = tan_proportion,
                             fill_hollow_embayments = fill_hollow_embayments,
                             hollow_area_threshold = hollow_area_threshold,
                             hollow_max_hole = hollow_max_hole,
                             hollow_min_patch = hollow_min_patch,
                             gorge_gradient = gorge_gradient,
                             gorge_max_hole = gorge_max_hole,
                             gorge_min_patch = gorge_min_patch,
                             slope_thresholds = slope_thresholds,
                             curve_thresholds = curve_thresholds,
                             edge_smoothing_iterations = edge_smoothing_iterations,
                             road_shapefile = road_shapefile,
                             road_buffer = road_buffer,
                             in_closest_node = in_closest_node,
                             in_dist_to_channel = in_dist_to_channel,
                             in_drainage_wing = in_drainage_wing,
                             in_valley_floor = in_valley_floor,
                             in_upgrad = in_upgrad,
                             in_downgrad = in_downgrad,
                             in_grad = in_grad,
                             in_tangential = in_tangential,
                             in_plan = in_plan,
                             in_prof = in_prof,
                             out_closest_node = out_closest_node,
                             out_dist_to_channel = out_dist_to_channel,
                             out_drainage_wing = out_drainage_wing,
                             out_valley_floor = out_valley_floor,
                             out_upgrad = out_upgrad,
                             out_downgrad = out_downgrad,
                             out_grad = out_grad,
                             out_tangential = out_tangential,
                             out_plan = out_plan,
                             out_prof = out_prof,
                             out_nodes = out_nodes,
                             out_zero_order = out_zero_order,
                             attribute_list = attribute_list,
                             use_ltd = use_ltd,
                             debug = debug,
                             overwrite = overwrite,
                             executable_dir = executable_dir)

  message("RIL finished. Output raster written to: ", out_RIL)
} else {
  message("Output raster already exists, skipping RIL: ", out_RIL_flt)
}
