## lshunter_lsstats_pipeline.R
##
## Driver script: runs LShunter.R, then LSstats.R, each as its own `Rscript <script> <params.txt>`
## subprocess (same convention ril_pipeline.R and run_align_pipeline.R use). LShunter finds
## candidate landslide patches on a DEM of difference and writes them as a patch raster; LSstats
## then summarizes a list of rasters over those patches (and over RIL's landform classes).
##
## Unlike ril_pipeline.R, which writes each stage a two-line mode-1 parameter file, this driver
## hands each stage a parameter file you've already written for it -- LShunter.R's and LSstats.R's
## own "keyword: value" files, completely unmodified -- so either stage may run in mode 1 (its
## parameter file sets input_file: to an existing Fortran input file) or mode 2 (its parameter file
## carries the individual parameters). See LShunter_params_template.txt /
## LShunter_params_mode1_example.txt and LSstats_params_template.txt /
## LSstats_params_mode1_example.txt.
##
## The one link between the two stages is the patch raster: LShunter's output patch raster
## (out_patch, or OUTPUT PATCH RASTER in a mode-1 input file) must be the landslide polygon raster
## LSstats reads (polygon_raster, or LANDSLIDE POLYGON RASTER in a mode-1 input file). Neither
## program knows about the other, so a mismatch would silently summarize some other, stale patch
## raster -- this driver reads both paths back out of the two parameter (or input) files before
## running anything, and stops if they don't name the same raster. Likewise, when LSstats reads
## LShunter's patch table for polygon volumes (patch_table, or LSHUNTER PATCH TABLE), it must be
## LShunter's out_table (OUTPUT PATCH TABLE).
##
## Resumability: each stage is skipped if its final output already exists -- the patch raster
## (.flt) for LShunter, and the last raster's <out_csv>_landforms_<n>.csv for LSstats (the last
## file LSstats.f90 writes). LShunter.R and LSstats.R already do this themselves in mode 2; this
## driver adds the same check for mode 1, where neither script has one. Delete a stage's output to
## force it to rerun. Because a rerun of LShunter would leave LSstats's existing csv files
## describing the old patches, this driver also stops if the LSstats csv is older than the patch
## raster, rather than skipping LSstats and silently keeping the stale statistics -- delete the
## LSstats csv files (<out_csv>_*.csv) to rerun it.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment),
## and LShunter.R/LSstats.R present alongside this script (same R/ folder).
##
## ---- Pipeline config file --------------------------------------------------------------------
##
## One plain-text config file (see lshunter_lsstats_pipeline_params_template.txt for a
## ready-to-copy example), read from pipeline_config_path below or a command-line argument:
##   Rscript lshunter_lsstats_pipeline.R path/to/your_pipeline_config.txt
##
## "keyword: value" lines (same rules as the per-stage scripts' own parameter files -- order
## doesn't matter, blank/"#"-only lines are ignored, a trailing "# comment" is stripped, only the
## first colon on a line splits keyword from value):
##
##   lshunter_params:  LShunter.R parameter file, mode 1 or mode 2 (required)
##   lsstats_params:   LSstats.R parameter file, mode 1 or mode 2 (required)

library(TWutils)

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# sibling per-stage scripts and the fallback config path below are found by this script's location
# on disk rather than by the caller's working directory. With no "--file=" argument (not run via
# Rscript), it tries, in order: the path source() was given (RStudio's Source button, or
# source("R/...") from the console), then the active RStudio editor document (running lines/chunks
# with Ctrl+Enter), and only then falls back to "." (assume cwd).
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) return(dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)))
  for (frame in rev(sys.frames())) {
    if (exists("ofile", envir = frame, inherits = FALSE) && is.character(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile, winslash = "/", mustWork = FALSE)))
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    editor_path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(editor_path)) return(dirname(normalizePath(editor_path, winslash = "/", mustWork = FALSE)))
  }
  "."
}

pipeline_config_path <- file.path(get_script_dir(),
                                  "lshunter_lsstats_pipeline_params_template.txt")  # fallback; overridden by a command-line argument below

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) pipeline_config_path <- cli_args[1]  # Rscript ... pipeline_config.txt takes precedence over the hardcoded fallback above

if (!file.exists(pipeline_config_path)) {
  stop("Pipeline config file not found: ", pipeline_config_path, ". Set pipeline_config_path in ",
       "this script, or pass its path as a command-line argument (Rscript ... pipeline_config.txt).")
}

# Sibling scripts, same folder as this one -- each is run unmodified, via its own Rscript
# invocation, exactly as if a user had typed `Rscript <script> <params.txt>` by hand.
lshunter_script <- file.path(get_script_dir(), "LShunter.R")
lsstats_script  <- file.path(get_script_dir(), "LSstats.R")

for (script_path in c(lshunter_script, lsstats_script)) {
  if (!file.exists(script_path)) {
    stop("Cannot find ", basename(script_path), " next to this pipeline script (expected at ",
         script_path, ").")
  }
}

## ---- Parse the pipeline config file and each stage's parameter file ---------------------------

is_nofile_text <- function(x) is.null(x) || tolower(trimws(x)) %in% c("nofile", "none", "na", "")

# Reads a "keyword: value" file -- the pipeline config, or LShunter.R's/LSstats.R's own parameter
# file, which use the same grammar -- into a named list of raw (character) values.
read_config <- function(path) {
  raw_lines <- trimws(readLines(path, warn = FALSE))
  raw_lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]  # drop blank lines and comment-only lines

  params <- list()
  for (line in raw_lines) {
    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the FIRST colon splits keyword from value -- keeps values like "c:\..." intact
    if (colon_pos < 0) {
      stop("Malformed line in ", path, " (expected \"keyword: value\"): ", line)
    }
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(substr(line, colon_pos + 1, nchar(line)))
    val <- trimws(sub("\\s+#.*$", "", val))  # strip a trailing " # comment", if any, from the value

    if (key %in% names(params)) stop("Duplicate keyword '", key, "' in ", path)
    params[[key]] <- val
  }
  params
}

required_keys <- c("lshunter_params", "lsstats_params")

config <- read_config(pipeline_config_path)

unknown_keys <- setdiff(names(config), required_keys)
if (length(unknown_keys) > 0) {
  stop("Unrecognized keyword(s) in ", pipeline_config_path, ": ", paste(unknown_keys, collapse = ", "),
       ". Expected: ", paste(required_keys, collapse = ", "), ".")
}

for (key in required_keys) {
  if (is_nofile_text(config[[key]])) {
    stop("Pipeline config file must set '", key, ":' to an existing parameter file: ", pipeline_config_path)
  }
  if (!file.exists(config[[key]])) {
    stop(key, " not found: ", config[[key]])
  }
}

lshunter_params_path <- config[["lshunter_params"]]
lsstats_params_path  <- config[["lsstats_params"]]
lshunter_params <- read_config(lshunter_params_path)
lsstats_params  <- read_config(lsstats_params_path)

# A stage runs in mode 1 when its parameter file sets input_file: to anything but NOFILE (the same
# test LShunter.R and LSstats.R themselves apply).
is_mode1 <- function(stage_params) !is_nofile_text(stage_params[["input_file"]])

# Value of a stage parameter-file keyword that the stage requires in mode 2.
require_param <- function(stage_params, key, path) {
  if (is_nofile_text(stage_params[[key]])) {
    stop("'", key, ":' is not set in ", path, " -- this pipeline needs it to link LShunter to LSstats.")
  }
  stage_params[[key]]
}

## ---- Read keywords back out of a mode-1 Fortran input file ------------------------------------
## Keywords are matched case-insensitively and anchored, so e.g. "OUTPUT CSV" isn't mistaken for
## some longer keyword containing it. Comment lines are skipped by TWutils::get_keyword().

read_keyword <- function(input_file, keyword) {
  pattern <- paste0("^", gsub(" ", "\\\\s+", keyword), "$")
  lines <- readLines(input_file, warn = FALSE)
  keywords <- vapply(seq_along(lines), function(i) TWutils::get_keyword(lines, i), character(1))
  matches <- which(!is.na(keywords) & grepl(pattern, keywords, ignore.case = TRUE))
  # get_args() keeps a trailing " # comment", which the Fortran grammar allows -- strip it first.
  vapply(matches, function(i) TWutils::parse_arg(sub("\\s+#.*$", "", TWutils::get_args(lines, i)))$Value,
         character(1))
}

read_single_keyword <- function(input_file, keyword) {
  values <- read_keyword(input_file, keyword)
  if (length(values) != 1) {
    stop("Expected exactly one '", keyword, ":' line in ", input_file, ", found ", length(values), ".")
  }
  values
}

# Number of entries in an LSstats input file's RASTER LIST block: every keyword line between
# "RASTER LIST:" and "END LIST:" (LSstats.f90 counts each one as a raster, whatever its keyword).
count_raster_list <- function(input_file) {
  lines <- readLines(input_file, warn = FALSE)
  keywords <- toupper(vapply(seq_along(lines), function(i) TWutils::get_keyword(lines, i), character(1)))
  start <- which(grepl("^RASTER\\s+LIST$", keywords))
  if (length(start) != 1) {
    stop("Expected exactly one 'RASTER LIST:' line in ", input_file, ", found ", length(start), ".")
  }
  end <- which(grepl("END\\s+LIST", keywords) & seq_along(keywords) > start)
  if (length(end) == 0) stop("RASTER LIST in ", input_file, " has no 'END LIST:' line.")
  in_list <- keywords[seq_len(end[1] - 1)]
  in_list <- in_list[-seq_len(start)]
  sum(!is.na(in_list))
}

## ---- Resolve each stage's linking paths --------------------------------------------------------

# LShunter's output patch raster.
if (is_mode1(lshunter_params)) {
  lshunter_input_file <- lshunter_params[["input_file"]]
  if (!file.exists(lshunter_input_file)) stop("LShunter input_file not found: ", lshunter_input_file)
  out_patch <- read_single_keyword(lshunter_input_file, "OUTPUT PATCH RASTER")
} else {
  out_patch <- require_param(lshunter_params, "out_patch", lshunter_params_path)
}

# LSstats's polygon raster, output csv base name, and number of rasters summarized.
if (is_mode1(lsstats_params)) {
  lsstats_input_file <- lsstats_params[["input_file"]]
  if (!file.exists(lsstats_input_file)) stop("LSstats input_file not found: ", lsstats_input_file)
  # LSstats.f90 accepts LANDSLIDE PATCH RASTER as an alias for LANDSLIDE POLYGON RASTER.
  polygon_raster <- read_single_keyword(lsstats_input_file, "LANDSLIDE (POLYGON|PATCH) RASTER")
  out_csv        <- read_single_keyword(lsstats_input_file, "OUTPUT CSV")
  num_rasters    <- count_raster_list(lsstats_input_file)
} else {
  polygon_raster <- require_param(lsstats_params, "polygon_raster", lsstats_params_path)
  out_csv        <- require_param(lsstats_params, "out_csv", lsstats_params_path)
  raster_files   <- trimws(strsplit(require_param(lsstats_params, "rasters", lsstats_params_path),
                                    ",", fixed = TRUE)[[1]])
  num_rasters    <- sum(nzchar(raster_files))
}
if (num_rasters == 0) stop("LSstats has no rasters to summarize -- check ", lsstats_params_path, ".")

# Comparable form of a raster path: extension-free, absolute, and lowercase (Windows paths are
# case-insensitive), so "c:\Data\patch_x" and "C:/data/patch_x.flt" compare equal.
raster_key <- function(path) {
  tolower(normalizePath(sub("\\.(flt|tif)$", "", path, ignore.case = TRUE), winslash = "/", mustWork = FALSE))
}

if (raster_key(out_patch) != raster_key(polygon_raster)) {
  stop("LShunter's output patch raster and LSstats's landslide polygon raster are not the same ",
       "file, so LSstats would not summarize the patches this pipeline's LShunter run produces:\n",
       "  LShunter output patch raster:     ", out_patch, "\n",
       "  LSstats landslide polygon raster: ", polygon_raster)
}

# LSstats's patch table (volumes) must be the one this LShunter run writes. Optional in both.
optional_value <- function(stage_params, mode1_file, key, keyword) {
  if (is_mode1(stage_params)) {
    v <- read_keyword(mode1_file, keyword)
    if (length(v)) v[1] else "NOFILE"
  } else if (is.null(stage_params[[key]])) "NOFILE" else stage_params[[key]]
}
lsstats_table <- optional_value(lsstats_params, if (is_mode1(lsstats_params)) lsstats_input_file,
                                "patch_table", "LSHUNTER PATCH TABLE")
if (!is_nofile_text(lsstats_table)) {
  lshunter_table <- optional_value(lshunter_params,
                                   if (is_mode1(lshunter_params)) lshunter_input_file,
                                   "out_table", "OUTPUT PATCH TABLE")
  csv_key <- function(p) tolower(normalizePath(paste0(sub("\\.csv$", "", p, ignore.case = TRUE), ".csv"),
                                               winslash = "/", mustWork = FALSE))
  if (is_nofile_text(lshunter_table) || csv_key(lsstats_table) != csv_key(lshunter_table)) {
    stop("LSstats's patch_table is not LShunter's out_table, so its volumes would describe other ",
         "patches:\n  LSstats patch_table: ", lsstats_table, "\n  LShunter out_table:  ", lshunter_table)
  }
}

out_patch_flt <- paste0(sub("\\.flt$", "", out_patch, ignore.case = TRUE), ".flt")
out_csv_base  <- sub("\\.csv$", "", out_csv, ignore.case = TRUE)
last_csv      <- paste0(out_csv_base, "_landforms_", num_rasters, ".csv")

## ---- Run one stage as its own Rscript subprocess ----------------------------------------------

run_stage <- function(label, script_path, param_path) {
  message("\n==== ", label, " (", if (is_mode1(read_config(param_path))) "mode 1" else "mode 2", ") ====")
  start_time <- Sys.time()
  status <- system2("Rscript", args = shQuote(c(script_path, param_path)))
  elapsed <- Sys.time() - start_time

  if (status != 0) {
    stop(label, " failed (Rscript exit status ", status, ") -- see its own output above. ",
         "Parameter file used for this run: ", param_path)
  }
  message(label, " finished in ", format(unclass(elapsed), digits = 4), " ", units(elapsed))
}

## ---- Run the pipeline ---------------------------------------------------------------------------

message("Patch raster linking LShunter -> LSstats: ", out_patch_flt)

if (file.exists(out_patch_flt)) {
  message("Patch raster already exists (", out_patch_flt, ") -- skipping LShunter.")
} else {
  run_stage("LShunter", lshunter_script, lshunter_params_path)
  if (!file.exists(out_patch_flt)) {
    stop("LShunter finished but its patch raster was not written: ", out_patch_flt)
  }
}

if (file.exists(last_csv)) {
  if (file.mtime(last_csv) < file.mtime(out_patch_flt)) {
    stop("LSstats output (", last_csv, ") is older than the patch raster (", out_patch_flt, "), so ",
         "it describes an earlier set of patches. Delete the LSstats csv files (", out_csv_base,
         "_*.csv) and rerun this pipeline to recompute them.")
  }
  message("LSstats output already exists (", last_csv, ") -- skipping LSstats.")
} else {
  run_stage("LSstats", lsstats_script, lsstats_params_path)
}

message("\nPipeline complete. Patch raster: ", out_patch_flt,
        "\nLSstats csv files: ", out_csv_base, "_<i>.csv and ", out_csv_base,
        "_landforms_<i>.csv, i = 1..", num_rasters)
