## ril_pipeline.R
##
## Driver script: runs bldgrds_enforce.R, then valleyfloor.R, then RIL.R, each in mode 1 (an
## existing, user-prepared Fortran input file) and each as its own `Rscript <script> <params.txt>`
## subprocess (same convention run_align_pipeline.R uses for ground_returns.R).
## bldgrds_enforce.R/valleyfloor.R/RIL.R each remain fully usable standalone, completely
## unmodified -- this driver only ever writes each of them a two-line mode-1 parameter file
## (input_file: plus, optionally, executable_dir:) and invokes them exactly as a user would by hand.
##
## Why chaining these three works at all (confirmed directly against the current
## ChannelUtilities/GridUtilities Fortran source, not just inferred from the R wrappers):
##   - bldgrds (enforce mode) traces the channel network from its channel mask into the DEM and
##     writes it out as a node-list database, NodeNet_<ID>.dat / NodeAttributes_<ID>.dat(.hdr),
##     unconditionally, next to the DEM -- <ID> is derived from the DEM's own filename (everything
##     after its first "_").
##   - ValleyFloor reads that same NodeNet_<ID>.dat back in (via the shared ChannelNode_Module --
##     there is no dedicated "node file" keyword; it's found purely by DEM path + DEM-derived ID)
##     to get its channel definition, then writes valleyfloor_<ID>.dat (per-channel
##     height-/depth-above-channel data), again unconditionally, next to the DEM.
##   - RIL reads BOTH of those back in: unless "INPUT VALLEY FLOOR RASTER" is supplied, RIL.f90
##     opens <DEM's folder>valleyfloor_<ID>.dat directly (`datafile =
##     TRIM(DEM%path)//'valleyfloor_'//TRIM(DEM%dataID)//'.dat'`) and *aborts* ("Error opening
##     binary data file") if that file is missing. So ValleyFloor must have already run, against
##     the SAME dem, before RIL does.
##
## The one thing that makes this chain work is that all three input files name the exact same
## DEM. None of the three programs takes an explicit "node file" or "previous stage's output"
## keyword; each derives its working directory/data ID purely from its own DEM: line. Point one of
## the three input files at a different DEM and the chain silently breaks -- RIL either aborts
## outright ("Cannot find nodenet file" / "Error opening binary data file") or, worse, silently
## reads a stale database left over from an unrelated previous run in the same folder. Since this
## driver no longer writes those input files itself, it instead reads the DEM: line back out of
## each one before running anything, and stops if they don't all match. It likewise stops if the
## ValleyFloor input file sets a DATA ID (on its ALL CHANNELS line) that differs from the
## DEM-derived ID, since RIL would then look for a valleyfloor_<ID>.dat that was never written.
##
## Resumability: none of bldgrds_enforce.R, valleyfloor.R, or RIL.R has a "skip if already done"
## check of its own in mode 1 (RIL.R's out_RIL check is mode-2 only -- see its header comment), so
## this driver adds one for each stage, keyed off the exact file each Fortran program writes:
## NodeNet_<ID>.dat for bldgrds, valleyfloor_<ID>.dat for ValleyFloor, and the RIL input file's own
## OUTPUT RIL RASTER (.flt) for RIL. Delete that file to force a stage to rerun.
##
## NOTE ON BUILD FRESHNESS: the ValleyFloor -> RIL file handoff described above (RIL.f90 reading
## valleyfloor_<ID>.dat) is very recent, active work in ChannelUtilities -- confirm
## valleyfloor.exe/RIL.exe are freshly rebuilt from current source before trusting this pipeline's
## output; see ChannelUtilities's own CLAUDE.md and this repo's own conventions note on
## executable_dir/program_name pointing at a sibling repo's build output.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment),
## and bldgrds_enforce.R/valleyfloor.R/RIL.R present alongside this script (same R/ folder).
##
## ---- Pipeline config file --------------------------------------------------------------------
##
## One plain-text config file (see ril_pipeline_params_template.txt for a ready-to-copy example),
## read from pipeline_config_path below or a command-line argument:
##   Rscript ril_pipeline.R path/to/your_pipeline_config.txt
##
## "keyword: value" lines (same rules as the three per-stage scripts' own parameter files -- order
## doesn't matter, blank/"#"-only lines are ignored, a trailing "# comment" is stripped, only the
## first colon on a line splits keyword from value):
##
##   bldgrds_input_file:          existing bldgrds "enforce" input file (required)
##   valleyfloor_input_file:      existing ValleyFloor input file (required)
##   ril_input_file:              existing RIL input file (required)
##   bldgrds_executable_dir:      folder containing bldgrds.exe (optional, default NOFILE)
##   valleyfloor_executable_dir:  folder containing ValleyFloor.exe (optional, default NOFILE)
##   ril_executable_dir:          folder containing RIL.exe (optional, default NOFILE)
##
## Each *_executable_dir is only a fallback, exactly as in each stage's own mode 1: an input file's
## own working "EXECUTABLE DIR:" line wins (see TWutils::resolve_executable_dir()). Set it when the
## input file has no such line, or the stage won't find its executable.
##
## The input files themselves use the Fortran programs' own `KEYWORD: value[, SUBFIELD=value]`
## grammar -- e.g. one written by a prior mode-2 run of bldgrds_enforce.R/valleyfloor.R/RIL.R (see
## scratch_dir), or a hand-written/edited one.

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
                                  "ril_pipeline_params_template.txt")  # fallback; overridden by a command-line argument below

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) pipeline_config_path <- cli_args[1]  # Rscript ... pipeline_config.txt takes precedence over the hardcoded fallback above

if (!file.exists(pipeline_config_path)) {
  stop("Pipeline config file not found: ", pipeline_config_path, ". Set pipeline_config_path in ",
       "this script, or pass its path as a command-line argument (Rscript ... pipeline_config.txt).")
}

# Sibling scripts, same folder as this one -- each is run unmodified, via its own Rscript
# invocation, exactly as if a user had typed `Rscript <script> <params.txt>` by hand.
bldgrds_enforce_script <- file.path(get_script_dir(), "bldgrds_enforce.R")
valleyfloor_script     <- file.path(get_script_dir(), "valleyfloor.R")
ril_script             <- file.path(get_script_dir(), "RIL.R")

for (script_path in c(bldgrds_enforce_script, valleyfloor_script, ril_script)) {
  if (!file.exists(script_path)) {
    stop("Cannot find ", basename(script_path), " next to this pipeline script (expected at ",
         script_path, ").")
  }
}

# Generated per-stage parameter files are written here (one per stage per run).
param_dir <- file.path(tempdir(), "ril_pipeline_params")
dir.create(param_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Parse the pipeline config file -----------------------------------------------------------

is_nofile_text <- function(x) is.null(x) || !nzchar(x) || toupper(x) == "NOFILE"

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

stage_names   <- c("bldgrds", "valleyfloor", "ril")
required_keys <- paste0(stage_names, "_input_file")
optional_keys <- paste0(stage_names, "_executable_dir")

config <- read_config(pipeline_config_path)

unknown_keys <- setdiff(names(config), c(required_keys, optional_keys))
if (length(unknown_keys) > 0) {
  stop("Unrecognized keyword(s) in ", pipeline_config_path, ": ", paste(unknown_keys, collapse = ", "),
       ". Expected: ", paste(c(required_keys, optional_keys), collapse = ", "), ".")
}

for (key in required_keys) {
  if (is_nofile_text(config[[key]])) {
    stop("Pipeline config file must set '", key, ":' to an existing input file: ", pipeline_config_path)
  }
  if (!file.exists(config[[key]])) {
    stop(key, " not found: ", config[[key]])
  }
}

input_files     <- setNames(vapply(required_keys, function(k) config[[k]], character(1)), stage_names)
executable_dirs <- setNames(vapply(optional_keys,
                                   function(k) if (is_nofile_text(config[[k]])) "NOFILE" else config[[k]],
                                   character(1)),
                            stage_names)

## ---- Read the DEM (and other linking keywords) back out of each input file ---------------------
## Keywords are matched case-insensitively and anchored, so e.g. ValleyFloor's "FILL DEM HOLES:"
## isn't mistaken for "DEM:" (TWutils::get_dem()'s own pattern is unanchored). Comment lines are
## skipped by TWutils::get_keyword().

read_keyword <- function(input_file, keyword) {
  pattern <- paste0("^", gsub(" ", "\\\\s+", keyword), "$")
  lines <- readLines(input_file, warn = FALSE)
  keywords <- vapply(seq_along(lines), function(i) TWutils::get_keyword(lines, i), character(1))
  matches <- which(!is.na(keywords) & grepl(pattern, keywords, ignore.case = TRUE))
  # get_args() keeps a trailing " # comment", which the Fortran grammar allows -- strip it first.
  vapply(matches, function(i) TWutils::parse_arg(sub("\\s+#.*$", "", TWutils::get_args(lines, i)))$Value,
         character(1))
}

read_single_keyword <- function(stage, keyword) {
  values <- read_keyword(input_files[[stage]], keyword)
  if (length(values) != 1) {
    stop("Expected exactly one '", keyword, ":' line in the ", stage, " input file (",
         input_files[[stage]], "), found ", length(values), ".")
  }
  values
}

# Comparable form of a DEM path: extension-free, absolute, and lowercase (Windows paths are
# case-insensitive), so "c:\Data\elev_x" and "C:/data/elev_x.flt" compare equal.
dem_key <- function(dem) {
  tolower(normalizePath(sub("\\.flt$", "", dem, ignore.case = TRUE), winslash = "/", mustWork = FALSE))
}

stage_dems <- vapply(stage_names, function(stage) read_single_keyword(stage, "DEM"), character(1))

if (length(unique(vapply(stage_dems, dem_key, character(1)))) != 1) {
  stop("The three input files do not all name the same DEM -- every stage of this pipeline must ",
       "use the same DEM so its DEM-derived data ID links bldgrds -> ValleyFloor -> RIL:\n",
       paste0("  ", stage_names, ": ", stage_dems, collapse = "\n"))
}
dem_path <- stage_dems[["bldgrds"]]

# RIL always opens valleyfloor_<DEM-derived ID>.dat, so a ValleyFloor DATA ID (a subfield of its
# ALL CHANNELS line) that differs from that ID would write a file RIL never finds.
valleyfloor_path <- TWutils::valleyfloor_dat_file(dem_path)
vf_lines <- readLines(input_files[["valleyfloor"]], warn = FALSE)
vf_lines <- vf_lines[!grepl("^\\s*#", vf_lines)]
data_id_match <- regmatches(vf_lines, regexec("(?i)DATA\\s+ID\\s*=\\s*([^,#]+)", vf_lines, perl = TRUE))
data_ids <- unique(trimws(unlist(lapply(data_id_match, function(m) if (length(m) == 2) m[2]))))
if (length(data_ids) > 0) {
  vf_with_id <- TWutils::valleyfloor_dat_file(dem_path, data_id = data_ids[1])
  if (length(data_ids) > 1 || tolower(vf_with_id) != tolower(valleyfloor_path)) {
    stop("The ValleyFloor input file (", input_files[["valleyfloor"]], ") sets DATA ID = ",
         paste(data_ids, collapse = ", "), ", so it would write ", vf_with_id, " -- but RIL reads ",
         valleyfloor_path, " (from the DEM-derived ID). Remove the DATA ID subfield for this pipeline.")
  }
}

out_RIL <- read_single_keyword("ril", "OUTPUT RIL RASTER")
out_RIL_flt <- if (grepl("\\.flt$", out_RIL, ignore.case = TRUE)) out_RIL else paste0(out_RIL, ".flt")

## ---- Run one stage as its own Rscript subprocess ----------------------------------------------
## Reuses that stage's own script completely unmodified, in mode 1: the parameter file written
## here holds only input_file: and (if set) executable_dir:.

run_stage <- function(label, script_path, stage) {
  param_path <- file.path(param_dir, paste0("params_", label, ".txt"))
  lines <- paste0("input_file: ", input_files[[stage]])
  if (!is_nofile_text(executable_dirs[[stage]])) {
    lines <- c(lines, paste0("executable_dir: ", executable_dirs[[stage]]))
  }
  writeLines(lines, param_path)

  message("\n==== ", label, " ====")
  start_time <- Sys.time()
  status <- system2("Rscript", args = shQuote(c(script_path, param_path)))
  elapsed <- Sys.time() - start_time

  if (status != 0) {
    stop(label, " failed (Rscript exit status ", status, ") -- see its own output above. ",
         "Parameter file used for this run: ", param_path)
  }
  message(label, " finished in ", format(unclass(elapsed), digits = 4), " ", units(elapsed))
}

## ---- Pre-flight: skip a stage whose output already exists -------------------------------------

nodenet_dat_file <- function(dem_path) {
  # NodeNet_<ID>.dat is written by bldgrds (bldGrds2.f90: TRIM(DEM%path)//'NodeNet_'//
  # TRIM(DEM%DEMID)//'.dat') using the same DEM-derived ID as ValleyFloor's own
  # valleyfloor_<ID>.dat. Reuse TWutils::valleyfloor_dat_file()'s directory/ID resolution
  # (extension-stripping, path normalization, resolveDEMname() logic) rather than duplicating it
  # here, and just swap the file-name prefix.
  sub("valleyfloor_([^\\\\/]+)\\.dat$", "NodeNet_\\1.dat", TWutils::valleyfloor_dat_file(dem_path))
}

nodenet_path <- nodenet_dat_file(dem_path)

## ---- Run the pipeline ---------------------------------------------------------------------------

message("Pipeline DEM (named by all three input files, and the source of the DEM-derived data ID ",
        "each stage's node/channel database is keyed by): ", dem_path)

if (file.exists(nodenet_path)) {
  message("NodeNet file already exists (", nodenet_path, ") -- skipping bldgrds_enforce.")
} else {
  run_stage("bldgrds_enforce", bldgrds_enforce_script, "bldgrds")
}

if (file.exists(valleyfloor_path)) {
  message("valleyfloor file already exists (", valleyfloor_path, ") -- skipping valleyfloor.")
} else {
  run_stage("valleyfloor", valleyfloor_script, "valleyfloor")
}

if (file.exists(out_RIL_flt)) {
  message("RIL output raster already exists (", out_RIL_flt, ") -- skipping RIL.")
} else {
  run_stage("RIL", ril_script, "ril")
}

message("\nPipeline complete. Final RIL landform raster: ", out_RIL_flt)
