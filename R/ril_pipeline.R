## ril_pipeline.R
##
## Driver script: runs bldgrds_enforce.R, then valleyfloor.R, then RIL.R, each as its own
## `Rscript <script> <params.txt>` subprocess (same convention run_align_pipeline.R uses for
## ground_returns.R) -- NOT by importing/duplicating their own (large) parameter-file grammars.
## bldgrds_enforce.R/valleyfloor.R/RIL.R each remain fully usable standalone, completely
## unmodified -- this driver only ever writes them a parameter file and invokes them exactly as a
## user would by hand, so it can never jeopardize (or drift out of sync with) any of their own
## independent behavior.
##
## Why chaining these three works at all (confirmed directly against the current
## ChannelUtilities/GridUtilities Fortran source, not just inferred from the R wrappers):
##   - bldgrds (enforce mode) traces the channel network from channel_mask into the DEM and writes
##     it out as a node-list database, NodeNet_<ID>.dat / NodeAttributes_<ID>.dat(.hdr),
##     unconditionally, next to the DEM -- <ID> is derived from the DEM's own filename (everything
##     after its first "_"), not from any keyword bldgrds_enforce.R exposes.
##   - ValleyFloor (valleyfloor.R) reads that same NodeNet_<ID>.dat back in (via the shared
##     ChannelNode_Module -- there is no dedicated "node file" keyword; it's found purely by DEM
##     path + DEM-derived ID) to get its channel definition, then writes valleyfloor_<ID>.dat
##     (per-channel height-/depth-above-channel data), again unconditionally, next to the DEM.
##   - RIL (RIL.R) reads BOTH of those back in, and unconditionally so: unless in_valley_floor /
##     "INPUT VALLEY FLOOR RASTER" is explicitly supplied, RIL.f90 opens
##     <DEM's folder>valleyfloor_<ID>.dat directly (`datafile =
##     TRIM(DEM%path)//'valleyfloor_'//TRIM(DEM%dataID)//'.dat'`) and reads it via
##     `readValleyHeader()`/`readChanRecord()` -- the same valleyFloor_Module routines
##     ValleyFloor.f90 itself uses to write it -- and *aborts* ("Error opening binary data file")
##     if that file is missing. So valleyfloor.R must have already run, against the SAME dem,
##     before RIL.R does.
##
## The one thing that makes this whole chain work -- and the one thing this driver exists to
## guarantee -- is that all three stages are pointed at the exact same DEM file. None of
## bldgrds_enforce()/valleyfloor()/RIL() (the TWutils wrappers) takes an explicit "node file",
## "valleyfloor file", or "previous stage's output" argument; each derives its own working
## directory/data-ID purely from its own `dem` argument's path and filename. Get that DEM path
## wrong (or use a different one) in any of the three per-stage sections below and the chain
## silently breaks -- RIL either aborts outright ("Cannot find nodenet file" /
## "Error opening binary data file") or, worse, silently reads a stale database left over from an
## unrelated previous run in the same folder. So this driver takes exactly ONE dem: value
## ([dem] section, below) and forces it onto all three stages -- a `dem:` (or `input_file:`, which
## would let a stage skip building from parameters -- and so skip using that shared dem entirely --
## or `data_id:`, which would break the DEM-name-derived ID every stage otherwise shares
## automatically) line inside [bldgrds_enforce]/[valleyfloor]/[ril] is dropped, with a warning,
## rather than honored.
##
## Resumability: RIL.R decides for itself whether to skip its work (it skips if out_RIL's .flt
## already exists -- see its own header comment for why), and since it runs as that exact script
## via Rscript, that behavior carries over here unchanged. bldgrds_enforce.R and valleyfloor.R have
## no such check of their own (by their own design -- again see their own header comments), so
## THIS driver adds one for each instead, checking for their DEM-derived output file
## (NodeNet_<ID>.dat / valleyfloor_<ID>.dat, respectively -- see the pre-flight check just above
## the "Run the pipeline" section below) before invoking that stage at all, and noting the skip if
## found.
##
## NOTE ON BUILD FRESHNESS: the ValleyFloor -> RIL file handoff described above (RIL.f90 reading
## valleyfloor_<ID>.dat) is very recent, active work in ChannelUtilities (its most recent commits
## touch exactly this) -- confirm valleyfloor.exe/RIL.exe at the executable_dir values below are
## freshly rebuilt from current source before trusting this pipeline's output; see
## ChannelUtilities's own CLAUDE.md and this repo's own conventions note on executable_dir/
## program_name pointing at a sibling repo's build output, not anything in this repo.
##
## Requires: TWutils (local package -- must already be installed/available in the R environment),
## and bldgrds_enforce.R/valleyfloor.R/RIL.R present alongside this script (same R/ folder).
##
## ---- Pipeline config file --------------------------------------------------------------------
##
## One plain-text config file (see ril_pipeline_params_template.txt for a
## ready-to-copy example), read from pipeline_config_path below or a command-line argument:
##   Rscript ril_pipeline.R path/to/your_pipeline_config.txt
##
## "keyword: value" lines (same rules as the three per-stage scripts' own parameter files -- order
## doesn't matter, blank/"#"-only lines are ignored, a trailing "# comment" is stripped, only the
## first colon on a line splits keyword from value) grouped under four section headers:
##
##   [dem]              -- exactly one line, dem: <path>. Used for every stage below.
##   [bldgrds_enforce]   -- every bldgrds_enforce.R parameter EXCEPT dem/input_file (see
##                          param_specs in bldgrds_enforce.R for the full keyword list).
##   [valleyfloor]       -- every valleyfloor.R parameter except dem/input_file/data_id (see
##                          param_specs in valleyfloor.R). write_data: TRUE and method: 4 (both
##                          valleyfloor.R's own defaults) matter here specifically: they're what
##                          make ValleyFloor actually compute and write the depth-above-channel
##                          data RIL reads back out of valleyfloor_<ID>.dat.
##   [ril]               -- every RIL.R parameter except dem/input_file (see param_specs in
##                          RIL.R). out_RIL (below) is this pipeline's final product.
##
## Each section's parameters are written straight through, verbatim, into that stage's own
## "keyword: value" parameter file (with one dem: line injected at the top) -- so anything
## bldgrds_enforce.R/valleyfloor.R/RIL.R's own param_specs recognizes can be set here, and any
## keyword neither recognizes is rejected by that script itself (with the same error it would give
## running standalone), not silently accepted by this driver.

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
## Same [section]-headed grammar/parser as run_align_pipeline.R's own read_sectioned_config():
## reads a [section]-headed config file into a named list of sections, each itself a named list of
## raw (still character) "keyword: value" pairs. Every section name here is expected exactly once.

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

    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the FIRST colon splits keyword from value -- keeps values like "c:\..." intact
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

dem_section         <- config[["dem"]]
bldgrds_section     <- config[["bldgrds_enforce"]]
valleyfloor_section <- config[["valleyfloor"]]
ril_section         <- config[["ril"]]

if (is.null(dem_section) || is.null(dem_section[["dem"]])) {
  stop("Pipeline config file must contain a [dem] section with a 'dem:' line: ", pipeline_config_path)
}
if (is.null(bldgrds_section) || is.null(valleyfloor_section) || is.null(ril_section)) {
  stop("Pipeline config file must contain [bldgrds_enforce], [valleyfloor], and [ril] sections: ",
       pipeline_config_path)
}

dem_path <- dem_section[["dem"]]

extra_dem_keys <- setdiff(names(dem_section), "dem")
if (length(extra_dem_keys) > 0) {
  warning("Ignoring extra keyword(s) in [dem] of ", pipeline_config_path, ": ",
          paste(extra_dem_keys, collapse = ", "), " -- [dem] only takes a single 'dem:' line.")
}

## ---- Write one stage's own "keyword: value" parameter file, forcing the shared dem: value ----
## Drops (with a warning) any dem/input_file/data_id line a section tried to set on its own --
## every stage must use the one shared dem_path above (input_file would let a stage skip building
## from parameters -- and so skip using dem_path -- entirely; a per-stage data_id would break the
## DEM-name-derived ID every stage otherwise shares automatically). Everything else is passed
## through untouched, so any keyword bldgrds_enforce.R/valleyfloor.R/RIL.R's own param_specs
## recognizes can be set here.

write_stage_param_file <- function(stage_name, section_params, path) {
  for (forced_out in c("dem", "input_file", "data_id")) {
    if (!is.null(section_params[[forced_out]])) {
      warning("[", stage_name, "]'s own '", forced_out, ":' line is ignored -- every stage in ",
              "this pipeline always uses the single dem: from [dem], with no per-stage override, ",
              "so the DEM-derived data ID stays identical across all three stages.")
      section_params[[forced_out]] <- NULL
    }
  }

  lines <- c(paste0("dem: ", dem_path),
             vapply(names(section_params), function(nm) paste0(nm, ": ", section_params[[nm]]),
                    character(1)))
  writeLines(lines, path)
  path
}

## ---- Run one stage as its own Rscript subprocess ----------------------------------------------
## Reuses that stage's own script completely unmodified -- including its own parameter parsing,
## defaults, validation, and (for RIL.R) its own "skip if output already exists" resumability
## check -- so this driver can never drift out of sync with whichever of the three is most
## recently updated, and each one stays exactly as usable standalone as it always was.

run_stage <- function(label, script_path, section_params) {
  param_path <- file.path(param_dir, paste0("params_", label, ".txt"))
  write_stage_param_file(label, section_params, param_path)

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

## ---- Pre-flight: skip a stage whose DEM-derived output already exists -------------------------
## bldgrds (enforce mode) writes NodeNet_<ID>.dat/NodeAttributes_<ID>.dat unconditionally next to
## the DEM, and valleyfloor writes valleyfloor_<ID>.dat unconditionally next to the DEM (see the
## header comment above) -- unlike RIL.R's own out_RIL check, neither bldgrds_enforce.R nor
## valleyfloor.R has a "skip if already done" resumability check of its own, so this driver adds
## one for each, keyed off the exact file each Fortran program itself reads/writes. <ID> is the
## DEM-derived data ID (DEM_module's resolveDEMname(): everything after the first underscore in
## the DEM's extension-free base name, or the whole base name if there's none) -- the same ID
## TWutils::valleyfloor_dat_file() already resolves, and the one every stage here is guaranteed to
## share since this pipeline never lets a per-stage data_id override it (see
## write_stage_param_file() above).

nodenet_dat_file <- function(dem_path) {
  # NodeNet_<ID>.dat is written by bldgrds (bldGrds2.f90: TRIM(DEM%path)//'NodeNet_'//
  # TRIM(DEM%DEMID)//'.dat') using the same DEM-derived ID as ValleyFloor's own
  # valleyfloor_<ID>.dat. Reuse TWutils::valleyfloor_dat_file()'s directory/ID resolution
  # (extension-stripping, path normalization, resolveDEMname() logic) rather than duplicating it
  # here, and just swap the file-name prefix.
  valleyfloor_path <- TWutils::valleyfloor_dat_file(dem_path)
  sub("valleyfloor_([^\\\\/]+)\\.dat$", "NodeNet_\\1.dat", valleyfloor_path)
}

nodenet_path     <- nodenet_dat_file(dem_path)
valleyfloor_path <- TWutils::valleyfloor_dat_file(dem_path)

## ---- Run the pipeline ---------------------------------------------------------------------------

message("Pipeline DEM (shared by all three stages, and the source of the DEM-derived data ID ",
        "each stage's node/channel database is keyed by): ", dem_path)

if (file.exists(nodenet_path)) {
  message("NodeNet file already exists (", nodenet_path, ") -- skipping bldgrds_enforce.")
} else {
  run_stage("bldgrds_enforce", bldgrds_enforce_script, bldgrds_section)
}

if (file.exists(valleyfloor_path)) {
  message("valleyfloor file already exists (", valleyfloor_path, ") -- skipping valleyfloor.")
} else {
  run_stage("valleyfloor", valleyfloor_script, valleyfloor_section)
}

run_stage("RIL", ril_script, ril_section)

## ---- Report the final output -------------------------------------------------------------------
## RIL.R requires out_RIL in its own [ril] section (same as running RIL.R by hand) -- reaching
## this point means that run already succeeded, so out_RIL is guaranteed to be set. The landform
## raster it names is this whole pipeline's final product.

out_RIL <- ril_section[["out_RIL"]]
out_RIL_flt <- if (grepl("\\.flt$", out_RIL, ignore.case = TRUE)) out_RIL else paste0(out_RIL, ".flt")

message("\nPipeline complete. Final RIL landform raster: ", out_RIL_flt)
