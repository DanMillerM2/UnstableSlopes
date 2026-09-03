## run_ground_returns_batches.R
##
## Driver script: runs ground_returns.R sequentially for two or more sets of LAZ tiles (e.g.
## different sites, or different acquisition dates for the same site), producing a separate
## ground-point-density raster for each set.
##
## ground_returns.R is meant to be invoked from the command line as:
##   Rscript ground_returns.R path/to/params.txt
## (see the "Parameter file" section at the top of ground_returns.R). This driver reads its own
## batch config file (see "Batch config file" below) describing one or more sets, generates one
## ground_returns.R parameter file per set, then calls that command once per set via system2(),
## waiting for each call to finish before starting the next. Runs are NOT parallelized against
## each other -- only the batches *within* one run are (via that run's own n_workers, inside
## ground_returns.R itself) -- because ground_returns.R sets up its own
## future::plan(multisession, ...) and holds an entire site's worth of batch/raster state in the
## global environment of the process it runs in; running each set as its own Rscript subprocess
## keeps that state fully isolated between sets instead of trying to source() the script
## repeatedly in one shared R session.
##
## ---- Batch config file --------------------------------------------------------------------
##
## Path given by batch_config_path below, overridable by a command-line argument (same
## convention as config_path in ground_returns.R):
##   Rscript run_ground_returns_batches.R path/to/your_batch_config.txt
##
## See run_ground_returns_batches_params_template.txt for a ready-to-copy example. Format:
## "keyword: value" lines (same rules as ground_returns.R's parameter file -- order within a
## section doesn't matter, blank/"#"-only lines are ignored, a trailing "# comment" is stripped,
## only the first colon on a line splits keyword from value) grouped under section headers:
##
##   [defaults]
##   n_workers: 6
##   batch_size_tiles: 50
##   ...
##
##   [set]
##   input_dir: c:\work\data\site1\laz
##   output_dir: c:\work\data\site1\grnd_den
##
##   [set]
##   input_dir: c:\work\data\site2\laz
##   output_dir: c:\work\data\site2\grnd_den
##   n_workers: 4                            # overrides [defaults] for this set only
##
## [defaults] is optional and, if present, must come before any [set] sections; anything set
## there applies to every [set] below unless that set repeats the same keyword itself. Each
## [set] section becomes one ground_returns.R run and must include input_dir and output_dir
## (no default -- same requirement ground_returns.R itself places on its own parameter file).
## Anything not set in either [defaults] or a [set] falls back to ground_returns.R's own
## built-in default (see param_specs in ground_returns.R).

batch_config_path <- "run_ground_returns_batches_params_template.txt"  # fallback; overridden by a command-line argument below

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) batch_config_path <- cli_args[1]  # Rscript ... batch_config.txt takes precedence over the hardcoded fallback above

if (!file.exists(batch_config_path)) {
  stop("Batch config file not found: ", batch_config_path, ". Set batch_config_path in this ",
       "script, or pass its path as a command-line argument (Rscript ... batch_config.txt).")
}

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so
# sibling scripts can be found by their location on disk rather than by the caller's working
# directory -- e.g. `Rscript R/run_ground_returns_batches.R ...` from the repo root and
# `Rscript run_ground_returns_batches.R ...` from inside R/ both correctly find
# R/ground_returns.R. Falls back to "." (assume cwd) when there's no "--file=" argument to read,
# e.g. when this script is source()'d from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

ground_returns_script <- file.path(get_script_dir(), "ground_returns.R")  # sibling script, same folder as this one
param_dir <- file.path(tempdir(), "ground_returns_params")  # per-set parameter files generated below are written here
dir.create(param_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Parse the batch config file -------------------------------------------------------------

# Reads a [defaults]/[set] batch config file into list(defaults = list(...), runs = list(...)),
# `runs` holding one named list of "keyword: value" pairs (all still raw character strings --
# ground_returns.R's own build_params() does the actual type conversion, once these are written
# out as that script's parameter file) per [set] section encountered, in file order.
read_batch_config <- function(path) {
  raw_lines <- readLines(path, warn = FALSE)
  raw_lines <- trimws(raw_lines)
  raw_lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]  # drop blank lines and comment-only lines

  defaults <- list()
  runs <- list()
  active <- NULL              # NULL until the first section header; then "defaults" or "set"
  current_keys <- character(0)  # keys seen so far in the section currently being filled, to catch duplicates within it

  for (line in raw_lines) {
    section_match <- regmatches(line, regexec("^\\[(.+)\\]$", line))[[1]]
    if (length(section_match) == 2) {
      section_name <- tolower(trimws(section_match[2]))
      if (section_name == "defaults") {
        if (length(runs) > 0) {
          stop("[defaults] must come before any [set] section in ", path)
        }
        active <- "defaults"
      } else if (section_name == "set") {
        runs[[length(runs) + 1]] <- list()
        active <- "set"
      } else {
        stop("Unrecognized section header in ", path, ": [", section_match[2],
             "] -- expected [defaults] or [set].")
      }
      current_keys <- character(0)  # reset duplicate-tracking for the new section
      next
    }

    if (is.null(active)) {
      stop("Parameter line found before any [defaults] or [set] section header in ", path,
           ": ", line)
    }

    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the FIRST colon splits keyword from value -- keeps values like "c:\..." or "EPSG:2856" intact
    if (colon_pos < 0) {
      stop("Malformed line in ", path, " (expected \"keyword: value\" or a [section] header): ",
           line)
    }
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(substr(line, colon_pos + 1, nchar(line)))
    val <- trimws(sub("\\s+#.*$", "", val))  # strip a trailing " # comment", if any, from the value

    if (key %in% current_keys) {
      stop("Duplicate keyword '", key, "' within one section of ", path)
    }
    current_keys <- c(current_keys, key)

    if (active == "defaults") {
      defaults[[key]] <- val
    } else {
      runs[[length(runs)]][[key]] <- val
    }
  }

  if (length(runs) == 0) {
    stop("No [set] sections found in ", path, " -- need at least one set of LAZ tiles to process.")
  }

  list(defaults = defaults, runs = runs)
}

batch_config <- read_batch_config(batch_config_path)
shared_defaults <- batch_config$defaults
runs <- batch_config$runs

message("Loaded ", length(runs), " set(s) from batch config: ", batch_config_path)

## ---- Helpers -----------------------------------------------------------------------------

# Merges one run's own parameters over shared_defaults (run-specific values win), writes the
# result out as a "keyword: value" parameter file ground_returns.R can read, and returns that
# file's path.
write_param_file <- function(run, defaults, path) {
  params <- modifyList(defaults, run)
  if (is.null(params$input_dir) || is.null(params$output_dir)) {
    stop("Each [set] must set input_dir and output_dir.")
  }
  lines <- vapply(names(params), function(nm) paste0(nm, ": ", params[[nm]]), character(1))
  writeLines(lines, path)
  path
}

## ---- Run ground_returns.R once per set, sequentially ------------------------------------------
## Each iteration blocks (system2() waits for the subprocess to exit) before starting the next
## set, so sets never compete for CPU/memory with each other -- only a set's own batches are
## parallelized, inside that set's ground_returns.R run.

results <- vector("list", length(runs))

for (i in seq_along(runs)) {
  run <- runs[[i]]
  message("\n==== Set ", i, " of ", length(runs), ": ", run$input_dir, " ====")

  param_path <- file.path(param_dir, sprintf("params_set%02d.txt", i))
  write_param_file(run, shared_defaults, param_path)

  start_time <- Sys.time()
  status <- system2("Rscript", args = shQuote(c(ground_returns_script, param_path)))
  elapsed <- Sys.time() - start_time

  if (status != 0) {
    stop("ground_returns.R failed for set ", i, " (input_dir: ", run$input_dir,
         ") -- exit status ", status, ". Stopping before remaining sets.")
  }

  message("Set ", i, " finished in ", format(unclass(elapsed), digits = 4), " ", units(elapsed))
  results[[i]] <- list(run          = run,
                       param_path   = param_path,
                       output_dir   = run$output_dir,
                       density_tif  = file.path(run$output_dir, "ground_density_epsg2856.tif"),
                       elapsed      = elapsed)
}

message("\nAll ", length(runs), " set(s) completed. Density rasters written to:")
for (i in seq_along(results)) {
  message(i, ": ", results[[i]]$density_tif)
}
