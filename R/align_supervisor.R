## align_supervisor.R
##
## Runs the Fortran program Align once on one reference/align DTM pair (TWutils::align_input()
## writes its input file), follows its progress, checks the executable supports every option
## asked for, and writes tables of the run for the pipeline report. Run on its own
## (Rscript align_supervisor.R params.txt) or from run_align_pipeline.R, whose Stage 3 always
## runs it.
##
## Needs Align at GridUtilities 8a08fa6 or later (alignSlope2 with align_modules\, built from
## Projects\align). The settings echo is checked a few seconds into the run: if Align reports
## unrecognized keywords, or doesn't echo an option set here (the convergence test on the change
## of the aligned DEM, cross-validated tile selection, central derivatives, the gradient raster,
## the density LoD), the executable is too old and the run is stopped.
##
## Tiles. With tileNx/tileNy: auto (the default), Align chooses the tile count and the prior
## strength itself by block cross-validation (TILES FIT = JOINT, SELECT = AUTO): square tiles from
## max_tile_size down to min_tile_size, scored on held-out 1 km blocks with a 100 m buffer (Align's
## defaults; cv_block, cv_buffer). min_tile_size is a judgement, not a fit: a shift field finer
## than the scale of real change (a large slow landslide, widespread deposition) can absorb it.
## With numbers, those tiles are used (tile_select NONE, or REPORT/CV to cross-validate around
## them).
##
## Iterations. Each loop (global, then tiles) runs at most `iterations` iterations and stops
## once the aligned DEM changes by less than `convergence` m: the 99th percentile of the change on
## steep ground and the predicted remaining change (Align's CONVERGENCE). 8 and 0.005 make the
## results independent of the iteration count (see Align's CLAUDE.md).
##
## Outputs. Align writes to <dirname(outOutlier)>/align_run/ (its scratch log to
## <scratch_dir>/align_run/), with its console output in align_run/align.log, streamed to this
## session's console as it runs. On success every output is moved to its final path (by base
## name: <outDif>_0, <outLod>_gain/_centre/_iqr, <outbins>_fit0.csv, ... go with their base;
## outOutlier's .flt last, since its existence marks the run complete). The gradient raster
## LShunter needs is written by Align next to the reference DTM, <refDTM>_grad<radius>m. Tables
## in align_run/, for the pipeline report:
##   align_run.csv        -- one row: version, chosen tiles and prior scale, iterations and
##                           convergence of each loop, starting/ending RMSE, run time;
##   align_iterations.csv -- every iteration of both loops: in-fence mean and RMSE, the change of
##                           the aligned DEM (P50/P99/max, all and steep points), its ratio to the
##                           previous P99 and the predicted remaining change; global shifts;
##                           tile-node increments and nodes held by their neighbours;
##   align_select.csv     -- the cross-validation table: grid, tile size, unknowns, prior scale,
##                           held-out gain (with and without the buffer), held-out RMS and
##                           block-mean spread, best/chosen;
##   align_nodes.csv      -- the joint fit's node totals relative to the global solution (mm),
##                           with each node's map position.
## Resumability: if outOutlier already exists nothing runs; a finished run in align_run/ made from
## the same input file is not rerun, only its outputs moved.
##
## Requires: TWutils, terra, processx.

library(TWutils)

## ---- Parameter file ------------------------------------------------------------------------
## "keyword: value" lines, as in align_dtms.R (order doesn't matter, "#" comments, only the
## first colon splits keyword from value). See param_specs below for the keywords and defaults;
## run_align_pipeline.R writes this file from its [align] section.

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

config_path <- file.path(get_script_dir(), "align_supervisor_params_template.txt")  # fallback; overridden by a command-line argument
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) config_path <- cli_args[1]
if (!file.exists(config_path)) {
  stop("Parameter file not found: ", config_path, ". Pass its path as a command-line argument ",
       "(Rscript align_supervisor.R params.txt).")
}

read_config <- function(path) {
  raw_lines <- trimws(readLines(path, warn = FALSE))
  raw_lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]
  params <- list()
  for (line in raw_lines) {
    colon_pos <- regexpr(":", line, fixed = TRUE)
    if (colon_pos < 0) stop("Malformed line in ", path, " (expected \"keyword: value\"): ", line)
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(sub("\\s+#.*$", "", trimws(substr(line, colon_pos + 1, nchar(line)))))
    if (key %in% names(params)) stop("Duplicate keyword '", key, "' in ", path)
    params[[key]] <- val
  }
  params
}

# A keyword with no default (and not required) is left out when not given, so Align's own
# default applies.
param_specs <- list(
  refDTM    = list(type = "character", required = TRUE),
  alignDTM  = list(type = "character", required = TRUE),
  refGrnd   = list(type = "character", default = "NOFILE"),
  alignGrnd = list(type = "character", default = "NOFILE"),
  refDSM    = list(type = "character", default = "NOFILE"),
  alignDSM  = list(type = "character", default = "NOFILE"),

  iterations  = list(type = "integer", default = 8),      # maximum iterations per loop
  convergence = list(type = "numeric", default = 0.005),  # m; 0 runs every iteration
  k           = list(type = "numeric", default = 1.5),
  dampener    = list(type = "numeric", default = 1),
  derivatives = list(type = "character", default = "CENTRAL"),  # SIDED or CENTRAL
  applySlope  = list(type = "character"),                       # YES (Align's default), NO, STABLE
  weights     = list(type = "character"),                       # NONE (Align's default) or IQR
  groundMaskRadius   = list(type = "numeric", default = 0),
  waterGradThreshold = list(type = "numeric"),
  waterCurvThreshold = list(type = "numeric"),

  tileNx  = list(type = "character", default = "auto"),  # "auto" or a number of tiles
  tileNy  = list(type = "character", default = "auto"),
  overlap = list(type = "numeric", default = 0.5),
  tile_fit        = list(type = "character", default = "JOINT"),  # JOINT or SEPARATE
  tile_select     = list(type = "character"),  # default AUTO with auto tiles, else NONE
  cv_block        = list(type = "numeric"),    # m (Align's default 1000)
  cv_buffer       = list(type = "numeric"),    # m (100)
  min_tile_size   = list(type = "numeric"),    # m (500)
  max_tile_size   = list(type = "numeric"),    # m (8000)
  prior_points    = list(type = "numeric"),
  prior_fraction  = list(type = "numeric"),
  prior           = list(type = "character"),  # NEIGHBOURS or GLOBAL
  min_tile_points = list(type = "integer", default = 100000),  # TILES MIN POINTS

  radius   = list(type = "numeric", default = 15),
  nslope   = list(type = "integer", default = 10),
  maxSlope = list(type = "numeric", default = 1.5),
  nAzimuth = list(type = "integer", default = 8),

  densityLod           = list(type = "character", default = "auto"),  # auto = with both ground rasters
  localOffsetWindow    = list(type = "numeric"),
  localOffsetMinPoints = list(type = "integer"),

  outDTM        = list(type = "character", required = TRUE),
  outbins       = list(type = "character", required = TRUE),
  outDif        = list(type = "character", required = TRUE),
  outOutlier    = list(type = "character", required = TRUE),
  outExceed     = list(type = "character", default = "NOFILE"),
  outLod        = list(type = "character", default = "NOFILE"),
  outCenDif     = list(type = "character", default = "NOFILE"),
  outDenOutlier = list(type = "character", default = "NOFILE"),
  outCorDif     = list(type = "character", default = "NOFILE"),
  outDX         = list(type = "character", default = "NOFILE"),
  outDY         = list(type = "character", default = "NOFILE"),
  outDZ         = list(type = "character", default = "NOFILE"),

  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),
  program_name   = list(type = "character", default = "align"),
  progress_interval_s = list(type = "numeric", default = 30)
)

build_params <- function(raw_params, specs, path) {
  unknown <- setdiff(names(raw_params), names(specs))
  if (length(unknown) > 0) {
    warning("Ignoring unrecognized parameter(s) in ", path, ": ", paste(unknown, collapse = ", "))
  }
  out <- list()
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    if (nm %in% names(raw_params) && nzchar(raw_params[[nm]])) {
      raw_val <- raw_params[[nm]]
      converted <- switch(spec$type,
        character = as.character(raw_val),
        numeric   = suppressWarnings(as.numeric(raw_val)),
        integer   = suppressWarnings(as.integer(raw_val)),
        logical   = suppressWarnings(as.logical(raw_val)))
      if (spec$type != "character" && is.na(converted)) {
        stop("Could not parse parameter '", nm, "' (value \"", raw_val, "\") as ", spec$type,
             " -- check ", path, ".")
      }
      out[[nm]] <- converted
    } else if (!is.null(spec$default)) {
      out[[nm]] <- spec$default
    } else if (isTRUE(spec$required)) {
      stop("Required parameter '", nm, "' is missing from ", path)
    }
  }
  out
}

# An environment, not a list: `$` on a list matches partially, so p$prior would return
# prior_points when prior isn't set.
p <- list2env(build_params(read_config(config_path), param_specs, config_path))

is_nofile <- function(x) is.null(x) || !nzchar(x) || toupper(x) == "NOFILE"

## ---- Tiles, density LoD -----------------------------------------------------------------------

auto_tiles <- tolower(p$tileNx) == "auto" || tolower(p$tileNy) == "auto"
if (auto_tiles) {
  if (tolower(p$tileNx) != tolower(p$tileNy)) {
    stop("tileNx and tileNy must both be \"auto\" or both be whole numbers (got ", p$tileNx,
         ", ", p$tileNy, ").")
  }
  tile_nx <- 0L; tile_ny <- 0L
  if (is.null(p$tile_select)) p$tile_select <- "AUTO"
  if (!startsWith(toupper(p$tile_select), "AUTO")) {
    stop("With tileNx/tileNy: auto, tile_select must be AUTO or AUTO REPORT (got ",
         p$tile_select, ").")
  }
} else {
  tile_nx <- suppressWarnings(as.integer(p$tileNx))
  tile_ny <- suppressWarnings(as.integer(p$tileNy))
  if (is.na(tile_nx) || is.na(tile_ny) || tile_nx < 1 || tile_ny < 1) {
    stop("tileNx and tileNy must both be \"auto\" or both be whole numbers of at least 1 ",
         "(got ", p$tileNx, ", ", p$tileNy, ").")
  }
  if (is.null(p$tile_select)) p$tile_select <- "NONE"
}
p$tile_select <- toupper(p$tile_select)

both_grnd <- !is_nofile(p$refGrnd) && !is_nofile(p$alignGrnd)
density_lod <- switch(tolower(p$densityLod),
  auto = both_grnd, true = TRUE, yes = TRUE, false = FALSE, no = FALSE,
  stop("densityLod must be auto, TRUE or FALSE (got ", p$densityLod, ")."))
lod_keys <- c("outExceed", "outLod", "outCenDif", "outDenOutlier")
if (!density_lod) {
  given <- lod_keys[!vapply(mget(lod_keys, envir = p), is_nofile, logical(1))]
  if (length(given)) {
    warning(paste(given, collapse = "/"), " ignored: they need the density LoD ",
            "(densityLod, with both ground-return rasters).")
    for (nm in given) p[[nm]] <- "NOFILE"
  }
}

## ---- Paths ----------------------------------------------------------------------------------

flt_path <- function(base) paste0(sub("\\.flt$", "", base, ignore.case = TRUE), ".flt")
out_keys <- c("outDTM", "outDif", "outOutlier", "outbins", "outExceed", "outLod", "outCenDif",
              "outDenOutlier", "outCorDif", "outDX", "outDY", "outDZ")
out_bases <- mget(out_keys, envir = p)
out_bases <- out_bases[!vapply(out_bases, is_nofile, logical(1))]
out_bases <- lapply(out_bases, function(b) sub("\\.(flt|csv)$", "", b, ignore.case = TRUE))
dup <- duplicated(tolower(basename(unlist(out_bases))))
if (any(dup)) {
  stop("Align writes every output into one run folder, so their file names must differ: ",
       paste(names(out_bases)[dup], collapse = ", "), " repeat an earlier output's name.")
}

final_outlier <- flt_path(out_bases$outOutlier)
if (file.exists(final_outlier)) {
  message("outOutlier already exists -- Align has already been run, nothing to do: ", final_outlier)
  quit(save = "no", status = 0)
}

executable <- file.path(p$executable_dir, paste0(p$program_name, ".exe"))
if (!file.exists(executable)) stop("Align executable not found: ", executable)
for (nm in c("refDTM", "alignDTM")) {
  if (!file.exists(flt_path(p[[nm]]))) stop(nm, " not found: ", p[[nm]])
}

run_dir     <- file.path(dirname(out_bases$outOutlier), "align_run")
run_scratch <- file.path(p$scratch_dir, "align_run")
run_log     <- file.path(run_dir, "align.log")
input_copy  <- file.path(run_dir, "input_align.txt")
run_out     <- lapply(out_bases, function(b) file.path(run_dir, basename(b)))
run_path    <- function(key) if (is.null(run_out[[key]])) "NOFILE" else run_out[[key]]
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(run_scratch, showWarnings = FALSE, recursive = TRUE)

input_file <- TWutils::align_input(
  p$refDTM, p$alignDTM, p$refGrnd, p$alignGrnd, p$refDSM, p$alignDSM,
  p$iterations, p$k, p$dampener, run_path("outDTM"),
  tile_nx, tile_ny, p$overlap, p$radius, p$nslope, p$maxSlope, p$nAzimuth,
  run_path("outbins"), run_path("outDif"), run_path("outOutlier"), run_scratch,
  convergence = p$convergence, minPoints = p$min_tile_points,
  derivatives = p$derivatives, applySlope = p$applySlope, weights = p$weights,
  groundMaskRadius = p$groundMaskRadius,
  waterGradThreshold = p$waterGradThreshold, waterCurvThreshold = p$waterCurvThreshold,
  tileFit = p$tile_fit, tileSelect = p$tile_select,
  cvBlock = p$cv_block, cvBuffer = p$cv_buffer,
  minTileSize = p$min_tile_size, maxTileSize = p$max_tile_size,
  priorPoints = p$prior_points, priorFraction = p$prior_fraction, prior = p$prior,
  outGradient = TRUE, densityLod = density_lod,
  outDenOutlier = run_path("outDenOutlier"), outLod = run_path("outLod"),
  outCenDif = run_path("outCenDif"), outExceed = run_path("outExceed"),
  outCorDif = run_path("outCorDif"),
  localOffsetWindow = p$localOffsetWindow, localOffsetMinPoints = p$localOffsetMinPoints,
  outDX = run_path("outDX"), outDY = run_path("outDY"), outDZ = run_path("outDZ"))

## ---- Run Align ------------------------------------------------------------------------------

# Messages Align echoes, before it reads any data, for each option asked for; one missing means
# the executable predates the option (Align drops unknown TILES sub-fields without a warning).
expected_echo <- c(
  if (p$convergence > 0) c(convergence = "Stop iterating when the aligned DEM changes"),
  if (p$tile_select %in% c("AUTO", "CV")) c(tile_select = "(SELECT = CV)"),
  if (p$tile_select %in% c("AUTO REPORT", "REPORT")) c(tile_select = "(SELECT = REPORT)"),
  if (toupper(p$tile_fit) == "JOINT") c(tile_fit = "(FIT = JOINT)"),
  if (toupper(p$derivatives) == "CENTRAL") c(derivatives = "x and y derivatives: central differences"),
  c(outGradient = "The reference gradient will be written next to the reference DEM"),
  if (density_lod) c(densityLod = "Level of detection from ground-return density"))

check_echo <- function(lines) {
  lines <- trimws(lines)
  unknown_at <- grep("^WARNING: unrecognized keywords", lines)
  if (length(unknown_at)) {
    following <- lines[seq(unknown_at[1] + 1, min(length(lines), unknown_at[1] + 5))]
    return(paste0("it reports unrecognized keywords (", paste(following, collapse = "; "),
                  " ...)"))
  }
  missing <- names(expected_echo)[!vapply(expected_echo, function(s)
    any(grepl(s, lines, fixed = TRUE)), logical(1))]
  if (length(missing)) return(paste0("it didn't acknowledge ", paste(missing, collapse = ", ")))
  NULL
}

# Echoes the log lines worth following while Align runs.
progress_pattern <- paste0("^(Global shift|Iterating over tiles|Fitting the tile shifts|",
                           "Iteration [0-9]+$|RMSE =|Change in aligned DEM|Ratio to the previous|",
                           ".*converged|Cross-validation of the tile count|Selected:|",
                           "The choice is the finest|Nodes held|Level of detection|",
                           "Ending RMSE|Starting RMSE|ERROR|Error|WARNING|No tiles)")

finished_log <- function() {
  file.exists(run_log) && any(grepl("Ending RMSE =", readLines(run_log, warn = FALSE),
                                    fixed = TRUE))
}
same_input <- file.exists(input_copy) &&
  identical(readLines(input_copy, warn = FALSE), readLines(input_file, warn = FALSE))
expected_flts <- vapply(run_out[setdiff(names(run_out), "outbins")], flt_path, character(1))

elapsed <- NA_real_; status <- 0L
if (same_input && finished_log() && all(file.exists(expected_flts))) {
  message("\n==== Align: a finished run with the same input is in ", run_dir,
          " -- not rerun ====")
} else {
  unlink(list.files(run_dir, full.names = TRUE))
  file.copy(input_file, input_copy, overwrite = TRUE)
  message("\n==== Align: ",
          if (auto_tiles) "tiles chosen by cross-validation (SELECT = AUTO)" else
            paste0(tile_nx, " x ", tile_ny, " tiles"),
          ", up to ", p$iterations, " iterations per loop, convergence ", p$convergence,
          " m ====\n  log: ", run_log)
  start_time <- Sys.time()
  proc <- processx::process$new(executable, input_file, stdout = run_log, stderr = "2>&1",
                                cleanup = TRUE)
  shown <- 0L; echo_checked <- FALSE
  show_new_lines <- function() {
    if (!file.exists(run_log)) return(invisible())
    lines <- readLines(run_log, warn = FALSE)
    # The settings echo ends before Align reads the elevations.
    if (!echo_checked && (any(grepl("Reading elevations", lines, fixed = TRUE)) ||
                          !proc$is_alive())) {
      echo_checked <<- TRUE
      problem <- check_echo(lines)
      if (!is.null(problem)) {
        proc$kill()
        stop("Align (", executable, ") is older than this pipeline: ", problem, ". Rebuild ",
             "Projects\\align from GridUtilities 8a08fa6 or later and rerun. Log: ", run_log)
      }
    }
    if (length(lines) > shown) {
      new <- trimws(lines[(shown + 1):length(lines)])
      new <- new[grepl(progress_pattern, new)]
      if (length(new)) message(paste0("  [", format(Sys.time(), "%H:%M:%S"), "] ", new,
                                      collapse = "\n"))
      shown <<- length(lines)
    }
  }
  proc$wait(timeout = 5000)  # long enough for the settings echo checked above
  while (proc$is_alive()) {
    show_new_lines()
    proc$wait(timeout = p$progress_interval_s * 1000)
  }
  show_new_lines()
  status <- proc$get_exit_status()
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message("Align finished in ", round(elapsed, 1), " min, exit status ", status)
  if (status != 0 || !finished_log()) {
    stop("Align failed (exit status ", status, "); see ", run_log)
  }
  missing <- expected_flts[!file.exists(expected_flts)]
  if (length(missing)) {
    stop("Align finished but didn't write ", paste(basename(missing), collapse = ", "),
         "; see ", run_log)
  }
}

## ---- Read Align's log -------------------------------------------------------------------------
## Align's console output is wrapped at 80 columns, but each output raster's .log copy isn't, and
## the copy written last (the longest) holds the whole run up to the final results.

logs <- list.files(run_dir, pattern = "\\.log$", full.names = TRUE)
logs <- setdiff(logs, run_log)
log_lines <- if (length(logs)) {
  sizes <- vapply(logs, function(f) length(readLines(f, warn = FALSE)), integer(1))
  readLines(logs[which.max(sizes)], warn = FALSE)
} else readLines(run_log, warn = FALSE)

num <- function(x) suppressWarnings(as.numeric(sub("[Dd]", "E", trimws(x))))
num_after <- function(line, sep = "=") num(sub(paste0("^.*", sep), "", line))
nums_in <- function(text) num(regmatches(text, gregexpr("-?[0-9]*\\.?[0-9]+([EeDd][-+]?[0-9]+)?",
                                                        text))[[1]])

parse_align_log <- function(lines) {
  lines <- trimws(lines)
  out <- list()
  phase <- "setup"; iter <- NA_integer_
  iters <- list()
  cur <- function() {
    key <- paste(phase, iter)
    if (is.null(iters[[key]])) iters[[key]] <<- list(phase = phase, iteration = iter)
    key
  }
  set <- function(...) { k <- cur(); v <- list(...); iters[[k]][names(v)] <<- v }
  cv_rows <- list(); node_rows <- list(); in_cv <- FALSE; in_nodes <- FALSE

  for (ln in lines) {
    # "Version: x" on the console, "Version x" in the log copies
    if (is.null(out$version) && grepl("^Version:? ", ln)) out$version <- sub("^Version:? +", "", ln)
    if (is.null(out$compiled) && grepl("^Compiled:? ", ln)) out$compiled <- sub("^Compiled:? +", "", ln)
    if (ln == "Global shift") { phase <- "global"; iter <- NA; next }
    if (startsWith(ln, "Cross-validation of the tile count")) { in_cv <- TRUE; next }
    if (startsWith(ln, "Fitting the tile shifts") || ln == "Iterating over tiles") {
      phase <- "tiles"; iter <- NA; in_cv <- FALSE
    }
    if (startsWith(ln, "Node totals")) { in_nodes <- TRUE; next }

    if (in_cv) {
      m <- regmatches(ln, regexec(paste0(
        "^([0-9]+) x +([0-9]+) +([0-9.]+) +([0-9]+) +([0-9.]+) +(-?[0-9.]+) \\+- *([0-9.]+)",
        "(?: +(-?[0-9.]+) \\+- *([0-9.]+))? +([0-9.]+) +([0-9.]+)(.*)$"), ln, perl = TRUE))[[1]]
      if (length(m)) {
        cv_rows[[length(cv_rows) + 1]] <- data.frame(
          tiles_x = as.integer(m[2]), tiles_y = as.integer(m[3]), tile_size_m = num(m[4]),
          unknowns = as.integer(m[5]), prior_scale = num(m[6]),
          gain_mm2 = num(m[7]), gain_se = num(m[8]),
          gain_no_buffer_mm2 = num(m[9]), gain_no_buffer_se = num(m[10]),
          heldout_rms_mm = num(m[11]), block_mean_spread_mm = num(m[12]),
          best = grepl("<- best", m[13], fixed = TRUE),
          chosen = grepl("<- chosen", m[13], fixed = TRUE))
        next
      }
      if (startsWith(ln, "Folds:")) out$cv_folds <- ln
      if (startsWith(ln, "Training cells per fold:")) out$cv_training <- ln
      if (startsWith(ln, "The choice is the finest grid tried")) out$finest_chosen <- TRUE
      if (startsWith(ln, "Selected: tiles")) {
        v <- nums_in(sub("Selected: tiles", "", ln))
        out$tiles_x <- as.integer(v[1]); out$tiles_y <- as.integer(v[2]); out$prior_scale <- v[3]
      }
      if (startsWith(ln, "No candidate could be fitted")) out$select_failed <- TRUE
    }

    if (in_nodes) {
      if (grepl("^i +j +points", ln)) next
      v <- nums_in(ln)
      if (length(v) == 8 && grepl("^[0-9]", ln)) {
        node_rows[[length(node_rows) + 1]] <- data.frame(
          i = as.integer(v[1]), j = as.integer(v[2]), points_k = v[3], x_mm = v[4],
          y_mm = v[5], z_mm = v[6], theta1_mm = v[7], theta2_mm = v[8])
        next
      }
      in_nodes <- FALSE
    }

    if (startsWith(ln, "Tiles X =")) out$tiles_x <- as.integer(num_after(ln))
    if (startsWith(ln, "Tiles Y =")) out$tiles_y <- as.integer(num_after(ln))
    if (startsWith(ln, "In limits RMSE =") && phase == "setup") out$starting_rmse <- num_after(ln)
    if (grepl("^Iteration [0-9]+$", ln)) { iter <- as.integer(sub("Iteration ", "", ln)); next }
    if (phase %in% c("global", "tiles") && !is.na(iter)) {
      if (startsWith(ln, "zmean =")) set(zmean = num_after(ln))
      if (startsWith(ln, "RMSE =")) set(rmse = num_after(ln))
      if (startsWith(ln, "Points included")) set(points = num_after(ln, "[:=]"))
      if (startsWith(ln, "Change in aligned DEM")) {
        # "... max: all  p50 p99 max;  g >= 0.5  p50 p99 max"
        parts <- strsplit(sub("^.*max: *all", "", ln), "g >= 0.5", fixed = TRUE)[[1]]
        a <- nums_in(parts[1]); s <- if (length(parts) > 1) nums_in(parts[2]) else rep(NA, 3)
        set(change_p50 = a[1], change_p99 = a[2], change_max = a[3],
            steep_change_p50 = s[1], steep_change_p99 = s[2], steep_change_max = s[3])
      }
      if (startsWith(ln, "Ratio to the previous P99")) {
        v <- nums_in(sub("^.*P99 =", "", ln))
        set(p99_ratio = v[1], predicted_remaining = v[length(v)])
      }
      if (phase == "global") {
        if (startsWith(ln, "delta x =")) set(delta_x = num_after(ln))
        if (startsWith(ln, "delta y =")) set(delta_y = num_after(ln))
        if (startsWith(ln, "beta 0 ="))  set(beta0 = num_after(ln))
        if (startsWith(ln, "beta 1 ="))  set(beta1 = num_after(ln))
        if (startsWith(ln, "beta 2 ="))  set(beta2 = num_after(ln))
      } else {
        if (startsWith(ln, "Nodes held by their neighbours")) set(nodes_held = num_after(ln))
        if (startsWith(ln, "Node increments, RMS")) {
          v <- nums_in(sub("^.*=", "", ln))
          set(inc_rms_x_mm = v[1], inc_rms_y_mm = v[2], inc_rms_z_mm = v[3])
        }
        if (startsWith(ln, "Node increments, max")) {
          v <- nums_in(sub("^.*=", "", ln))
          set(inc_max_x_mm = v[1], inc_max_y_mm = v[2], inc_max_z_mm = v[3])
        }
      }
    }
    m <- regmatches(ln, regexec("^(Global shift|Tile shifts) (converged after iteration|not converged, iterations =) *([0-9]+)", ln))[[1]]
    if (length(m)) {
      loop <- if (m[2] == "Global shift") "global" else "tiles"
      out[[paste0(loop, "_converged")]] <- startsWith(m[3], "converged")
      out[[paste0(loop, "_iterations")]] <- as.integer(m[4])
    }
    if (startsWith(ln, "Starting RMSE =")) out$starting_rmse <- num_after(ln)
    if (startsWith(ln, "Ending RMSE ="))   out$ending_rmse <- num_after(ln)
    if (startsWith(ln, "Starting mean =")) out$starting_mean <- num_after(ln)
    if (startsWith(ln, "Ending mean ="))   out$ending_mean <- num_after(ln)
  }

  all_cols <- unique(unlist(lapply(iters, names)))
  if (length(iters)) {
    out$iterations <- do.call(rbind, lapply(iters, function(r) {
      r[setdiff(all_cols, names(r))] <- NA
      as.data.frame(r[all_cols])
    }))
    rownames(out$iterations) <- NULL
  }
  if (length(cv_rows)) out$select <- do.call(rbind, cv_rows)
  if (length(node_rows)) out$nodes <- do.call(rbind, node_rows)
  out
}

res <- parse_align_log(log_lines)

## ---- Move the outputs into place ------------------------------------------------------------
## Each file in the run folder goes to the output whose base name it starts with (longest match
## first, so "elev_dif.flt" goes with outDif "elev_dif" rather than outDTM "elev"); outOutlier's
## .flt moves last.

move_file <- function(from, to) {
  if (file.exists(to)) unlink(to)
  if (!file.rename(from, to)) {
    if (!file.copy(from, to, overwrite = TRUE)) stop("Could not move ", from, " to ", to)
    unlink(from)
  }
}

run_tables <- c("align_run.csv", "align_iterations.csv", "align_select.csv", "align_nodes.csv")
files <- setdiff(list.files(run_dir), c("align.log", "input_align.txt", run_tables))
bases <- vapply(run_out, basename, character(1))
bases <- bases[order(-nchar(bases))]
moves <- list()
for (f in files) {
  hit <- names(bases)[vapply(bases, function(b) startsWith(tolower(f), tolower(b)), logical(1))]
  if (length(hit)) moves[[f]] <- file.path(dirname(out_bases[[hit[1]]]), f)
}
outlier_flt <- basename(flt_path(run_out$outOutlier))
for (f in c(setdiff(names(moves), outlier_flt), intersect(names(moves), outlier_flt))) {
  dir.create(dirname(moves[[f]]), showWarnings = FALSE, recursive = TRUE)
  move_file(file.path(run_dir, f), moves[[f]])
}
message("Moved Align's outputs to their final paths (", length(moves), " files).")

## ---- Tables for the report --------------------------------------------------------------------

iters <- res$iterations
n_phase <- function(ph) if (is.null(iters)) 0L else sum(iters$phase == ph)
final_nodes_held <- if (!is.null(iters) && "nodes_held" %in% names(iters) &&
                        any(iters$phase == "tiles")) {
  tail(iters$nodes_held[iters$phase == "tiles"], 1)
} else NA

run_row <- data.frame(
  version = res$version %||% NA, compiled = res$compiled %||% NA,
  tile_select = p$tile_select, tile_fit = toupper(p$tile_fit),
  tiles_x = res$tiles_x %||% NA, tiles_y = res$tiles_y %||% NA,
  prior_scale = res$prior_scale %||% NA, finest_chosen = isTRUE(res$finest_chosen),
  nodes_held = final_nodes_held,
  global_iterations = res$global_iterations %||% n_phase("global"),
  global_converged = res$global_converged %||% NA,
  tile_iterations = res$tiles_iterations %||% n_phase("tiles"),
  tile_converged = res$tiles_converged %||% NA,
  convergence_m = p$convergence, max_iterations = p$iterations,
  starting_rmse = res$starting_rmse %||% NA, ending_rmse = res$ending_rmse %||% NA,
  starting_mean = res$starting_mean %||% NA, ending_mean = res$ending_mean %||% NA,
  density_lod = density_lod, elapsed_min = elapsed, exit_status = status,
  gradient_raster = flt_path(TWutils::align_gradient_raster(p$refDTM, p$radius)),
  log = run_log)
write.csv(run_row, file.path(run_dir, "align_run.csv"), row.names = FALSE)
if (!is.null(iters)) write.csv(iters, file.path(run_dir, "align_iterations.csv"), row.names = FALSE)
if (!is.null(res[["select"]])) {
  write.csv(res[["select"]], file.path(run_dir, "align_select.csv"), row.names = FALSE)
}

# Node positions: the joint fit's nodes are the tiles' midpoints (getTiles() in align_fit.f90),
# placed as raster points are -- (col, row) at xllcorner + (col-1)*cellsize,
# yllcorner + (nrows-row)*cellsize.
nint <- function(x) as.integer(sign(x) * floor(abs(x) + 0.5))  # Fortran NINT
tile_midpoints <- function(ncols, nrows, nx, ny, overlap) {
  dcol <- nint(ncols / nx); drow <- nint(nrows / ny)
  edge_end_c <- min(dcol + nint(dcol * overlap), ncols)
  edge_end_r <- min(drow + nint(drow * overlap), nrows)
  edge_start_c <- max(ncols - (dcol + nint(dcol * overlap)), 1L)
  edge_start_r <- max(nrows - (drow + nint(drow * overlap)), 1L)
  half_c <- nint(dcol * overlap * 0.5); half_r <- nint(drow * overlap * 0.5)
  span <- function(k, n, d, half, end1, startn, len) {
    if (k == 1) c(1L, end1) else if (k == n) c(startn, len)
    else c(max(d * (k - 1) - half, 1L), min(d * k + half, len))
  }
  grid <- expand.grid(i = seq_len(nx), j = seq_len(ny))
  grid$midcol <- vapply(grid$i, function(i) {
    s <- span(i, nx, dcol, half_c, edge_end_c, edge_start_c, ncols); nint(sum(s) * 0.5)
  }, integer(1))
  grid$midrow <- vapply(grid$j, function(j) {
    s <- span(j, ny, drow, half_r, edge_end_r, edge_start_r, nrows); nint(sum(s) * 0.5)
  }, integer(1))
  grid
}

if (!is.null(res$nodes) && !is.null(res$tiles_x)) {
  ref <- terra::rast(flt_path(p$refDTM))
  cs <- terra::res(ref)[1]
  e <- terra::ext(ref)
  overlap_used <- p$overlap
  mid <- tile_midpoints(terra::ncol(ref), terra::nrow(ref), res$tiles_x, res$tiles_y, overlap_used)
  mid$x <- e$xmin + (mid$midcol - 1) * cs
  mid$y <- e$ymin + (terra::nrow(ref) - mid$midrow) * cs
  nodes <- merge(res$nodes, mid[c("i", "j", "x", "y")], by = c("i", "j"), all.x = TRUE)
  nodes <- nodes[order(nodes$j, nodes$i), ]
  write.csv(nodes, file.path(run_dir, "align_nodes.csv"), row.names = FALSE)
}

message("\nAlign: ", res$tiles_x %||% "?", " x ", res$tiles_y %||% "?", " tiles",
        if (!is.null(res$prior_scale)) paste0(", prior x ", res$prior_scale),
        if (isTRUE(res$finest_chosen)) " (the finest grid tried -- consider a smaller min_tile_size)",
        "; ending RMSE ", signif(res$ending_rmse %||% NA, 4), " m.")
message("Run tables: ", run_dir)
