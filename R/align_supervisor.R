## align_supervisor.R
##
## Runs the Fortran program Align (TWutils::align_input() writes its input file) up to
## max_attempts times on one reference/align DTM pair, choosing the number of tiles for each
## attempt, and keeps the attempt with the best RMSE. Run on its own (Rscript align_supervisor.R
## params.txt) or from run_align_pipeline.R, whose Stage 3 runs it when [align] sets
## tileNx: auto.
##
## Iterations. Each attempt runs at most `iterations` iterations, and Align itself stops early
## once an iteration changes the RMSE by less than `convergence` meters (Align's CONVERGENCE
## keyword -- alignSlope2 1.0.6 or later; an older executable ignores the keyword, so this
## script stops if the log shows it wasn't applied). The global loop and the tiled loop are
## tested separately: the global loop starts from the unaligned RMSE, the tiled loop from the
## last global RMSE.
##
## Tiles. Align lays its tiles over the reference DTM's whole grid (getTiles() in
## alignSlope2.f90), not the study-area polygon, and solves a tile only if it holds at least
## MIN POINTS regression points (100000 by default) and a quarter of its cells with data; the
## rest are filled from their neighbors. So the useful number of tiles depends on the polygon's
## shape, not just its area. This script reproduces getTiles() and counts the cells with data in
## each tile (reference DTM, align DTM if it shares the reference grid, and both ground-return
## rasters > 0) to predict which tiles Align will solve:
##   attempt 1 -- the tileNx x tileNy layout (tile aspect ratio following the DTM's) whose
##                predicted number of solved tiles is closest to (area with data)/tile_size_m^2;
##                or tileNx/tileNy themselves, if given as numbers;
##   attempt n -- the layout whose predicted number of solved tiles is closest to tile_growth
##                times the previous attempt's.
## An attempt is kept only if its ending RMSE is at least min_rmse_reduction (a fraction) below
## the best attempt so far; once one isn't, or no finer layout adds solved tiles, or
## max_attempts is reached, the search stops and the best attempt is used. More tiles let the
## shift -- and so the LoD -- vary across the study area, but each tile has fewer points, so a
## finer layout must earn its place with a lower RMSE.
##
## Outputs. Each attempt writes to <dirname(outOutlier)>/align_attempts/attempt_<n>/ (Align's
## scratch files to <scratch_dir>/align_attempt_<n>/), with Align's console output in
## align.log there, streamed to this session's console as it runs. The chosen attempt's outputs
## are then moved to outDTM/outDif/outOutlier/outbins (outOutlier last, since its existence
## marks the run complete), and with keep_attempt_rasters FALSE the other attempts' rasters are
## deleted. Three tables go to the align_attempts folder, for the pipeline report:
##   align_attempts.csv   -- one row per attempt: tiles, iterations used, convergence,
##                           starting/global/ending RMSE, run time, and the decision taken;
##   align_iterations.csv -- RMSE and mean difference for every iteration of every attempt;
##   align_tiles.csv      -- each attempt's tiles: extent, center, and last-iteration shifts,
##                           points, and standard errors (or insufficient = TRUE).
## Resumability: an attempt whose outputs and log are complete, with the same tiles, is read
## back instead of rerun; if outOutlier already exists, nothing runs.
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

param_specs <- list(
  refDTM    = list(type = "character", required = TRUE),
  alignDTM  = list(type = "character", required = TRUE),
  refGrnd   = list(type = "character", default = "NOFILE"),
  alignGrnd = list(type = "character", default = "NOFILE"),
  refDSM    = list(type = "character", default = "NOFILE"),
  alignDSM  = list(type = "character", default = "NOFILE"),

  iterations  = list(type = "integer", default = 10),     # maximum iterations per loop
  convergence = list(type = "numeric", default = 0.001),  # m; stop once RMSE changes by less
  k           = list(type = "numeric", default = 1.5),
  dampener    = list(type = "numeric", default = 1),

  tileNx  = list(type = "character", default = "auto"),  # "auto" or a number of tiles
  tileNy  = list(type = "character", default = "auto"),
  overlap = list(type = "numeric", default = 0.5),

  radius   = list(type = "numeric", default = 15),
  nslope   = list(type = "integer", default = 7),
  maxSlope = list(type = "numeric", default = 1),
  nAzimuth = list(type = "integer", default = 8),

  outDTM     = list(type = "character", required = TRUE),
  outbins    = list(type = "character", required = TRUE),
  outDif     = list(type = "character", required = TRUE),
  outOutlier = list(type = "character", required = TRUE),

  scratch_dir    = list(type = "character", required = TRUE),
  executable_dir = list(type = "character", required = TRUE),
  program_name   = list(type = "character", default = "align"),

  # Tile search
  tile_size_m          = list(type = "numeric", default = 2000),  # attempt 1's target tile spacing
  max_attempts         = list(type = "integer", default = 3),
  tile_growth          = list(type = "numeric", default = 2),     # solved-tile multiplier per attempt
  min_rmse_reduction   = list(type = "numeric", default = 0.02),  # fraction of the best RMSE
  min_tile_points      = list(type = "integer", default = 100000),  # Align's TILES MIN POINTS default
  keep_attempt_rasters = list(type = "logical", default = FALSE),
  progress_interval_s  = list(type = "numeric", default = 30)
)

build_params <- function(raw_params, specs, path) {
  unknown <- setdiff(names(raw_params), names(specs))
  if (length(unknown) > 0) {
    warning("Ignoring unrecognized parameter(s) in ", path, ": ", paste(unknown, collapse = ", "))
  }
  out <- list()
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    if (nm %in% names(raw_params)) {
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

p <- build_params(read_config(config_path), param_specs, config_path)

auto_tiles <- tolower(p$tileNx) == "auto" || tolower(p$tileNy) == "auto"
if (!auto_tiles) {
  fixed_nx <- suppressWarnings(as.integer(p$tileNx))
  fixed_ny <- suppressWarnings(as.integer(p$tileNy))
  if (is.na(fixed_nx) || is.na(fixed_ny) || fixed_nx < 1 || fixed_ny < 1) {
    stop("tileNx and tileNy must both be \"auto\" or both be whole numbers of at least 1 ",
         "(got ", p$tileNx, ", ", p$tileNy, ").")
  }
}
if (p$max_attempts < 1) stop("max_attempts must be at least 1.")
if (p$tile_growth <= 1) stop("tile_growth must be greater than 1.")

flt_path  <- function(base) paste0(sub("\\.flt$", "", base, ignore.case = TRUE), ".flt")
out_bases <- list(outDTM = p$outDTM, outDif = p$outDif, outOutlier = p$outOutlier,
                  outbins = p$outbins)
out_bases <- lapply(out_bases, function(b) sub("\\.(flt|csv)$", "", b, ignore.case = TRUE))

final_outlier <- flt_path(out_bases$outOutlier)
if (file.exists(final_outlier)) {
  message("outOutlier already exists -- Align has already been run, nothing to do: ", final_outlier)
  quit(save = "no", status = 0)
}

executable <- file.path(p$executable_dir, paste0(p$program_name, ".exe"))
if (!file.exists(executable)) stop("Align executable not found: ", executable)
for (nm in c("refDTM", "alignDTM")) {
  if (!file.exists(p[[nm]])) stop(nm, " not found: ", p[[nm]])
}

attempts_dir <- file.path(dirname(out_bases$outOutlier), "align_attempts")
dir.create(attempts_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Tile layout, as Align's getTiles() lays it out ------------------------------------------

# Fortran NINT: rounds halves away from zero (R's round() rounds them to even).
nint <- function(x) as.integer(sign(x) * floor(abs(x) + 0.5))

# Same assignments, in the same order, as getTiles() -- with one tile in x or y, later
# assignments overwrite earlier ones there too.
tile_layout <- function(ncols, nrows, nx, ny, overlap) {
  sc <- ec <- mc <- sr <- er <- mr <- matrix(NA_integer_, nx, ny)
  dcol <- nint(ncols / nx)
  drow <- nint(nrows / ny)
  colov <- nint(dcol * overlap)
  rowov <- nint(drow * overlap)

  sc[1, 1] <- 1L; sr[1, 1] <- 1L
  ec[1, 1] <- min(dcol + colov, ncols); er[1, 1] <- min(drow + rowov, nrows)
  mc[1, 1] <- nint((1 + ec[1, 1]) * 0.5); mr[1, 1] <- nint((1 + er[1, 1]) * 0.5)

  sc[1, ny] <- 1L; ec[1, ny] <- ec[1, 1]
  sr[1, ny] <- max(nrows - (drow + rowov), 1L); er[1, ny] <- nrows
  mc[1, ny] <- mc[1, 1]; mr[1, ny] <- nint((sr[1, ny] + er[1, ny]) * 0.5)

  sc[nx, 1] <- max(ncols - (dcol + colov), 1L); ec[nx, 1] <- ncols
  sr[nx, 1] <- 1L; er[nx, 1] <- er[1, 1]
  mc[nx, 1] <- nint((sc[nx, 1] + ec[nx, 1]) * 0.5); mr[nx, 1] <- nint((1 + er[nx, 1]) * 0.5)

  sc[nx, ny] <- sc[nx, 1]; ec[nx, ny] <- ncols
  sr[nx, ny] <- sr[1, ny]; er[nx, ny] <- nrows
  mc[nx, ny] <- nint((sc[nx, ny] + ec[nx, ny]) * 0.5)
  mr[nx, ny] <- nint((sr[nx, ny] + er[nx, ny]) * 0.5)

  colov <- nint(dcol * overlap * 0.5)
  rowov <- nint(drow * overlap * 0.5)
  inner_cols <- function(i) c(max(dcol * (i - 1) - colov, 1L), min(dcol * i + colov, ncols))
  inner_rows <- function(j) c(max(drow * (j - 1) - rowov, 1L), min(drow * j + rowov, nrows))

  for (i in seq_len(max(nx - 2, 0)) + 1) {
    for (j in c(1, ny)) {
      sc[i, j] <- inner_cols(i)[1]; ec[i, j] <- inner_cols(i)[2]
      sr[i, j] <- if (j == 1) 1L else sr[1, ny]
      er[i, j] <- if (j == 1) er[1, 1] else nrows
      mc[i, j] <- nint((sc[i, j] + ec[i, j]) * 0.5); mr[i, j] <- nint((sr[i, j] + er[i, j]) * 0.5)
    }
  }
  for (j in seq_len(max(ny - 2, 0)) + 1) {
    for (i in c(1, nx)) {
      sc[i, j] <- if (i == 1) 1L else sc[nx, 1]
      ec[i, j] <- if (i == 1) ec[1, 1] else ec[nx, 1]
      sr[i, j] <- inner_rows(j)[1]; er[i, j] <- inner_rows(j)[2]
      mc[i, j] <- nint((sc[i, j] + ec[i, j]) * 0.5); mr[i, j] <- nint((sr[i, j] + er[i, j]) * 0.5)
    }
  }
  for (j in seq_len(max(ny - 2, 0)) + 1) {
    for (i in seq_len(max(nx - 2, 0)) + 1) {
      sc[i, j] <- inner_cols(i)[1]; ec[i, j] <- inner_cols(i)[2]
      sr[i, j] <- inner_rows(j)[1]; er[i, j] <- inner_rows(j)[2]
      mc[i, j] <- nint((sc[i, j] + ec[i, j]) * 0.5); mr[i, j] <- nint((sr[i, j] + er[i, j]) * 0.5)
    }
  }

  data.frame(i = as.vector(row(sc)), j = as.vector(col(sc)),
             startcol = as.vector(sc), endcol = as.vector(ec), midcol = as.vector(mc),
             startrow = as.vector(sr), endrow = as.vector(er), midrow = as.vector(mr))
}

## ---- Cells with data, in blocks ---------------------------------------------------------------
## Cells Align can use: reference DTM, align DTM (when it shares the reference grid -- Align
## resamples it otherwise), and both ground-return rasters > 0. Summed into blocks of
## block x block cells so a tile's count is a sum over a small matrix.

ref <- terra::rast(flt_path(p$refDTM))
dtm_ncols <- terra::ncol(ref)
dtm_nrows <- terra::nrow(ref)
cellsize  <- terra::res(ref)[1]
dtm_ext   <- terra::ext(ref)

valid <- !is.na(ref)
aln <- terra::rast(flt_path(p$alignDTM))
if (terra::compareGeom(ref, aln, stopOnError = FALSE)) {
  valid <- valid & !is.na(aln)
} else {
  message("alignDTM is not on the reference DTM's grid -- predicting tiles from the reference ",
          "DTM alone.")
}
for (nm in c("refGrnd", "alignGrnd")) {
  if (toupper(p[[nm]]) != "NOFILE") {
    g <- terra::rast(flt_path(p[[nm]]))
    if (!terra::compareGeom(ref, g, stopOnError = FALSE)) {
      stop(nm, " is not on the reference DTM's grid: ", p[[nm]])
    }
    valid <- valid & !is.na(g) & g > 0
  }
}

block <- max(1L, ceiling(max(dtm_ncols, dtm_nrows) / 1500))
counts <- if (block > 1) terra::aggregate(valid, fact = block, fun = "sum", na.rm = TRUE) else valid
counts <- terra::as.matrix(counts, wide = TRUE)  # rows north to south, as Align's rows
counts[is.na(counts)] <- 0
data_cells <- sum(counts)
if (data_cells == 0) stop("The reference DTM has no cells with data (after the ground-return masks).")

# Cells with data in each tile (to block resolution), and whether Align should solve it:
# Align needs MIN POINTS regression points and a quarter of the cells with data; outlier and
# edge cells drop out of the regression, so include_fraction of the cells with data is assumed
# to make it in.
include_fraction <- 0.85
predict_tiles <- function(nx, ny) {
  tl <- tile_layout(dtm_ncols, dtm_nrows, nx, ny, p$overlap)
  bc1 <- (tl$startcol - 1) %/% block + 1; bc2 <- pmin((tl$endcol - 1) %/% block + 1, ncol(counts))
  br1 <- (tl$startrow - 1) %/% block + 1; br2 <- pmin((tl$endrow - 1) %/% block + 1, nrow(counts))
  tl$data_cells <- mapply(function(a, b, c, d) sum(counts[c:d, a:b]), bc1, bc2, br1, br2)
  tl$predicted_ok <- tl$data_cells * include_fraction >= p$min_tile_points
  tl
}

# Candidate layouts keep the tile aspect ratio near the DTM's (the short side's count rounded
# both down and up). Picks the one predicting more than more_than solved tiles whose count is
# closest to target (fewer tiles on a tie); NULL if none predicts more than more_than.
width_m  <- dtm_ncols * cellsize
height_m <- dtm_nrows * cellsize
choose_layout <- function(target, more_than = -Inf) {
  max_n <- max(1, ceiling(max(width_m, height_m) / max(cellsize * 100, 50)))
  ratio <- min(width_m, height_m) / max(width_m, height_m)
  cand <- do.call(rbind, lapply(seq_len(min(max_n, 60)), function(n) {
    short <- unique(pmax(1, c(floor(n * ratio), ceiling(n * ratio))))
    if (width_m >= height_m) cbind(n, short) else cbind(short, n)
  }))
  cand <- unique(cand)
  scored <- do.call(rbind, lapply(seq_len(nrow(cand)), function(r) {
    tl <- predict_tiles(cand[r, 1], cand[r, 2])
    data.frame(nx = cand[r, 1], ny = cand[r, 2], solved = sum(tl$predicted_ok))
  }))
  scored <- scored[scored$solved > more_than, ]
  if (nrow(scored) == 0) return(NULL)
  scored <- scored[order(abs(scored$solved - target), scored$nx * scored$ny), ]
  scored[1, ]
}

## ---- Reading Align's log -------------------------------------------------------------------

num_after <- function(line, label) {
  as.numeric(sub("[Dd]", "E", trimws(sub(paste0(".*", label), "", line))))
}

parse_align_log <- function(lines) {
  lines <- trimws(lines)
  out <- list(iterations = NULL, tiles = NULL, insufficient = NULL)
  phase <- "setup"; iter <- NA; zmean <- NA; ti <- NA; tj <- NA; rec <- NULL
  iter_rows <- list(); tile_rows <- list(); insuff_rows <- list(); global_rows <- list()
  flush_tile <- function() {
    if (!is.null(rec)) tile_rows[[length(tile_rows) + 1]] <<- as.data.frame(rec)
    rec <<- NULL
  }
  pending_insufficient <- FALSE
  for (ln in lines) {
    if (ln == "Global shift") { phase <- "global"; next }
    if (ln == "Iterating over tiles") { flush_tile(); phase <- "tiles"; next }
    if (startsWith(ln, "Version:")) out$version <- trimws(sub("Version:", "", ln))
    if (startsWith(ln, "Stop iterating when RMSE changes")) out$convergence_applied <- TRUE
    if (startsWith(ln, "In limits RMSE =")) out$starting_rmse <- num_after(ln, "=")
    if (grepl("^Iteration [0-9]+$", ln)) { iter <- as.integer(sub("Iteration ", "", ln)); next }
    if (startsWith(ln, "zmean =")) { zmean <- num_after(ln, "="); next }
    if (startsWith(ln, "RMSE =")) {
      flush_tile()
      iter_rows[[length(iter_rows) + 1]] <- data.frame(phase = phase, iteration = iter,
                                                       zmean = zmean, rmse = num_after(ln, "="))
      next
    }
    if (grepl("^Global shift converged after iteration", ln)) out$global_converged <- TRUE
    if (grepl("^Global shift not converged", ln)) out$global_converged <- FALSE
    if (grepl("^Tile shifts converged after iteration", ln)) out$tiles_converged <- TRUE
    if (grepl("^Tile shifts not converged", ln)) out$tiles_converged <- FALSE
    if (phase == "global") {
      if (startsWith(ln, "delta x =")) global_rows[[length(global_rows) + 1]] <- list(dx = num_after(ln, "="))
      if (startsWith(ln, "delta y =")) global_rows[[length(global_rows)]]$dy <- num_after(ln, "=")
      if (startsWith(ln, "beta 0 ="))  global_rows[[length(global_rows)]]$beta0 <- num_after(ln, "=")
    }
    if (phase == "tiles") {
      if (startsWith(ln, "X tile =")) {
        flush_tile(); rec <- list(iteration = iter, i = as.integer(num_after(ln, "=")))
        next
      }
      if (!is.null(rec)) {
        if (startsWith(ln, "Y tile ="))         rec$j <- as.integer(num_after(ln, "="))
        if (startsWith(ln, "Points included:")) rec$points <- num_after(ln, ":")
        if (startsWith(ln, "delta x ="))        rec$dx <- num_after(ln, "=")
        if (startsWith(ln, "delta y ="))        rec$dy <- num_after(ln, "=")
        if (startsWith(ln, "beta 0 ="))         rec$beta0 <- num_after(ln, "=")
        if (startsWith(ln, "SE delta x (m) =")) rec$se_dx <- num_after(ln, "=")
        if (startsWith(ln, "SE delta y (m) =")) rec$se_dy <- num_after(ln, "=")
        if (startsWith(ln, "SE beta 0 (m) ="))  rec$se_beta0 <- num_after(ln, "=")
      }
      if (ln == "Insufficient points for tile") { flush_tile(); pending_insufficient <- TRUE; ti <- NA; next }
      if (pending_insufficient && startsWith(ln, "Tile x =")) { ti <- as.integer(num_after(ln, "=")); next }
      if (pending_insufficient && startsWith(ln, "Tile y =")) {
        tj <- as.integer(num_after(ln, "="))
        insuff_rows[[length(insuff_rows) + 1]] <- data.frame(iteration = iter, i = ti, j = tj)
        pending_insufficient <- FALSE
        next
      }
    }
    if (startsWith(ln, "Starting RMSE ="))  out$starting_rmse <- num_after(ln, "=")
    if (startsWith(ln, "Ending RMSE ="))    out$ending_rmse <- num_after(ln, "=")
    if (startsWith(ln, "Ending mean ="))    out$ending_mean <- num_after(ln, "=")
    if (startsWith(ln, "Total mean ="))     out$total_mean <- num_after(ln, "=")
    if (grepl("No tiles had enough points", ln, fixed = TRUE)) out$no_tiles <- TRUE
  }
  flush_tile()
  if (length(iter_rows))   out$iterations   <- do.call(rbind, iter_rows)
  if (length(tile_rows))   out$tiles        <- do.call(rbind, lapply(tile_rows, function(r) {
    for (nm in c("j", "points", "dx", "dy", "beta0", "se_dx", "se_dy", "se_beta0")) {
      if (is.null(r[[nm]])) r[[nm]] <- NA
    }
    r[c("iteration", "i", "j", "points", "dx", "dy", "beta0", "se_dx", "se_dy", "se_beta0")]
  }))
  if (length(insuff_rows)) out$insufficient <- do.call(rbind, insuff_rows)
  # Each iteration's shift is applied on top of the last, so the global shift is their sum.
  out$global_shift <- c(dx = sum(vapply(global_rows, function(r) r$dx %||% 0, 0)),
                        dy = sum(vapply(global_rows, function(r) r$dy %||% 0, 0)),
                        beta0 = sum(vapply(global_rows, function(r) r$beta0 %||% 0, 0)))
  out
}

## ---- One Align attempt ----------------------------------------------------------------------

attempt_paths <- function(n) {
  dir <- file.path(attempts_dir, sprintf("attempt_%d", n))
  list(dir = dir,
       scratch = file.path(p$scratch_dir, sprintf("align_attempt_%d", n)),
       log = file.path(dir, "align.log"),
       out = lapply(out_bases, function(b) file.path(dir, basename(b))))
}

# Echoes the log lines worth following while Align runs: phases, RMSE, convergence, problems.
progress_pattern <- paste0("^(Global shift$|Iterating over tiles|RMSE =|In limits RMSE|",
                           "Iteration [0-9]+$|.*converged|Insufficient points|Ending RMSE|",
                           "Starting RMSE|ERROR|Error|No tiles)")

run_attempt <- function(n, nx, ny) {
  ap <- attempt_paths(n)
  dir.create(ap$dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(ap$scratch, showWarnings = FALSE, recursive = TRUE)

  # A finished attempt with these tiles is read back rather than rerun.
  input_copy <- file.path(ap$dir, "input_align.txt")
  if (file.exists(flt_path(ap$out$outOutlier)) && file.exists(ap$log) && file.exists(input_copy)) {
    tiles_line <- grep("^TILES:", readLines(input_copy, warn = FALSE), value = TRUE)
    log_lines <- readLines(ap$log, warn = FALSE)
    if (length(tiles_line) == 1 &&
        startsWith(tiles_line, sprintf("TILES: X = %d, Y = %d, OVERLAP = %s", nx, ny,
                                       format(p$overlap))) &&
        any(grepl("Ending RMSE =", log_lines, fixed = TRUE))) {
      message("\n==== Align attempt ", n, " (", nx, " x ", ny, " tiles): already complete, ",
              "reading ", ap$log, " ====")
      return(list(log = log_lines, elapsed_min = NA_real_, status = 0L))
    }
  }
  unlink(list.files(ap$dir, full.names = TRUE))

  input_file <- TWutils::align_input(p$refDTM, p$alignDTM, p$refGrnd, p$alignGrnd,
                                     p$refDSM, p$alignDSM,
                                     p$iterations, p$k, p$dampener, ap$out$outDTM,
                                     nx, ny, p$overlap, p$radius,
                                     p$nslope, p$maxSlope, p$nAzimuth,
                                     ap$out$outbins, ap$out$outDif, ap$out$outOutlier,
                                     ap$scratch, convergence = p$convergence,
                                     minPoints = p$min_tile_points)

  message("\n==== Align attempt ", n, ": ", nx, " x ", ny, " tiles, up to ", p$iterations,
          " iterations, convergence ", p$convergence, " m ====\n  log: ", ap$log)
  start_time <- Sys.time()
  proc <- processx::process$new(executable, input_file, stdout = ap$log, stderr = "2>&1",
                                cleanup = TRUE)
  shown <- 0L
  show_new_lines <- function() {
    if (!file.exists(ap$log)) return(invisible())
    lines <- readLines(ap$log, warn = FALSE)
    # Align echoes its settings before it starts work; an executable older than alignSlope2
    # 1.0.6 silently ignores CONVERGENCE, so stop it now rather than after a full run.
    if (p$convergence > 0 && any(startsWith(trimws(lines), "Number of iterations =")) &&
        !any(startsWith(trimws(lines), "Stop iterating when RMSE changes"))) {
      proc$kill()
      stop("Align (", executable, ") ignored CONVERGENCE -- it predates alignSlope2 1.0.6. ",
           "Rebuild Projects\\align and rerun.")
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
  file.copy(input_file, input_copy, overwrite = TRUE)
  message("Align attempt ", n, " finished in ", round(elapsed, 1), " min, exit status ", status)
  list(log = readLines(ap$log, warn = FALSE), elapsed_min = elapsed, status = status)
}

## ---- The search ---------------------------------------------------------------------------

attempt_rows <- list(); iteration_rows <- list(); tile_rows <- list()
best <- NA_integer_
# target: solved tiles wanted from the next attempt; must_exceed: a finer layout has to
# predict more solved tiles than this (the last successful attempt's), or the search stops.
target <- if (auto_tiles) max(1, round(data_cells * cellsize^2 / p$tile_size_m^2)) else NA
must_exceed <- -Inf
stop_reason <- "max_attempts reached"

message("Reference DTM: ", dtm_ncols, " x ", dtm_nrows, " cells of ", signif(cellsize, 6),
        " m; ", signif(data_cells * cellsize^2 / 1e6, 4), " km^2 with data.")

for (n in seq_len(p$max_attempts)) {

  # Choose this attempt's tiles.
  if (n == 1 && !auto_tiles) {
    layout <- data.frame(nx = fixed_nx, ny = fixed_ny,
                         solved = sum(predict_tiles(fixed_nx, fixed_ny)$predicted_ok))
  } else {
    layout <- choose_layout(target, more_than = must_exceed)
    if (is.null(layout)) {
      stop_reason <- paste0("no finer layout is predicted to solve more than ", must_exceed,
                            " tiles")
      break
    }
  }
  tl <- predict_tiles(layout$nx, layout$ny)
  message("\nAttempt ", n, ": ", layout$nx, " x ", layout$ny, " tiles, ", layout$solved,
          " predicted to have enough points",
          if (!is.na(target)) paste0(" (target ", round(target), ")"))

  run <- run_attempt(n, layout$nx, layout$ny)
  res <- parse_align_log(run$log)

  # Per-tile table: geometry from the layout; points and standard errors from the last tiled
  # iteration; dx/dy/beta0 (the last iteration's correction, which should be small once
  # converged); and the total shift at the tile center -- the global shift plus every tiled
  # iteration's, NA if the tile was ever filled from its neighbors rather than solved.
  last_iter <- if (!is.null(res$iterations) && any(res$iterations$phase == "tiles")) {
    max(res$iterations$iteration[res$iterations$phase == "tiles"])
  } else NA
  if (!is.null(res$tiles)) {
    n_tile_iters <- sum(res$iterations$phase == "tiles")
    totals <- aggregate(cbind(dx, dy, beta0) ~ i + j, data = res$tiles, FUN = sum)
    solved_iters <- aggregate(iteration ~ i + j, data = res$tiles, FUN = length)
    totals <- merge(totals, solved_iters, by = c("i", "j"))
    always <- totals$iteration == n_tile_iters
    totals$dx_total    <- ifelse(always, res$global_shift[["dx"]] + totals$dx, NA)
    totals$dy_total    <- ifelse(always, res$global_shift[["dy"]] + totals$dy, NA)
    totals$beta0_total <- ifelse(always, res$global_shift[["beta0"]] + totals$beta0, NA)
    totals <- totals[c("i", "j", "dx_total", "dy_total", "beta0_total")]
  }
  tiles_out <- tl
  tiles_out$x_min <- dtm_ext$xmin + (tl$startcol - 1) * cellsize
  tiles_out$x_max <- dtm_ext$xmin + tl$endcol * cellsize
  tiles_out$y_max <- dtm_ext$ymax - (tl$startrow - 1) * cellsize
  tiles_out$y_min <- dtm_ext$ymax - tl$endrow * cellsize
  tiles_out$x_center <- dtm_ext$xmin + (tl$midcol - 0.5) * cellsize
  tiles_out$y_center <- dtm_ext$ymax - (tl$midrow - 0.5) * cellsize
  if (!is.null(res$tiles) && !is.na(last_iter)) {
    solved <- res$tiles[res$tiles$iteration == last_iter, c("i", "j", "points", "dx", "dy",
                                                            "beta0", "se_dx", "se_dy",
                                                            "se_beta0")]
    tiles_out <- merge(tiles_out, solved, by = c("i", "j"), all.x = TRUE)
    tiles_out <- merge(tiles_out, totals, by = c("i", "j"), all.x = TRUE)
  }
  tiles_out$global_dx <- res$global_shift[["dx"]]
  tiles_out$global_dy <- res$global_shift[["dy"]]
  tiles_out$global_beta0 <- res$global_shift[["beta0"]]
  insufficient_ij <- if (!is.null(res$insufficient) && !is.na(last_iter)) {
    res$insufficient[res$insufficient$iteration == last_iter, c("i", "j")]
  } else data.frame(i = integer(0), j = integer(0))
  tiles_out$insufficient <- paste(tiles_out$i, tiles_out$j) %in%
    paste(insufficient_ij$i, insufficient_ij$j)
  tiles_out$attempt <- n
  tile_rows[[n]] <- tiles_out

  if (!is.null(res$iterations)) iteration_rows[[n]] <- cbind(attempt = n, res$iterations)

  global_rmse <- if (!is.null(res$iterations) && any(res$iterations$phase == "global")) {
    tail(res$iterations$rmse[res$iterations$phase == "global"], 1)
  } else NA_real_
  ok <- run$status == 0 && !is.null(res$ending_rmse) && is.finite(res$ending_rmse)

  row <- data.frame(
    attempt = n, tileNx = layout$nx, tileNy = layout$ny, tiles = layout$nx * layout$ny,
    tiles_predicted_ok = layout$solved, tiles_insufficient = nrow(insufficient_ij),
    global_iterations = sum(res$iterations$phase == "global"),
    global_converged = res$global_converged %||% NA,
    tile_iterations = sum(res$iterations$phase == "tiles"),
    tile_converged = res$tiles_converged %||% NA,
    starting_rmse = res$starting_rmse %||% NA, global_rmse = global_rmse,
    ending_rmse = if (ok) res$ending_rmse else NA, ending_mean = res$ending_mean %||% NA,
    rmse_reduction = NA_real_, elapsed_min = run$elapsed_min, exit_status = run$status,
    selected = FALSE, decision = "", log = attempt_paths(n)$log)

  if (!ok) {
    row$decision <- if (isTRUE(res$no_tiles)) "failed: no tile had enough points" else "failed"
    attempt_rows[[n]] <- row
    # Too many tiles for the data and nothing better yet: try a coarser layout.
    if (is.na(best) && isTRUE(res$no_tiles) && n < p$max_attempts && layout$nx * layout$ny > 1) {
      target <- max(1, floor(layout$solved / p$tile_growth))
      next
    }
    stop_reason <- paste0("attempt ", n, " failed")
    break
  }

  if (is.na(best)) {
    row$decision <- "first successful attempt"
    best <- n
  } else {
    best_rmse <- attempt_rows[[best]]$ending_rmse
    row$rmse_reduction <- (best_rmse - row$ending_rmse) / best_rmse
    if (row$rmse_reduction >= p$min_rmse_reduction) {
      row$decision <- sprintf("kept: RMSE %.1f%% below attempt %d", 100 * row$rmse_reduction, best)
      best <- n
    } else {
      row$decision <- sprintf("not kept: RMSE %.1f%% below attempt %d, less than %.1f%%",
                              100 * row$rmse_reduction, best, 100 * p$min_rmse_reduction)
      attempt_rows[[n]] <- row
      stop_reason <- "a finer layout no longer reduced RMSE enough"
      break
    }
  }
  attempt_rows[[n]] <- row
  must_exceed <- layout$solved
  target <- max(layout$solved, 1) * p$tile_growth
}

attempts <- do.call(rbind, attempt_rows)
if (is.na(best)) {
  write.csv(attempts, file.path(attempts_dir, "align_attempts.csv"), row.names = FALSE)
  stop("No Align attempt succeeded (", stop_reason, "); see the logs in ", attempts_dir)
}
attempts$selected <- attempts$attempt == best
message("\nTile search stopped: ", stop_reason, ". Using attempt ", best, " (",
        attempts$tileNx[best], " x ", attempts$tileNy[best], " tiles, ending RMSE ",
        signif(attempts$ending_rmse[best], 4), " m).")

write.csv(attempts, file.path(attempts_dir, "align_attempts.csv"), row.names = FALSE)
write.csv(do.call(rbind, iteration_rows), file.path(attempts_dir, "align_iterations.csv"),
          row.names = FALSE)
all_cols <- unique(unlist(lapply(tile_rows, names)))
write.csv(do.call(rbind, lapply(tile_rows, function(t) {
  t[setdiff(all_cols, names(t))] <- NA
  t[all_cols]
})), file.path(attempts_dir, "align_tiles.csv"), row.names = FALSE)

## ---- Move the chosen attempt's outputs into place -------------------------------------------
## Each file in the attempt folder goes to the output whose base name it starts with (longest match
## first, so "elev_dif.flt" goes with outDif "elev_dif" rather than outDTM "elev") -- including
## Align's extra outputs named from a base, such as <outOutlier>DSM (the canopy-based outlier
## raster) and <outDif>_0 (the unaligned DoD); outOutlier's .flt moves last.

move_file <- function(from, to) {
  if (file.exists(to)) unlink(to)
  if (!file.rename(from, to)) {
    if (!file.copy(from, to, overwrite = TRUE)) stop("Could not move ", from, " to ", to)
    unlink(from)
  }
}

chosen <- attempt_paths(best)
files <- setdiff(list.files(chosen$dir), c("align.log", "input_align.txt"))
bases <- vapply(chosen$out, basename, character(1))
bases <- bases[order(-nchar(bases))]
moves <- list()
for (f in files) {
  hit <- names(bases)[vapply(bases, function(b) startsWith(f, b), logical(1))]
  if (length(hit)) {
    moves[[f]] <- file.path(dirname(out_bases[[hit[1]]]), f)
  }
}
outlier_flt <- basename(flt_path(chosen$out$outOutlier))
for (f in c(setdiff(names(moves), outlier_flt), intersect(names(moves), outlier_flt))) {
  dir.create(dirname(moves[[f]]), showWarnings = FALSE, recursive = TRUE)
  move_file(file.path(chosen$dir, f), moves[[f]])
}
message("Moved attempt ", best, "'s outputs to outDTM/outDif/outOutlier/outbins.")

if (!isTRUE(p$keep_attempt_rasters)) {
  for (n in setdiff(attempts$attempt, best)) {
    rasters <- list.files(attempt_paths(n)$dir, pattern = "\\.(flt|hdr|prj)$", full.names = TRUE)
    unlink(rasters)
  }
}

message("Attempt summaries: ", attempts_dir)
