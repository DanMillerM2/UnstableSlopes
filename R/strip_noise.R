## strip_noise.R
##
## Noise of one lidar survey measured inside the survey itself: where two flight passes overlap,
## each pass's ground points are triangulated on their own and the two surfaces differenced. The
## passes were flown minutes apart, so their difference holds no real change -- only what makes two
## surveys' DTMs differ for reasons other than the ground moving: look angle, which returns reach
## the ground under canopy, ground classification, strip georeferencing. It is an independent
## reference for the level of detection, tabulated by canopy-height and gradient class.
##
## 1. Flight passes. Ground points per tile, second of GPS time and PointSourceID come from
##    ground_returns.R's pass_tally.csv (pass_tally), or are tallied here for the sampled tiles.
##    noise_method chooses how passes are told apart; auto tries them in this order and records why
##    each one it skips failed:
##      pointsource  PointSourceID, checked: not all 0, more than one, not one ID per tile when GPS
##                   time shows several passes there. An ID flown again later (a gap in its GPS
##                   time longer than reflight_gap_min) is split into one pass per flight; IDs flown at the
##                   same time (e.g. one per scanner channel) are merged into one pass.
##      gpstime      GPS time alone: a new pass starts after a gap longer than pass_gap_s (a turn
##                   between flight lines). Needs a point format with GPS time (not 0 or 2).
##      none         no pass information: nothing is measured, and the LoD rests on Align's
##                   stable-ground estimate alone (the summary says so).
##    Each method also needs at least min_overlap_tiles tiles where two passes each have
##    min_pass_points ground points. Passes lasting longer than max_pass_min, or covering tiles in
##    no line (elongation), are flagged in the summary.
## 2. A sample of the tiles with overlapping passes (tile_fraction of them, between min_tiles and
##    max_tiles; seed) is read with a buffer_m margin from neighbouring tiles. Each pass's ground
##    points (ground_classes, withheld points left out) are triangulated (edges up to max_edge_m)
##    and rasterized on a res_m grid in the point cloud's own CRS.
## 3. Every pair of passes overlapping on at least min_pair_cells cells is differenced (m). Each
##    cell is classed by canopy height (dsm - dtm of this survey; canopy_breaks, m) and gradient
##    (gradient_raster, or the gradient of template_raster smoothed over slope_radius_m;
##    slope_breaks), read at the cell centre after projecting it into template_raster's CRS, and by
##    ground-return density (density_breaks, returns/m2): each pass's ground points are counted per
##    cell and averaged over a square window about density_window_m wide, and a pair's cell takes the
##    sparser pass's density; the class also records the mean of 1/d_a + 1/d_b over its cells, which
##    the pair's variance scales with (IQR^2 = 2a + b (1/d_a + 1/d_b) gives one pass's a + b/d).
##    Noise against single-pass density lets it be read at a full DTM's density, which combines
##    the passes.
##
## Outputs in output_dir, named <label>_noise_*:
##   _summary.txt  method used, why others were rejected, pass and overlap checks
##   _passes.csv   one row per pass: PointSourceIDs, GPS time span, ground points, tiles, elongation
##   _pairs.csv    one row per tile and pass pair: cells, median offset, IQR (m)
##   _hist.csv     differences by class as a histogram (bin_m bins over +-range_m, plus tails)
##   _classes.csv  per canopy and gradient class: cells, median, IQR, 1st/99th percentiles, share
##                 beyond 3/5/10/20 IQRs
##   _density_classes.csv  the same per canopy, gradient and density class, with the mean density
##                 and mean 1/d_a + 1/d_b
##   tiles/        one .rds per sampled tile (resume: existing ones are reused)
## The .lax index files lasR writes beside the LAZ tiles it reads are removed afterwards (those
## that existed before the run are kept).
## The run is skipped if _classes.csv exists; delete it (and tiles/) after changing a parameter.
##
## Requires: lasR, terra, sf, lidR (LAS headers only), future, future.apply
##
##   Rscript strip_noise.R path/to/params.txt        (see strip_noise_params_template.txt)

library(terra)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}
config_path <- file.path(get_script_dir(), "strip_noise_params_template.txt")
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) config_path <- cli_args[1]
if (!file.exists(config_path)) stop("Parameter file not found: ", config_path)

read_param_file <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  pos <- regexpr(":", lines, fixed = TRUE)
  if (any(pos < 0)) stop("Malformed line in parameter file (expected \"keyword: value\"): ", lines[pos < 0][1])
  keys <- trimws(substr(lines, 1, pos - 1))
  vals <- trimws(sub("\\s+#.*$", "", trimws(substr(lines, pos + 1, nchar(lines)))))
  if (anyDuplicated(keys)) stop("Duplicate keyword(s) in ", path, ": ", paste(unique(keys[duplicated(keys)]), collapse = ", "))
  setNames(as.list(vals), keys)
}

param_specs <- list(
  input_dir  = list(type = "character", required = TRUE),   # folder of LAZ/LAS tiles (one survey)
  output_dir = list(type = "character", required = TRUE),
  label      = list(type = "character", required = TRUE),   # survey label, prefixes the outputs
  # Grid the classes are read on (the reference DTM, in the target CRS), and this survey's DTM and
  # DSM in the same CRS (canopy height = dsm - dtm)
  template_raster = list(type = "character", required = TRUE),
  dtm = list(type = "character", required = TRUE),
  dsm = list(type = "character", required = TRUE),
  gradient_raster = list(type = "character", default = "NOFILE"),  # e.g. Align's <ref DTM>_grad15m
  slope_radius_m  = list(type = "numeric", default = 15),         # used when gradient_raster is NOFILE
  pass_tally = list(type = "character", default = "NOFILE"),       # ground_returns.R's pass_tally.csv
  noise_method = list(type = "character", default = "auto"),       # auto, pointsource, gpstime, none
  pass_gap_s = list(type = "numeric", default = 20),
  reflight_gap_min = list(type = "numeric", default = 10),         # an ID flown again after a gap this long is a new pass
  max_pass_min = list(type = "numeric", default = 60),
  min_pass_points = list(type = "integer", default = 2000),        # ground points of a pass in a tile
  min_overlap_tiles = list(type = "integer", default = 3),
  min_pair_cells = list(type = "integer", default = 1000),
  tile_fraction = list(type = "numeric", default = 0.15),
  min_tiles = list(type = "integer", default = 10),
  max_tiles = list(type = "integer", default = 60),
  seed = list(type = "integer", default = 1),
  res_m = list(type = "numeric", default = 1.0),
  max_edge_m = list(type = "numeric", default = 10),
  buffer_m = list(type = "numeric", default = 20),
  ground_classes = list(type = "character", default = "2"),
  z_units = list(type = "character", default = "auto"),            # auto (= horizontal), m, ft, ftUS
  canopy_breaks = list(type = "character", default = "2 10 20"),
  slope_breaks = list(type = "character", default = "0.3 0.5"),
  density_breaks = list(type = "character", default = "0.25 0.75 2 4 8 16"),  # returns/m2; Align's classes, split above 4
  density_window_m = list(type = "numeric", default = 5),
  bin_m = list(type = "numeric", default = 0.002),
  range_m = list(type = "numeric", default = 5),
  clip_polygon = list(type = "character", default = "NOFILE"),     # only tiles intersecting it are sampled
  file_pattern = list(type = "character", default = "\\.la[sz]$"),
  n_workers = list(type = "integer", default = 8),
  threads_per_worker = list(type = "integer", default = 1)
)

build_params <- function(raw, specs, path) {
  unknown <- setdiff(names(raw), names(specs))
  if (length(unknown)) stop("Unrecognized parameter(s) in ", path, ": ", paste(unknown, collapse = ", "))
  out <- list()
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    if (nm %in% names(raw)) {
      v <- switch(sp$type, character = raw[[nm]], numeric = suppressWarnings(as.numeric(raw[[nm]])),
                  integer = suppressWarnings(as.integer(raw[[nm]])), logical = suppressWarnings(as.logical(raw[[nm]])))
      if (sp$type != "character" && is.na(v)) stop("Could not parse '", nm, "' (\"", raw[[nm]], "\") as ", sp$type)
      out[[nm]] <- v
    } else if (isTRUE(sp$required)) {
      stop("Required parameter '", nm, "' is missing from ", path)
    } else {
      out[[nm]] <- sp$default
    }
  }
  out
}
p <- build_params(read_param_file(config_path), param_specs, config_path)
message("Loaded parameters from: ", config_path)
is_set <- function(x) nzchar(x) && toupper(x) != "NOFILE"
nums <- function(s) as.numeric(strsplit(trimws(s), "[ ,]+")[[1]])
ground_cls <- as.integer(nums(p$ground_classes))
p$noise_method <- tolower(p$noise_method)
if (!p$noise_method %in% c("auto", "pointsource", "gpstime", "none")) stop("noise_method must be auto, pointsource, gpstime or none")
for (f in c(p$template_raster, p$dtm, p$dsm, if (is_set(p$gradient_raster)) p$gradient_raster)) {
  if (!file.exists(f)) stop("Not found: ", f)
}

out_file <- function(what) file.path(p$output_dir, paste0(p$label, "_noise_", what))
classes_csv <- out_file("classes.csv")
if (file.exists(classes_csv)) {
  message("Already done (delete it to rerun): ", classes_csv)
  quit(save = "no", status = 0)
}
tile_dir <- file.path(p$output_dir, "tiles")
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)
summary_lines <- character(0)
note <- function(...) {
  s <- paste0(...)
  summary_lines <<- c(summary_lines, s)
  message(s)
}

## ---- Tiles and their CRS --------------------------------------------------------------------
files <- list.files(p$input_dir, pattern = p$file_pattern, full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No files matching ", p$file_pattern, " in ", p$input_dir)
## lasR indexes the files it reads (.lax beside each LAZ); those this run creates are removed at the end
lax_before <- list.files(p$input_dir, "\\.lax$", full.names = TRUE, ignore.case = TRUE)
remove_new_lax <- function() {
  new <- setdiff(list.files(p$input_dir, "\\.lax$", full.names = TRUE, ignore.case = TRUE), lax_before)
  if (length(new)) message("Removed ", sum(file.remove(new)), " .lax index files this run created in ", p$input_dir)
}
ctg <- lidR::readLAScatalog(files, progress = FALSE)
hd <- ctg@data
tiles <- data.frame(file = normalizePath(hd$filename, winslash = "/"), tile = basename(hd$filename),
                    xmin = hd$Min.X, xmax = hd$Max.X, ymin = hd$Min.Y, ymax = hd$Max.Y)
tiles$xc <- (tiles$xmin + tiles$xmax) / 2
tiles$yc <- (tiles$ymin + tiles$ymax) / 2
crs_las <- sf::st_crs(ctg)
unit_tab <- c(`US survey foot` = 1200 / 3937, foot = 0.3048, metre = 1)
unit_m <- unit_tab[crs_las$units_gdal]
if (is.na(unit_m)) stop("Unrecognized horizontal unit of the point cloud's CRS: ", crs_las$units_gdal)
z_m <- switch(tolower(p$z_units), auto = unit_m, m = 1, ft = 0.3048, ftus = 1200 / 3937,
              stop("z_units must be auto, m, ft or ftUS"))
res <- p$res_m / unit_m
buf <- p$buffer_m / unit_m
max_edge <- p$max_edge_m / unit_m
template_crs <- crs(rast(p$template_raster))
fmt <- unique(hd$Point.Data.Format.ID)

## Tiles that may be sampled: those intersecting clip_polygon, else the template's extent
cand <- tiles
area <- if (is_set(p$clip_polygon)) sf::st_union(sf::st_read(p$clip_polygon, quiet = TRUE)) else
  sf::st_as_sfc(sf::st_bbox(as.vector(ext(rast(p$template_raster)))[c(1, 3, 2, 4)] |>
                              setNames(c("xmin", "ymin", "xmax", "ymax")), crs = sf::st_crs(template_crs)))
area <- sf::st_transform(area, crs_las)
boxes <- do.call(c, lapply(seq_len(nrow(tiles)), function(i)
  sf::st_as_sfc(sf::st_bbox(c(xmin = tiles$xmin[i], ymin = tiles$ymin[i], xmax = tiles$xmax[i], ymax = tiles$ymax[i]), crs = crs_las))))
cand <- tiles[lengths(sf::st_intersects(boxes, area)) > 0, ]
note(sprintf("Survey %s: %d tiles (%d in the study area), point format %s, CRS unit %s, z unit factor %.6f m",
             p$label, nrow(tiles), nrow(cand), paste(fmt, collapse = "/"), crs_las$units_gdal, z_m))

no_reference <- function() {
  note("Method: none -- no independent noise reference for this survey; the level of detection rests on Align's stable-ground estimate alone")
  writeLines(summary_lines, out_file("summary.txt"))
  write.csv(data.frame(), classes_csv, row.names = FALSE)
  quit(save = "no", status = 0)
}
if (p$noise_method == "none") no_reference()

## ---- 1. Pass tally ----------------------------------------------------------------------------
tally_tiles <- function(files_to_tally) {
  future::plan(future::multisession, workers = min(p$n_workers, length(files_to_tally)))
  out <- future.apply::future_lapply(files_to_tally, function(f) {
    suppressMessages(library(lasR))
    acc <- list()
    cb <- function(d) {
      if (!nrow(d)) return(0L)
      sec <- if (is.null(d$gpstime)) rep(-1, nrow(d)) else floor(d$gpstime)
      key <- sec * 65536 + d$PointSourceID
      u <- unique(key)
      acc[[length(acc) + 1]] <<- data.frame(second = floor(u / 65536), psid = u %% 65536, n = tabulate(match(key, u)))
      0L
    }
    exec(reader(filter = keep_class(ground_cls)) + callback(cb, expose = "tp", no_las_update = TRUE),
         on = f, ncores = p$threads_per_worker)
    a <- do.call(rbind, acc)
    if (is.null(a)) return(NULL)
    a <- aggregate(n ~ second + psid, data = a, FUN = sum)
    cbind(tile = basename(f), a)
  }, future.seed = TRUE, future.globals = c("ground_cls", "p"))
  future::plan(future::sequential)
  do.call(rbind, out)
}

tally <- if (is_set(p$pass_tally)) read.csv(p$pass_tally) else NULL
if (!is.null(tally)) {
  tally <- tally[tally$tile %in% cand$tile, ]
  note(sprintf("Pass tally from %s: %d of %d study-area tiles", p$pass_tally, length(unique(tally$tile)), nrow(cand)))
}
missing_tiles <- setdiff(cand$tile, tally$tile)
if (length(missing_tiles)) {
  ## Not tallied by ground_returns.R: tally a random sample here (enough to choose from)
  set.seed(p$seed)
  n_need <- min(length(missing_tiles), 2 * max(p$min_tiles, min(p$max_tiles, round(p$tile_fraction * nrow(cand)))))
  pick <- sample(missing_tiles, n_need)
  t0 <- Sys.time()
  tally <- rbind(tally, tally_tiles(cand$file[match(pick, cand$tile)]))
  note(sprintf("Tallied %d tiles here (%.1f min)", length(pick), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
if (is.null(tally) || !nrow(tally)) stop("No ground points found in the study-area tiles")

## ---- Pass definitions -----------------------------------------------------------------------
## Each method returns list(ok, reason, passes, unit_of): passes is a data.frame (pass, psids, t0,
## t1), and unit_of maps tally rows to passes.
has_time <- any(tally$second > 0) && length(unique(tally$second)) > 1
week_time <- has_time && max(tally$second) < 604800
segment <- function(secs, gap = p$pass_gap_s) {   # seconds -> segment number, a new one after each gap > gap
  u <- sort(unique(secs))
  seg <- cumsum(c(1, diff(u) > gap))
  seg[match(secs, u)]
}

passes_gpstime <- function() {
  if (!has_time) return(list(ok = FALSE, reason = if (all(fmt %in% c(0, 2))) "point format has no GPS time" else "GPS time is missing or constant"))
  seg <- segment(tally$second)
  ps <- do.call(rbind, lapply(split(seq_len(nrow(tally)), seg), function(i)
    data.frame(t0 = min(tally$second[i]), t1 = max(tally$second[i]),
               psids = paste(sort(unique(tally$psid[i])), collapse = " "))))
  ps$pass <- seq_len(nrow(ps))
  list(ok = TRUE, passes = ps, unit_of = seg, by_psid = FALSE)
}

passes_pointsource <- function() {
  ids <- unique(tally$psid)
  if (all(ids == 0)) return(list(ok = FALSE, reason = "PointSourceID is 0 for every ground point"))
  if (length(ids) == 1) return(list(ok = FALSE, reason = paste0("only one PointSourceID (", ids, ")")))
  ids_per_tile <- tapply(tally$psid, tally$tile, function(x) length(unique(x)))
  if (all(ids_per_tile == 1) && has_time) {
    multi <- mean(tapply(tally$second, tally$tile, function(s) max(segment(s))) > 1)
    if (multi > 0.2) return(list(ok = FALSE, reason = sprintf(
      "one PointSourceID per tile, but GPS time shows several passes in %.0f%% of tiles (IDs look like tile numbers)", 100 * multi)))
  }
  ## Units: an ID's flights (split at GPS-time gaps); units flown at the same time are one pass
  key <- if (has_time) paste(tally$psid, ave(tally$second, tally$psid, FUN = function(x) segment(x, 60 * p$reflight_gap_min))) else as.character(tally$psid)
  un <- do.call(rbind, lapply(split(seq_len(nrow(tally)), key), function(i)
    data.frame(key = key[i[1]], psid = tally$psid[i[1]], t0 = min(tally$second[i]), t1 = max(tally$second[i]))))
  un <- un[order(un$t0), ]
  grp <- seq_len(nrow(un))
  if (has_time) for (a in seq_len(nrow(un))) for (b in seq_len(nrow(un))) {
    if (b <= a) next
    ov <- min(un$t1[a], un$t1[b]) - max(un$t0[a], un$t0[b]) + 1
    if (ov > 0.5 * (min(un$t1[a] - un$t0[a], un$t1[b] - un$t0[b]) + 1)) grp[grp == grp[b]] <- grp[a]
  }
  un$pass <- match(grp, unique(grp))
  n_merged <- sum(duplicated(un$pass))
  ps <- do.call(rbind, lapply(split(un, un$pass), function(u)
    data.frame(pass = u$pass[1], t0 = min(u$t0), t1 = max(u$t1), psids = paste(sort(unique(u$psid)), collapse = " "))))
  n_split <- length(unique(un$psid[duplicated(un$psid)]))
  list(ok = TRUE, passes = ps, unit_of = un$pass[match(key, un$key)], by_psid = TRUE,
       info = sprintf("%d PointSourceIDs -> %d passes (%d IDs flown more than once, split by GPS time; %d merged with an ID flown at the same time)",
                      length(ids), nrow(ps), n_split, n_merged))
}

## Tiles where at least two passes have min_pass_points ground points
overlap_tiles <- function(def) {
  per <- aggregate(n ~ tile + pass, data = data.frame(tile = tally$tile, pass = def$unit_of, n = tally$n), FUN = sum)
  per <- per[per$n >= p$min_pass_points, ]
  k <- table(per$tile)
  list(per = per, tiles = names(k)[k >= 2])
}

methods <- if (p$noise_method == "auto") c("pointsource", "gpstime") else setdiff(p$noise_method, "none")
def <- NULL
for (m in methods) {
  d <- if (m == "pointsource") passes_pointsource() else passes_gpstime()
  if (d$ok) {
    ov <- overlap_tiles(d)
    if (length(ov$tiles) < p$min_overlap_tiles) {
      d <- list(ok = FALSE, reason = sprintf("only %d tiles with two passes of at least %d ground points (need %d)",
                                            length(ov$tiles), p$min_pass_points, p$min_overlap_tiles))
    }
  }
  if (d$ok) { def <- d; def$method <- m; def$ov <- ov; break }
  note("Method ", m, " rejected: ", d$reason)
  if (p$noise_method != "auto") stop("noise_method ", m, " is not usable for this survey: ", d$reason)
}
if (is.null(def)) no_reference()
note("Method: ", def$method, if (!is.null(def$info)) paste0(" -- ", def$info))
if (week_time) note("Warning: GPS time is GPS week time (all below 604,800 s); flights on different weeks can share times")

## Pass diagnostics: duration, ground points, tiles, elongation (line-like spread of tile centres)
ps <- def$passes
pt <- aggregate(n ~ tile + pass, data = data.frame(tile = tally$tile, pass = def$unit_of, n = tally$n), FUN = sum)
ps$ground_pts <- tapply(pt$n, pt$pass, sum)[as.character(ps$pass)]
big <- pt[pt$n >= p$min_pass_points, ]
ps$tiles <- as.integer(table(factor(big$pass, levels = ps$pass)))
ps$elongation <- vapply(ps$pass, function(k) {
  tl <- big$tile[big$pass == k]
  if (length(tl) < 4) return(NA_real_)
  ev <- eigen(cov(cbind(tiles$xc, tiles$yc)[match(tl, tiles$tile), ]), only.values = TRUE)$values
  sqrt(ev[1] / max(ev[2], 1e-9 * ev[1]))
}, numeric(1))
ps$duration_min <- if (has_time) (ps$t1 - ps$t0 + 1) / 60 else NA_real_
write.csv(ps[, c("pass", "psids", "t0", "t1", "duration_min", "ground_pts", "tiles", "elongation")],
          out_file("passes.csv"), row.names = FALSE)
note(sprintf("%d passes; %d of %d tallied tiles have two or more passes with at least %d ground points",
             nrow(ps), length(def$ov$tiles), length(unique(tally$tile)), p$min_pass_points))
if (has_time) {
  note(sprintf("Pass duration (min): median %.1f, max %.1f", median(ps$duration_min), max(ps$duration_min)))
  long <- ps$duration_min > p$max_pass_min
  if (any(long)) note(sprintf("Warning: %d passes last longer than %g min -- they may hold several flight lines", sum(long), p$max_pass_min))
}
if (any(!is.na(ps$elongation))) {
  note(sprintf("Pass elongation (tile-centre spread, major/minor axis; passes over 4+ tiles): median %.1f, %d of %d below 2",
               median(ps$elongation, na.rm = TRUE), sum(ps$elongation < 2, na.rm = TRUE), sum(!is.na(ps$elongation))))
}

## ---- 2. Sample tiles ----------------------------------------------------------------------------
set.seed(p$seed)
pool <- def$ov$tiles
n_sample <- min(length(pool), max(p$min_tiles, min(p$max_tiles, round(p$tile_fraction * nrow(cand)))))
sample_tiles <- sort(sample(pool, n_sample))
note(sprintf("Sampled %d of %d tiles with overlapping passes (seed %d)", n_sample, length(pool), p$seed))

## Filter for one pass: ground classes, not withheld, its PointSourceIDs and GPS-time span
pass_filter <- function(k) {
  r <- ps[ps$pass == k, ]
  f <- c(keep_class_str = paste("Classification %in%", paste(ground_cls, collapse = " ")), "Withheld == 0")
  if (def$by_psid) f <- c(f, paste("PointSourceID %in%", r$psids))
  if (has_time) f <- c(f, sprintf("gpstime >= %.0f", r$t0), sprintf("gpstime < %.0f", r$t1 + 1))
  unname(f)
}

canopy_breaks <- nums(p$canopy_breaks)
slope_breaks <- nums(p$slope_breaks)
lab <- function(b, pre) c(paste0(pre, "<", b[1]), if (length(b) > 1) paste0(pre, b[-length(b)], "-", b[-1]), paste0(pre, ">=", b[length(b)]))
canopy_labels <- lab(canopy_breaks, "")
slope_labels <- lab(slope_breaks, "g")
density_breaks <- nums(p$density_breaks)
density_labels <- lab(density_breaks, "d")
dens_win <- 2 * floor(p$density_window_m / p$res_m / 2) + 1     # cells, odd
nb <- round(2 * p$range_m / p$bin_m)
n_class <- length(canopy_labels) * length(slope_labels) * length(density_labels)

## ---- 3. One tile ------------------------------------------------------------------------------
noise_tile <- function(tile_name) {
  rds <- file.path(tile_dir, paste0(p$label, "_", sub("\\.la[sz]$", "", tile_name, ignore.case = TRUE), ".rds"))
  if (file.exists(rds)) {
    old <- readRDS(rds)
    if (identical(nrow(old$hist), as.integer(n_class)) && !is.null(old$isum)) return(rds)
  }
  suppressMessages({library(lasR); library(terra)})
  t0 <- Sys.time()
  t <- tiles[tiles$tile == tile_name, ]
  core <- c(floor(t$xmin / res) * res, ceiling(t$xmax / res) * res, floor(t$ymin / res) * res, ceiling(t$ymax / res) * res)
  rect <- core + c(-buf, buf, -buf, buf)
  src <- tiles$file[tiles$xmax >= rect[1] & tiles$xmin <= rect[2] & tiles$ymax >= rect[3] & tiles$ymin <= rect[4]]
  kp <- def$ov$per$pass[def$ov$per$tile == tile_name]
  tmp <- file.path(tempdir(), paste0("sn_", sub("\\.la[sz]$", "", tile_name, ignore.case = TRUE), "_p", kp, ".tif"))
  tmp_n <- sub("[.]tif$", "_n.tif", tmp)
  pipe <- reader_rectangles(rect[1], rect[3], rect[2], rect[4])
  for (j in seq_along(kp)) {
    tri <- triangulate(max_edge, filter = pass_filter(kp[j]))
    pipe <- pipe + tri + lasR::rasterize(res, tri, ofile = tmp[j]) +
      lasR::rasterize(res, "count", filter = pass_filter(kp[j]), ofile = tmp_n[j])
  }
  exec(pipe, on = src, ncores = p$threads_per_worker)
  ok <- file.exists(tmp) & file.exists(tmp_n)
  r <- crop(rast(tmp[ok]), ext(core), snap = "near")
  ## Ground returns per m2 of each pass, averaged over the density window (on the buffered grid)
  dn <- lapply(tmp_n[ok], function(f) {
    n <- subst(rast(f), NA, 0)
    terra::focal(n, dens_win, fun = "mean", na.rm = TRUE) / p$res_m^2
  })
  dn <- values(resample(rast(dn), r, method = "near"))
  kp <- kp[ok]
  v <- values(r) * z_m
  unlink(c(tmp, tmp_n))
  have <- rowSums(!is.na(v)) >= 2
  cells <- which(have)
  v <- v[cells, , drop = FALSE]
  dn <- dn[cells, , drop = FALSE]
  ## Classes at the cell centres, in the template's CRS
  xy <- xyFromCell(r, cells)
  xy_t <- sf::sf_project(crs_las$wkt, template_crs, xy)
  cd <- extract(c(rast(p$dsm), rast(p$dtm)), xy_t, method = "simple")
  canopy <- cd[[1]] - cd[[2]]
  if (is_set(p$gradient_raster)) {
    g <- extract(rast(p$gradient_raster), xy_t, method = "simple")[[1]]
  } else {
    tp <- rast(p$template_raster)
    e <- ext(min(xy_t[, 1]), max(xy_t[, 1]), min(xy_t[, 2]), max(xy_t[, 2])) + 2 * p$slope_radius_m
    d <- crop(tp, e, snap = "out")
    w <- focalMat(d, p$slope_radius_m, "circle")
    w[w > 0] <- 1; w[w == 0] <- NA
    g <- extract(tan(terrain(focal(d, w, fun = "mean", na.rm = TRUE), "slope", unit = "radians")), xy_t, method = "simple")[[1]]
  }
  cs <- findInterval(canopy, canopy_breaks) * length(slope_labels) + findInterval(g, slope_breaks)
  cs[is.na(canopy) | is.na(g)] <- NA
  hist <- matrix(0, n_class, nb + 2)    # columns: below -range, nb bins, above +range
  dsum <- numeric(n_class)              # sum of density, for each class's mean
  isum <- numeric(n_class)              # sum of 1/d_a + 1/d_b (the pair's variance scales with it)
  pairs <- list()
  if (length(kp) >= 2) for (a in 1:(length(kp) - 1)) for (b in (a + 1):length(kp)) {
    i <- which(!is.na(v[, a]) & !is.na(v[, b]))
    if (length(i) < p$min_pair_cells) next
    dz <- v[i, a] - v[i, b]
    pairs[[length(pairs) + 1]] <- data.frame(tile = tile_name, pass_a = kp[a], pass_b = kp[b], cells = length(i),
                                             median_m = median(dz), iqr_m = IQR(dz))
    d <- pmin(dn[i, a], dn[i, b])       # the sparser pass
    cls <- cs[i] * length(density_labels) + findInterval(d, density_breaks) + 1
    k <- !is.na(cls)
    bin <- pmin(pmax(floor((dz[k] + p$range_m) / p$bin_m) + 2, 1), nb + 2)
    hist <- hist + matrix(tabulate((bin - 1) * n_class + cls[k], nbins = n_class * (nb + 2)), n_class)
    fk <- factor(cls[k], levels = seq_len(n_class))
    dsum <- dsum + vapply(split(d[k], fk), sum, numeric(1))
    inv <- 1 / pmax(dn[i, a], 0.02) + 1 / pmax(dn[i, b], 0.02)
    isum <- isum + vapply(split(inv[k], fk), sum, numeric(1))
  }
  saveRDS(list(tile = tile_name, passes = length(kp), pairs = do.call(rbind, pairs), hist = hist, dsum = dsum, isum = isum,
               seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))), rds)
  rds
}

t_start <- Sys.time()
future::plan(future::multisession, workers = min(p$n_workers, length(sample_tiles)))
shared <- c("tiles", "def", "ps", "p", "res", "buf", "max_edge", "z_m", "crs_las", "template_crs", "tile_dir",
            "ground_cls", "has_time", "pass_filter", "canopy_breaks", "slope_breaks", "slope_labels",
            "density_breaks", "density_labels", "dens_win", "nb", "n_class", "is_set")
rds_files <- tryCatch(future.apply::future_lapply(sample_tiles, noise_tile, future.seed = TRUE, future.globals = shared),
                      finally = { future::plan(future::sequential); remove_new_lax() })
res_tiles <- lapply(unlist(rds_files), readRDS)
note(sprintf("Built and differenced %d tiles in %.1f min (median %.0f s per tile, one worker)",
             length(res_tiles), as.numeric(difftime(Sys.time(), t_start, units = "mins")),
             median(vapply(res_tiles, `[[`, numeric(1), "seconds"))))

## ---- Outputs ---------------------------------------------------------------------------------
pairs <- do.call(rbind, lapply(res_tiles, `[[`, "pairs"))
write.csv(pairs, out_file("pairs.csv"), row.names = FALSE)
H <- Reduce(`+`, lapply(res_tiles, `[[`, "hist"))
D <- Reduce(`+`, lapply(res_tiles, `[[`, "dsum"))
I <- Reduce(`+`, lapply(res_tiles, `[[`, "isum"))
edges <- -p$range_m + p$bin_m * (0:nb)
grid_cls <- expand.grid(density = density_labels, slope = slope_labels, canopy = canopy_labels,
                        stringsAsFactors = FALSE)[, c("canopy", "slope", "density")]
hist_long <- data.frame(grid_cls[rep(seq_len(n_class), nb + 2), ],
                        bin_lo = rep(c(-Inf, edges[1:nb], p$range_m), each = n_class),
                        bin_hi = rep(c(-p$range_m, edges[2:(nb + 1)], Inf), each = n_class), n = as.vector(H))
write.csv(hist_long[hist_long$n > 0, ], out_file("hist.csv"), row.names = FALSE)

## Class statistics from the histograms (to bin_m)
hq <- function(h, q) {           # quantile from a histogram row, interpolated within the bin
  cs <- cumsum(h) / sum(h)
  j <- which(cs >= q)[1]
  lo <- c(-p$range_m - p$bin_m, edges[1:nb], p$range_m)[j]
  prev <- if (j == 1) 0 else cs[j - 1]
  lo + p$bin_m * (q - prev) / max(cs[j] - prev, 1e-12)
}
mid <- c(-p$range_m - p$bin_m, (edges[1:nb] + edges[2:(nb + 1)]) / 2, p$range_m + p$bin_m)
stats_of <- function(H, cl) do.call(rbind, lapply(seq_len(nrow(H)), function(k) {
  h <- H[k, ]
  if (sum(h) == 0) return(NULL)
  med <- hq(h, 0.5); iqr <- hq(h, 0.75) - hq(h, 0.25)
  z <- abs(mid - med) / iqr
  data.frame(cl[k, , drop = FALSE], cells = sum(h), median_m = med, iqr_m = iqr, p01_m = hq(h, 0.01), p99_m = hq(h, 0.99),
             gt3 = sum(h[z > 3]) / sum(h), gt5 = sum(h[z > 5]) / sum(h), gt10 = sum(h[z > 10]) / sum(h),
             gt20 = sum(h[z > 20]) / sum(h))
}))
## Canopy and gradient classes (density pooled), and the full classes with their mean density
cs_key <- paste(grid_cls$canopy, grid_cls$slope)
H_cs <- rowsum(H, factor(cs_key, levels = unique(cs_key)), reorder = FALSE)
cls_stats <- stats_of(H_cs, grid_cls[!duplicated(cs_key), c("canopy", "slope")])
den_cls <- cbind(grid_cls, density_mean = D / pmax(rowSums(H), 1), inv_density_sum = I / pmax(rowSums(H), 1))
den_stats <- stats_of(H, den_cls)
write.csv(den_stats, out_file("density_classes.csv"), row.names = FALSE)
off <- abs(pairs$median_m)
note(sprintf("%d pass pairs over %.1f M cells; |median offset| (m): median %.3f, 90th percentile %.3f, max %.3f; pair IQR median %.3f m",
             nrow(pairs), sum(pairs$cells) / 1e6, median(off), quantile(off, 0.9), max(off), median(pairs$iqr_m)))
writeLines(summary_lines, out_file("summary.txt"))
write.csv(cls_stats, classes_csv, row.names = FALSE)   # last: marks the run complete
print(cls_stats, digits = 3, row.names = FALSE)
