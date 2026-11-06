## scar_review.R
##
## Builds a scar-review page for one run: a stratified random sample of LShunter's candidate
## patches, an 8-panel lidar image of each, optionally NAIP aerial photos of each, and the labelling
## page (scar_review_page.html) that shows them one at a time for a verdict -- landslide scar, other
## real change, artifact, can't tell -- with the outline's fit and the reasons for a "not a scar".
## The labels are the ground truth for the false-positive rate and for testing which patch
## attributes separate scars from the rest.
##
## Writes, in out_dir:
##   index.html          the page (scar_review_page.html with the title and the sample filled in)
##   patches.json        the sample: every patch's attributes, its stratum, the stratum
##                       populations (for weighting estimates back to the whole run), the strata
##                       definitions, and whether it has a NAIP image
##   img/p<label>.jpg    lidar panels: both surveys' hillshades (the align survey's after Align),
##                       centred elevation change -1..+1 m, a 700 m context hillshade, both
##                       surveys' canopy height 0-50 m and ground returns per m2
##   img/p<label>_naip.jpg  with naip: TRUE (naip_images.R), plus img/naip_images.csv
## The page saves labels only when published as a claude.ai artifact with the db and user
## capabilities (the landslide-pipelines skill does this): collection "labels", one document per
## patch and reviewer, id p<label>__<user id>, fields verdict (scar/change/artifact/unsure), outline,
## reasons[], note, stratum. Opened from disk, the page shows the images but can't save.
##
## Inputs come from the run folder's configs (run_site.R writes them): the patch shapefile
## (patch_polygons_params.txt), patch raster, centred difference and both surveys' DTMs and DSMs
## (lshunter_params.txt), the aligned DTM, ground-return count rasters and survey labels
## (align_pipeline.txt). Any of them can be given directly instead (see the template).
##
## Sampling: each patch is put in one level of every stratum factor (factor_<name>: level =
## condition; ...; the first true condition wins, "otherwise" catches the rest, and a condition that
## is NA counts as false), the strata are their combinations, and per_stratum patches are drawn from
## each (all of them when there are fewer), with seed, then shuffled. Conditions are R expressions
## over the shapefile's fields. With the same shapefile, strata and seed the sample is the same, so
## a rebuilt page keeps its labels. Lidar images newer than the shapefile are kept.
##
## Requires: terra, sf, jsonlite (and future, future.apply for NAIP).
##
##   Rscript scar_review.R path/to/params.txt        (see scar_review_params_template.txt)

suppressMessages({ library(terra); library(sf); library(jsonlite) })
terraOptions(progress = 0)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}
script_dir <- get_script_dir()
config_path <- file.path(script_dir, "scar_review_params_template.txt")
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) config_path <- cli_args[1]
if (!file.exists(config_path)) stop("Parameter file not found: ", config_path)

## "keyword: value" files; [section] lines (run_align_pipeline.R's config) start a section
read_kv <- function(path, sections = FALSE) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  out <- list(); sec <- ""
  for (line in lines) {
    if (sections && grepl("^\\[.*\\]$", line)) { sec <- sub("^\\[(.*)\\]$", "\\1", line); next }
    pos <- regexpr(":", line, fixed = TRUE)
    if (pos < 0) stop("Malformed line in ", path, " (expected \"keyword: value\"): ", line)
    key <- trimws(substr(line, 1, pos - 1))
    if (sections) key <- paste0(sec, ".", key)
    if (key %in% names(out)) stop("Duplicate keyword '", key, "' in ", path)
    out[[key]] <- trimws(sub("\\s+#.*$", "", trimws(substr(line, pos + 1, nchar(line)))))
  }
  out
}
is_set <- function(x) !is.null(x) && !is.na(x) && nzchar(x) && !toupper(x) %in% c("NOFILE", "NONE")

input_keys <- c("patch_shp", "patch_raster", "cendif", "ref_dtm", "aligned_dtm", "align_dtm", "ref_dsm",
                "align_dsm", "ref_ground_count", "align_ground_count", "ref_label", "align_label")
p <- read_kv(config_path)
mask_keys <- c("water_mask", "ref_water_count", "align_water_count", "coverage_min", "coverage_block",
               "water_min_density", "water_window")
known <- c("run_dir", "out_dir", input_keys, mask_keys, "title", "subtitle", "filter", "per_stratum", "seed",
           "naip", "naip_years", "naip_survey_years", "naip_n_workers")
unknown <- setdiff(names(p)[!startsWith(names(p), "factor_")], known)
if (length(unknown)) stop("Unrecognized parameter(s) in ", config_path, ": ", paste(unknown, collapse = ", "))
if (!is_set(p$out_dir)) stop("'out_dir' is required in ", config_path)
message("Loaded parameters from: ", config_path)

## ---- Inputs: from the run folder's configs, unless given -----------------------------------------
from_run <- list()
if (is_set(p$run_dir)) {
  cfg <- file.path(p$run_dir, "configs")
  if (!dir.exists(cfg)) stop("No configs folder in run_dir (written by run_site.R): ", cfg)
  lsh <- read_kv(file.path(cfg, "lshunter_params.txt"))
  alp <- read_kv(file.path(cfg, "align_pipeline.txt"), sections = TRUE)
  ply <- read_kv(file.path(cfg, "patch_polygons_params.txt"))
  from_run <- list(patch_shp = ply$out_shp, patch_raster = lsh$out_patch, cendif = lsh$cen_dif,
    ref_dtm = lsh$canopy_ref_dem, align_dtm = lsh$canopy_align_dem, ref_dsm = lsh$canopy_ref_dsm,
    align_dsm = lsh$canopy_align_dsm, aligned_dtm = alp[["align.outDTM"]],
    ref_ground_count = if (is_set(alp[["ref_set.output_dir"]])) file.path(alp[["ref_set.output_dir"]], "ground_count_ref_grid.flt"),
    align_ground_count = if (is_set(alp[["align_set.output_dir"]])) file.path(alp[["align_set.output_dir"]], "ground_count_ref_grid.flt"),
    ref_label = alp[["strip_noise.ref_label"]], align_label = alp[["strip_noise.align_label"]],
    water_mask = alp[["align.outWaterMask"]], coverage_min = alp[["align.coverageMin"]],
    coverage_block = alp[["align.coverageBlock"]], water_min_density = alp[["align.waterMinDensity"]],
    water_window = alp[["align.waterWindow"]])
}
inp <- list()
for (k in input_keys) inp[[k]] <- if (is_set(p[[k]])) p[[k]] else from_run[[k]]
missing <- input_keys[!vapply(inp, is_set, TRUE)]
if (length(missing)) stop("Not found in run_dir's configs and not given: ", paste(missing, collapse = ", "))
flt <- function(x) paste0(sub("\\.flt$", "", x, ignore.case = TRUE), ".flt")
for (k in setdiff(input_keys, c("patch_shp", "ref_label", "align_label"))) {
  inp[[k]] <- flt(inp[[k]])
  if (!file.exists(inp[[k]])) stop(k, " not found: ", inp[[k]])
}
inp$patch_shp <- paste0(sub("\\.shp$", "", inp$patch_shp, ignore.case = TRUE), ".shp")
if (!file.exists(inp$patch_shp)) stop("patch_shp not found: ", inp$patch_shp)
rl <- inp$ref_label; al <- inp$align_label

## Masks, as Align draws them (Align 1.1.0+ leaves these cells out of everything): coverage gaps --
## blocks of coverage_block m (counted from the reference DTM's top-left corner) where a survey has
## fewer than coverage_min ground returns per m2 (missing or snow-covered tiles) -- per survey; and
## water, from Align's water mask (codes 1-3: polygons, water returns, flat open water) and from each
## survey's water-return counts (water_count_ref_grid beside the ground counts; at least
## water_min_density returns/m2 over a (2 * water_window + 1)-cell square). Run first with an older
## Align, the coverage gaps and the water counts still show.
opt <- function(k, default = NULL) { v <- if (is_set(p[[k]])) p[[k]] else from_run[[k]]; if (is_set(v)) v else default }
cov_min <- as.numeric(opt("coverage_min", 0.02)); cov_block <- as.numeric(opt("coverage_block", 50))
wat_min <- as.numeric(opt("water_min_density", 0.1)); wat_win <- as.integer(opt("water_window", 1))
wmask_f <- opt("water_mask")
if (!is.null(wmask_f)) { wmask_f <- flt(wmask_f); if (!file.exists(wmask_f)) { message("Water mask not found, not used: ", wmask_f); wmask_f <- NULL } }
wcount_f <- list()
for (sv in c("ref", "align")) {
  f <- opt(paste0(sv, "_water_count"), file.path(dirname(inp[[paste0(sv, "_ground_count")]]), "water_count_ref_grid.flt"))
  f <- flt(f); if (file.exists(f)) wcount_f[[sv]] <- f
}
message("Masks: coverage gaps (< ", cov_min, " returns/m2 per ", cov_block, " m block)",
        if (!is.null(wmask_f)) paste0("; Align's water mask ", wmask_f),
        if (length(wcount_f)) paste0("; water returns (", paste(names(wcount_f), collapse = ", "), ")"))

out_dir <- p$out_dir
img_dir <- file.path(out_dir, "img")
dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Strata and sample -----------------------------------------------------------------------------
poly <- st_read(inp$patch_shp, quiet = TRUE)
if (!"label" %in% names(poly)) stop(inp$patch_shp, " has no 'label' field (expected patch_polygons.R's shapefile)")
# With LShunter's DEPOSITS the shapefile also holds deposits (type "deposit"), whose labels can
# repeat the scars' (each numbers its own raster) -- images and label ids use the label, so the
# sample, its populations and the dashed neighbours are scars only
if ("type" %in% names(poly)) poly <- poly[!is.na(poly$type) & poly$type == "scar", ]
for (nm in names(poly)) if (is.numeric(poly[[nm]])) poly[[nm]][!is.na(poly[[nm]]) & poly[[nm]] <= -9998] <- NA
dat <- st_drop_geometry(poly)
cond <- function(expr, what) {
  v <- tryCatch(eval(parse(text = expr), dat, baseenv()), error = function(e)
    stop("Can't evaluate ", what, " '", expr, "': ", conditionMessage(e), ". Fields: ", paste(names(dat), collapse = ", ")))
  if (length(v) == 1) v <- rep(v, nrow(dat))
  !is.na(v) & as.logical(v)
}
if (is_set(p$filter)) {
  keep <- cond(p$filter, "filter")
  poly <- poly[keep, ]; dat <- dat[keep, , drop = FALSE]
  message("filter ", p$filter, ": ", nrow(poly), " patches")
}
fdefs <- p[startsWith(names(p), "factor_")]
if (!length(fdefs)) {   # the Site 2 design: sign-flip FDR x canopy change. (sim_fdr, the null
  # simulation's rate, is empty unless the run set lshunter_null_realizations.)
  fdefs <- list(factor_fdr = "sig = flip_fdr <= 0.05; nonsig = otherwise",
                factor_canopy = "loss = canopy_chg < -10; growth = canopy_chg > 5; mid = otherwise")
}
strata_desc <- character(0)
lev <- list()
for (fn in names(fdefs)) {
  parts <- trimws(strsplit(fdefs[[fn]], ";", fixed = TRUE)[[1]]); parts <- parts[nzchar(parts)]
  # "name = condition": the name is a plain word, so "<=" or "==" in the condition isn't the split
  m <- regmatches(parts, regexec("^([A-Za-z0-9_.-]+)\\s*=\\s*(.+)$", parts))
  if (any(lengths(m) != 3)) stop(fn, ": expected 'level = condition; ...; level = otherwise': ", fdefs[[fn]])
  nms <- vapply(m, `[`, "", 2); exprs <- trimws(vapply(m, `[`, "", 3))
  if (anyDuplicated(nms)) stop(fn, ": duplicate level names")
  v <- rep(NA_character_, nrow(dat))
  for (i in seq_along(nms)) {
    hit <- if (tolower(exprs[i]) == "otherwise") rep(TRUE, nrow(dat)) else cond(exprs[i], fn)
    v[is.na(v) & hit] <- nms[i]
  }
  if (anyNA(v)) stop(fn, ": ", sum(is.na(v)), " patches match no level; end it with 'name = otherwise'")
  lev[[sub("^factor_", "", fn)]] <- v
  strata_desc <- c(strata_desc, paste0(sub("^factor_", "", fn), " (", paste0(nms, ": ", exprs, collapse = "; "), ")"))
}
poly$stratum <- do.call(paste, c(lev, sep = "_"))
pop <- table(poly$stratum); print(pop)
per <- as.integer(if (is_set(p$per_stratum)) p$per_stratum else 40)
set.seed(as.integer(if (is_set(p$seed)) p$seed else 2026))
pick <- unlist(lapply(names(pop), function(s) { i <- which(poly$stratum == s); if (length(i) <= per) i else sample(i, per) }))
pick <- sample(pick)   # shuffled presentation order
message("sampled ", length(pick), " of ", nrow(poly), " patches, up to ", per, " per stratum")

## ---- Lidar panels ------------------------------------------------------------------------------------
e1 <- rast(inp$ref_dtm); e2a <- rast(inp$aligned_dtm); s1 <- rast(inp$ref_dsm)
e2 <- rast(inp$align_dtm); s2 <- rast(inp$align_dsm)
cd <- rast(inp$cendif); pr <- rast(inp$patch_raster)
g1 <- rast(inp$ref_ground_count); g2 <- rast(inp$align_ground_count)
wm <- if (!is.null(wmask_f)) rast(wmask_f)
wc <- lapply(wcount_f, rast)
pv <- vect(poly)

## Coverage gap and water masks over a window (1 = masked, 0 = not, NA off the DTM), on tmpl's grid
cs <- res(e1)[1]; bsz <- max(1L, as.integer(round(cov_block / cs))); bw <- bsz * cs; oe <- ext(e1)
grab <- function(r, x) extend(crop(r, x), x)
gap_mask <- function(g, w, tmpl) {
  bx <- ext(xmin(oe) + floor((w$xmin - xmin(oe)) / bw) * bw, xmin(oe) + ceiling((w$xmax - xmin(oe)) / bw) * bw,
            ymax(oe) - ceiling((ymax(oe) - w$ymin) / bw) * bw, ymax(oe) - floor((ymax(oe) - w$ymax) / bw) * bw)
  bx <- intersect(bx, oe)
  valid <- !is.na(grab(e1, bx)); n <- grab(g, bx)
  s <- aggregate(ifel(valid & !is.na(n), n, 0), bsz, sum); nv <- aggregate(valid, bsz, sum)
  low <- ifel(nv > 0 & s / (nv * cs^2) < cov_min, 1, 0)
  m <- resample(low, tmpl, "near"); m[is.na(tmpl)] <- NA
  # A block straddling a tile edge can pass on the returns in its covered part, so the gap grows
  # into the cells with no returns over a 5 x 5 window that it touches (a missing tile's edge is
  # straight; dense canopy next to a gap, also without returns, can't be told from it either)
  if (global(m, "max", na.rm = TRUE)[[1]] == 1) {
    n <- on_tmpl(g, tmpl); n[is.na(n)] <- 0
    empty <- focal(n, 5, sum, na.rm = TRUE) == 0
    pa <- patches(ifel(m == 1 | empty, 1, NA), directions = 8)
    keep <- unique(as.vector(values(mask(pa, ifel(m == 1, 1, NA)), na.rm = TRUE)))
    m <- ifel(!is.na(pa) & pa %in% keep, 1, m)
  }
  m
}
# r on tmpl's grid (a window at the DTM's edge is cut to it, so crop(r, w) can be larger)
on_tmpl <- function(r, tmpl) { x <- grab(r, ext(tmpl))
  if (!compareGeom(x, tmpl, stopOnError = FALSE)) x <- resample(x, tmpl, "near"); x }
water_mask <- function(w, tmpl) {
  m <- ifel(is.na(tmpl), NA, 0)
  if (!is.null(wm)) { a <- on_tmpl(wm, tmpl); m <- ifel(!is.na(a) & a >= 1 & a <= 3, 1, m) }
  k <- 2L * wat_win + 1L
  for (r in wc) {
    n <- grab(r, extend(ext(tmpl), k * cs)); n[is.na(n)] <- 0
    d <- on_tmpl(focal(n, k, sum) / (k * k * cs^2), tmpl)
    m <- ifel(!is.na(d) & d >= wat_min, 1, m)
  }
  m[is.na(tmpl)] <- NA; m
}
masked <- function(m, col) { x <- ifel(m == 1, 1, NA); if (!all(is.na(values(x)))) plot(x, col = col, legend = FALSE, axes = FALSE, add = TRUE) }

hshade <- function(x) { sl <- terrain(x, "slope", unit = "radians"); as <- terrain(x, "aspect", unit = "radians")
  mean(shade(sl, as, angle = 40, direction = c(225, 270, 315, 360))) }
grey <- gray.colors(100, 0.05, 0.98)
dodpal <- colorRampPalette(c("#8c1d18", "#d7604d", "#f4d6cc", "#f7f7f7", "#cfe0ee", "#5e9ac7", "#173f73"))(101)
canpal <- colorRampPalette(c("#f3efe0", "#c9d99a", "#6fa457", "#2f6b3a", "#123b24"))(51)
grnpal <- c("#b2182b", "#f1a340", "#fee8a6", "#c7e3c0", "#4dac6b")
st <- function(x) { q <- as.numeric(quantile(values(x), c(0.02, 0.98), na.rm = TRUE)); clamp(x, q[1], q[2], values = TRUE) }
dens <- function(n) focal(n, 3, sum, na.rm = TRUE) / (9 * prod(res(n)))   # ground returns per m2 over 3 x 3 cells
gb <- c(-0.01, 0.01, 0.5, 1, 3, 1e9)
num_or_na <- function(x) if (is.null(x)) NA else x

rows <- list(); t0 <- Sys.time(); made <- 0
for (k in seq_along(pick)) {
  pp <- poly[pick[k], ]; lab <- pp$label; bb <- st_bbox(pp)
  cx <- unname((bb$xmin + bb$xmax) / 2); cy <- unname((bb$ymin + bb$ymax) / 2)
  half <- max(45, 0.75 * max(bb$xmax - bb$xmin, bb$ymax - bb$ymin) + 20)
  w <- ext(cx - half, cx + half, cy - half, cy + half); wb <- ext(cx - half - 4, cx + half + 4, cy - half - 4, cy + half + 4)
  tmpl <- crop(e1, w)
  n1 <- crop(g1, w); n2 <- crop(g2, w); n1[is.na(n1) & !is.na(tmpl)] <- 0; n2[is.na(n2) & !is.na(tmpl)] <- 0
  # coverage: share of patch cells with >= 1 ground return in each survey
  pm <- crop(pr, w) == lab
  npm <- global(pm, "sum", na.rm = TRUE)[[1]]
  cov1 <- global(n1 >= 1 & pm, "sum", na.rm = TRUE)[[1]] / npm
  cov2 <- global(n2 >= 1 & pm, "sum", na.rm = TRUE)[[1]] / npm
  # masks: each survey's coverage gaps and water (both surveys'), and their share of the patch
  gm1 <- gap_mask(g1, w, tmpl); gm2 <- gap_mask(g2, w, tmpl); wtr <- water_mask(w, tmpl)
  share <- function(m) global(m == 1 & pm, "sum", na.rm = TRUE)[[1]] / npm
  gap1 <- share(gm1); gap2 <- share(gm2); wshare <- share(wtr)

  f <- file.path(img_dir, sprintf("p%d.jpg", lab))
  # redrawn when older than the shapefile or than this script (a change to the panels)
  if (!(file.exists(f) && file.mtime(f) > max(file.mtime(inp$patch_shp), file.mtime(file.path(script_dir, "scar_review.R"))))) {
    hs1 <- crop(hshade(crop(e1, wb)), w); hs2 <- crop(hshade(crop(e2a, wb)), w)
    dz <- clamp(crop(cd, w), -1, 1, values = TRUE)
    c1 <- clamp(crop(s1, w) - tmpl, 0, 50, values = TRUE)
    c2 <- resample(crop(s2, extend(w, 2)), tmpl, "bilinear") - resample(crop(e2, extend(w, 2)), tmpl, "bilinear")
    c2 <- clamp(c2, 0, 50, values = TRUE)
    d1 <- dens(n1); d2 <- dens(n2); d1[is.na(tmpl)] <- NA; d2[is.na(tmpl)] <- NA
    me <- pv[pv$label == lab]; others <- crop(pv[pv$label != lab], w)
    cw <- ext(cx - 350, cx + 350, cy - 350, cy + 350)
    ctx <- crop(hshade(aggregate(crop(e1, extend(cw, 10)), 3, mean)), cw)
    out <- function() {
      if (length(others)) lines(others, col = "#ffd400", lwd = 1.2, lty = 2)
      lines(me, col = "#ff00c8", lwd = 1.8)
      sbar(20, xy = "bottomleft", type = "bar", divs = 2, below = "m", cex = 1)
    }
    tmp <- sub("\\.jpg$", "_partial.jpg", f)
    jpeg(tmp, width = 1120, height = 600, quality = 80)
    par(mfrow = c(2, 4), oma = c(0, 0, 0, 0), bg = "white")
    m <- c(0.4, 0.4, 1.6, 0.4)
    plot(st(hs1), col = grey, legend = FALSE, axes = FALSE, mar = m, main = paste(rl, "hillshade")); out()
    plot(st(hs2), col = grey, legend = FALSE, axes = FALSE, mar = m, main = paste(al, "hillshade (aligned)")); out()
    plot(st(hs1), col = grey, legend = FALSE, axes = FALSE, mar = m, main = "Elevation change, -1 to +1 m")
    plot(dz, col = dodpal, range = c(-1, 1), alpha = 0.85, legend = FALSE, axes = FALSE, add = TRUE); out()
    plot(st(ctx), col = grey, legend = FALSE, axes = FALSE, mar = m, main = "Context, 700 m")
    rect(w$xmin, w$ymin, w$xmax, w$ymax, border = "#ff00c8", lwd = 2)
    sbar(200, xy = "bottomleft", type = "bar", divs = 2, below = "m", cex = 1)
    plot(c1, col = canpal, range = c(0, 50), legend = FALSE, axes = FALSE, mar = m, main = paste(rl, "canopy height, 0-50 m")); out()
    plot(c2, col = canpal, range = c(0, 50), legend = FALSE, axes = FALSE, mar = m, main = paste(al, "canopy height, 0-50 m")); out()
    # coverage gaps grey and water blue over the density: no returns there isn't thin ground cover
    plot(d1, col = grnpal, breaks = gb, legend = FALSE, axes = FALSE, mar = m, main = paste(rl, "ground returns per m2"))
    masked(gm1, "#9a9a9a"); masked(wtr, "#4a90d9"); out()
    plot(d2, col = grnpal, breaks = gb, legend = FALSE, axes = FALSE, mar = m, main = paste(al, "ground returns per m2"))
    masked(gm2, "#9a9a9a"); masked(wtr, "#4a90d9"); out()
    dev.off()
    file.rename(tmp, f); made <- made + 1
  }

  ll <- st_coordinates(st_transform(st_sfc(st_point(c(cx, cy)), crs = st_crs(poly)), 4326))
  rows[[k]] <- data.frame(id = sprintf("p%d", lab), label = lab, order = k, stratum = pp$stratum,
    x = round(cx, 1), y = round(cy, 1), lat = round(ll[2], 6), lon = round(ll[1], 6),
    area_m2 = num_or_na(pp$area_m2), length_m = round(num_or_na(pp$length_m), 1), width_m = round(num_or_na(pp$width_m), 1),
    mean_dz = num_or_na(pp$mean_dz), min_dz = num_or_na(pp$min_dz), vol_m3 = num_or_na(pp$vol_m3),
    sd_vol_m3 = num_or_na(pp$sd_vol_m3), grad_mean = num_or_na(pp$grad_mean), sim_fdr = num_or_na(pp$sim_fdr),
    flip_fdr = num_or_na(pp$flip_fdr), canopy_ref = num_or_na(pp$canopy_ref), canopy_aln = num_or_na(pp$canopy_aln),
    canopy_chg = num_or_na(pp$canopy_chg), cov_ref = round(cov1, 2), cov_aln = round(cov2, 2),
    gap_ref = round(gap1, 2), gap_aln = round(gap2, 2), water = round(wshare, 2),
    landform = num_or_na(pp$landform), road_flag = num_or_na(pp$road_flag), dist_chan = num_or_na(pp$dist_chan))
  if (k %% 20 == 0) message(sprintf("%d / %d patches (%.0f s)", k, length(pick), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
x <- do.call(rbind, rows)
message(made, " lidar images made, ", length(pick) - made, " kept")

## ---- NAIP (optional) -------------------------------------------------------------------------------
x$naip <- FALSE
naip_on <- is_set(p$naip) && isTRUE(as.logical(p$naip))
naip_status <- NULL
if (naip_on) {
  lab_file <- file.path(out_dir, "sample_labels.txt")
  writeLines(as.character(x$label), lab_file)
  sy <- if (is_set(p$naip_survey_years)) p$naip_survey_years else
    if (all(grepl("^(19|20)[0-9]{2}$", c(rl, al)))) paste(rl, al, sep = ", ") else NULL
  yrs <- if (is_set(p$naip_years)) p$naip_years else "auto"
  if (tolower(yrs) == "auto" && is.null(sy)) stop("naip needs naip_survey_years or naip_years: the survey labels aren't years")
  np <- c(paste0("patch_shp: ", inp$patch_shp), paste0("out_dir: ", img_dir), paste0("labels: ", lab_file),
          paste0("years: ", yrs), if (!is.null(sy)) paste0("survey_years: ", sy),
          if (is_set(p$naip_n_workers)) paste0("n_workers: ", p$naip_n_workers))
  np_file <- file.path(out_dir, "naip_images_params.txt")
  writeLines(c("# naip_images.R parameters, written by scar_review.R", np), np_file)
  message("\n==== NAIP images ====")
  naip_status <- system2("Rscript", shQuote(c(file.path(script_dir, "naip_images.R"), np_file)))
  if (naip_status != 0) warning("naip_images.R exit status ", naip_status, ": some patches have no NAIP image; ",
                                "rerun scar_review.R to retry them (the page hides the missing ones)")
}
x$naip <- file.exists(file.path(img_dir, paste0(x$id, "_naip.jpg")))

## ---- patches.json and the page -------------------------------------------------------------------------
title <- if (is_set(p$title)) p$title else paste0("Scar review, ", rl, "–", al)
subtitle <- if (is_set(p$subtitle)) p$subtitle else
  paste0("Loss candidates from ", basename(dirname(inp$patch_shp)), ". Is each one a landslide scar?")
meta <- list(title = title, subtitle = subtitle, ref_label = rl, align_label = al,
             patch_shp = normalizePath(inp$patch_shp, winslash = "/"), filter = if (is_set(p$filter)) p$filter else NULL,
             per_stratum = per, strata = strata_desc, population = as.list(pop), n_total = nrow(poly),
             naip = naip_on, patches = x)
meta_json <- toJSON(meta, auto_unbox = TRUE, digits = NA, na = "null", null = "null")
writeLines(meta_json, file.path(out_dir, "patches.json"), useBytes = TRUE)

tpl <- file.path(script_dir, "scar_review_page.html")
if (!file.exists(tpl)) stop("Page template not found: ", tpl)
html <- readLines(tpl, warn = FALSE, encoding = "UTF-8")
esc <- function(s) gsub("<", "&lt;", gsub("&", "&amp;", s, fixed = TRUE), fixed = TRUE)
html <- gsub("{{TITLE}}", esc(title), html, fixed = TRUE)
# The sample also goes in the page, so it loads from disk (browsers block fetch() on file://);
# "</" is escaped so a value can't close the script element
html <- gsub("{{PATCHES}}", gsub("</", "<\\/", meta_json, fixed = TRUE), html, fixed = TRUE)
con <- file(file.path(out_dir, "index.html"), open = "w", encoding = "UTF-8"); writeLines(html, con); close(con)

message(sprintf("Scar review: %d patches in %s (index.html, patches.json, img/)%s", nrow(x), out_dir,
                if (naip_on) sprintf("; NAIP for %d", sum(x$naip)) else ""))
message("Publish the folder as a claude.ai artifact with the db and user capabilities to collect labels.")
if (!is.null(naip_status) && naip_status != 0) quit(save = "no", status = 2)
