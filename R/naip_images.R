## naip_images.R
##
## NAIP aerial photos of each candidate landslide, one jpg per patch: for every chosen NAIP year a
## close-up (the patch, as the scar-review images frame it) and a wide view (wide_m across), the
## patch outlined in magenta and its neighbours dashed yellow, the flight date in each title.
## Leaf-on 0.6-1 m imagery shows canopy, not the ground under it: it helps most where the canopy
## opened or was never closed (harvest, roads, fresh bare soil, debris tracks).
##
## Imagery: Microsoft Planetary Computer's NAIP collection (cloud-optimized GeoTIFFs, no account),
## read window by window over the network -- nothing is downloaded whole. The script searches the
## catalog over the patches' extent itself and writes the tiles it found to <out_dir>/naip_items.csv.
##
## years: auto (the default) picks, for each survey year in survey_years, the latest NAIP year
## flown in or before it (or the earliest after it, if none was). E.g. lidar 2017 and 2025 with
## NAIP 2017 and 2023 available: 2017 and 2023 -- change after the 2023 flight won't show.
##
## Outputs, in out_dir:
##   p<label>_naip.jpg   one per patch (2 x years panels, panel_px square each)
##   naip_images.csv     label, years, flight dates per panel, status (ok / failed: reason)
##   naip_items.csv      the NAIP tiles found (year, date, CRS, footprint, url)
## Resumable: an image newer than patch_shp is skipped; a patch whose imagery failed to read (the
## network, an expired token) gets no image, so a rerun retries it. Exits with status 2 if any
## patch failed.
##
## Requires: terra, sf, jsonlite, future, future.apply; internet access.
##
##   Rscript naip_images.R path/to/params.txt        (see naip_images_params_template.txt)

suppressMessages({ library(terra); library(sf); library(jsonlite) })
terraOptions(progress = 0)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}
config_path <- file.path(get_script_dir(), "naip_images_params_template.txt")
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
defaults <- list(patch_shp = NULL, out_dir = NULL, years = "auto", survey_years = NULL, labels = "ALL",
                 near_min_half_m = "45", wide_m = "350", panel_px = "280", n_workers = "4",
                 stac_url = "https://planetarycomputer.microsoft.com/api/stac/v1",
                 collection = "naip",
                 token_url = "https://planetarycomputer.microsoft.com/api/sas/v1/token/naip")
p <- read_param_file(config_path)
unknown <- setdiff(names(p), names(defaults))
if (length(unknown)) stop("Unrecognized parameter(s) in ", config_path, ": ", paste(unknown, collapse = ", "))
is_set <- function(x) !is.null(x) && nzchar(x) && !toupper(x) %in% c("NOFILE", "NONE")
for (k in names(defaults)) if (!is_set(p[[k]]) && !is.null(defaults[[k]])) p[[k]] <- defaults[[k]]
for (k in c("patch_shp", "out_dir")) if (!is_set(p[[k]])) stop("'", k, "' is required in ", config_path)
num <- function(k) { v <- suppressWarnings(as.numeric(p[[k]])); if (is.na(v) || v <= 0) stop(k, " must be a positive number"); v }
near_min_half <- num("near_min_half_m"); wide_half <- num("wide_m") / 2
panel_px <- as.integer(num("panel_px")); n_workers <- as.integer(num("n_workers"))
split_list <- function(x) { x <- trimws(strsplit(x, ",", fixed = TRUE)[[1]]); x[nzchar(x)] }
message("Loaded parameters from: ", config_path)

shp <- paste0(sub("\\.shp$", "", p$patch_shp, ignore.case = TRUE), ".shp")
if (!file.exists(shp)) stop("patch_shp not found: ", shp)
out_dir <- p$out_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Patches -----------------------------------------------------------------------------------
poly <- st_read(shp, quiet = TRUE)
if (!"label" %in% names(poly)) stop(shp, " has no 'label' field (expected patch_polygons.R's shapefile)")
# With LShunter's DEPOSITS the shapefile also holds deposits (type "deposit"), whose labels can
# repeat the scars' (each numbers its own raster) -- images are named by label, so scars only
if ("type" %in% names(poly)) poly <- poly[!is.na(poly$type) & poly$type == "scar", ]
if (is.na(st_crs(poly))) stop(shp, " has no CRS (.prj)")
sel <- if (toupper(p$labels) == "ALL") {
  poly$label
} else if (file.exists(p$labels)) {
  scan(p$labels, quiet = TRUE)
} else as.numeric(split_list(p$labels))
bad <- setdiff(sel, poly$label)
if (length(bad)) stop(length(bad), " label(s) not in ", shp, ", e.g. ", bad[1])
sel <- unique(as.integer(sel))
message(length(sel), " of ", nrow(poly), " patches selected")

## ---- NAIP catalog over the selected patches ------------------------------------------------------
bb <- st_bbox(st_transform(st_as_sfc(st_bbox(st_buffer(poly[poly$label %in% sel, ], wide_half * 1.5))), 4326))
search <- paste0(p$stac_url, "/search?collections=", p$collection, "&limit=1000&bbox=",
                 paste(sprintf("%.6f", bb[c("xmin", "ymin", "xmax", "ymax")]), collapse = ","))
feats <- list(); url <- search
while (!is.null(url)) {
  pg <- fromJSON(url, simplifyVector = FALSE)
  feats <- c(feats, pg$features)
  nxt <- Filter(function(l) identical(l$rel, "next"), pg$links)
  url <- if (length(nxt)) nxt[[1]]$href else NULL
}
if (!length(feats)) stop("No NAIP imagery found over ", shp, " (", search, ")")
items <- do.call(rbind, lapply(feats, function(f) {
  pr <- f$properties; pb <- unlist(pr[["proj:bbox"]])
  epsg <- if (!is.null(pr[["proj:epsg"]])) pr[["proj:epsg"]] else as.integer(sub("^EPSG:", "", pr[["proj:code"]]))
  data.frame(year = as.integer(if (!is.null(pr[["naip:year"]])) pr[["naip:year"]] else substr(pr$datetime, 1, 4)),
             date = substr(pr$datetime, 1, 10), gsd = if (is.null(pr$gsd)) NA else pr$gsd, epsg = epsg,
             xmin = pb[1], ymin = pb[2], xmax = pb[3], ymax = pb[4], href = f$assets$image$href)
}))
items <- items[order(items$year, items$date), ]
utils::write.csv(items, file.path(out_dir, "naip_items.csv"), row.names = FALSE)
avail <- sort(unique(items$year))
message("NAIP years over the patches: ", paste(avail, collapse = ", "))

if (tolower(p$years) == "auto") {
  if (!is_set(p$survey_years)) stop("years: auto needs survey_years (e.g. 2017, 2025)")
  sy <- suppressWarnings(as.integer(split_list(p$survey_years)))
  if (anyNA(sy)) stop("survey_years must be years, e.g. 2017, 2025: ", p$survey_years)
  years <- vapply(sy, function(y) if (any(avail <= y)) max(avail[avail <= y]) else min(avail), integer(1))
  for (i in seq_along(sy)) message("Survey ", sy[i], ": NAIP ", years[i])
  if (anyDuplicated(years)) message("Note: surveys share a NAIP year; showing it once")
  years <- unique(years)
} else {
  years <- as.integer(split_list(p$years))
  miss <- setdiff(years, avail)
  if (length(miss)) stop("No NAIP for year(s) ", paste(miss, collapse = ", "), " here; available: ", paste(avail, collapse = ", "))
}
items <- items[items$year %in% years, ]

## ---- Work list (skip images newer than the shapefile) ----------------------------------------------
img_of <- function(lab) file.path(out_dir, sprintf("p%d_naip.jpg", lab))
done <- vapply(sel, function(l) { f <- img_of(l); file.exists(f) && file.mtime(f) > file.mtime(shp) }, TRUE)
todo <- sel[!done]
message(sum(done), " images already up to date; ", length(todo), " to make")

## ---- One patch -----------------------------------------------------------------------------------
make_chunk <- function(labs) {
  suppressMessages({ library(terra); library(sf); library(jsonlite) }); terraOptions(progress = 0)
  setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
  setGDALconfig("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
  setGDALconfig("GDAL_HTTP_MAX_RETRY", "5"); setGDALconfig("GDAL_HTTP_RETRY_DELAY", "2")
  tk <- new.env(); tk$time <- as.POSIXct(0); tk$tok <- ""
  token <- function(force = FALSE) {   # Planetary Computer SAS tokens last ~15 minutes
    if (force || difftime(Sys.time(), tk$time, units = "mins") > 8) {
      tk$tok <- fromJSON(p$token_url)$token; tk$time <- Sys.time() }
    tk$tok }
  read_win <- function(href, w) {
    for (a in 1:4) {
      r <- tryCatch({
        x <- rast(paste0("/vsicurl/", href, "?", token(a > 1)))[[1:3]]
        e <- intersect(ext(x), w)
        if (is.null(e) || (e$xmax - e$xmin) < 2 || (e$ymax - e$ymin) < 2) return(NULL)
        crop(x, e) * 1 }, error = function(err) err)   # * 1 reads the window now, inside the retry
      if (!inherits(r, "error")) return(r)
      Sys.sleep(3 * a)
    }
    stop(conditionMessage(r))
  }
  # Square window around (cx, cy) in the patch CRS, as an extent in EPSG:epsg
  win <- function(cx, cy, half, epsg) {
    sq <- st_as_sfc(st_bbox(c(xmin = cx - half, ymin = cy - half, xmax = cx + half, ymax = cy + half), crs = st_crs(poly)))
    b <- st_bbox(st_transform(sq, epsg)); ext(b$xmin, b$xmax, b$ymin, b$ymax)
  }
  # RGB for one year over a window: tiles touching it, in the CRS of the tile at the centre
  naip <- function(yr, cx, cy, half) {
    it <- items[items$year == yr, ]
    hit <- NULL
    for (e in unique(it$epsg)) {
      w <- win(cx, cy, half, e)
      s <- it[it$epsg == e & it$xmin < w$xmax & it$xmax > w$xmin & it$ymin < w$ymax & it$ymax > w$ymin, ]
      if (nrow(s)) hit <- rbind(hit, s)
    }
    if (is.null(hit)) return(NULL)
    epsg <- as.integer(names(which.max(table(hit$epsg))))
    w <- win(cx, cy, half, epsg)
    r <- list()
    for (i in seq_len(nrow(hit))) {
      wi <- if (hit$epsg[i] == epsg) w else win(cx, cy, half * 1.1, hit$epsg[i])
      x <- read_win(hit$href[i], wi)
      if (is.null(x)) next
      if (hit$epsg[i] != epsg) x <- project(x, paste0("EPSG:", epsg))
      r[[length(r) + 1]] <- x
    }
    if (!length(r)) return(NULL)
    x <- if (length(r) == 1) r[[1]] else do.call(merge, r)
    list(r = x, w = w, epsg = epsg, dates = paste(sort(unique(hit$date)), collapse = ", "))
  }

  rows <- list()
  for (lab in labs) {
    me <- poly[poly$label == lab, ]; bx <- st_bbox(me)
    cx <- unname((bx$xmin + bx$xmax) / 2); cy <- unname((bx$ymin + bx$ymax) / 2)
    half <- max(near_min_half, 0.75 * max(bx$xmax - bx$xmin, bx$ymax - bx$ymin) + 20)   # as the scar-review images
    panels <- list(); err <- NULL
    for (sc in c("near", "wide")) for (yr in years) {
      h <- if (sc == "near") half else wide_half
      pn <- tryCatch(naip(yr, cx, cy, h), error = function(e) { err <<- conditionMessage(e); NULL })
      panels[[length(panels) + 1]] <- list(yr = yr, sc = sc, half = h, img = pn)
    }
    if (!is.null(err)) {
      rows[[length(rows) + 1]] <- data.frame(label = lab, status = paste("failed:", substr(err, 1, 200)))
      next
    }
    f <- img_of(lab); tmp <- sub("\\.jpg$", "_partial.jpg", f)
    jpeg(tmp, width = panel_px * length(panels), height = panel_px + 20, quality = 82)
    dev <- dev.cur()
    drawn <- tryCatch({ draw(lab, me, panels); TRUE }, error = function(e) { err <<- conditionMessage(e); FALSE })
    dev.off(dev)
    if (!drawn) {
      unlink(tmp)
      rows[[length(rows) + 1]] <- data.frame(label = lab, status = paste("failed: drawing:", substr(err, 1, 200)))
      next
    }
    file.rename(tmp, f)
    d <- vapply(panels, function(pn) if (is.null(pn$img)) "none" else pn$img$dates, "")
    rows[[length(rows) + 1]] <- data.frame(label = lab, status = "ok",
      t(setNames(d, paste0(vapply(panels, `[[`, "", "sc"), "_", vapply(panels, function(pn) as.character(pn$yr), "")))))
  }
  do.call(rbind, lapply(rows, function(r) { r[setdiff(c("label", "status"), names(r))] <- NA; r }))
}

draw <- function(lab, me, panels) {
    par(mfrow = c(1, length(panels)), mar = c(0.3, 0.3, 1.6, 0.3), bg = "white")
    for (pn in panels) {
      ttl <- paste0("NAIP ", pn$yr, if (pn$sc == "wide") paste0(", ", round(2 * pn$half), " m"))
      if (is.null(pn$img)) { plot.new(); title(paste(ttl, "- no imagery")); next }
      plotRGB(pn$img$r, stretch = "lin", ext = pn$img$w, axes = FALSE, mar = c(0.3, 0.3, 1.6, 0.3))
      title(paste0(ttl, " (", pn$img$dates, ")"), cex.main = 1.05, line = 0.3)
      pv <- vect(st_transform(st_geometry(me), pn$img$epsg))
      if (pn$sc == "near") {
        oth <- st_transform(st_geometry(poly[poly$label != lab, ]), pn$img$epsg)
        e <- as.vector(pn$img$w)   # xmin, xmax, ymin, ymax
        box <- st_as_sfc(st_bbox(c(xmin = e[[1]], ymin = e[[3]], xmax = e[[2]], ymax = e[[4]]), crs = st_crs(pn$img$epsg)))
        oth <- oth[st_intersects(oth, box, sparse = FALSE)[, 1]]
        if (length(oth)) lines(vect(oth), col = "#ffd400", lwd = 1.2, lty = 2)
        lines(pv, col = "#ff00c8", lwd = 1.8)
        sbar(20, xy = "bottomleft", type = "bar", divs = 2, below = "m", cex = 1, col = "white")
      } else {
        lines(pv, col = "#ff00c8", lwd = 1.5)
        sbar(100, xy = "bottomleft", type = "bar", divs = 2, below = "m", cex = 1, col = "white")
      }
    }
}

## ---- Run -------------------------------------------------------------------------------------------
t0 <- Sys.time()
if (length(todo)) {
  chunks <- split(todo, ceiling(seq_along(todo) / 10))
  if (n_workers > 1 && length(chunks) > 1) {
    future::plan(future::multisession, workers = min(n_workers, length(chunks)))
    res <- future.apply::future_lapply(chunks, make_chunk, future.seed = NULL)
    future::plan(future::sequential)
  } else {
    res <- list()
    for (i in seq_along(chunks)) {
      res[[i]] <- make_chunk(chunks[[i]])
      message(sprintf("%d / %d patches", min(i * 10, length(todo)), length(todo)))
    }
  }
  cols <- unique(unlist(lapply(res, names)))
  res <- do.call(rbind, lapply(res, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
} else res <- NULL

## Index: this run's rows replace earlier rows for the same labels
idx_f <- file.path(out_dir, "naip_images.csv")
if (file.exists(idx_f) && !is.null(res)) {
  old <- utils::read.csv(idx_f, check.names = FALSE)
  old <- old[!old$label %in% res$label, , drop = FALSE]
  if (nrow(old)) {
    cols <- union(names(old), names(res))
    old[setdiff(cols, names(old))] <- NA; res[setdiff(cols, names(res))] <- NA
    res <- rbind(old[cols], res[cols])
  }
}
if (!is.null(res)) utils::write.csv(res[order(res$label), ], idx_f, row.names = FALSE)

nfail <- if (is.null(res)) 0 else sum(res$label %in% todo & startsWith(res$status, "failed"))
message(sprintf("NAIP images: %d made, %d failed, %d already up to date, in %s (%.0f min); years %s",
                length(todo) - nfail, nfail, sum(done), out_dir,
                as.numeric(difftime(Sys.time(), t0, units = "mins")), paste(years, collapse = ", ")))
if (nfail) {
  message("Failed patches have no image; rerun to retry them. Reasons in ", idx_f)
  quit(save = "no", status = 2)
}
