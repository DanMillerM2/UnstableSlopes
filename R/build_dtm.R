## build_dtm.R
##
## Builds DTM and DSM tiles from a folder of lidar LAZ/LAS tiles with one interpolation scheme, so
## two surveys built with the same settings can be differenced like for like (instead of using
## each vendor's own DTM, which may come from a different ground classification or interpolation
## than the points delivered with it). Writes GeoTIFF tiles in the point cloud's own CRS; mosaic_dtm.R
## then mosaics them (as for downloaded tiles), converts to meters and relabels or reprojects to
## the target CRS -- point a site's ref_dtm_tile_dir / ref_dsm_tile_dir (and align_...) at
## <output_dir>/dtm and <output_dir>/dsm.
##
## For each tile (one future worker each):
##   1. The tile's extent, snapped outward to the res grid (cells are multiples of res in the CRS's
##      units, so all tiles -- and both surveys, at the same res -- share one grid), is read with a
##      buffer of buffer_m from every tile overlapping it (lasR reader_rectangles: points in the
##      rectangle, each read once even where tiles' boxes overlap).
##   2. DTM: a Delaunay triangulation of the ground returns (ground_classes; withheld points
##      dropped), rasterized by linear interpolation at the cell centres. Triangles with an edge
##      longer than max_edge_m are dropped, so where ground returns are too sparse (dense canopy,
##      water) the DTM is nodata rather than invented.
##   3. DSM: the highest return in each cell, noise classes (dsm_drop_classes) left out; cells
##      without a return are nodata.
##   4. Both are cropped back to the tile's extent and written to <output_dir>/dtm/<tile>.tif and
##      <output_dir>/dsm/<tile>.tif, under a _partial name renamed into place when complete.
## A tile whose two outputs exist is skipped, so an interrupted run resumes. Delete the outputs
## after changing a parameter.
##
## Requires: lasR, terra, lidR (LAS headers only), sf, future, future.apply
##
##   Rscript build_dtm.R path/to/params.txt        (see build_dtm_params_template.txt)

library(terra)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}
config_path <- file.path(get_script_dir(), "build_dtm_params_template.txt")
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
  input_dir  = list(type = "character", required = TRUE),   # folder of LAZ/LAS tiles
  output_dir = list(type = "character", required = TRUE),   # dtm/ and dsm/ subfolders are written here
  # Cell size in meters (converted to the point cloud's horizontal units). Build both surveys at
  # the same res: their grids then coincide.
  res_m      = list(type = "numeric", default = 1.0),
  # Longest triangle edge kept in the ground TIN (m); longer triangles leave nodata.
  max_edge_m = list(type = "numeric", default = 10),
  # Buffer read around each tile (m), so the TIN and the DSM are seamless across tile edges.
  buffer_m   = list(type = "numeric", default = 20),
  # LAS classes triangulated for the DTM, and left out of the DSM (space- or comma-separated).
  ground_classes   = list(type = "character", default = "2"),
  dsm_drop_classes = list(type = "character", default = "7 18"),
  build_dsm  = list(type = "logical", default = TRUE),
  # Only tiles whose LAS-header box intersects this polygon (any CRS) are built; NOFILE builds all.
  # Neighbouring tiles outside it are still read for the buffer.
  clip_polygon = list(type = "character", default = "NOFILE"),
  file_pattern = list(type = "character", default = "\\.la[sz]$"),
  n_workers    = list(type = "integer", default = 4),
  threads_per_worker = list(type = "integer", default = 2)
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
classes <- function(s) as.integer(strsplit(trimws(s), "[ ,]+")[[1]])
ground_cls <- classes(p$ground_classes)
drop_cls <- classes(p$dsm_drop_classes)

## ---- Tiles and their CRS --------------------------------------------------------------------
files <- list.files(p$input_dir, pattern = p$file_pattern, full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No files matching ", p$file_pattern, " in ", p$input_dir)
ctg <- lidR::readLAScatalog(files, progress = FALSE)
hd <- ctg@data
tiles <- data.frame(file = normalizePath(hd$filename, winslash = "/"), xmin = hd$Min.X, xmax = hd$Max.X,
                    ymin = hd$Min.Y, ymax = hd$Max.Y)
crs_las <- sf::st_crs(ctg)
unit_m <- c(`US survey foot` = 1200 / 3937, foot = 0.3048, metre = 1)[crs_las$units_gdal]
if (is.na(unit_m)) stop("Unrecognized horizontal unit of the point cloud's CRS: ", crs_las$units_gdal)
res <- p$res_m / unit_m
buf <- p$buffer_m / unit_m
max_edge <- p$max_edge_m / unit_m
message(sprintf("%d tiles; CRS unit %s: res %.6f, buffer %.2f, max edge %.2f (CRS units)",
                nrow(tiles), crs_las$units_gdal, res, buf, max_edge))

build <- tiles
if (is_set(p$clip_polygon)) {
  poly <- sf::st_transform(sf::st_union(sf::st_read(p$clip_polygon, quiet = TRUE)), crs_las)
  boxes <- do.call(c, lapply(seq_len(nrow(tiles)), function(i)
    sf::st_as_sfc(sf::st_bbox(c(xmin = tiles$xmin[i], ymin = tiles$ymin[i], xmax = tiles$xmax[i], ymax = tiles$ymax[i]), crs = crs_las))))
  build <- tiles[lengths(sf::st_intersects(boxes, poly)) > 0, ]
  message(nrow(build), " of ", nrow(tiles), " tiles intersect ", p$clip_polygon)
}
dir.create(file.path(p$output_dir, "dtm"), recursive = TRUE, showWarnings = FALSE)
if (p$build_dsm) dir.create(file.path(p$output_dir, "dsm"), recursive = TRUE, showWarnings = FALSE)
out_of <- function(f, kind) file.path(p$output_dir, kind, sub("\\.la[sz]$", ".tif", basename(f), ignore.case = TRUE))
done <- file.exists(out_of(build$file, "dtm")) & (!p$build_dsm | file.exists(out_of(build$file, "dsm")))
todo <- build[!done, ]
message(sum(done), " tiles already built; ", nrow(todo), " to build")

## ---- One tile -------------------------------------------------------------------------------
build_tile <- function(i) {
  suppressMessages({library(lasR); library(terra)})
  t <- todo[i, ]
  core <- c(floor(t$xmin / res) * res, ceiling(t$xmax / res) * res, floor(t$ymin / res) * res, ceiling(t$ymax / res) * res)
  rect <- core + c(-buf, buf, -buf, buf)
  src <- tiles$file[tiles$xmax >= rect[1] & tiles$xmin <= rect[2] & tiles$ymax >= rect[3] & tiles$ymin <= rect[4]]
  tmp <- file.path(p$output_dir, paste0("tmp_", sub("\\.la[sz]$", "", basename(t$file), ignore.case = TRUE)))
  tri <- triangulate(max_edge, filter = keep_class(ground_cls))
  pipe <- reader_rectangles(rect[1], rect[3], rect[2], rect[4], filter = "Withheld == 0") + tri +
    lasR::rasterize(res, tri, ofile = paste0(tmp, "_dtm.tif"))
  if (p$build_dsm) pipe <- pipe + lasR::rasterize(res, "zmax", filter = drop_class(drop_cls), ofile = paste0(tmp, "_dsm.tif"))
  exec(pipe, on = src, ncores = p$threads_per_worker)
  for (kind in c("dtm", if (p$build_dsm) "dsm")) {
    r <- crop(rast(paste0(tmp, "_", kind, ".tif")), ext(core), snap = "near")
    names(r) <- kind
    fin <- out_of(t$file, kind); part <- sub("\\.tif$", "_partial.tif", fin)
    writeRaster(r, part, overwrite = TRUE, datatype = "FLT4S", NAflag = -9999, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    rm(r); gc()
    file.remove(paste0(tmp, "_", kind, ".tif"))
    if (!file.rename(part, fin)) stop("Could not rename ", part, " to ", fin)
  }
  basename(t$file)
}

## ---- Run ------------------------------------------------------------------------------------
t0 <- Sys.time()
if (nrow(todo)) {
  future::plan(future::multisession, workers = min(p$n_workers, nrow(todo)))
  shared <- c("todo", "tiles", "res", "buf", "max_edge", "ground_cls", "drop_cls", "p", "out_of")
  res_list <- future.apply::future_lapply(seq_len(nrow(todo)), build_tile, future.seed = TRUE, future.globals = shared)
  future::plan(future::sequential)
}
message(sprintf("Built %d tiles in %.1f min; outputs in %s", nrow(todo), as.numeric(difftime(Sys.time(), t0, units = "mins")),
                normalizePath(p$output_dir, winslash = "/")))
