## patch_polygons.R
##
## LShunter's candidate landslides -- the scars (loss patches) and, with LShunter's DEPOSITS, the
## deposits (gain patches linked to a scar) -- as one polygon shapefile, one feature per patch,
## with everything known about each patch in its attribute table:
##   - type: "scar" or "deposit"; label is the patch's value in its own raster (a scar and a deposit
##     can share a label), and a deposit's scar field is the label of the scar it came from;
##   - LShunter's patch table (out_table) for the scars and its deposit table
##     (<out_table stem>_deposits.csv) for the deposits: points, area, centroid, mean and most
##     extreme centred difference, volume and its SD and z, gradient, most extreme S (or outlier
##     value), road, valley-floor, hollow and inner-gorge shares, distance to channel,
##     false-discovery rates, canopy height in both surveys, seed type; for deposits also the scar,
##     travel distance, drop and reach angle -- whichever columns this LShunter version wrote;
##   - the polygon's own shape: perimeter, and the length and width of its minimum rotated
##     rectangle;
##   - LSstats's per-polygon csvs (<out_csv>_<i>.csv), if given, for the scars: the landform the
##     patch lies in, distance to the floodplain edge, and min / median / max of every raster
##     LSstats summarized (null for deposits, which LSstats doesn't summarize).
## Shapefile field names are at most 10 characters, so <out_shp>_fields.csv lists every field with
## its source column, units and meaning. -9999 (no data) is written as null. A field shared by
## scars and deposits (mean_dz, vol_m3, ...) is signed: negative for loss, positive for gain.
##
## A patch raster's cell values are its table's labels (patches dropped by MINIMUM PATCH
## SIGNIFICANCE have label 0 in the table and aren't in the raster). Patches are dissolved by label,
## so a patch whose cells touch only diagonally becomes one multipart feature.
##
## Skipped if out_shp is newer than the rasters and every table it reads; otherwise rewritten.
##
## Requires: terra, sf
##
##   Rscript patch_polygons.R path/to/params.txt        (see patch_polygons_params_template.txt)

library(terra)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}
config_path <- file.path(get_script_dir(), "patch_polygons_params_template.txt")
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
param_specs <- c("patch_raster", "patch_table", "lsstats_csv", "out_shp", "deposit_raster", "deposit_table")
p <- read_param_file(config_path)
unknown <- setdiff(names(p), param_specs)
if (length(unknown)) stop("Unrecognized parameter(s) in ", config_path, ": ", paste(unknown, collapse = ", "))
is_set <- function(x) !is.null(x) && nzchar(x) && !toupper(x) %in% c("NOFILE", "NONE")
for (k in c("patch_raster", "out_shp")) if (!is_set(p[[k]])) stop("'", k, "' is required in ", config_path)
message("Loaded parameters from: ", config_path)

flt <- paste0(sub("\\.flt$", "", p$patch_raster, ignore.case = TRUE), ".flt")
out_shp <- paste0(sub("\\.shp$", "", p$out_shp, ignore.case = TRUE), ".shp")
fields_csv <- sub("\\.shp$", "_fields.csv", out_shp)
table_csv <- if (is_set(p$patch_table)) paste0(sub("\\.csv$", "", p$patch_table, ignore.case = TRUE), ".csv")
if (!file.exists(flt)) stop("Patch raster not found: ", flt)
## A shapefile writes <stem>.prj, which is also a GridFloat raster's projection file: never share a
## stem with a raster (e.g. out_shp patch_x.shp beside patch_x.flt would replace patch_x.prj)
stem <- sub("\\.shp$", "", out_shp)
if (file.exists(paste0(stem, ".prj")) && !file.exists(out_shp)) {
  stop("out_shp would overwrite ", stem, ".prj, which belongs to another file (",
       paste(basename(setdiff(Sys.glob(paste0(stem, ".*")), paste0(stem, ".prj"))), collapse = ", "),
       "); choose another name, e.g. ", basename(stem), "_polygons.shp")
}
if (!is.null(table_csv) && !file.exists(table_csv)) stop("Patch table not found: ", table_csv)
## Deposits (LShunter DEPOSITS): their raster and table, added as type "deposit"
dep_flt <- if (is_set(p$deposit_raster)) paste0(sub("\\.flt$", "", p$deposit_raster, ignore.case = TRUE), ".flt")
dep_csv <- if (is_set(p$deposit_table)) paste0(sub("\\.csv$", "", p$deposit_table, ignore.case = TRUE), ".csv")
if (!is.null(dep_csv) && is.null(dep_flt)) stop("deposit_table needs deposit_raster in ", config_path)
if (!is.null(dep_flt) && !file.exists(dep_flt)) stop("Deposit raster not found: ", dep_flt)
if (!is.null(dep_csv) && !file.exists(dep_csv)) stop("Deposit table not found: ", dep_csv)

## LSstats's per-polygon csvs: <base>_1.csv, <base>_2.csv, ... (not the _landforms_ ones)
ls_files <- character(0)
if (is_set(p$lsstats_csv)) {
  base <- sub("\\.csv$", "", p$lsstats_csv, ignore.case = TRUE)
  i <- 1
  while (file.exists(f <- paste0(base, "_", i, ".csv"))) { ls_files <- c(ls_files, f); i <- i + 1 }
  if (!length(ls_files)) stop("No LSstats csv files ", base, "_<i>.csv")
}
inputs <- c(flt, table_csv, ls_files, dep_flt, dep_csv)
if (file.exists(out_shp) && file.exists(fields_csv) && all(file.mtime(out_shp) > file.mtime(inputs))) {
  message("Already up to date (delete it to rebuild): ", out_shp)
  quit(save = "no", status = 0)
}

## ---- Polygons ---------------------------------------------------------------------------------
polygonize <- function(f, type) {
  t0 <- Sys.time()
  r <- trim(classify(rast(f), cbind(-Inf, 0.5, NA)))   # labels are 1..n; 0 and nodata are background
  pv <- as.polygons(r, dissolve = TRUE, na.rm = TRUE)
  names(pv) <- "label"
  pv$label <- as.integer(round(pv$label))
  message(sprintf("%d %ss polygonized (%.0f s)", nrow(pv), type,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  pv
}
v <- polygonize(flt, "scar")
n_scar <- nrow(v)
if (!is.null(dep_flt)) {
  vd <- polygonize(dep_flt, "deposit")
  if (!same.crs(v, vd)) stop("The deposit raster's CRS differs from the patch raster's: ", dep_flt)
  v <- rbind(v, vd)
}
type <- rep(c("scar", "deposit"), c(n_scar, nrow(v) - n_scar))
values(v) <- data.frame(type = type, label = v$label)

fields <- data.frame(field = c("type", "label"), source = c("patch rasters", "patch raster"), units = "",
                     description = c("scar (loss patch) or deposit (gain patch linked to a scar)",
                                     "LShunter patch label (the value in the scar or deposit raster)"))
add_field <- function(name, source, units, description) {
  if (name %in% fields$field) {  # a field both tables fill (e.g. mean_dz): note both sources
    i <- match(name, fields$field)
    if (!grepl(source, fields$source[i], fixed = TRUE)) {
      fields$source[i] <<- paste0(fields$source[i], "; ", source)
      if (nzchar(description) && !grepl(description, fields$description[i], fixed = TRUE)) {
        fields$description[i] <<- paste0(fields$description[i], "; ", description)
      }
    }
    return(invisible())
  }
  fields <<- rbind(fields, data.frame(field = name, source = source, units = units, description = description))
}

## Shape: perimeter, and length and width of the minimum rotated rectangle
g <- sf::st_as_sf(v)
add_field("perim_m", "polygon", "m", "perimeter")
v$perim_m <- round(perim(v), 1)
mrr <- sf::st_minimum_rotated_rectangle(sf::st_geometry(g))
sides <- t(vapply(mrr, function(pg) {
  xy <- sf::st_coordinates(pg)[1:3, 1:2]
  d <- sqrt(rowSums(diff(xy)^2))
  c(max(d), min(d))
}, numeric(2)))
v$length_m <- round(sides[, 1], 1)
v$width_m <- round(sides[, 2], 1)
add_field("length_m", "polygon", "m", "length of the minimum rotated rectangle (longer side)")
add_field("width_m", "polygon", "m", "width of the minimum rotated rectangle (shorter side)")

na9999 <- function(x) { if (is.numeric(x)) x[!is.na(x) & x <= -9998] <- NA; x }

## ---- LShunter's patch table -------------------------------------------------------------------
## Known columns: short name, units, meaning. Columns not listed keep a name cut to 10 characters.
known <- list(
  candidate = c("cand", "", "candidate number before patches were dropped by significance"),
  points = c("points", "cells", "cells in the patch"),
  area_m2 = c("area_m2", "m2", "area (cells x cell area)"),
  x = c("cent_x", "m", "centroid x"),
  y = c("cent_y", "m", "centroid y"),
  mean_depth_m = c("mean_dz", "m", "mean centred difference (negative: loss)"),
  min_depth_m = c("min_dz", "m", "most negative centred difference (deepest point)"),
  volume_m3 = c("vol_m3", "m3", "volume = area x mean centred difference (negative: loss, positive: gain)"),
  sd_volume_m3 = c("sd_vol_m3", "m3", "SD of the volume from spatially correlated noise"),
  z_volume = c("z_vol", "", "volume / SD of volume"),
  n_eff = c("n_eff", "", "effective number of independent points"),
  mean_grad_tan = c("grad_mean", "tan", "mean gradient"),
  max_grad_tan = c("grad_max", "tan", "maximum gradient"),
  min_S = c("min_S", "", "most extreme exceedance S (negative: loss; S = -2 is p = 1%)"),
  min_outlier = c("min_outl", "", "most extreme outlier value"),
  road_frac = c("road_frac", "", "share of cells in the road buffer"),
  valley_frac = c("valley_frc", "", "share of cells on valley floor (RIL)"),
  hollow_frac = c("hollow_frc", "", "share of cells in bedrock hollows (RIL)"),
  gorge_frac = c("gorge_frc", "", "share of cells in inner gorges (RIL)"),
  min_dist_channel_m = c("dist_chan", "m", "minimum distance to a channel"),
  n_at_least = c("n_atleast", "", "null patches at least this extreme"),
  sim_expected = c("sim_exp", "", "expected count at least this extreme in the null simulation"),
  sim_fdr = c("sim_fdr", "", "false-discovery rate from the null simulation"),
  flip_count = c("flip_cnt", "", "gain patches at least this extreme (sign-flip control)"),
  flip_fdr = c("flip_fdr", "", "false-discovery rate from the sign-flip control"),
  ref_canopy_m = c("canopy_ref", "m", "mean canopy height, reference survey"),
  align_canopy_m = c("canopy_aln", "m", "mean canopy height, aligned survey"),
  canopy_change_m = c("canopy_chg", "m", "canopy height change (aligned - reference)"),
  seed_type = c("seed_type", "", "1: has a threshold-1 seed zone; 2: area seeds only"),
  min_area_z = c("min_area_z", "", "lowest area-evidence z (area seeds; positive for gain patches)"),
  # Deposit and unlinked gain-patch tables (LShunter 1.4.0 DEPOSITS)
  mean_thickness_m = c("mean_dz", "m", "mean centred difference (positive: gain)"),
  max_thickness_m = c("max_dz", "m", "largest centred difference (thickest point)"),
  max_S = c("max_S", "", "most extreme exceedance S (positive: gain; S = 2 is p = 1%)"),
  channel_frac = c("chan_frac", "", "share of cells in channels (accumulation >= the deposit seed limit)"),
  scar_label = c("scar", "", "label of the scar (loss patch) the deposit is linked to"),
  link_points = c("link_pts", "cells", "seed cells in the scar's runout corridor that link it"),
  travel_m = c("travel_m", "m", "least runout-path distance of the linking cells beyond the last coherent loss"),
  drop_m = c("drop_m", "m", "elevation drop, scar's highest point to the deposit's lowest (reference DTM)"),
  reach_m = c("reach_m", "m", "horizontal distance, scar's highest point to the deposit's lowest"),
  reach_angle_deg = c("reach_ang", "degrees", "reach angle atan(drop / reach)"))
## One table's columns onto the features of one type (rows = that type's features, matched by
## label); other features get null there
add_table <- function(csv, what, rows) {
  tb <- utils::read.csv(csv, check.names = FALSE, strip.white = TRUE)
  names(tb) <- trimws(names(tb))
  tb <- tb[tb$label > 0, , drop = FALSE]
  labs <- v$label[rows]
  missing <- setdiff(labs, tb$label)
  if (length(missing)) stop(length(missing), " ", what, " label(s) in the raster are not in its table (",
                            csv, "), e.g. ", missing[1], ": they don't belong to the same LShunter run.")
  if (nrow(tb) != length(labs)) warning(nrow(tb) - length(labs), " ", what, "-table rows have no cells in the raster.")
  tb <- tb[match(labs, tb$label), setdiff(names(tb), "label"), drop = FALSE]
  own <- character(0)  # this table's field names, so two of its columns never share one
  for (col in names(tb)) {
    k <- known[[col]]
    nm <- if (is.null(k)) substr(gsub("[^A-Za-z0-9_]", "_", col), 1, 10) else k[1]
    if (nm %in% c("type", "label", own) || (is.null(k) && nm %in% fields$field)) {
      nm <- substr(paste0(substr(nm, 1, 8), "_", length(fields$field)), 1, 10)
    }
    own <- c(own, nm)
    vals <- na9999(tb[[col]])
    att <- values(v)
    if (!nm %in% names(att)) att[[nm]] <- if (is.numeric(vals)) NA_real_ else NA_character_
    att[[nm]][rows] <- vals
    values(v) <<- att
    add_field(nm, paste0(what, " table: ", col), if (is.null(k)) "" else k[2], if (is.null(k)) "" else k[3])
  }
}
if (!is.null(table_csv)) add_table(table_csv, "patch", which(v$type == "scar"))
if (!is.null(dep_csv)) add_table(dep_csv, "deposit", which(v$type == "deposit"))

## ---- LSstats per-polygon csvs -----------------------------------------------------------------
## Each starts with "Raster file: <path>", then one row per polygon ID (= label): landform, area,
## volume (already in the patch table), quantiles of that raster, and distance to the floodplain.
## Field prefix and units for a raster, from its name (Align's outputs and the gradient raster)
short_of <- function(path) {
  b <- tolower(basename(sub("\\.flt$", "", path, ignore.case = TRUE)))
  pre <- list(cendif = c("cdif", "m"), exceed = c("S", ""), grad = c("grd", "tan"), outlier = c("outl", ""),
              lod = c("lod", "m"), dif = c("dif", "m"), elev = c("elev", "m"))
  hit <- names(pre)[vapply(names(pre), function(k) startsWith(b, k) || grepl(paste0("_", k), b), TRUE)]
  if (length(hit)) pre[[hit[1]]] else c(substr(gsub("[^a-z0-9]", "", b), 1, 5), "")
}
for (i in seq_along(ls_files)) {
  lines <- readLines(ls_files[i], warn = FALSE)
  rpath <- trimws(sub("^\\s*Raster file:\\s*", "", lines[1]))
  d <- utils::read.csv(text = lines[-1], check.names = FALSE, strip.white = TRUE)
  names(d) <- trimws(names(d))
  # LSstats summarizes the scars only: deposits match no row and get null
  d <- d[match(ifelse(v$type == "scar", v$label, NA), d$ID), , drop = FALSE]
  if (i == 1) {
    v$landform <- trimws(d[["Landform Name"]])
    v$lf_code <- d[["Landform Code"]]
    add_field("landform", "LSstats", "", "landform (RIL class) the patch lies in")
    add_field("lf_code", "LSstats", "", "landform code (1-10)")
    if ("DistTo" %in% names(d)) {
      dist <- na9999(d$DistTo)
      dist[!is.na(dist) & dist >= 1e5] <- NA     # LSstats's "beyond the search radius"
      v$dist_fp <- dist
      add_field("dist_fp", "LSstats", "m", "distance to the floodplain edge (null beyond the search radius)")
    }
  }
  pu <- short_of(rpath)
  pre <- pu[1]
  if (any(startsWith(fields$field, paste0(pre, "_")))) pre <- paste0(substr(pre, 1, 3), i)
  stat <- c(Min = "minimum", Q50 = "median", Max = "maximum")
  for (col in names(stat)) {
    nm <- paste0(pre, "_", tolower(col))
    v[[nm]] <- na9999(d[[col]])
    add_field(nm, paste0("LSstats: ", basename(rpath)), pu[2], paste0(stat[[col]], " of ", rpath, " over the patch"))
  }
}

## ---- Write --------------------------------------------------------------------------------------
dir.create(dirname(out_shp), recursive = TRUE, showWarnings = FALSE)
tmp <- sub("\\.shp$", "_partial.shp", out_shp)
writeVector(v, tmp, filetype = "ESRI Shapefile", overwrite = TRUE)
for (ext in c(".shp", ".shx", ".dbf", ".prj", ".cpg")) {
  from <- paste0(sub("\\.shp$", "", tmp), ext)
  if (file.exists(from)) file.rename(from, paste0(stem, ext))
}
utils::write.csv(fields, fields_csv, row.names = FALSE)
message(sprintf("Wrote %s: %d scars and %d deposits, %d fields (described in %s); total area %.0f m2",
                out_shp, sum(v$type == "scar"), sum(v$type == "deposit"), ncol(v), basename(fields_csv),
                sum(expanse(v, transform = FALSE))))
