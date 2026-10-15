## mosaic_dtm.R
##
## Mosaics a folder of DTM tiles (GeoTIFFs, e.g. as downloaded from a lidar portal) into one
## contiguous DTM covering an input polygon, converts elevations to meters, and reprojects it to
## target_crs if the tiles aren't already in it. The result is written as an ESRI GridFloat
## (.flt) raster -- the format the Fortran programs downstream (Align, bldgrds, ...) read.
## Nothing here is DTM-specific: DSM (surface-model) tiles are mosaicked the same way, and
## run_align_pipeline.R uses it for its [ref_dsm]/[align_dsm] sections too. With clip_polygon
## NOFILE, every tile in tile_dir is mosaicked, unmasked.
##
## 1. List the tiles in tile_dir matching file_pattern and read each one's extent/CRS (header
##    metadata only -- no pixel values are read at this point)
## 2. Keep only the tiles that intersect clip_polygon (reprojected to the tiles' own CRS)
## 3. Mosaic the kept tiles, then crop to the polygon plus a margin of clip_buffer_cells cells,
##    so the reprojection in step 5 has real data to interpolate from right up to the polygon's
##    edge instead of leaving a rim of NA cells along it
## 4. Convert elevations to meters (see z_units below)
## 5. If the tiles' CRS differs from target_crs, reproject (bilinear) at the tiles' own pixel
##    size expressed in target_crs units -- e.g. 3 US-survey-foot tiles come out at
##    3 * 1200/3937 = 0.9144018 m pixels -- so the output keeps the input's resolution
## 6. Mask to the exact polygon and write <output_dtm>.flt
##
## Resumability: if <output_dtm>.flt already exists, the script does nothing -- same
## exists-then-skip convention used throughout this repo. Delete it to rebuild (e.g. after
## changing clip_polygon or any other parameter below).
##
## Requires: terra
##
## All inputs are read from a plain-text "keyword: value" parameter file (same rules as
## ground_returns.R's and align_dtms.R's own parameter files) -- see mosaic_dtm_params_template.txt
## for a ready-to-copy example. Pass its path as a command-line argument:
##   Rscript mosaic_dtm.R path/to/your_params.txt
## or, when sourcing from R/RStudio, set config_path below before running the script.
## run_align_pipeline.R runs this script once each for its [ref_dtm] and [align_dtm] sections.

library(terra)

## ---- Parameter file ------------------------------------------------------------------------

# Resolves the folder this script itself lives in (from Rscript's "--file=" argument), so the
# fallback config_path below is found by this script's location on disk rather than by the
# caller's working directory. Falls back to "." (assume cwd) when there's no "--file="
# argument to read, e.g. when this script is source()'d from RStudio instead of run via Rscript.
get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}

#config_path <- "path/to/params.txt"   # template -- point this at your own project's parameter file
config_path <- file.path(get_script_dir(), "mosaic_dtm_params_template.txt")

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) config_path <- cli_args[1]  # Rscript ... params.txt takes precedence over the hardcoded fallback above

if (!file.exists(config_path)) {
  stop("Parameter file not found: ", config_path, ". Set config_path in this script, or pass ",
       "the parameter file's path as a command-line argument (Rscript ... params.txt).")
}

# Reads a "keyword: value" parameter file into a named list of raw (still character) values.
# Same as read_param_file() in ground_returns.R.
read_param_file <- function(path) {
  raw_lines <- readLines(path, warn = FALSE)
  raw_lines <- trimws(raw_lines)
  raw_lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]  # drop blank lines and comment-only lines

  keys <- character(0)
  vals <- character(0)
  for (line in raw_lines) {
    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the FIRST colon splits keyword from value -- keeps values like "c:\..." or "EPSG:2856" intact
    if (colon_pos < 0) {
      stop("Malformed line in parameter file (expected \"keyword: value\"): ", line)
    }
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(substr(line, colon_pos + 1, nchar(line)))
    val <- trimws(sub("\\s+#.*$", "", val))  # strip a trailing " # comment", if any, from the value
    keys <- c(keys, key)
    vals <- c(vals, val)
  }

  if (anyDuplicated(keys)) {
    stop("Duplicate keyword(s) in parameter file ", path, ": ",
         paste(unique(keys[duplicated(keys)]), collapse = ", "))
  }

  setNames(as.list(vals), keys)
}

# The full set of recognized parameter-file keywords: expected type, and default value used
# when a keyword is omitted from the file (required = TRUE means there is no default).
param_specs <- list(
  # Folder containing the DTM tiles. No default -- must be set.
  tile_dir = list(type = "character", required = TRUE),

  # Output DTM path, without extension -- ".flt" (plus GDAL's ".hdr"/".prj" sidecars) is
  # written there. A trailing ".flt" is accepted and ignored. Keep the Fortran programs'
  # "elev_<id>" naming (e.g. c:\work\data\site1\elev_2017): Align takes the DEM's ID from the
  # text after the first "_" in the file name. No default -- must be set.
  output_dtm = list(type = "character", required = TRUE),

  # Polygon (any vector format terra/GDAL reads, e.g. a .shp) delimiting the area to build the
  # DTM for; it may be in any CRS. Only tiles intersecting it are mosaicked, and the output is
  # masked to it. NOFILE mosaics every tile in tile_dir with no masking.
  clip_polygon = list(type = "character", default = "NOFILE"),

  # Regular expression matched (case-insensitively) against file names in tile_dir.
  file_pattern = list(type = "character", default = "\\.tiff?$"),

  # CRS the output DTM must be in. The tiles are reprojected to it only if they aren't in it
  # already. Must be a projected CRS.
  target_crs = list(type = "character", default = "EPSG:2856"),  # NAD83(HARN) / Washington South (Meters)

  # Units of the tiles' elevation values: "m", "ft" (international foot, 0.3048 m), "ftUS"
  # (US survey foot, 1200/3937 m), or "auto". "auto" uses the vertical CRS's unit if the tiles
  # carry a compound (horizontal + vertical) CRS, and otherwise assumes elevations share the
  # horizontal CRS's unit -- with a warning, since that's an assumption, not something the
  # file actually states. Set this explicitly whenever you know it.
  z_units = list(type = "character", default = "auto"),

  # terra::mosaic() function applied where tiles overlap: "first", "last", "mean", "median",
  # "min", "max", or "sum".
  mosaic_fun = list(type = "character", default = "first"),

  # Margin (in input cells) kept around the polygon before reprojecting -- see step 3 above.
  clip_buffer_cells = list(type = "integer", default = 5)
)

# Same as build_params() in ground_returns.R: converts each recognized keyword to its declared
# type, falls back to its default when omitted, and stops on a missing required keyword or an
# unparsable value. Keywords not in param_specs are ignored with a warning.
build_params <- function(raw_params, specs, path) {
  unknown <- setdiff(names(raw_params), names(specs))
  if (length(unknown) > 0) {
    warning("Ignoring unrecognized parameter(s) in ", path, ": ", paste(unknown, collapse = ", "))
  }

  out <- list()
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    if (nm %in% names(raw_params)) {
      raw_val   <- raw_params[[nm]]
      converted <- switch(spec$type,
        character = as.character(raw_val),
        numeric   = suppressWarnings(as.numeric(raw_val)),
        integer   = suppressWarnings(as.integer(raw_val)),
        logical   = suppressWarnings(as.logical(raw_val)),
        stop("Internal error: unknown parameter type '", spec$type, "' for '", nm, "'")
      )
      if (spec$type != "character" && is.na(converted)) {
        stop("Could not parse parameter '", nm, "' (value \"", raw_val, "\") as ", spec$type,
             " -- check ", path, ".")
      }
      out[[nm]] <- converted
    } else if (!is.null(spec$default)) {
      out[[nm]] <- spec$default
    } else if (isTRUE(spec$required)) {
      stop("Required parameter '", nm, "' is missing from parameter file: ", path)
    }
  }
  out
}

params <- build_params(read_param_file(config_path), param_specs, config_path)
invisible(list2env(params, envir = globalenv()))  # makes tile_dir, output_dtm, ... ordinary top-level variables

message("Loaded parameters from: ", config_path)

z_unit_factors <- c(m = 1, ft = 0.3048, ftUS = 1200 / 3937)  # meters per unit
if (!z_units %in% c("auto", names(z_unit_factors))) {
  stop("z_units must be one of auto, m, ft, ftUS (got \"", z_units, "\").")
}
mosaic_funs <- c("first", "last", "mean", "median", "min", "max", "sum")
if (!mosaic_fun %in% mosaic_funs) {
  stop("mosaic_fun must be one of ", paste(mosaic_funs, collapse = ", "),
       " (got \"", mosaic_fun, "\").")
}

use_clip <- nzchar(clip_polygon) && toupper(clip_polygon) != "NOFILE"

output_flt <- paste0(sub("\\.flt$", "", output_dtm, ignore.case = TRUE), ".flt")

## ---- Skip if the output already exists ------------------------------------------------------

if (file.exists(output_flt)) {
  message("Output DTM already exists, nothing to do: ", output_flt)
} else {

  ## ---- Discover tiles and select those intersecting the polygon ----------------------------

  tile_files <- list.files(tile_dir, pattern = file_pattern, full.names = TRUE,
                           ignore.case = TRUE)
  if (length(tile_files) == 0) {
    stop("No DTM tiles matching \"", file_pattern, "\" found in tile_dir: ", tile_dir)
  }

  tiles <- lapply(tile_files, function(f) rast(f)[[1]])  # header only; first band, in case a tile carries more than one

  # Every tile must share one CRS and pixel size -- terra::mosaic() needs a common resolution,
  # and mixed-CRS tiles slipping into tile_dir would otherwise be silently misplaced.
  ref_tile <- tiles[[1]]
  if (crs(ref_tile) == "") {
    stop("DTM tile has no CRS assigned: ", tile_files[1])
  }
  crs_ok <- vapply(tiles, function(t) same.crs(t, ref_tile), logical(1))
  if (!all(crs_ok)) {
    stop("DTM tiles do not share one CRS -- e.g. ", tile_files[!crs_ok][1], " differs from ",
         tile_files[1], ".")
  }
  res_all <- t(vapply(tiles, res, numeric(2)))
  if (max(abs(sweep(res_all, 2, res(ref_tile)))) > 1e-6 * max(res(ref_tile))) {
    stop("DTM tiles do not share one pixel size -- check tile_dir for tiles from a different ",
         "product.")
  }

  native_crs     <- crs(ref_tile)
  native_m_per_u <- linearUnits(ref_tile)  # meters per horizontal CRS unit; 0 for a geographic (lat/long) CRS
  if (native_m_per_u == 0) {
    stop("DTM tiles are in a geographic (lat/long) CRS -- a projected CRS is required to ",
         "express their pixel size in meters.")
  }

  if (use_clip) {
    if (!file.exists(clip_polygon)) stop("clip_polygon not found: ", clip_polygon)
    clip_vect <- aggregate(vect(clip_polygon))  # dissolve every feature into one (multi)polygon
    if (crs(clip_vect) == "") stop("clip_polygon has no CRS assigned: ", clip_polygon)

    clip_native     <- project(clip_vect, native_crs)
    clip_native_buf <- buffer(clip_native, width = clip_buffer_cells * max(res(ref_tile)))

    in_poly <- vapply(tiles, function(t) {
      is.related(as.polygons(ext(t), crs = native_crs), clip_native_buf, "intersects")
    }, logical(1))

    if (!any(in_poly)) {
      stop("None of the ", length(tiles), " DTM tiles in ", tile_dir, " intersect clip_polygon ",
           clip_polygon, ".")
    }
    message(sum(in_poly), " of ", length(tiles), " DTM tiles intersect the clip polygon.")
    tiles <- tiles[in_poly]
  } else {
    message("No clip_polygon -- mosaicking all ", length(tiles), " DTM tiles.")
  }

  ## ---- Mosaic --------------------------------------------------------------------------------

  message("Mosaicking ", length(tiles), " tile(s)...")
  dtm <- if (length(tiles) == 1) tiles[[1]] else mosaic(sprc(tiles), fun = mosaic_fun)
  if (use_clip) dtm <- crop(dtm, clip_native_buf)

  ## ---- Vertical units -> meters --------------------------------------------------------------

  # Meters per vertical unit from a compound CRS's VERTCRS (WKT2) / VERT_CS (WKT1) block, or NA
  # when the CRS has no vertical component. The first unit following that block's opening is
  # the vertical axis's own unit, since a compound CRS lists its horizontal part first.
  vertical_m_per_unit <- function(wkt) {
    vert_start <- regexpr("VERT(CRS|_CS)\\[", wkt)
    if (vert_start < 0) return(NA_real_)
    vert_wkt   <- substr(wkt, vert_start, nchar(wkt))
    unit_match <- regmatches(vert_wkt,
                             regexec("UNIT\\[\"([^\"]+)\",\\s*([0-9.eE+-]+)", vert_wkt))[[1]]
    if (length(unit_match) != 3) return(NA_real_)
    as.numeric(unit_match[3])
  }

  if (z_units == "auto") {
    z_factor <- vertical_m_per_unit(native_crs)
    if (!is.na(z_factor)) {
      message("Elevation units taken from the tiles' vertical CRS: ", signif(z_factor, 10),
              " m per unit.")
    } else {
      z_factor <- native_m_per_u
      warning("DTM tiles carry no vertical CRS; assuming elevations are in the horizontal CRS's ",
              "unit (", signif(z_factor, 10), " m per unit). Set z_units explicitly if that's ",
              "wrong.")
    }
  } else {
    z_factor <- z_unit_factors[[z_units]]
  }

  if (abs(z_factor - 1) > 1e-12) {
    message("Converting elevations to meters (x ", signif(z_factor, 10), ").")
    dtm <- dtm * z_factor
  }

  ## ---- Reproject to target_crs if needed -----------------------------------------------------

  if (same.crs(dtm, target_crs)) {
    message("DTM tiles are already in ", target_crs, " -- no reprojection needed.")
  } else {
    target_m_per_u <- linearUnits(rast(crs = target_crs))
    if (target_m_per_u == 0) {
      stop("target_crs must be a projected CRS: ", target_crs)
    }
    # Keep the input's pixel size: native units -> meters -> target_crs units.
    target_res <- res(dtm) * native_m_per_u / target_m_per_u
    message("Reprojecting to ", target_crs, " at ", paste(signif(target_res, 8), collapse = " x "),
            " (target CRS units) pixels.")
    dtm <- project(dtm, target_crs, method = "bilinear", res = target_res)
  }

  if (use_clip) {
    clip_out <- project(clip_vect, crs(dtm))
    dtm <- mask(crop(dtm, clip_out), clip_out)
  }

  ## ---- Write ---------------------------------------------------------------------------------

  # The Fortran grid reader (..\modules\Grid_Module.f90) only understands ESRI GridFloat header
  # keywords (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value, byteorder), but GDAL's
  # EHdr driver writes its own BIL-style header (ULXMAP/ULYMAP/XDIM/YDIM), so the .hdr is
  # rewritten below from terra's own extent. Not via TWutils' convert_hdr(): that treats
  # ULXMAP/ULYMAP as the upper-left cell's corner, but GDAL defines them as its center, so its
  # output is shifted by half a cell.
  cell_size <- res(dtm)
  if (abs(cell_size[1] - cell_size[2]) > 1e-6 * cell_size[1]) {
    stop("Output DTM pixels are not square (", paste(cell_size, collapse = " x "), ") -- the ",
         ".flt header supports only one cellsize.")
  }

  names(dtm) <- "elevation"
  dir.create(dirname(output_flt), showWarnings = FALSE, recursive = TRUE)
  writeRaster(dtm, output_flt, filetype = "EHdr", datatype = "FLT4S", NAflag = -9999,
              gdal = "BYTEORDER=LSB", overwrite = TRUE)

  hdr_path <- sub("\\.flt$", ".hdr", output_flt, ignore.case = TRUE)
  dtm_ext  <- ext(dtm)
  hdr_lines <- sprintf("%-14s %s",
                       c("ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "NODATA_value",
                         "byteorder"),
                       c(ncol(dtm), nrow(dtm), sprintf("%.6f", xmin(dtm_ext)),
                         sprintf("%.6f", ymin(dtm_ext)), sprintf("%.12f", cell_size[1]), "-9999",
                         "LSBFIRST"))
  writeLines(hdr_lines, hdr_path)
  message("DTM written to: ", output_flt)
  print(dtm)
}
