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
## Instead of tile_dir, input_file can name one already-mosaicked raster (a GeoTIFF, or a .flt
## with its .hdr/.prj, ...). It's treated as a single tile: steps 2-6 below still apply, so it is
## clipped to the polygon, converted to meters, and reprojected exactly like a mosaic would be.
##
## 1. List the tiles in tile_dir matching file_pattern (or take input_file as the only tile) and
##    read each one's extent/CRS (header metadata only -- no pixel values are read at this point)
## 2. Keep only the tiles that intersect clip_polygon (reprojected to the tiles' own CRS)
## 3. Mosaic the kept tiles -- as a GDAL virtual raster (VRT) when mosaic_fun is "first"/"last"
##    and the tiles share one cell grid, so nothing is written yet (otherwise terra::mosaic()).
##    If step 5 must reproject, it first crops to the polygon plus a margin of
##    clip_buffer_cells cells, so the reprojection has real data to interpolate from right up to
##    the polygon's edge instead of leaving a rim of NA cells along it
## 4. Convert elevations to meters (see z_units below) -- as scale metadata on the mosaic, applied
##    when the pixels are read in step 6
## 5. If the tiles' CRS differs from target_crs, reproject (bilinear) at the tiles' own pixel
##    size expressed in target_crs units -- e.g. 3 US-survey-foot tiles come out at
##    3 * 1200/3937 = 0.9144018 m pixels -- so the output keeps the input's resolution. When
##    the two CRSs are the same projection in different units (e.g. EPSG:2927 ftUS vs
##    EPSG:2856 m), the coordinates are just rescaled instead: no resampling, same cells
## 6. Mask to the exact polygon and write <output_dtm>.flt, in one pass over the pixels: written
##    as <output_dtm>_partial.* and renamed into place when complete, so an interrupted run never
##    leaves a truncated .flt that a later run would skip as finished
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
  # Folder containing the DTM tiles, or (input_file) a single raster to process as the only
  # tile. Exactly one of the two must be set.
  tile_dir   = list(type = "character", default = "NOFILE"),
  input_file = list(type = "character", default = "NOFILE"),

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

is_set <- function(x) nzchar(x) && toupper(x) != "NOFILE"

use_clip <- is_set(clip_polygon)

output_flt <- paste0(sub("\\.flt$", "", output_dtm, ignore.case = TRUE), ".flt")

if (is_set(tile_dir) == is_set(input_file)) {
  stop("Set exactly one of tile_dir (a folder of tiles to mosaic) or input_file (a single ",
       "raster) in ", config_path, ".")
}
if (is_set(input_file)) {
  if (!file.exists(input_file)) stop("input_file not found: ", input_file)
  # Writing the output over the input would destroy it -- and once it existed, the skip below
  # would hand back the unprocessed input as if it were the output.
  same_file <- function(a, b) {
    tolower(normalizePath(a, winslash = "/", mustWork = FALSE)) ==
      tolower(normalizePath(b, winslash = "/", mustWork = FALSE))
  }
  if (same_file(input_file, output_flt)) {
    stop("output_dtm must differ from input_file (", input_file, ").")
  }
}

## ---- Skip if the output already exists ------------------------------------------------------

if (file.exists(output_flt)) {
  message("Output DTM already exists, nothing to do: ", output_flt)
} else {

  ## ---- Discover tiles and select those intersecting the polygon ----------------------------

  if (is_set(input_file)) {
    tile_files <- input_file
    source_desc <- input_file
  } else {
    tile_files <- list.files(tile_dir, pattern = file_pattern, full.names = TRUE,
                             ignore.case = TRUE)
    if (length(tile_files) == 0) {
      stop("No DTM tiles matching \"", file_pattern, "\" found in tile_dir: ", tile_dir)
    }
    source_desc <- tile_dir
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
      stop("None of the ", length(tiles), " DTM tile(s) in ", source_desc,
           " intersect clip_polygon ", clip_polygon, ".")
    }
    message(sum(in_poly), " of ", length(tiles), " DTM tiles intersect the clip polygon.")
    tiles      <- tiles[in_poly]
    tile_files <- tile_files[in_poly]
  } else {
    message("No clip_polygon -- mosaicking all ", length(tiles), " DTM tiles.")
  }

  ## ---- Mosaic (virtually) --------------------------------------------------------------------
  ## The mosaic is never written out on its own: with mosaic_fun "first" or "last" and tiles on
  ## one common grid (the usual case), it's a GDAL virtual raster (VRT) over the tiles, and the
  ## unit conversion and rescaling below are metadata changes on it, so the tiles' pixels are
  ## read exactly once -- in the final masked write. Other mosaic_fun values (or tiles that
  ## aren't cell-aligned) fall back to terra::mosaic(), which writes a temporary mosaic first.

  # TRUE when every tile's cell edges line up with the first tile's (origins equal modulo the
  # cell size) -- required for the VRT, which would otherwise resample silently.
  half <- res(ref_tile) / 2
  aligned <- all(vapply(tiles, function(t) {
    max(abs((origin(t) - origin(ref_tile) + half) %% res(ref_tile) - half)) <
      1e-6 * max(res(ref_tile))
  }, logical(1)))

  if (length(tiles) == 1) {
    dtm <- tiles[[1]]
  } else if (mosaic_fun %in% c("first", "last") && aligned) {
    message("Building a virtual mosaic (VRT) of ", length(tiles), " tile(s)...")
    # GDAL paints later VRT sources over earlier ones, so "first" (the first tile wins where
    # tiles overlap, as in terra::mosaic()) needs the tiles listed in reverse.
    dtm <- vrt(if (mosaic_fun == "first") rev(tile_files) else tile_files, options = c("-b", "1"))
  } else {
    message("Mosaicking ", length(tiles), " tile(s)",
            if (!aligned) " (not all on one cell grid)", "...")
    dtm <- mosaic(sprc(tiles), fun = mosaic_fun)
  }

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
    # As scale metadata (composed with any scale/offset the source already carries), so values
    # are converted as they're read in the final write instead of in a pass of their own. An
    # in-memory raster has no file to attach it to, so it's multiplied directly instead.
    dtm <- tryCatch({
      sc <- scoff(dtm)
      scoff(dtm) <- cbind(sc[, 1] * z_factor, sc[, 2] * z_factor)
      dtm
    }, error = function(e) dtm * z_factor)
  }

  ## ---- Reproject to target_crs if needed -----------------------------------------------------

  # TRUE when native_crs and target_crs are the same projection in different linear units (e.g.
  # EPSG:2927, WA South in US survey feet, and EPSG:2856, WA South in meters): then every
  # coordinate converts by the one factor unit_scale, checked here by projecting the extent's
  # corners and center and comparing against plain scaling. Tolerance is 1 mm (in target_crs
  # units, assumed meters-like).
  is_unit_rescale <- function(r, target_crs, unit_scale) {
    e  <- ext(r)
    xy <- cbind(c(xmin(e), xmin(e), xmax(e), xmax(e), mean(c(xmin(e), xmax(e)))),
                c(ymin(e), ymax(e), ymin(e), ymax(e), mean(c(ymin(e), ymax(e)))))
    projected <- tryCatch(project(xy, from = crs(r), to = target_crs), error = function(e) NULL)
    !is.null(projected) && max(abs(projected - xy * unit_scale)) < 1e-3
  }

  if (same.crs(dtm, target_crs)) {
    message("DTM tiles are already in ", target_crs, " -- no reprojection needed.")
  } else {
    target_m_per_u <- linearUnits(rast(crs = target_crs))
    if (target_m_per_u == 0) {
      stop("target_crs must be a projected CRS: ", target_crs)
    }
    unit_scale <- native_m_per_u / target_m_per_u  # target_crs units per native unit
    if (is_unit_rescale(dtm, target_crs, unit_scale)) {
      # Same projection, different units: relabel instead of resampling. Scaling the extent
      # rescales the cells with it, so every cell keeps its own value (no bilinear smoothing)
      # and the grid is exactly the input grid expressed in target_crs units.
      message("DTM tiles' CRS differs from ", target_crs, " only in linear units -- rescaling ",
              "coordinates by ", signif(unit_scale, 10), " instead of reprojecting.")
      ext(dtm) <- as.vector(ext(dtm)) * unit_scale
      crs(dtm) <- target_crs
    } else {
      # Keep the input's pixel size: native units -> meters -> target_crs units. Cropped first
      # to the polygon plus its margin, so only that much is resampled (see step 3 above).
      if (use_clip) dtm <- crop(dtm, clip_native_buf)
      target_res <- res(dtm) * unit_scale
      message("Reprojecting to ", target_crs, " at ", paste(signif(target_res, 8), collapse = " x "),
              " (target CRS units) pixels.")
      dtm <- project(dtm, target_crs, method = "bilinear", res = target_res)
    }
  }

  cell_size <- res(dtm)
  if (abs(cell_size[1] - cell_size[2]) > 1e-6 * cell_size[1]) {
    stop("Output DTM pixels are not square (", paste(cell_size, collapse = " x "), ") -- the ",
         ".flt header supports only one cellsize.")
  }

  ## ---- Clip and write --------------------------------------------------------------------------
  ## One pass: window() limits reading to the polygon's extent -- snapped to the nearest cell
  ## edges, as crop() does by default -- and mask() masks to the polygon while writing the
  ## output. Written under a "_partial" name and renamed into place only once complete, so an
  ## interrupted run never leaves a truncated <output_dtm>.flt that the skip-if-exists check
  ## (here, and in run_align_pipeline.R) would then accept as finished.

  # The Fortran grid reader (..\modules\Grid_Module.f90) only understands ESRI GridFloat header
  # keywords (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value, byteorder), but GDAL's
  # EHdr driver writes its own BIL-style header (ULXMAP/ULYMAP/XDIM/YDIM), so the .hdr is
  # rewritten below from terra's own extent. Not via TWutils' convert_hdr(): that treats
  # ULXMAP/ULYMAP as the upper-left cell's corner, but GDAL defines them as its center, so its
  # output is shifted by half a cell.
  dir.create(dirname(output_flt), showWarnings = FALSE, recursive = TRUE)
  out_base     <- sub("\\.flt$", "", output_flt, ignore.case = TRUE)
  partial_base <- paste0(out_base, "_partial")
  partial_flt  <- paste0(partial_base, ".flt")
  # The files making up one GridFloat raster <base>: .flt, .hdr, .prj, and GDAL's .stx /
  # .flt.aux.xml, if present.
  raster_exts <- c(".flt", ".hdr", ".prj", ".stx", ".flt.aux.xml")
  sidecars <- function(base) {
    candidates <- paste0(base, raster_exts)
    candidates[file.exists(candidates)]
  }
  unlink(sidecars(partial_base))  # leftovers from an earlier interrupted run

  wopt <- list(filetype = "EHdr", datatype = "FLT4S", NAflag = -9999, gdal = "BYTEORDER=LSB",
               names = "elevation")
  if (use_clip) {
    clip_out  <- project(clip_vect, crs(dtm))
    clip_ext  <- intersect(ext(clip_out), ext(dtm))
    if (is.null(clip_ext)) stop("clip_polygon does not overlap the mosaicked tiles.")
    clip_ext <- align(clip_ext, dtm, snap = "near")
    if (all(nzchar(sources(dtm)))) {
      window(dtm) <- clip_ext       # file-backed (VRT, tile, or rescaled): limit reading, no copy
    } else {
      dtm <- crop(dtm, clip_ext)    # in memory (small mosaic() / project() result): window() needs a file
    }
    dtm <- mask(dtm, clip_out, filename = partial_flt, overwrite = TRUE, wopt = wopt)
  } else {
    dtm <- writeRaster(dtm, partial_flt, overwrite = TRUE, filetype = wopt$filetype,
                       datatype = wopt$datatype, NAflag = wopt$NAflag, gdal = wopt$gdal,
                       names = wopt$names)
  }

  dtm_ext  <- ext(dtm)
  hdr_lines <- sprintf("%-14s %s",
                       c("ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "NODATA_value",
                         "byteorder"),
                       c(ncol(dtm), nrow(dtm), sprintf("%.6f", xmin(dtm_ext)),
                         sprintf("%.6f", ymin(dtm_ext)), sprintf("%.12f", cell_size[1]), "-9999",
                         "LSBFIRST"))
  writeLines(hdr_lines, paste0(partial_base, ".hdr"))

  # Rename every file of the finished raster into place -- the .flt last, since its existence is
  # what marks the output complete. Stale sidecars of an earlier output are removed first.
  rm(dtm); invisible(gc())  # release terra's handle on the partial file before renaming it
  unlink(sidecars(out_base))  # output_flt itself can't exist here (see the skip above)
  for (e in c(setdiff(raster_exts, ".flt"), ".flt")) {
    if (!file.exists(paste0(partial_base, e))) next
    if (!file.rename(paste0(partial_base, e), paste0(out_base, e))) {
      stop("Could not rename ", partial_base, e, " to ", out_base, e)
    }
  }
  message("DTM written to: ", output_flt)
  print(rast(output_flt))
}
