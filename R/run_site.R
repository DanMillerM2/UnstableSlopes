## run_site.R
##
## One site file in, the whole landslide-mapping method out. Reads a single "keyword: value" site
## file (see site_template.txt) holding every input for one study area and one pair of lidar
## surveys, writes the configuration files the pipelines already use into
## <site_dir>\<site_run_name>\configs\, and runs run_landslide_pipelines.R on them:
##
##   align_pipeline.txt            run_align_pipeline.R (ground returns, mosaics, Align)
##   ril_pipeline.txt              ril_pipeline.R (bldgrds, ValleyFloor, RIL)
##   lshunter_params.txt           LShunter.R
##   lsstats_params.txt            LSstats.R
##   lshunter_lsstats_pipeline.txt lshunter_lsstats_pipeline.R
##   landslide_pipelines.txt       run_landslide_pipelines.R
##
## Every link between programs is derived here, never typed twice: Align's output names come from
## the run folder and the survey labels; LShunter reads Align's exceedance raster, gradient,
## centred difference, IQR raster and noise correlogram, and RIL's landform raster; LSstats reads
## LShunter's patch raster and patch table and RIL's landforms. In lsstats_rasters, @cendif,
## @exceed, @dif, @dif0, @outlier, @lod, @gradient and @aligned stand for this run's Align outputs.
## The RIL chain's method parameters stay in its three Fortran input files, which the site file
## names.
##
## The configs are rewritten on every run, so the site file is the one place to change anything.
## Outputs are not: every pipeline skips the outputs it already has, so after changing a parameter
## delete the outputs it affects (each pipeline's header says which).
##
##   Rscript R/run_site.R path\to\site.txt                 write the configs, then run everything
##   Rscript R/run_site.R path\to\site.txt --configs-only  write and check the configs only
##
## Requires: TWutils, plus whatever the pipelines require.

library(TWutils)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}
script_dir <- get_script_dir()

cli_args <- commandArgs(trailingOnly = TRUE)
configs_only <- "--configs-only" %in% cli_args
cli_args <- setdiff(cli_args, "--configs-only")
site_path <- if (length(cli_args) >= 1) cli_args[1] else file.path(script_dir, "site_template.txt")
if (!file.exists(site_path)) {
  stop("Site file not found: ", site_path, ". Pass its path: Rscript run_site.R site.txt")
}

## ---- Read the site file ------------------------------------------------------------------------

read_site_file <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  out <- list()
  for (line in lines) {
    colon_pos <- regexpr(":", line, fixed = TRUE)  # only the first colon: keeps "c:\..." intact
    if (colon_pos < 0) stop("Malformed line in ", path, " (expected \"keyword: value\"): ", line)
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(sub("\\s+#.*$", "", trimws(substr(line, colon_pos + 1, nchar(line)))))
    if (key %in% names(out)) stop("Duplicate keyword '", key, "' in ", path)
    out[[key]] <- val
  }
  out
}

s <- read_site_file(site_path)
is_set <- function(x) !is.null(x) && !is.na(x) && nzchar(x) && toupper(x) != "NOFILE"
get <- function(key, default = NULL) if (is_set(s[[key]])) s[[key]] else default

## ---- Recognized keywords ---------------------------------------------------------------------
## Survey keywords are a fixed set; every other align_ keyword is an Align setting, passed to
## align_supervisor.R without the prefix. ground_ keywords go to ground_returns.R the same way.

survey_keys <- c("label", "laz_dir", "ground_dir", "dtm", "dtm_tile_dir", "dtm_input", "dsm",
                 "dsm_tile_dir", "dsm_input", "z_units", "file_pattern", "clip_buffer_cells")
site_keys <- c("site_name", "site_dir", "site_run_name", "site_polygon", "site_crs", "scratch_dir",
               "report_dir", "narrative_md", "render_report", "progress_interval_s")
align_settings <- c("iterations", "convergence", "k", "dampener", "derivatives", "applySlope",
                    "weights", "groundMaskRadius", "waterGradThreshold", "waterCurvThreshold",
                    "tileNx", "tileNy", "overlap", "tile_fit", "tile_select", "cv_block",
                    "cv_buffer", "min_tile_size", "max_tile_size", "prior_points",
                    "prior_fraction", "prior", "min_tile_points", "radius", "nslope", "maxSlope",
                    "nAzimuth", "densityLod", "localOffsetWindow", "localOffsetMinPoints",
                    "executable_dir", "program_name", "progress_interval_s",
                    # optional extra outputs, TRUE to write them into the run folder
                    "outCorDif", "outDenOutlier", "outDX", "outDY", "outDZ")
ril_keys <- c("ril_bldgrds_input", "ril_valleyfloor_input", "ril_input",
              "ril_bldgrds_executable_dir", "ril_valleyfloor_executable_dir", "ril_executable_dir")
lshunter_settings <- c("threshold1", "threshold2", "min1", "min2", "accum", "max_accum1",
                       "max_accum2", "min_size", "road_shapefile", "road_buffer", "noise_cover",
                       "dist_channel", "min_patch_significance", "null_realizations",
                       "null_tile", "null_seed", "seed_den_min", "seed_den_window", "sign_flip",
                       "executable_dir")
lsstats_settings <- c("rasters", "distance_raster", "radius", "out_distance_raster",
                      "executable_dir")

known <- c(site_keys, paste0("ref_", survey_keys), paste0("align_", survey_keys),
           paste0("align_", align_settings), ril_keys, paste0("lshunter_", lshunter_settings),
           paste0("lsstats_", lsstats_settings))
unknown <- setdiff(names(s)[!startsWith(names(s), "ground_")], known)
if (length(unknown)) {
  stop("Unrecognized keyword(s) in ", site_path, ": ", paste(unknown, collapse = ", "),
       ". See site_template.txt for the keywords.")
}

required <- c("site_dir", "site_run_name", "scratch_dir",
              "ref_label", "ref_laz_dir", "ref_ground_dir", "ref_dtm",
              "align_label", "align_laz_dir", "align_ground_dir", "align_dtm",
              "align_executable_dir", ril_keys[1:3], "lshunter_accum", "lshunter_executable_dir",
              "lsstats_rasters", "lsstats_executable_dir")
missing <- required[!vapply(required, function(k) is_set(s[[k]]), logical(1))]
if (length(missing)) stop("Required keyword(s) missing from ", site_path, ": ", paste(missing, collapse = ", "))

## ---- Checks on the inputs ----------------------------------------------------------------------

flt_path <- function(base) paste0(sub("\\.flt$", "", base, ignore.case = TRUE), ".flt")
no_ext   <- function(path) sub("\\.(flt|csv)$", "", path, ignore.case = TRUE)
problems <- character(0)
need <- function(ok, msg) if (!ok) problems <<- c(problems, msg)

need(dir.exists(get("scratch_dir")), paste0("scratch_dir not found: ", get("scratch_dir")))
if (is_set(s$site_polygon)) need(file.exists(s$site_polygon), paste0("site_polygon not found: ", s$site_polygon))
for (sv in c("ref", "align")) {
  need(dir.exists(get(paste0(sv, "_laz_dir"))),
       paste0(sv, "_laz_dir not found: ", get(paste0(sv, "_laz_dir"))))
  for (kind in c("dtm", "dsm")) {
    out <- get(paste0(sv, "_", kind))
    src <- c(tile_dir = get(paste0(sv, "_", kind, "_tile_dir")), input = get(paste0(sv, "_", kind, "_input")))
    if (length(src) > 1) problems <- c(problems, paste0("set only one of ", sv, "_", kind, "_tile_dir and ", sv, "_", kind, "_input"))
    if (length(src) && names(src)[1] == "tile_dir") need(dir.exists(src[1]), paste0(sv, "_", kind, "_tile_dir not found: ", src[1]))
    if (length(src) && names(src)[1] == "input") need(file.exists(src[1]), paste0(sv, "_", kind, "_input not found: ", src[1]))
    if (length(src) && is.null(out)) problems <- c(problems, paste0(sv, "_", kind, " (the output name) is not set"))
    if (!length(src) && !is.null(out)) need(file.exists(flt_path(out)),
      paste0(sv, "_", kind, " doesn't exist and no ", sv, "_", kind, "_tile_dir or _input is given to build it: ", out))
  }
}
if (is.null(get("ref_dsm")) != is.null(get("align_dsm"))) {
  problems <- c(problems, "give both ref_dsm and align_dsm, or neither (Align needs both for canopy terms)")
}
for (k in ril_keys[1:3]) need(file.exists(s[[k]]), paste0(k, " not found: ", s[[k]]))
need(file.exists(flt_path(s$lshunter_accum)), paste0("lshunter_accum not found: ", s$lshunter_accum))
exe_dirs <- c(align = "align_executable_dir", LShunter = "lshunter_executable_dir",
              LSstats = "lsstats_executable_dir")
for (prog in names(exe_dirs)) {
  exe <- file.path(s[[exe_dirs[[prog]]]], paste0(if (prog == "align") get("align_program_name", "align") else prog, ".exe"))
  need(file.exists(exe), paste0(prog, " executable not found: ", exe))
}
if (length(problems)) stop("Problems in ", site_path, ":\n  ", paste(problems, collapse = "\n  "))

## ---- Derived paths -------------------------------------------------------------------------------

run_dir    <- file.path(s$site_dir, s$site_run_name)
config_dir <- file.path(run_dir, "configs")
dir.create(config_dir, showWarnings = FALSE, recursive = TRUE)
win  <- function(path) gsub("/", "\\\\", path)  # the Fortran programs and configs use backslashes
in_run <- function(name) win(file.path(run_dir, name))
tag <- paste0(s$ref_label, "_", s$align_label)

align_out <- list(
  outDTM     = in_run(paste0("elev_", s$align_label, "aligned")),
  outDif     = in_run(paste0("dif_", tag)),
  outOutlier = in_run(paste0("outlier_", tag)),
  outbins    = in_run(paste0("outbins_", tag)),
  outExceed  = in_run(paste0("exceed_", tag)),
  outLod     = in_run(paste0("lod_", tag)),
  outCenDif  = in_run(paste0("cendif_", tag)))
extra_out <- c(outCorDif = "difcor_", outDenOutlier = "denoutlier_", outDX = "dx_", outDY = "dy_",
               outDZ = "dz_")
for (k in names(extra_out)) {
  v <- get(paste0("align_", k))
  if (!is.null(v) && isTRUE(as.logical(v))) align_out[[k]] <- in_run(paste0(extra_out[[k]], tag))
}
radius <- as.numeric(get("align_radius", "15"))
gradient <- win(TWutils::align_gradient_raster(flt_path(s$ref_dtm), radius))
count_raster <- function(sv) win(file.path(s[[paste0(sv, "_ground_dir")]], "ground_count_ref_grid.flt"))

# RIL's landform raster: the RIL input file's OUTPUT RIL RASTER.
read_fortran_value <- function(input_file, keyword) {
  pattern <- paste0("^", gsub(" ", "\\\\s+", keyword), "$")
  lines <- readLines(input_file, warn = FALSE)
  kw <- vapply(seq_along(lines), function(i) TWutils::get_keyword(lines, i), character(1))
  hit <- which(!is.na(kw) & grepl(pattern, kw, ignore.case = TRUE))
  if (length(hit) != 1) stop("Expected one '", keyword, ":' line in ", input_file, ", found ", length(hit))
  TWutils::parse_arg(sub("\\s+#.*$", "", TWutils::get_args(lines, hit)))$Value
}
ril_raster <- read_fortran_value(s$ril_input, "OUTPUT RIL RASTER")
ril_dem    <- read_fortran_value(s$ril_input, "DEM")
if (tolower(normalizePath(no_ext(ril_dem), winslash = "/", mustWork = FALSE)) !=
    tolower(normalizePath(no_ext(s$ref_dtm), winslash = "/", mustWork = FALSE))) {
  message("Note: the RIL chain uses DEM ", ril_dem, ", not the reference DTM ", s$ref_dtm,
          ". LShunter and LSstats read its landform raster, and LShunter its flow accumulation, ",
          "by map position; check that they cover the study area.")
}

out_patch <- in_run(paste0("patch_", tag))
out_table <- in_run(paste0("patch_", tag, "_table.csv"))

## ---- Writers -------------------------------------------------------------------------------------

write_lines <- function(lines, name, header) {
  path <- file.path(config_dir, name)
  writeLines(c(paste0("# ", header), paste0("# Written by run_site.R from ", normalizePath(site_path, winslash = "\\"),
                                         " on ", format(Sys.time(), "%Y-%m-%d %H:%M"),
                                         ". Edit the site file, not this file."), "", lines), path)
  win(path)
}
kv <- function(lst) {
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (!length(lst)) character(0) else paste0(names(lst), ": ", unlist(lst))
}
section <- function(name, lst) c(paste0("[", name, "]"), kv(lst), "")

survey_raster_section <- function(sv, kind) {
  out <- get(paste0(sv, "_", kind))
  if (is.null(out)) return(list(existing_dsm = "NOFILE"))
  tile <- get(paste0(sv, "_", kind, "_tile_dir")); input <- get(paste0(sv, "_", kind, "_input"))
  extra <- list(z_units = get(paste0(sv, "_z_units")), file_pattern = get(paste0(sv, "_file_pattern")),
                clip_buffer_cells = get(paste0(sv, "_clip_buffer_cells")))
  if (!is.null(tile)) {
    c(setNames(list(tile, out), c("tile_dir", paste0("output_", kind))), extra)
  } else if (!is.null(input)) {
    extra$file_pattern <- NULL
    c(setNames(list(input, out), c(paste0("input_", kind), paste0("output_", kind))), extra)
  } else {
    setNames(list(flt_path(out)), paste0("existing_", kind))
  }
}

## ---- 1. Align pipeline ----------------------------------------------------------------------------

ground <- s[startsWith(names(s), "ground_")]
names(ground) <- sub("^ground_", "", names(ground))
align_section <- list()
for (k in setdiff(align_settings, names(extra_out))) {
  v <- get(paste0("align_", k))
  if (!is.null(v)) align_section[[k]] <- v
}
align_section <- c(align_section, align_out, list(scratch_dir = s$scratch_dir))

align_cfg <- write_lines(c(
  section("study_area", list(clip_polygon = get("site_polygon", "NOFILE"), target_crs = get("site_crs", "EPSG:2856"))),
  section("defaults", ground),
  section("ref_set", list(input_dir = s$ref_laz_dir, output_dir = s$ref_ground_dir)),
  section("align_set", list(input_dir = s$align_laz_dir, output_dir = s$align_ground_dir)),
  section("ref_dtm", survey_raster_section("ref", "dtm")),
  section("align_dtm", survey_raster_section("align", "dtm")),
  section("ref_dsm", survey_raster_section("ref", "dsm")),
  section("align_dsm", survey_raster_section("align", "dsm")),
  section("align", align_section)),
  "align_pipeline.txt", "run_align_pipeline.R config")

## ---- 2. RIL pipeline ------------------------------------------------------------------------------

ril_cfg <- write_lines(kv(list(
  bldgrds_input_file = s$ril_bldgrds_input,
  valleyfloor_input_file = s$ril_valleyfloor_input,
  ril_input_file = s$ril_input,
  bldgrds_executable_dir = get("ril_bldgrds_executable_dir"),
  valleyfloor_executable_dir = get("ril_valleyfloor_executable_dir"),
  ril_executable_dir = get("ril_executable_dir"))),
  "ril_pipeline.txt", "ril_pipeline.R config")

## ---- 3. LShunter and LSstats ----------------------------------------------------------------------

lsh <- list(
  exceedance = align_out$outExceed,
  gradient = gradient,
  gradient_units = "SIN",
  accum = s$lshunter_accum,
  out_patch = out_patch,
  scratch_dir = s$scratch_dir,
  executable_dir = s$lshunter_executable_dir,
  threshold1 = get("lshunter_threshold1", "-4"),   # S = -log10(p)
  threshold2 = get("lshunter_threshold2", "-2"),
  cen_dif = align_out$outCenDif,
  iqr = paste0(align_out$outLod, "_iqr"),
  noise_corr = paste0(align_out$outbins, "_densityLodCorr.csv"),
  landform = win(no_ext(ril_raster)),
  out_table = out_table)
for (k in setdiff(lshunter_settings, c("threshold1", "threshold2", "accum", "executable_dir"))) {
  v <- get(paste0("lshunter_", k))
  if (!is.null(v)) lsh[[k]] <- v
}
if (!is.null(lsh$seed_den_min)) {  # seed density: the two count rasters Align uses
  lsh$seed_den_ref <- count_raster("ref")
  lsh$seed_den_align <- count_raster("align")
}
if (!is.null(lsh$road_shapefile) && is.null(lsh$road_buffer)) {
  stop("lshunter_road_shapefile needs lshunter_road_buffer (m).")
}
lsh_params <- write_lines(kv(lsh), "lshunter_params.txt", "LShunter.R parameters (mode 2)")

align_token <- c(`@cendif` = align_out$outCenDif, `@exceed` = align_out$outExceed,
                 `@dif` = align_out$outDif, `@dif0` = paste0(align_out$outDif, "_0"),
                 `@outlier` = align_out$outOutlier, `@lod` = align_out$outLod,
                 `@gradient` = gradient, `@aligned` = align_out$outDTM)
rasters <- trimws(strsplit(s$lsstats_rasters, ",", fixed = TRUE)[[1]])
rasters <- rasters[nzchar(rasters)]
bad_token <- rasters[startsWith(rasters, "@") & !rasters %in% names(align_token)]
if (length(bad_token)) {
  stop("Unknown lsstats_rasters token(s): ", paste(bad_token, collapse = ", "), ". Known: ",
       paste(names(align_token), collapse = ", "))
}
rasters <- ifelse(rasters %in% names(align_token), align_token[rasters], rasters)
for (r in rasters[!rasters %in% align_token]) {
  if (!file.exists(flt_path(r))) stop("lsstats_rasters: raster not found: ", r)
}

lss <- list(
  polygon_raster = out_patch,
  landform_raster = win(no_ext(ril_raster)),
  rasters = paste(rasters, collapse = ", "),
  out_csv = in_run(paste0("lsstats_", tag)),
  scratch_dir = s$scratch_dir,
  executable_dir = s$lsstats_executable_dir,
  patch_table = out_table)
if (is_set(s$lsstats_distance_raster)) {
  lss$distance_raster <- s$lsstats_distance_raster
} else {
  lss$radius <- get("lsstats_radius", "1000")
  lss$out_distance_raster <- get("lsstats_out_distance_raster", in_run("distance_floodplain"))
}
lss_params <- write_lines(kv(lss), "lsstats_params.txt", "LSstats.R parameters (mode 2)")

ls_cfg <- write_lines(kv(list(lshunter_params = lsh_params, lsstats_params = lss_params)),
                      "lshunter_lsstats_pipeline.txt", "lshunter_lsstats_pipeline.R config")

## ---- 4. Driver --------------------------------------------------------------------------------------

driver_cfg <- write_lines(kv(list(
  study_name = get("site_name", paste0(basename(s$site_dir), ", ", s$ref_label, "-", s$align_label)),
  align_pipeline_config = align_cfg,
  ril_pipeline_config = ril_cfg,
  lshunter_lsstats_config = ls_cfg,
  report_dir = get("report_dir", in_run("report")),
  narrative_md = get("narrative_md"),
  render_report = get("render_report", "TRUE"),
  progress_interval_s = get("progress_interval_s", "30"))),
  "landslide_pipelines.txt", "run_landslide_pipelines.R config")

message("Configs written to ", win(config_dir), ":\n  ",
        paste(basename(c(align_cfg, ril_cfg, lsh_params, lss_params, ls_cfg, driver_cfg)), collapse = "\n  "),
        "\nOutputs go to ", win(run_dir), " (tag ", tag, "); gradient raster ", gradient, ".")

if (configs_only) quit(save = "no", status = 0)

## ---- Run --------------------------------------------------------------------------------------------

status <- system2("Rscript", shQuote(c(file.path(script_dir, "run_landslide_pipelines.R"), driver_cfg)))
if (status != 0) stop("run_landslide_pipelines.R failed (exit status ", status, ").")
