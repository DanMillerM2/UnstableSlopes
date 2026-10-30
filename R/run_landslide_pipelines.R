## run_landslide_pipelines.R
##
## Driver script for the whole method: runs run_align_pipeline.R and ril_pipeline.R side by side
## (neither needs the other's output), then lshunter_lsstats_pipeline.R on their results, then
## renders PipelineReport.qmd (repo root) to one self-contained html file summarizing the run.
## Each pipeline is its own `Rscript <script> <config>` subprocess, handed its own config file
## unmodified, with its console output in <report_dir>/logs/<pipeline>.log; the align and RIL
## logs are echoed here as they grow, so Align's tile search can be followed from this console.
##
## Links checked before anything runs (each pipeline only knows its own files, so a mismatch
## would silently use some other, stale raster):
##   - LShunter's input raster must be Align's: its exceedance raster [align] outExceed, or its
##     outlier raster outOutlier (or outDenOutlier);
##   - LShunter's gradient must be the one Align writes, <reference DTM>_grad<radius>m;
##   - LShunter's optional patch-table inputs, when used, must be Align's (cen_dif = outCenDif,
##     iqr = <outLod>_iqr, noise_corr = <outbins>_densityLodCorr.csv) and the seed-density
##     rasters the two LAZ sets' count rasters;
##   - LSstats's landform raster is reported, and a warning given if it isn't RIL's output
##     raster (LSstats may use a landform raster that extends RIL's classes, so this isn't an
##     error).
##
## Resumability: every pipeline skips the outputs it already has (see each one's header), so
## rerunning this script after a failure, or just to re-render the report (e.g. after editing
## the narrative file), redoes only what's missing.
##
## The report reads <report_dir>/pipeline_manifest.txt, written here: every output path the
## report needs, plus the path of an optional narrative file (narrative_md, default
## <report_dir>/narrative.md) whose Markdown is included as the report's interpretation section.
##
## Requires: TWutils, processx, terra, and Quarto (the quarto CLI on the PATH) with the R
## packages PipelineReport.qmd loads; plus whatever the three pipelines require.
##
## ---- Config file --------------------------------------------------------------------------
##   Rscript run_landslide_pipelines.R path/to/landslide_pipelines_config.txt
## "keyword: value" lines (see run_landslide_pipelines_params_template.txt):
##   study_name:              name used in the report title (default: the config file's name)
##   align_pipeline_config:   run_align_pipeline.R config (required)
##   ril_pipeline_config:     ril_pipeline.R config (required)
##   lshunter_lsstats_config: lshunter_lsstats_pipeline.R config (required)
##   report_dir:              folder for the report, manifest, and logs (default: a "report"
##                            folder next to outOutlier)
##   narrative_md:            Markdown file included in the report (default
##                            <report_dir>/narrative.md; leaving it absent is fine)
##   render_report:           TRUE (default) or FALSE
##   progress_interval_s:     how often the align/RIL logs are echoed (default 30)

library(TWutils)

get_script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(file_arg) == 1) dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE)) else "."
}
script_dir <- get_script_dir()

config_path <- file.path(script_dir, "run_landslide_pipelines_params_template.txt")
cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) >= 1) config_path <- cli_args[1]
if (!file.exists(config_path)) stop("Config file not found: ", config_path)

## ---- Reading config and input files ----------------------------------------------------------

# "keyword: value" lines, optionally under [section] headers (sections are returned as named
# sub-lists; lines before any header go in "").
read_keyword_file <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  out <- list(); section <- ""
  for (line in lines) {
    m <- regmatches(line, regexec("^\\[(.+)\\]$", line))[[1]]
    if (length(m) == 2) { section <- tolower(trimws(m[2])); out[[section]] <- list(); next }
    colon_pos <- regexpr(":", line, fixed = TRUE)
    if (colon_pos < 0) stop("Malformed line in ", path, ": ", line)
    key <- trimws(substr(line, 1, colon_pos - 1))
    val <- trimws(sub("\\s+#.*$", "", trimws(substr(line, colon_pos + 1, nchar(line)))))
    if (section == "") out[[key]] <- val else out[[section]][[key]] <- val
  }
  out
}

# One keyword's value from a Fortran-grammar input file (KEYWORD: value[, SUBFIELD=value]), as
# lshunter_lsstats_pipeline.R reads them.
read_fortran_keyword <- function(input_file, keyword) {
  pattern <- paste0("^", gsub(" ", "\\\\s+", keyword), "$")
  lines <- readLines(input_file, warn = FALSE)
  keywords <- vapply(seq_along(lines), function(i) TWutils::get_keyword(lines, i), character(1))
  matches <- which(!is.na(keywords) & grepl(pattern, keywords, ignore.case = TRUE))
  vals <- vapply(matches, function(i) {
    TWutils::parse_arg(sub("\\s+#.*$", "", TWutils::get_args(lines, i)))$Value
  }, character(1))
  if (length(vals) != 1) {
    stop("Expected exactly one '", keyword, ":' line in ", input_file, ", found ", length(vals), ".")
  }
  vals
}

# One sub-field of a keyword (e.g. GRADIENT: FILE = ...), or with field NULL its first value;
# NA when the keyword is absent, which is fine for an optional input.
read_fortran_field <- function(input_file, keyword, field = NULL) {
  pattern <- paste0("^", gsub(" ", "\\\\s+", keyword), "$")
  lines <- readLines(input_file, warn = FALSE)
  keywords <- vapply(seq_along(lines), function(i) TWutils::get_keyword(lines, i), character(1))
  matches <- which(!is.na(keywords) & grepl(pattern, keywords, ignore.case = TRUE))
  if (length(matches) == 0) return(NA_character_)
  if (length(matches) > 1) {
    stop("Expected at most one '", keyword, ":' line in ", input_file, ", found ", length(matches), ".")
  }
  args <- TWutils::parse_args(sub("\\s+#.*$", "", TWutils::get_args(lines, matches)))
  if (is.null(field)) return(args$Value[[1]])
  hit <- which(toupper(args$Parameter) == field)
  if (length(hit)) args$Value[[hit[1]]] else NA_character_
}

is_set <- function(x) !is.null(x) && !is.na(x) && nzchar(x) && toupper(x) != "NOFILE"
raster_key <- function(path) {
  tolower(normalizePath(sub("\\.(flt|tif)$", "", path, ignore.case = TRUE), winslash = "/",
                        mustWork = FALSE))
}
flt_path <- function(base) paste0(sub("\\.flt$", "", base, ignore.case = TRUE), ".flt")

cfg <- read_keyword_file(config_path)
for (key in c("align_pipeline_config", "ril_pipeline_config", "lshunter_lsstats_config")) {
  if (!is_set(cfg[[key]])) stop("'", key, ":' must be set in ", config_path)
  if (!file.exists(cfg[[key]])) stop(key, " not found: ", cfg[[key]])
}
study_name <- if (is_set(cfg$study_name)) cfg$study_name else
  tools::file_path_sans_ext(basename(config_path))
render_report <- if (is.null(cfg$render_report)) TRUE else as.logical(cfg$render_report)
progress_interval_s <- if (is.null(cfg$progress_interval_s)) 30 else as.numeric(cfg$progress_interval_s)

## ---- Each pipeline's outputs ------------------------------------------------------------------

align_cfg <- read_keyword_file(cfg$align_pipeline_config)
al <- align_cfg[["align"]]
for (key in c("outDTM", "outDif", "outOutlier", "outbins")) {
  if (!is_set(al[[key]])) stop("[align] ", key, " is not set in ", cfg$align_pipeline_config)
}
section_raster <- function(section, kind) {
  if (is.null(section)) return("NOFILE")
  for (key in c(paste0("existing_", kind), paste0("output_", kind))) {
    if (is_set(section[[key]])) return(flt_path(section[[key]]))
  }
  "NOFILE"
}

ril_cfg <- read_keyword_file(cfg$ril_pipeline_config)
if (!is_set(ril_cfg$ril_input_file)) stop("ril_input_file is not set in ", cfg$ril_pipeline_config)
ril_raster <- read_fortran_keyword(ril_cfg$ril_input_file, "OUTPUT RIL RASTER")
ril_dem    <- read_fortran_keyword(ril_cfg$ril_input_file, "DEM")

ls_cfg <- read_keyword_file(cfg$lshunter_lsstats_config)
lshunter_p <- read_keyword_file(ls_cfg$lshunter_params)
lsstats_p  <- read_keyword_file(ls_cfg$lsstats_params)
# LShunter's inputs that Align writes: the raster it hunts on (exceedance or outlier), the
# gradient, and the optional patch-table inputs. In mode 1 they come from its input file.
if (is_set(lshunter_p$input_file)) {
  ls_in <- lshunter_p$input_file
  ls_exceedance <- read_fortran_field(ls_in, "EXCEEDANCE RASTER")
  ls_outlier    <- read_fortran_field(ls_in, "OUTLIER")
  ls_gradient   <- read_fortran_field(ls_in, "GRADIENT", "FILE")
  ls_cendif     <- read_fortran_field(ls_in, "CENTERED DIFFERENCE RASTER")
  ls_iqr        <- read_fortran_field(ls_in, "IQR RASTER")
  ls_corr       <- read_fortran_field(ls_in, "NOISE CORRELATION")
  ls_seed_ref   <- read_fortran_field(ls_in, "SEED DENSITY", "REFERENCE")
  ls_seed_align <- read_fortran_field(ls_in, "SEED DENSITY", "ALIGN")
  ls_table      <- read_fortran_field(ls_in, "OUTPUT PATCH TABLE")
  ls_patch      <- read_fortran_keyword(ls_in, "OUTPUT PATCH RASTER")
} else {
  ls_exceedance <- lshunter_p$exceedance
  ls_outlier    <- lshunter_p$outlier
  ls_gradient   <- lshunter_p$gradient
  ls_cendif     <- lshunter_p$cen_dif
  ls_iqr        <- lshunter_p$iqr
  ls_corr       <- lshunter_p$noise_corr
  ls_seed_ref   <- lshunter_p$seed_den_ref
  ls_seed_align <- lshunter_p$seed_den_align
  ls_table      <- lshunter_p$out_table
  ls_patch      <- lshunter_p$out_patch
}
if (is_set(lsstats_p$input_file)) {
  landform_raster <- read_fortran_keyword(lsstats_p$input_file, "LANDFORM RASTER")
  lsstats_csv     <- read_fortran_keyword(lsstats_p$input_file, "OUTPUT CSV")
  input_lines <- readLines(lsstats_p$input_file, warn = FALSE)
  kw <- toupper(vapply(seq_along(input_lines), function(i) TWutils::get_keyword(input_lines, i),
                       character(1)))
  start <- which(grepl("^RASTER\\s+LIST$", kw)); end <- which(grepl("END\\s+LIST", kw))
  in_list <- seq_along(kw) > start[1] & seq_along(kw) < end[end > start[1]][1] & !is.na(kw)
  lsstats_rasters <- vapply(which(in_list), function(i) {
    TWutils::parse_arg(sub("\\s+#.*$", "", TWutils::get_args(input_lines, i)))$Value
  }, character(1))
} else {
  landform_raster <- lsstats_p$landform_raster
  lsstats_csv     <- lsstats_p$out_csv
  lsstats_rasters <- trimws(strsplit(lsstats_p$rasters, ",", fixed = TRUE)[[1]])
}

# Stops unless LShunter's input `ls_path` is the align pipeline's output `al_path` (each
# extension-free); label names them in the message.
check_link <- function(ls_path, al_path, label) {
  if (!is_set(al_path) || raster_key(ls_path) != raster_key(al_path)) {
    stop("LShunter's ", label, " is not the align pipeline's, so LShunter would not use this ",
         "run's alignment:\n  LShunter: ", ls_path,
         "\n  align:    ", if (is_set(al_path)) al_path else "(not set in [align])", call. = FALSE)
  }
}
if (is_set(ls_exceedance)) {
  check_link(ls_exceedance, al$outExceed, "exceedance raster ([align] outExceed)")
  ls_input <- ls_exceedance
} else if (is_set(ls_outlier)) {
  # Align's fence outlier raster, or its density-LoD outlier raster.
  if (is_set(al$outDenOutlier) && raster_key(ls_outlier) == raster_key(al$outDenOutlier)) {
    check_link(ls_outlier, al$outDenOutlier, "outlier raster")
  } else {
    check_link(ls_outlier, al$outOutlier, "outlier raster ([align] outOutlier)")
  }
  ls_input <- ls_outlier
} else {
  stop("LShunter has neither an exceedance nor an outlier raster: ", ls_cfg$lshunter_params)
}

# The gradient must be the one Align writes, so patches use Align's slope and length scale.
align_radius <- if (is_set(al$radius)) as.numeric(al$radius) else 15
align_gradient <- TWutils::align_gradient_raster(section_raster(align_cfg$ref_dtm, "dtm"),
                                                 align_radius)
if (!is_set(ls_gradient)) stop("LShunter has no gradient raster: ", ls_cfg$lshunter_params)
check_link(ls_gradient, align_gradient,
           paste0("gradient raster (Align writes <reference DTM>_grad<radius>m, radius ",
                  align_radius, ")"))

# Optional patch-table inputs, each checked only when LShunter uses it.
al_bins <- sub("\\.csv$", "", al$outbins, ignore.case = TRUE)
if (is_set(ls_cendif)) check_link(ls_cendif, al$outCenDif, "centred difference raster ([align] outCenDif)")
if (is_set(ls_iqr)) {
  check_link(ls_iqr, if (is_set(al$outLod)) paste0(sub("\\.flt$", "", al$outLod, ignore.case = TRUE), "_iqr"),
             "IQR raster (<[align] outLod>_iqr)")
}
if (is_set(ls_corr) && tolower(normalizePath(ls_corr, winslash = "/", mustWork = FALSE)) !=
    tolower(normalizePath(paste0(al_bins, "_densityLodCorr.csv"), winslash = "/", mustWork = FALSE))) {
  stop("LShunter's noise correlation is not the align pipeline's correlogram:\n  LShunter: ",
       ls_corr, "\n  align:    ", al_bins, "_densityLodCorr.csv", call. = FALSE)
}
grnd_path <- function(set) {
  if (is.null(set) || !is_set(set$output_dir)) NA_character_ else
    file.path(set$output_dir, "ground_count_ref_grid.flt")
}
if (is_set(ls_seed_ref))   check_link(ls_seed_ref, grnd_path(align_cfg$ref_set), "seed density reference raster ([ref_set] count raster)")
if (is_set(ls_seed_align)) check_link(ls_seed_align, grnd_path(align_cfg$align_set), "seed density align raster ([align_set] count raster)")
if (raster_key(landform_raster) != raster_key(ril_raster)) {
  warning("LSstats's landform raster is not RIL's output raster -- the report describes the ",
          "one LSstats uses:\n  LSstats landform raster: ", landform_raster,
          "\n  RIL output raster:       ", ril_raster, call. = FALSE)
}

report_dir <- if (is_set(cfg$report_dir)) cfg$report_dir else
  file.path(dirname(al$outOutlier), "report")
log_dir <- file.path(report_dir, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
narrative_md <- if (is_set(cfg$narrative_md)) cfg$narrative_md else file.path(report_dir, "narrative.md")

## ---- Run the pipelines ------------------------------------------------------------------------

start_pipeline <- function(script, config, label) {
  log <- file.path(log_dir, paste0(label, ".log"))
  message("\n==== Starting ", label, " (", script, ")\n  log: ", log)
  proc <- processx::process$new("Rscript", c(file.path(script_dir, script), config),
                                stdout = log, stderr = "2>&1", cleanup = TRUE)
  list(proc = proc, log = log, label = label, shown = 0L, start = Sys.time())
}

# Echoes new log lines from each running pipeline, prefixed with its label.
echo_logs <- function(runs) {
  for (k in seq_along(runs)) {
    r <- runs[[k]]
    if (!file.exists(r$log)) next
    lines <- readLines(r$log, warn = FALSE)
    if (length(lines) > r$shown) {
      new <- lines[(r$shown + 1):length(lines)]
      new <- new[nzchar(trimws(new)) & !grepl("^\\s*\\|-+", new)]
      if (length(new)) message(paste0("[", r$label, "] ", new, collapse = "\n"))
      runs[[k]]$shown <- length(lines)
    }
  }
  runs
}

finish <- function(r) {
  status <- r$proc$get_exit_status()
  elapsed <- as.numeric(difftime(Sys.time(), r$start, units = "mins"))
  if (status != 0) {
    stop(r$label, " failed (exit status ", status, "); last lines of ", r$log, ":\n",
         paste(tail(readLines(r$log, warn = FALSE), 15), collapse = "\n"), call. = FALSE)
  }
  message("==== ", r$label, " finished in ", round(elapsed, 1), " min ====")
  elapsed
}

runs <- list(start_pipeline("run_align_pipeline.R", cfg$align_pipeline_config, "align_pipeline"),
             start_pipeline("ril_pipeline.R", cfg$ril_pipeline_config, "ril_pipeline"))
while (any(vapply(runs, function(r) r$proc$is_alive(), logical(1)))) {
  for (r in runs) if (r$proc$is_alive()) { r$proc$wait(timeout = progress_interval_s * 1000); break }
  runs <- echo_logs(runs)
}
runs <- echo_logs(runs)
elapsed <- c(align_pipeline = finish(runs[[1]]), ril_pipeline = finish(runs[[2]]))

ls_run <- start_pipeline("lshunter_lsstats_pipeline.R", cfg$lshunter_lsstats_config,
                         "lshunter_lsstats_pipeline")
while (ls_run$proc$is_alive()) {
  ls_run$proc$wait(timeout = progress_interval_s * 1000)
  ls_run <- echo_logs(list(ls_run))[[1]]
}
ls_run <- echo_logs(list(ls_run))[[1]]
elapsed["lshunter_lsstats_pipeline"] <- finish(ls_run)

## ---- Manifest and report ----------------------------------------------------------------------

manifest <- list(
  study_name         = study_name,
  created            = format(Sys.time(), "%Y-%m-%d %H:%M"),
  config             = normalizePath(config_path, winslash = "/"),
  clip_polygon       = if (is_set(align_cfg$study_area$clip_polygon)) align_cfg$study_area$clip_polygon else "NOFILE",
  target_crs         = if (is_set(align_cfg$study_area$target_crs)) align_cfg$study_area$target_crs else "EPSG:2856",
  ref_dtm            = section_raster(align_cfg$ref_dtm, "dtm"),
  align_dtm          = section_raster(align_cfg$align_dtm, "dtm"),
  ref_dsm            = section_raster(align_cfg$ref_dsm, "dsm"),
  align_dsm          = section_raster(align_cfg$align_dsm, "dsm"),
  out_dtm            = flt_path(al$outDTM),
  out_dif            = flt_path(al$outDif),
  out_dif_0          = flt_path(paste0(sub("\\.flt$", "", al$outDif, ignore.case = TRUE), "_0")),
  out_outlier        = flt_path(al$outOutlier),
  out_bins           = sub("\\.csv$", "", al$outbins, ignore.case = TRUE),
  out_exceed         = if (is_set(al$outExceed)) flt_path(al$outExceed) else "NOFILE",
  out_lod            = if (is_set(al$outLod)) flt_path(al$outLod) else "NOFILE",
  out_cendif         = if (is_set(al$outCenDif)) flt_path(al$outCenDif) else "NOFILE",
  out_gradient       = flt_path(TWutils::align_gradient_raster(
                         section_raster(align_cfg$ref_dtm, "dtm"),
                         if (is_set(al$radius)) as.numeric(al$radius) else 15)),
  k                  = if (is_set(al$k)) al$k else "1.5",
  align_run_dir      = file.path(dirname(al$outOutlier), "align_run"),
  ril_dem            = ril_dem,
  ril_raster         = flt_path(ril_raster),
  patch_raster       = flt_path(ls_patch),
  lshunter_input     = flt_path(ls_input),
  lshunter_table     = if (is_set(ls_table)) paste0(sub("\\.csv$", "", ls_table, ignore.case = TRUE), ".csv") else "NOFILE",
  landform_raster    = flt_path(landform_raster),
  lsstats_csv        = sub("\\.csv$", "", lsstats_csv, ignore.case = TRUE),
  lsstats_rasters    = paste(lsstats_rasters, collapse = "; "),
  narrative_md       = narrative_md,
  log_dir            = log_dir,
  elapsed_min        = paste(names(elapsed), round(elapsed, 1), sep = "=", collapse = "; ")
)
manifest_path <- file.path(report_dir, "pipeline_manifest.txt")
writeLines(paste0(names(manifest), ": ", unlist(manifest)), manifest_path)
message("\nManifest: ", manifest_path)

if (isTRUE(render_report)) {
  qmd <- normalizePath(file.path(script_dir, "..", "PipelineReport.qmd"), winslash = "/")
  out_name <- paste0(gsub("[^A-Za-z0-9_-]+", "_", study_name), "_report.html")
  message("\n==== Rendering ", basename(qmd), " -> ", file.path(report_dir, out_name), " ====")
  status <- system2("quarto", c("render", shQuote(qmd), "--to", "html",
                                "-P", shQuote(paste0("manifest:", normalizePath(manifest_path, winslash = "/"))),
                                "--output", shQuote(out_name)))
  if (status != 0) stop("quarto render failed (exit status ", status, ").")
  rendered <- file.path(dirname(qmd), out_name)
  if (!file.copy(rendered, file.path(report_dir, out_name), overwrite = TRUE)) {
    stop("Could not copy the rendered report to ", report_dir)
  }
  unlink(rendered)
  message("Report: ", file.path(report_dir, out_name))
}

message("\nAll pipelines complete.")
