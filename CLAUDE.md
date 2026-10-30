# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Research documentation for the "Unstable Slopes Criteria" project: detecting and characterizing
landslides from repeat airborne lidar. It is a set of Quarto (`.qmd`) documents that mix narrative
methodology (with LaTeX math and citations) and executable R code chunks — not a software package.
There is no build/lint/test tooling; "running" a document means rendering it with Quarto/knitr.

This repo is the *consumer* end of a dependency chain, not where any of the underlying tools are
built. Its `.qmd` chunks and the standalone scripts in `R/` (see below) drive named Fortran
programs — Align, HuntLS, LShunter, RIL, Quantiles — through the `TWutils` R package, which wraps
precompiled executables built from source in sibling repos `..\GridUtilities`, `..\ChannelUtilities`,
and `..\LandslideUtilities` (a shared `..\modules` library underlies all three). See each of those
repos' own `CLAUDE.md`, and `..\TWutils\CLAUDE.md` in particular — it cross-checks `TWutils`'s R
wrapper functions' keyword grammars against the current Fortran source and records several found
inconsistencies (e.g. `align()`'s build project actually compiles `alignSlope2.f90`, not `align.f90`)
that are directly relevant to any `.qmd` chunk or `R/` script calling `TWutils::align()`/`huntLS()`/
`LShunter()`/`RIL()`/`quantiles()`.

## Rendering documents

This is an RStudio project (`UnstableSlopes.Rproj`). Render a single document from R/RStudio:

```r
quarto::quarto_render("MappingLandslides.qmd")
```

or from a shell with the Quarto CLI:

```
quarto render MappingLandslides.qmd
```

Rendering executes every R chunk in the file, so it can be slow and I/O-heavy (lidar tiles, large
rasters/point clouds) — don't render as a way to "check syntax." To sanity-check R code without a
full render, source or step through the relevant chunk interactively in R instead.

Some documents render to multiple formats (`html`, `docx`, `pdf`, or `gfm`) per their YAML
frontmatter; `.docx`/`.pdf` outputs and `*_cache/` directories are gitignored, but rendered `.html`
and its `_files/` sidecar directory are frequently committed alongside the source `.qmd` in this repo.

## Document map and reading order

The `.qmd` files build on each other conceptually; read them in this order to understand the
overall method:

1. **Coregistration.qmd** — aligning two sequential lidar DTMs (least-squares horizontal/vertical
   shift, corrected for hillslope gradient and canopy-height change) and estimating the minimum
   level of detection (LoD) for elevation change via Tukey fences on the DEM-of-difference (DoD).
2. **Quantiles.qmd** — theoretical LoD estimation from DTM error-propagation statistics.
3. **MappingLandslides.qmd** — using DoD "k-values" (position relative to the Tukey-fence LoD) plus
   slope and accumulation-area thresholds to detect candidate landslide patches ("Program HuntLS"),
   match them to field-inventoried landslide points, and compare patch area/volume against the
   inventory.
4. **Landslide_Density.qmd** — density/statistics of mapped landslides.
5. **Mosaic_Reproject_DTM.qmd** / **SusNRun_PointCloudWorkflow_blank(_*).qmd** — upstream data-prep
   pipeline: mosaicking downloaded DTM tiles and reprojecting them, and a batched lidar
   point-cloud pre-processing workflow (ground-point filtering, per-batch merge, optional
   GeoPackage export, count/density rasterization, mosaic, reprojection to EPSG:2856). Optional
  `clip_polygon` keeps only tiles whose LAS-header bounding box intersects it (read via
  `lidR::readLAScatalog()`, called by namespace so lidR never masks lasR) and masks the final
  density raster to it; the points themselves aren't clipped. The
   `_blank` files are templates meant to be copied and filled in with a project's own
   `input_dir`/`output_dir` rather than edited in place.

This checkout (branch `align_pipeline`) only has **Coregistration.qmd**, **MappingLandslides.qmd**,
and **Quantiles.qmd** tracked at repo root — `Landslide_Density.qmd`, `Mosaic_Reproject_DTM.qmd`, and
the `SusNRun_PointCloudWorkflow_blank(_*).qmd` files above (items 4-5) currently exist only on the
`elise` branch (`git show elise:<file>` / `git ls-tree -r --name-only elise`), not on
`align_pipeline` or `main`. Don't assume they're present in the working tree; check `git branch -a`
before editing or referencing one of them. (The rendered `SusNRun_PointCloudWorkflow_blank_082626.html`
and `..._sequential.html` sitting untracked at repo root are leftover renders from that branch's
`.qmd`, not evidence the source is checked out here.)

`references.bib` holds shared citations used via `[@key]` across the `.qmd` files.

## Standalone R/ scripts

`R/` holds `Rscript`-able equivalents of code that also appears (or will appear) as `.qmd` chunks —
meant to be run standalone (e.g. as a batch/scheduled job) rather than read as narrative
documentation. Each reads its parameters from a plain-text `keyword: value` file — a different,
simpler grammar than the Fortran executables' own `KEYWORD: value[, SUBFIELD=value] # comment`
input files (see Conventions below) — passed as a command-line argument (`Rscript RIL.R
path/to/params.txt`) or set via a `config_path`-style variable when sourcing the script from
RStudio; each script has a matching `*_params_template.txt` (or, for RIL's separate attribute-list
block, `RIL_attributes_example.txt`) to copy as a starting point.

- **align_dtms.R** — assembled from Coregistration.qmd's first code chunk (`@sec-coregister`); runs
  `TWutils::align()` (Fortran program Align, see `..\GridUtilities`).
- **run_align_pipeline.R** — driver for one study-area polygon (`[study_area]` `clip_polygon`/
  `target_crs`, passed to every sub-run; `clip_polygon: NOFILE`, the default, uses every tile in
  each folder unmasked): runs `ground_returns.R` once each for a "reference" and
  an "align" set of LAZ tiles (ground-point count rasters → `refGrnd`/`alignGrnd`), runs
  `mosaic_dtm.R` once each for a reference and an align folder of DTM tiles (→ `refDTM`/
  `alignDTM`; or a section sets `input_dtm` to run one raster file through `mosaic_dtm.R` as a
  single tile — clipped, unit-converted, reprojected into `output_dtm` — or `existing_dtm` to use
  an already-mosaicked `.flt` as-is, which must already be in meters and `target_crs`; optional
  `[ref_dsm]`/`[align_dsm]` sections work the same way with `output_dsm`/`input_dsm`/`existing_dsm`
  and default to no DSM), then runs Align through `align_supervisor.R`. The reference DTM is built first;
  the align DTM and both DSMs are then built in a background R session (`future`, logging to
  `<output>_mosaic_dtm.log`) while the ground returns are counted, and Align waits for them — a
  counting failure waits for them too before stopping. Each `ground_returns.R` run gets
  `template_raster` = the reference DTM: ground points are
  counted per reference DTM cell (`ground_count_ref_grid.flt` in each set's `output_dir`, used as
  `refGrnd`/`alignGrnd`), so the counts share the DTM's extent, cell corners, and cellsize; the
  pipeline stops if a count raster is older than the reference DTM. Each sub-run is its own `Rscript`
  subprocess and is skipped if its output already exists — apart from that count-raster age check,
  existence is the only check, so outputs must be deleted by hand after changing the polygon or
  other parameters. Stage 3 always runs **align_supervisor.R**; `[align]` is passed to it as is.
- **align_supervisor.R** — runs Align once (needs GridUtilities 8a08fa6+). With `tileNx`/`tileNy:
  auto` (default) Align itself chooses the tile count and prior by block cross-validation
  (`TILES FIT = JOINT, SELECT = AUTO`; `min_tile_size`, `max_tile_size`, `cv_block`, `cv_buffer`);
  defaults follow Align's CLAUDE.md: `iterations` 8 with `convergence` 0.005 (change of the aligned
  DEM), `derivatives` CENTRAL, `nslope` 10 / `maxSlope` 1.5, `densityLod` auto (with both
  ground-return rasters), and always `OUTPUT GRADIENT RASTER` (`<refDTM>_grad<radius>m`, which
  LShunter needs). A few seconds in it checks Align's settings echo — unrecognized keywords, or an
  option not acknowledged (Align drops unknown `TILES` sub-fields silently) — and kills an old exe.
  It runs the exe via `processx` (input file from `TWutils::align_input()`) into
  `<outOutlier folder>/align_run/`, streaming `align.log`, then moves every output to its final
  path by base name (`_0`, `_gain`/`_centre`/`_iqr`, `_fit0.csv`, ... with their base; outOutlier
  last) and writes `align_run.csv`, `align_iterations.csv`, `align_select.csv` (the CV table) and
  `align_nodes.csv` (joint-fit node totals, mm, with map positions) there. It parses the longest
  `<output>.log` copy, because Align's console output wraps at 80 columns.
- **run_landslide_pipelines.R** — driver for the whole method: runs `run_align_pipeline.R` and
  `ril_pipeline.R` concurrently (processx, logs in `<report_dir>/logs/`), then
  `lshunter_lsstats_pipeline.R`, checking first that LShunter's inputs are this run's Align
  outputs: the exceedance (`outExceed`) or outlier (`outOutlier`/`outDenOutlier`) raster, the
  gradient `<ref DTM>_grad<radius>m`, and when used `cen_dif`/`iqr`/`noise_corr` and the
  seed-density count rasters; writes `pipeline_manifest.txt` and renders **PipelineReport.qmd** (repo root,
  `embed-resources: true`) to `<report_dir>/<study_name>_report.html`, including
  `narrative.md` as its Interpretation section. The report reads only the manifest's outputs, so
  a rerun (every pipeline skipping) just re-renders. `.claude/skills/landslide-pipelines/` is the
  Claude Code skill that drives this end to end.
- **mosaic_dtm.R** — mosaics the GeoTIFF DTM tiles that intersect a polygon (or, with
  `input_file` in place of `tile_dir`, takes one raster as the only tile) — as a GDAL VRT, with the
  unit conversion as scale metadata and a units-only CRS change (e.g. EPSG:2927 → 2856) as a
  relabel, so pixels are read once, in a single masked write to `<output>_partial.*` that's renamed
  into place when complete (no truncated `.flt` for a resumed run to skip) — converts elevations
  to meters (`z_units`: auto/m/ft/ftUS), reprojects to `target_crs` only if the tiles aren't
  already in it (keeping the input pixel size, e.g. 3 US ft → 0.9144018 m), masks to the polygon,
  and writes `<output_dtm>.flt`. It rewrites GDAL's BIL-style `.hdr` as an ESRI GridFloat header
  itself, because `..\modules\Grid_Module.f90` only reads ESRI keywords (`ncols`, `xllcorner`,
  `cellsize`, ...); TWutils' `convert_hdr()` does the same job but shifts the grid by half a cell
  (it treats GDAL's cell-center `ULXMAP`/`ULYMAP` as a corner), so don't swap it in.
- **ground_returns.R** — assembled from `SusNRun_PointCloudWorkflow_blank_082626.qmd` (see the
  branch caveat above); the batched lidar point-cloud pre-processing workflow described in the
  document map (ground filtering → optional GeoPackage → per-batch count raster → mosaic →
  reprojection to EPSG:2856). With `template_raster` set, lasR's `rasterize()` (which lays its own
  native-CRS grid) is bypassed: each batch's window of the template grid is fixed first (from the
  tiles' LAS headers), points are read straight from each tile — one `exec()` per tile, ground
  filtered at read time, no merged LAZ unless `export_ground_points_gpkg` needs one — projected into
  the template's CRS, and tallied per template cell in a lasR `callback()`; the output is the point
  count per template cell, `ground_count_ref_grid.flt`, instead of the density raster. It's masked
  with the polygon built exactly as `mosaic_dtm.R` builds it (terra `aggregate()` + `project()`, not
  sf), so it keeps the same boundary cells as the DTM, and — like the per-batch `.tif`s — it's
  written under a `_partial` name and renamed into place when complete.
- **run_ground_returns_batches.R** — driver: runs `ground_returns.R` sequentially, once per named
  set in a batch config file, each as its own `Rscript` subprocess (isolating each set's
  `future::plan` and global-environment state rather than `source()`-ing repeatedly in one session).
- **RIL.R** — runs `TWutils::RIL()` (Fortran program RIL, see `..\ChannelUtilities\RIL\RIL.f90`) to
  trace the channel network from a DEM and classify valley floor / hollow / inner-gorge landforms.
- **LShunter.R** — runs `TWutils::LShunter()` (Fortran program LShunter, see
  `..\LandslideUtilities\LShunter.f90`) to find candidate landslide patches from Align's exceedance
  raster (recommended, S thresholds such as −4/−2) or outlier raster, plus Align's gradient raster
  (required) and a flow-accumulation raster, in two rounds of paired thresholds. LShunter 1.1's
  options are keywords too: patch table (`cen_dif`, `iqr` + `noise_corr`, `landform`,
  `dist_channel`, `out_table`, `min_patch_significance`), `null_*` simulation, `seed_den_*` and
  `sign_flip`; with `out_table` set, the run counts as done only once the table exists. Supports
  both `LShunter()` modes like RIL.R (mode 1: existing input file via `input_file`). Note
  `LShunter()`'s current wrapper signature takes named, capitalized arguments (`Outlier`, `Accum`,
  `OutPatch`, `ScratchDir`, `Executable_dir`, plus `...` to `LShunterInput()`); MappingLandslides.qmd's
  LShunter chunk still uses an older positional call and would fail against it.
- **LSstats.R** — runs `TWutils::LSstats()` (Fortran program LSstats, see
  `..\LandslideUtilities\LSstats.f90`) to summarize a comma-separated list of rasters over
  LShunter's landslide polygons and RIL's landform classes (`<out_csv>_<i>.csv` /
  `<out_csv>_landforms_<i>.csv` per raster). Supports both modes like LShunter.R; distance to the
  floodplain edge comes from `distance_raster` or is computed out to `radius`. With `patch_table`
  (LShunter's `out_table`; LSstats 1.1.0+) each per-polygon csv also has Volume, SD_Volume and
  Z_Volume; LSstats stops if the table's labels or point counts don't match the patch raster.
  `TWutils::LSstats()` also checks that every csv was written, because a Fortran runtime error
  (`forrtl: severe`) still exits with code 0.
- **lshunter_lsstats_pipeline.R** — driver: runs `LShunter.R` then `LSstats.R` as `Rscript`
  subprocesses, each handed its own unmodified parameter file (so either stage may be mode 1 or
  2). Stops unless LShunter's output patch raster is LSstats's polygon raster; skips a stage whose
  final output exists, and stops (rather than skipping) if LSstats's csv is older than the patch
  raster.
- **bldgrds.R** — runs `TWutils::bldgrds()` (Fortran program bldgrds, see
  `..\GridUtilities\bldGrds2.f90`) to compute D-infinity flow direction and flow accumulation for a
  DEM and, unless `no_channels` is set, trace the channel network downstream into a node-list
  database. Unlike RIL.R/align_dtms.R's programs, bldgrds has no single output raster/keyword to
  read back — outputs are named from `path`/`dem_id` (`ang_<dem_id>.flt`, `accum_<dem_id>.flt`,
  `NodeAttributes_<dem_id>.dat`/`NodeNet_<dem_id>.dat`), and `path` must end with a trailing path
  separator since bldgrds concatenates it directly with those prefixes rather than joining it.

## Key R dependencies

Chunks rely on: `terra`, `sf`, `lasR` (lidar point-cloud processing), `data.table`, `future` /
`future.apply` (multi-process batch parallelism, `multisession` — required for Windows),
`ggplot2` / `ggExtra` / `patchwork` / `latex2exp` / `RColorBrewer` / `scales`, `stringr`, and a
project-specific package, `TWutils`, which is not on CRAN and must already be installed/available
in the R environment. `TWutils` is the current name; **`TerrainWorksUtils` is the same package under
its old name**, not a second dependency — the rename is incomplete (per `..\TWutils\CLAUDE.md`) and
`Quantiles.qmd` still has `library(TerrainWorksUtils)`/`TerrainWorksUtils::quantiles(...)` rather
than the `TWutils::` form used everywhere else in this repo.

## Conventions specific to this repo

- **lasR pipelines are lazy and namespace-sensitive**: `terra` is loaded after `lasR` in setup
  chunks and its `rasterize()` generic masks `lasR::rasterize()`; always call `lasR::rasterize()`
  explicitly inside a lasR pipeline, or it silently dispatches to `terra::rasterize()` and errors.
- **lasR double-reads overlapping tiles**: given several files, lasR's `exec()` treats them as one
  seamless point cloud and, by default, makes each file's bounding box a chunk that collects every
  file's points inside it — so where tiles' boxes overlap (e.g. quarter-quad tiles, whose boxes
  overlap at their skewed corners), those points are processed once per overlapping tile (+3.9% on
  an 8-tile site2 2017 block; no point is actually duplicated in the files). Either `exec()` one
  file at a time, or pass `with = list(chunk = <size>)`: regular-grid chunks don't overlap, and
  with one chunk covering everything the result is exact (smaller chunks can still double-read the
  odd point lying exactly on a chunk edge).
- **Batch-oriented processing**: large point-cloud/raster workflows split input tiles into batches
  (`batch_size_tiles`), namespace every intermediate file by batch number (e.g.
  `b003_ground_merged.laz`), and process batches in parallel worker processes via
  `future::plan(multisession, ...)`. This bounds memory/file size, keeps lasR's 32-bit point-index
  ceiling (~2.1B points) from overflowing, and lets an interrupted run resume without redoing
  batches whose output files already exist — preserve this resumability when editing these chunks.
- **`terra::SpatRaster` objects hold an external C++ pointer** and don't survive being returned from
  a `future` worker process — workers write rasters to disk and the main session re-reads them,
  rather than returning the raster object itself.
- Elevation-difference sign convention: negative DoD/k-values indicate elevation loss (erosion,
  landslide scars); positive indicate elevation gain (deposition, aggradation).
- **`executable_dir`/`program_name` point at a sibling repo's build output**, not anything in this
  repo: e.g. MappingLandslides.qmd's HuntLS/LShunter chunks hardcode
  `c:\work\sandbox\landslideutilities\projects\HuntLS\x64\release\` /
  `...\projects\LShunter\x64\release`, and `R/RIL.R`'s/`R/align_dtms.R`'s param templates default to
  `...\channelutilities\projects\ril\x64\release\` / `...\gridutilities\projects\align\x64\release\`
  — these are the same `Projects\<Name>\` build trees documented in `..\ChannelUtilities\CLAUDE.md`/
  `..\GridUtilities\CLAUDE.md`. If a `TWutils` call here starts failing or behaving unexpectedly,
  confirm the executable at that path is up to date with the Fortran source before assuming the bug
  is in this repo's R code.
- This repo's `align_pipeline` branch and `..\GridUtilities`'s current `align_Sep1` branch are the
  same in-progress alignment work seen from opposite ends of the pipeline (this repo drives
  `TWutils::align()`; `align_Sep1` is where the `align`/`alignSlope2` Fortran source itself is being
  changed) — when editing `Coregistration.qmd` or `R/align_dtms.R`/`R/run_align_pipeline.R`, check
  which branch of `GridUtilities` is actually checked out rather than assuming `master`.
