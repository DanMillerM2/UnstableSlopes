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
2. **Quantiles.md** — theoretical LoD estimation from DTM error-propagation statistics. Only the
   rendered `gfm` output is tracked; its source `Quantiles.qmd` was removed in 5fc726d (recover it
   with `git show 5fc726d^:Quantiles.qmd`).
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

Of the method documents above, only **Coregistration.qmd** and **MappingLandslides.qmd** are tracked
as `.qmd` on `main` (with `Quantiles.md`, above) — `Landslide_Density.qmd`, `Mosaic_Reproject_DTM.qmd`,
and the `SusNRun_PointCloudWorkflow_blank(_*).qmd` files (items 4-5) exist only on the `elise`
branch (`git show elise:<file>` / `git ls-tree -r --name-only elise`). Don't assume they're present
in the working tree; check `git branch -a` before editing or referencing one of them. The rendered
`SusNRun_PointCloudWorkflow_blank_082626.html` and `..._sequential.html` at repo root are tracked
(committed in 3940764) but are renders of that branch's `.qmd`, not evidence the source is here;
their code now lives on as `R/ground_returns.R`.

Three further documents describe the `R/` pipeline rather than the method's theory:

- **Align_LShunter_Review.qmd** — what Align and LShunter do, how they work together, how well they
  meet the two goals (accurate elevation change; a defensible LoD), and improvements tried (its
  rendered `.html` is gitignored).
- **LandslidePipelines_Guide.qmd** — user guide to running the whole method through the
  `landslide-pipelines` Claude Code skill, including setup on another computer (`.html` gitignored).
  **Keep it current, and keep the user's wording.** Whenever a change to the pipeline (the `R/`
  scripts, their keywords and defaults, the site file, the skill, the report, or the versions of
  the programs it needs) would alter anything the guide presents, update the guide in the same
  change. The guide's wording, structure and tone are the user's: update only the facts the
  change affects, edited in place, keeping the surrounding wording, headings and order, and add
  new material in the guide's existing style. Don't reorganize or rewrite sections unless asked,
  and say which passages changed.
- **PipelineReport.qmd** — the report template `run_landslide_pipelines.R` renders for every site;
  tracked (since 1f66259) so every run starts from the same template. Site-specific content goes
  in that site's `narrative.md`, not here; commit only edits meant for every site's report.

The `R/` pipeline, these documents and the skill were developed on branch `align_pipeline`
(Sep-Nov 2026, merged into `main` with `--no-ff`, tagged `align-pipeline-v1`; the branch itself has
since been deleted — the tag's second parent, `align-pipeline-v1^2`, is its last commit). It was
revived for the NAIP images and scar review page (a5cba7a, merged after LShunter 1.4.0's deposits,
which were developed on branch `deposits`), then deleted again. Its commit messages
name the matching Align/TWutils/LShunter versions (e.g. `(Align 2083edb, TWutils 9d80279)`); use
`git log main^1..align-pipeline-v1 -- <file>` or `git log -S <keyword>` to trace why a script or
keyword is the way it is.

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
  an already-mosaicked `.flt` as-is, which must already be in meters and `target_crs`; the
  required `[ref_dsm]`/`[align_dsm]` sections work the same way with `output_dsm`/`input_dsm`/
  `existing_dsm` — Align needs both DSMs), then runs Align through `align_supervisor.R`. The reference DTM is built first;
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
  Water: `[study_area]` `water_polygons` (with `water_layer`, and `water_where`, an SQL condition
  such as `FCode = 39000`) is rasterized onto the reference DTM grid as
  `water_mask_ref_grid.flt` in the run folder (rebuilt when its source changes) and becomes
  Align's `waterMask` — lakes only: rivers move between surveys, so river polygons match neither
  survey's channel (run_site.R asks Align for `outWaterMask`, `watermask_<tag>`, unless
  `align_outWaterMask: FALSE`); and, by default (`extra_counts: water=9` unless `[defaults]` or a set sets
  `extra_counts: NONE`), each set's `water_count_ref_grid.flt` becomes Align's
  `refWater`/`alignWater` (WATER POINTS), which is how rivers are found.
- **align_supervisor.R** — runs Align once (needs Align 1.2.0+, GridUtilities f568852+, which no
  longer requires the outlier raster: the method relies solely on the exceedance raster
  `outExceed`, and `outOutlier` is optional and read by nothing downstream; Align 1.3.0+ writes
  `outDif` and `<outDif>_0` wherever both DTMs have data, not only inside the ground-return mask,
  so they line up cell for cell with the exceedance and patch rasters). With `tileNx`/`tileNy:
  auto` (default) Align itself chooses the tile count and prior by block cross-validation
  (`TILES FIT = JOINT, SELECT = AUTO`; `min_tile_size`, `max_tile_size`, `cv_block`, `cv_buffer`);
  defaults follow Align's CLAUDE.md: `iterations` 8 with `convergence` 0.005 (change of the aligned
  DEM), `derivatives` CENTRAL, `nslope` 10 / `maxSlope` 1.5, `densityLod` auto (with both
  ground-return rasters), and always `OUTPUT GRADIENT RASTER` (`<refDTM>_grad<radius>m`, which
  LShunter needs). Align 1.1.0+ also takes `waterMask`/`waterMaskBuffer`, `refWater`/`alignWater`
  (`waterMinDensity`, `waterWindow`), `openWaterMinArea`, `coverage` (default TRUE: blocks where
  a survey has almost no ground returns -- missing or snow-covered tiles -- are left out;
  `coverageMin`/`coverageBlock` change the rule) and `outWaterMask`
  (why each cell was left out); all left-out cells are nodata in every output. `matchSmoothness`
  (LOCAL/GLOBAL, default NONE; `matchFit` CURVATURE/ROUGHNESS, `matchMaxSigma`, `outSigma`) smooths
  one DTM after the global fit so dz stops depending on surface curvature (2017 is rougher at cell
  scale); the gradient raster and aligned-DEM output keep the original surfaces. `suspectGround`
  (REFERENCE/ALIGN/BOTH, default NONE; `suspectHeight`, `suspectWindow`, `suspectCanopy`,
  `suspectBuffer`) replaces with the local median (code 5; nodata where too few neighbours remain)
  bumps under canopy in one DTM that the other lacks -- vegetation classified as ground, which
  read as loss (2017) or gain. A few seconds in it
  checks Align's settings echo — unrecognized keywords, an
  option not acknowledged (Align drops unknown `TILES` sub-fields silently), or an Align before
  1.2.0 stopping for want of OUTPUT OUTLIER RASTER — and kills an old exe; it also stops if the
  installed TWutils' `align_input()` still requires `outOutlier` (before ac87ce5).
  It runs the exe via `processx` (input file from `TWutils::align_input()`) into
  `<outDTM folder>/align_run/`, streaming `align.log`, then moves every output to its final
  path by base name (`_0`, `_iqr`, `_fit0.csv`, ... with their base; `lodGain`/`lodCentre`
  TRUE also ask Align 1.4.0+ for `<outLod>_gain`/`_centre`, which nothing reads; the aligned
  DEM, outDTM, last: its `.flt` existing means the run is done, here and in
  `run_align_pipeline.R`) and writes `align_run.csv`, `align_iterations.csv`, `align_select.csv` (the CV table) and
  `align_nodes.csv` (joint-fit node totals, mm, with map positions) there. It parses the longest
  `<output>.log` copy, because Align's console output wraps at 80 columns.
- **run_site.R** — front end for one site: reads a single flat `keyword: value` site file
  (`R/site_template.txt`; prefixes `site_`, `ref_`/`align_` survey keys, `ground_`, `align_`
  Align settings, `ril_`, `lshunter_`, `lsstats_`; unknown keywords stop it), writes the six
  configs below into `<site_dir>\<site_run_name>\configs\` and runs `run_landslide_pipelines.R`
  (`--configs-only` writes and checks them). It derives every link: Align's outputs are
  `<run folder>\elev_<align>aligned`, `dif_/exceed_/lod_/cendif_/outbins_<ref>_<align>` (no outlier raster);
  LShunter gets Align's exceedance, gradient, cendif, IQR and correlogram and RIL's landforms
  (`OUTPUT RIL RASTER` of `ril_input`), and with both DSMs the two surveys' DSMs and original DTMs
  (`canopy_*`), writing `patch_<tag>` and `patch_<tag>_table.csv`; LSstats
  gets the patch raster and table; `@cendif`, `@exceed`, ... in `lsstats_rasters` stand for Align
  outputs. With `lshunter_deposits: TRUE` (LShunter 1.4.0) LShunter also gets the reference DTM
  and `<run folder>\ang_<ref>` (its flow directions, which LShunter.R builds with RIL's bldgrds
  executable; no corridor raster is asked for, as nothing reads it); `lshunter_dep_*` keys pass through as `dep_*`, and the
  deposits join the scars in `patch_<tag>_polygons.shp` (type `deposit`, with their scar).
  The configs are regenerated every run; edit the site file, not them.
- **run_landslide_pipelines.R** — driver for the whole method: runs `run_align_pipeline.R` and
  `ril_pipeline.R` concurrently (processx, logs in `<report_dir>/logs/`), then
  `lshunter_lsstats_pipeline.R`, checking first that LShunter's inputs are this run's Align
  outputs: the exceedance raster (`outExceed`, required -- the only change raster used), the
  gradient `<ref DTM>_grad<radius>m`, and when used `cen_dif`/`iqr`/`noise_corr` and the
  seed-density count rasters; writes `pipeline_manifest.txt` and renders **PipelineReport.qmd** (repo root,
  `embed-resources: true`) to `<report_dir>/<study_name>_report.html`, including
  `narrative.md` as its Interpretation section and, with deposits, a Deposits section (map,
  deposit sizes and travel, the scars' volume balance, from the deposit, balance and unlinked
  tables). The report reads only the manifest's outputs, so
  a rerun (every pipeline skipping) just re-renders, keeping the previous run's logs and the
  manifest's `created`/`elapsed_min` (no output file created or changed; each pipeline logs to
  `<pipeline>.log.new`, renamed into place only when the run did something or failed). `.claude/skills/landslide-pipelines/` is the
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
  count per template cell, `ground_count_ref_grid.flt`, instead of the density raster;
  `extra_counts` (e.g. `water=9, snow=21`) counts other LAS classes in the same pass, as
  `<name>_count_ref_grid.flt`; the same pass tallies ground points per tile, second of GPS time and
  PointSourceID into `pass_tally.csv` (for strip_noise.R). It's masked
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
  raster (required, S thresholds −4/−2 by default; the `outlier` keyword is no longer accepted),
  plus Align's gradient raster
  (required) and a flow-accumulation raster, in two rounds of paired thresholds. LShunter 1.1's
  options are keywords too: patch table (`cen_dif`, `iqr` + `noise_corr`, `landform`,
  `dist_channel`, `out_table`, `min_patch_significance`), `null_*` simulation, `seed_den_*`,
  `sign_flip`, and (LShunter 1.2.0+) `canopy_ref_dsm`/`canopy_ref_dem`/`canopy_align_dsm`/
  `canopy_align_dem`, which add each patch's mean canopy height in both surveys and its change to
  the table (`canopy_align_dem` is the align survey's original DTM, not Align's shifted one); with
  `out_table` set, the run counts as done only once the table exists. LShunter 1.3.0+ adds area
  seeds (`area_windows`, e.g. `5 10 20`, default NONE; `area_threshold`, `area_min_fraction`,
  `out_area_z`): seeds from evidence pooled over windows, which need the exceedance raster and
  `noise_corr`. LShunter 1.4.0+ adds deposits (`dep_flow_dir` + `dep_ref_dem`, the reference DTM;
  `dep_*` settings, `out_deposits`/`out_unlinked`/`out_corridor`/`out_paths`): gain patches
  downslope of the scars, linked to them along runout paths traced on the reference DTM's
  D-infinity flow directions (distance counted from the last coherent loss, so scour along a
  debris-flow track doesn't use it up; `dep_max_distance` 500 m beyond it, `dep_spread` 100 m
  either side), seeded out of channels and off the valley floor but free to grow into them; it
  writes `<out_table stem>_deposits.csv` (with scar, travel, drop, reach angle), `_unlinked.csv`
  (gain patches with no scar above them) and `_balance.csv` (each scar's volume against its
  deposits', the net and their SDs), last. `dep_flow_dir` is built here with
  `TWutils::bldgrds_nochannels(out_ang = )` when missing or older than `dep_ref_dem`
  (`dep_bldgrds_executable_dir`), and the contributing area bldgrds always writes beside it
  (`<dep_flow_dir>_accum`, and an empty `BaseNode_<id>.dat`) is deleted; LShunter 1.5.0+ writes
  the unlinked raster only with `out_unlinked`; with deposits the run counts as done once `_balance.csv`
  exists. The scars' sign-flip control then leaves out gain patches in a scar's corridor, and
  `null_*` gives deposits false-discovery rates from simulated noise patches in the same
  corridors. Supports
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
  raster. With `polygons_params` (run_site.R always writes it) it then runs `patch_polygons.R`,
  whose `deposit_raster`, if set, must be LShunter's deposit raster. With deposits,
  LShunter is skipped only once its `_balance.csv` exists too.
- **patch_polygons.R** — LShunter's patch raster as a polygon shapefile (`patch_<tag>_polygons.shp`
  in run_site.R's runs), one feature per patch label (dissolved, so diagonal-only contacts make
  multipart features), with the patch table's row (every column LShunter wrote, -9999 as null),
  perimeter and minimum-rotated-rectangle length and width, and from LSstats's per-polygon csvs
  the landform, distance to the floodplain edge, and min/median/max of each summarized raster.
  With `deposit_raster`/`deposit_table` (LShunter's deposits) the deposits are features of the
  same shapefile: field `type` is `scar` or `deposit`, `label` is the value in that type's raster,
  a deposit's `scar` is its scar's label; fields both tables fill (`mean_dz`, `vol_m3`, ...) are
  signed, and LSstats fields are null for deposits.
  Shapefile field names are 10 characters at most; `<out_shp>_fields.csv` maps each to its
  source column, units and meaning. Rewritten when older than any input. It refuses an `out_shp`
  sharing a stem with another file's `.prj`: `patch_x.shp` beside `patch_x.flt` would replace the
  raster's `.prj`.
- **naip_images.R** — NAIP aerial photos of each scar in patch_polygons.R's shapefile (deposits are left out: their labels can repeat the scars')
  (`p<label>_naip.jpg` in `out_dir`: per NAIP year a close-up, framed like the scar-review images,
  and a `wide_m` view; `naip_images.csv` index with flight dates and status), read window by window
  from Microsoft Planetary Computer's COGs (no account; it searches the STAC catalog itself).
  `years: auto` takes the latest NAIP in or before each of `survey_years`. Planetary Computer SAS
  tokens last ~15 min and reads can come back truncated, so it refreshes tokens every 8 min and
  retries each window; a patch that still fails gets no image (rerun retries it) and the script
  exits 2. Images newer than the shapefile are skipped; `n_workers` patches run in parallel.
  lshunter_lsstats_pipeline.R runs it after the polygons when given `naip_params` (a failure only
  warns); run_site.R writes that config only with `naip_run: TRUE` (`naip_` keywords; off by
  default, ~5-10 s per patch per worker). Washington has no 2025 NAIP yet (2023 is the latest).
- **scar_review.R** — a labelling page for one run (`out_dir`: `index.html` from
  **scar_review_page.html** with `{{TITLE}}` filled in and the sample embedded as `{{PATCHES}}`, so
  it also opens from disk, where browsers block `fetch()`; `patches.json`, `img/p<label>.jpg`): a
  stratified random sample of patch_polygons.R's scars (deposits left out) (`factor_<name>: level = condition; ...;
  level = otherwise`, R expressions over the shapefile's fields, first true wins, NA = false;
  strata are the factors' combinations; `per_stratum`, `seed`, optional `filter`; default factors
  flip_fdr <= 0.05 x canopy_chg < -10 / > 5, which reproduce the Site 2 sample; `sim_fdr` is
  empty unless the run used the null simulation), each drawn as 8
  lidar panels (hillshades, centred difference, 700 m context, canopy height, ground returns/m2).
  The ground-return panels mask what Align leaves out: each survey's coverage gaps (Align's rule,
  `coverage_min`/`coverage_block`, grown into adjacent cells without returns because a 50 m block
  can straddle a missing tile's edge) in grey, and water in blue (Align's `outWaterMask` codes 1-3,
  and each survey's `water_count_ref_grid` at `water_min_density`); `gap_ref`/`gap_aln`/`water` in
  patches.json are each patch's share. Lidar images are redrawn when older than the shapefile or
  the script.
  Inputs are read from `run_dir`'s configs (patch_polygons_params, lshunter_params,
  align_pipeline) or given directly. `naip: TRUE` runs naip_images.R into `img/`; the page shows
  the NAIP strip only for patches with `naip: true` in patches.json. The page saves labels only as
  a claude.ai artifact with the `db` and `user` capabilities (collection `labels`, doc id
  `p<label>__<uid>`: verdict, outline, reasons[], note, stratum); patches.json's `population`
  weights a stratum's labels back to the whole run. Same shapefile, strata and seed give the same
  sample, so a rebuilt page keeps its labels.
- **bldgrds.R** — runs `TWutils::bldgrds()` (Fortran program bldgrds, see
  `..\GridUtilities\bldGrds2.f90`) to compute D-infinity flow direction and flow accumulation for a
  DEM and, unless `no_channels` is set, trace the channel network downstream into a node-list
  database. Unlike RIL.R/align_dtms.R's programs, bldgrds has no single output raster/keyword to
  read back — outputs are named from `path`/`dem_id` (`ang_<dem_id>.flt`, `accum_<dem_id>.flt`,
  `NodeAttributes_<dem_id>.dat`/`NodeNet_<dem_id>.dat`), and `path` must end with a trailing path
  separator since bldgrds concatenates it directly with those prefixes rather than joining it.
- **build_dtm.R** — builds DTM (TIN of `ground_classes`, edges up to `max_edge_m`) and DSM
  (highest return, `dsm_drop_classes` left out) GeoTIFF tiles from LAZ with lasR, one output tile
  per input tile on a `res_m` grid in the LAZ's own CRS and units, each with a `buffer_m` margin
  read from neighbouring tiles; optional `clip_polygon`; resumable (existing tiles are skipped,
  each is written as `_partial` and renamed). lasR leaves `.lax` index files beside the LAZ. A
  Site1 pilot found the rebuilt DTMs give the same patches as the vendor's; vendor DTMs remain the
  default, but the rebuilt ones exposed a half-cell offset in the vendor 2017 GeoTIFFs.
- **strip_noise.R** — one survey's noise measured inside the survey: where two flight passes
  overlap, each pass's ground points are triangulated alone and the surfaces differenced (no real
  change between passes flown minutes apart), on a sample of tiles; histograms and statistics by
  canopy and gradient class (`<label>_noise_classes.csv`, `_hist.csv`, `_pairs.csv`,
  `_passes.csv`, `_summary.txt`), and also by ground-return density of the sparser pass
  (`_density_classes.csv`; `density_breaks`, `density_window_m`). Passes come from `noise_method: auto`: PointSourceID if it passes
  its checks (not all 0, not one per tile; IDs reflown are split by GPS time, IDs flown together
  merged), else GPS-time gaps, else none (no independent reference; the summary says so). Pass
  tallies come from ground_returns.R's `pass_tally.csv` or are made for the sampled tiles. The
  pipeline runs it as run_align_pipeline.R's Stage 4 (`[strip_noise]` section; run_site.R's
  `noise_` keywords; an optional diagnostic, off unless `noise_run: TRUE`; outputs in
  `<run folder>\strip_noise\`), after
  Align and beside LShunter (`--no-strip-noise` / `--strip-noise-only`), and the report's "Noise
  reference" section shows each survey's checks, strip offsets, class table, and the pass-noise
  floor (each survey's single-pass noise read at its ground density, from `ref_ground_count`/
  `align_ground_count` in the manifest) beside the observed DoD IQR and Align's IQR by canopy and
  density. It is a diagnostic and doesn't feed Align. Pilots: within-pass splits (scanner
  channels) understate DoD noise 3-7x. Pooled over density, single-pass differences come within
  about 0.8-1.4x of observed DoD noise under forest (Site1, Site2) and 2.7-3.2x below it on open
  flat ground, but see below. Site1 2017 has strip offsets up to
  0.22 m. Pass noise falls steeply with density (Site1: IQR^2 = a + b/density, r2 >= 0.99), but read
  at the surveys' full density it predicts a DoD IQR 1.5-5x below the observed one, which
  falls less steeply with density: the DoD holds a between-survey term (about 0.1 m under forest at
  Site1) that passes of one survey don't share. The pooled ratio near 1 under forest is a
  coincidence of overlap cells being sparse; pass noise is a floor, not the whole LoD. Nor is the
  rest a constant: pass noise plus one between-survey term per canopy x gradient class, fitted on
  dense ground, underpredicts the DoD IQR at 0.25-2 returns/m2 by 1.2-1.6x under cover (Site1),
  where Align's own density LoD is within 1.0-1.2x; the between-survey term grows as ground
  returns thin (about 0.09 m at 2-4 returns/m2, 0.17 m at 0.25-0.75 under forest).

## Key R dependencies

Chunks rely on: `terra`, `sf`, `lasR` (lidar point-cloud processing), `data.table`, `future` /
`future.apply` (multi-process batch parallelism, `multisession` — required for Windows),
`ggplot2` / `ggExtra` / `patchwork` / `latex2exp` / `RColorBrewer` / `scales`, `stringr`, and a
project-specific package, `TWutils`, which is not on CRAN and must already be installed/available
in the R environment. `TWutils` is the current name; **`TerrainWorksUtils` is the same package under
its old name**, not a second dependency — the rename is incomplete (per `..\TWutils\CLAUDE.md`), and
the removed `Quantiles.qmd` (see the document map) used `library(TerrainWorksUtils)`/
`TerrainWorksUtils::quantiles(...)` rather than the `TWutils::` form used everywhere else here.

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
- The Align work here and in `..\GridUtilities` is the same alignment work seen from opposite ends
  of the pipeline (this repo drives `TWutils::align()`; GridUtilities holds the `align`/`alignSlope2`
  Fortran source; Align 1.1.0 was developed on `align_clean` and merged into `master`, tag
  `align-1.1.0`, but new work may again be on a feature branch) — when
  editing `Coregistration.qmd` or `R/align_dtms.R`/`R/run_align_pipeline.R`/`R/align_supervisor.R`,
  check which branch of `GridUtilities` is actually checked out rather than assuming `master`.
