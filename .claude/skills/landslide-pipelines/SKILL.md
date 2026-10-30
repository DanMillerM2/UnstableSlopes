---
name: landslide-pipelines
description: Run the full landslide-mapping method for one study area -- DTM coregistration with Align (tiles chosen by cross-validation, density-based level of detection), RIL landforms, LShunter candidate landslides and LSstats summaries -- then write the interpretation and render the shareable html report. Use when someone wants to process a new site or pair of lidar dates, rerun a site, or regenerate its report.
---

# Landslide pipelines

This repo chains three pipelines and a report:

| Step | Script | Produces |
|---|---|---|
| 1. Align (and its inputs) | `R/run_align_pipeline.R` -> `R/align_supervisor.R` | aligned DTM, DoD (`outDif`, `<outDif>_0`), outlier raster, exceedance raster, LoD rasters (`<outLod>`, `_gain`, `_centre`, `_iqr`), centred difference, gradient raster `<ref DTM>_grad<radius>m`, `<outbins>_*.csv`, `align_run\*.csv` |
| 2. RIL (independent of 1) | `R/ril_pipeline.R` | landform raster (codes 1-10 as LSstats names them: Floodplain, Bedrock_Hollow, Inner_Gorge, Steep_/Convergent, Divergent, Planar, Other) |
| 3. LShunter + LSstats | `R/lshunter_lsstats_pipeline.R` | patch raster, patch table (`out_table`, `_null.csv`), `<out_csv>_<i>.csv`, `<out_csv>_landforms_<i>.csv` |
| Report | `PipelineReport.qmd` (repo root, gitignored -- see section 6) | one self-contained html file |

`R/run_landslide_pipelines.R` runs 1 and 2 side by side, then 3, then renders the report. **Run
that script; don't reimplement its steps.** Every pipeline skips outputs that already exist, so
rerunning it is how you resume after a failure and how you re-render after editing the
narrative.

Align makes its own choices; you don't search for them. With `tileNx: auto` it chooses the tile
count and prior strength by block cross-validation (`TILES FIT = JOINT, SELECT = AUTO`: square
tiles from `max_tile_size` 8000 m down to `min_tile_size` 500 m, held-out 1 km blocks with a 100 m
buffer), and each loop stops once the aligned DEM changes by less than `convergence` (0.005 m:
the 99th percentile on steep ground and the predicted remaining change), at most `iterations` (8).
`min_tile_size` is a judgement, not a fit: a shift field finer than the scale of real change (a
large slow landslide, widespread deposition) can absorb that change. Change it only if the user
has a view on that scale. Read the header of `R/align_supervisor.R` once before your first run;
Align's own `CLAUDE.md` (`..\GridUtilities\Projects\align`) explains the method and the evidence
behind the defaults.

## 1. Collect the inputs

Ask for anything not already given, in one message:

- **Study name**, and an **output folder** for this site.
- **Study-area polygon** (any CRS) and the **target CRS** (default EPSG:2856).
- For each of the reference (older) and align (newer) lidar: the **DTM tile folder** (or one DTM
  file), the **LAZ tile folder**, and optionally the **DSM tile folder** (or file). Without both
  LAZ sets there is no density LoD and no exceedance raster; without both DSMs, no canopy terms.
- The **complete DTM for RIL** (ESRI GridFloat `.flt`, named `<prefix>_<ID>.flt` -- the programs
  derive their data ID from the text after the first `_`).
- **RIL-chain input files** to use as templates: the bldgrds (enforce), ValleyFloor, and RIL
  Fortran input files from a previous site. `ril_pipeline.R` runs in mode 1 only, so these carry
  the method's parameters; you will copy them and change only the paths.
- **LShunter and LSstats parameters**: their `*_params.txt` files from a previous site, or the
  values (S thresholds, flow-accumulation raster, rasters to summarize, landform raster). The
  gradient is always Align's; don't ask for one.
- The **executable folders** (`Projects\<Name>\x64\release\` under the sibling repos), if they
  differ from the templates' defaults.

## 2. Check prerequisites

Run these and fix anything that fails before writing configs:

```
Rscript -e "a <- names(formals(TWutils::align_input)); stopifnot(all(c('tileSelect','minTileSize','outGradient','densityLod','outExceed','outCenDif') %in% a)); l <- names(formals(TWutils::LShunterInput)); stopifnot(all(c('Exceedance','Gradient','OutTable','SignFlip') %in% l)); stopifnot(exists('align_gradient_raster', asNamespace('TWutils')), 'patch_table' %in% names(formals(TWutils::LSstats_input))); for (p in c('terra','sf','processx','ggplot2','knitr')) stopifnot(requireNamespace(p, quietly = TRUE))"
quarto --version
```

- If TWutils fails the check, reinstall it from `..\TWutils` (`R CMD INSTALL`; current `main`).
- **Align** must be built from GridUtilities `align_clean` 8a08fa6 or later (`Projects\align`,
  Release|x64). An older `align.exe` is caught a few seconds into the run -- the supervisor stops
  with "Align ... is older than this pipeline" and names what it didn't recognize. **LShunter**
  must be 1.1.0 or later (LandslideUtilities `LShunter_update` e37b98c; its log's first lines show
  the version), and so must **LSstats** (patch volumes; earlier versions crash on a raster path
  containing "DoD"). The fix for any of them is to rebuild in Visual Studio. Tell the user; don't try to
  work around it.
- **Memory**: Align holds whole DEMs in memory (about 48 GB on a 370 M-cell site). Run one site at
  a time, and don't start other Align jobs while it runs.

## 3. Write the configs

Put everything for the site in its output folder. Start each file from the repo template and
change only what the user's inputs require; keep the templates' comments.

1. **Align** -- copy `R/run_align_pipeline_params_template.txt`. Fill `[study_area]`, `[ref_set]`/
   `[align_set]` (LAZ `input_dir`, `output_dir`), `[ref_dtm]`/`[align_dtm]` (`tile_dir` +
   `output_dtm`, or `input_dtm` + `output_dtm`), `[ref_dsm]`/`[align_dsm]` if DSMs were given,
   and the `[align]` outputs: `outDTM`, `outbins`, `outDif`, `outOutlier`, and with both LAZ sets
   `outExceed`, `outLod`, `outCenDif`. Every output's file name must differ (Align writes them all
   into one run folder first). Keep the template's settings (`tileNx: auto`, `iterations: 8`,
   `convergence: 0.005`, `derivatives: CENTRAL`, `nslope: 10`, `maxSlope: 1.5`, `densityLod:
   auto`) unless the user says otherwise.
2. **RIL** -- copy the three template input files, and in each change the `DEM:` (bldgrds:
   `DEM FILE:` too) line to the complete DTM, `SCRATCH DIRECTORY:` to the site's scratch folder,
   and every output path (e.g. `OUTPUT RIL RASTER:`) to the output folder. All three must name
   the same DEM. Then copy `R/ril_pipeline_params_template.txt` and point it at them.
3. **LShunter/LSstats** -- start from the user's params files (or `R/LShunter_params_template.txt`
   / `R/LSstats_params_template.txt`). The driver stops unless LShunter's inputs are this run's
   Align outputs, so set them from the `[align]` section:
   - `exceedance` = `outExceed` (recommended; thresholds in S = -log10 p, e.g. `threshold1: -4`,
     `threshold2: -2`), or `outlier` = `outOutlier` with k thresholds (e.g. -5 / -1.5);
   - `gradient` = `<ref DTM>_grad<radius>m` (`TWutils::align_gradient_raster(refDTM, radius)`
     gives the path; radius 15 unless `[align]` sets it), `gradient_units: SIN`;
   - for the patch table: `cen_dif` = `outCenDif`, `iqr` = `<outLod>_iqr`, `noise_corr` =
     `<outbins>_densityLodCorr.csv`, and an `out_table`; `sign_flip: TRUE` gives the gain patches
     as a false-alarm control. `null_realizations` (simulated noise) and `seed_den_*` (the two LAZ
     sets' `ground_count_ref_grid.flt`) are optional -- use them if the previous site did.
   - LSstats `polygon_raster` must equal LShunter `out_patch`. With an LShunter `out_table`, set
     LSstats `patch_table` to it: every per-polygon csv then gets Volume, SD_Volume and Z_Volume
     (the driver stops if they differ, and LSstats stops if the table doesn't match the patches).
   Then copy `R/lshunter_lsstats_pipeline_params_template.txt`.
4. **Driver** -- copy `R/run_landslide_pipelines_params_template.txt` and set `study_name` and
   the three config paths.

Show the user the `[align]` section, LShunter's inputs and thresholds, and the paths you changed
before running.

## 4. Run and monitor

Start the driver in the background (Align alone takes about 25 min on a 370 M-cell site; the
ground-return counts and mosaics can take hours):

```
Rscript R/run_landslide_pipelines.R <site folder>\landslide_pipelines.txt
```

Logs go to `<report_dir>\logs\` (`align_pipeline.log`, `ril_pipeline.log`,
`lshunter_lsstats_pipeline.log`); Align's own console output is in
`<outOutlier folder>\align_run\align.log`. Don't poll in a tight loop -- wait for the completion
notification, checking the logs at long intervals. Give the user one line at each Align milestone:
the global loop converged (iterations, RMSE), the tiles chosen (`Selected: tiles nx x ny, prior
scaled by s`, and whether it says the choice is the finest grid tried), the tile loop converged.

On failure, read the end of the failing log and report the cause. Rerunning resumes from what
exists. **Ask before deleting any output** -- after a parameter change, the outputs it affects
must be deleted for the change to take effect (each pipeline's header says which; for Align,
deleting the outlier raster `.flt` reruns it).

## 5. Write the interpretation and render

When the driver finishes, read these, then write `<report_dir>\narrative.md` (Markdown, no title
-- the report supplies the heading):

- `align_run\align_run.csv` (tiles, prior scale, `finest_chosen`, iterations and convergence of
  each loop, starting/ending RMSE), `align_iterations.csv` (per iteration: RMSE, change of the
  aligned DEM P50/P99/max overall and steep, P99 ratio, predicted remaining change, nodes held),
  `align_select.csv` (held-out gain ± SE per grid and prior scale; `best`, `chosen`),
  `align_nodes.csv` (node shifts in mm relative to the global solution, with x/y);
- `<outbins>_fitFinal.csv` and `_finalCoef.csv` (the fence), `<outbins>_densityLod.csv` (density
  LoD bins and fits) and `_densityLodCorr.csv` (noise correlogram, SDz);
- LShunter's patch table and `_null.csv`, and its log in `<report_dir>\logs\` (candidate counts,
  sign-flip and simulated false-discovery counts);
- the LSstats csv files.

Cover, briefly:

- the tiles and prior Align chose and how clear the choice was (the chosen gain against the best,
  in SE); if `finest_chosen`, say that finer tiles weren't tried and that `min_tile_size` set the
  limit; whether each loop converged (one that hit `iterations` deserves a sentence);
- RMSE before and after (it is measured inside each run's own fence, so it can't be compared
  with other runs), and the level of detection: how the loss limit varies with slope, canopy and
  ground-return density;
- the spatial pattern of node shifts, and the nodes held by their neighbours (too few points);
- how many candidate landslides LShunter found, their sizes, depths and volumes ± SD, how many
  gain patches the sign-flip control found (false-alarm level), and which landforms hold them;
- caveats a reader needs (e.g. LSstats used a landform raster other than RIL's output).

Every number must come from those files -- quote, don't estimate. Then rerun the driver: it skips
every pipeline and re-renders the report, which includes the narrative.

Give the user the report path (`<report_dir>\<study_name>_report.html`). It is self-contained
and can be emailed as is; offer to publish it as a private claude.ai page if they want a link.

## 6. If PipelineReport.qmd is missing

The report template is gitignored, so a fresh checkout has none, and the driver's render then
fails. Before the first run on a machine, check for `PipelineReport.qmd` at the repo root; if it
isn't there, write it, then show it to the user:

- Quarto html, `embed-resources: true`, one parameter `manifest` (the driver passes
  `-P manifest:<report_dir>/pipeline_manifest.txt`). Read the manifest's `key: value` lines
  (split at the first colon) and use only the paths it lists -- it's written by
  `R/run_landslide_pipelines.R`, whose `manifest <- list(...)` is the authoritative key list.
- Sections: study area and inputs (`ref_dtm`, `align_dtm`, DSMs, `clip_polygon`, `target_crs`);
  alignment (the `align_run_dir` tables: run summary, convergence by iteration, the
  cross-validation table, a map of node shifts); difference maps (`out_dif_0` and `out_dif`,
  same colour scale, plus `out_exceed` if set); level of detection (`out_lod` and its `_iqr`, and
  the `<out_bins>_densityLod.csv` fits); candidate landslides (`patch_raster` over a hillshade of
  `ref_dtm`, the `lshunter_table` summary); landforms (`landform_raster`, the LSstats csvs);
  Interpretation (include `narrative_md` if it exists, else a note that it hasn't been written);
  run details (`created`, `config`, `elapsed_min`, `log_dir`).
- Each optional input (`NOFILE` in the manifest, or a file that doesn't exist) gets a one-line
  note instead of its figure or table, so the report renders for every configuration.
- Downsample rasters for display (`terra::aggregate`) so the html stays emailable.
