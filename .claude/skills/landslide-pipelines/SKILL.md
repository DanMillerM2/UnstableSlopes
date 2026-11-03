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
| 3. LShunter + LSstats | `R/lshunter_lsstats_pipeline.R` | patch raster, patch table (`out_table`, `_null.csv`), `<out_csv>_<i>.csv`, `<out_csv>_landforms_<i>.csv`, and (`R/patch_polygons.R`) `patch_<tag>_polygons.shp` with `_fields.csv`; with `lshunter_deposits: TRUE`, also `ang_<ref>` (bldgrds), `corridor_<tag>`, `patch_<tag>_deposits`/`_unlinked` rasters, `patch_<tag>_table_deposits.csv`/`_unlinked.csv`/`_balance.csv`, and the deposits in `patch_<tag>_polygons.shp` (field `type`: scar or deposit; a deposit's `scar` field links it) |
| Report | `PipelineReport.qmd` (repo root, gitignored -- see section 6) | one self-contained html file |

**Everything for a site goes in one site file** (`keyword: value` lines; template
`R/site_template.txt`). `R/run_site.R <site file>` writes each pipeline's config into
`<site_dir>\<site_run_name>\configs\`, deriving every link between the programs (Align's
output names, LShunter's inputs from Align and RIL, LSstats's from LShunter), then runs
`R/run_landslide_pipelines.R`, which runs 1 and 2 side by side, then 3, then renders the report.
**Run `run_site.R`; don't write the pipeline configs by hand or reimplement their steps.** Every
pipeline skips outputs that already exist, so rerunning it is how you resume after a failure and
how you re-render after editing the narrative.

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

If the user gives a site file, use it: check it with `--configs-only` (section 3) and go on. If
there's a site file from a previous run of this site, start from that. Otherwise ask for anything
not already given, in one message:

- **Study name**, a **site folder**, and a **run name** (the outputs' folder under it).
- **Study-area polygon** (any CRS) and the **target CRS** (default EPSG:2856).
- For each of the reference (older) and align (newer) lidar: a **label** (e.g. the year), the
  **DTM tile folder** (or one DTM file), the **LAZ tile folder**, and the **DSM tile folder** (or
  file). Both LAZ sets are required (they give the density LoD and the exceedance raster), and so
  are both DSMs (Align uses canopy height for its canopy classes, stable ground and the LoD).
- The **complete DTM for RIL** (ESRI GridFloat `.flt`, named `<prefix>_<ID>.flt` -- the programs
  derive their data ID from the text after the first `_`).
- **RIL-chain input files** to use as templates: the bldgrds (enforce), ValleyFloor, and RIL
  Fortran input files from a previous site. `ril_pipeline.R` runs in mode 1 only, so these carry
  the method's parameters; you will copy them and change only the paths.
- **LShunter and LSstats settings**: S thresholds, minimum gradients and size, the
  flow-accumulation raster, optional road layer and distance-to-channel raster, the null
  simulation and sign-flip control, and the rasters LSstats summarizes. Don't ask for the
  gradient, landform, exceedance or patch inputs: `run_site.R` fills those in.
- The **executable folders** (`Projects\<Name>\x64\release\` under the sibling repos), if they
  differ from the templates' defaults.

## 2. Check prerequisites

Run these and fix anything that fails before writing configs:

```
Rscript -e "a <- names(formals(TWutils::align_input)); stopifnot(all(c('tileSelect','minTileSize','outGradient','densityLod','outExceed','outCenDif','matchSmoothness','suspectGround','exceedanceTrim') %in% a)); l <- names(formals(TWutils::LShunterInput)); stopifnot(all(c('Exceedance','Gradient','OutTable','SignFlip','CanopyRefDSM','AreaWindows') %in% l)); stopifnot(exists('align_gradient_raster', asNamespace('TWutils')), 'patch_table' %in% names(formals(TWutils::LSstats_input))); for (p in c('terra','sf','processx','ggplot2','knitr')) stopifnot(requireNamespace(p, quietly = TRUE))"
quarto --version
```

- If TWutils fails the check, reinstall it from `..\TWutils` (`R CMD INSTALL`; current `main`,
  c2290a4 or later). An older TWutils is also caught when Align starts: the supervisor stops with
  "The installed TWutils is older than this pipeline" and names the missing arguments.
- If `quarto` isn't found, the report can't render. RStudio has a copy that isn't on the PATH
  (`C:\Program Files\RStudio\resources\app\bin\quarto\bin`); ask the user to add it to the PATH
  (or install Quarto), or set `render_report: FALSE` in the site file until then.
- **Align** must be built from GridUtilities `master`, 8a08fa6 or later (`Projects\align`,
  Release|x64); the options `align_matchSmoothness`, `align_suspectGround` and
  `align_exceedanceTrim` need 2083edb, 78b4911 and ab662ea or later respectively. An older
  `align.exe` is caught a few seconds into the run -- the supervisor stops with "Align ... is
  older than this pipeline" and names what it didn't recognize. **LShunter**
  must be 1.1.0 or later (LandslideUtilities `master`, e37b98c or later; its log's first lines show
  the version), and 1.4.0 or later with `lshunter_deposits: TRUE` (LandslideUtilities branch
  `deposits` until merged), and so must **LSstats** (patch volumes; earlier versions crash on a raster path
  containing "DoD"). The fix for any of them is to rebuild in Visual Studio. Tell the user; don't try to
  work around it.
- **Memory**: Align holds whole DEMs in memory (about 48 GB on a 370 M-cell site). Run one site at
  a time, and don't start other Align jobs while it runs.

## 3. Write the site file

1. **RIL-chain input files.** For a new site, copy the three template input files into the site
   folder, and in each change the `DEM:` (bldgrds: `DEM FILE:` too) line to the complete DTM,
   `SCRATCH DIRECTORY:` to the scratch folder, and every output path (e.g. `OUTPUT RIL RASTER:`)
   to the site folder. All three must name the same DEM. The method's parameters stay as they are.
2. **Site file.** Copy `R/site_template.txt` (or the site's previous file) to
   `<site_dir>\<site name>_<ref label>_<align label>.txt` and fill it in; keep its comments.
   Keywords are prefixed by what they configure: `site_`, `ref_`/`align_` for the surveys,
   `ground_` (ground_returns.R), `align_` (any `align_supervisor.R` setting), `ril_`,
   `lshunter_` and `lsstats_`. Leave out anything at its default: Align's defaults (`tileNx`
   auto, 8 iterations, convergence 0.005, `CENTRAL`, 10 slope bins to 1.5, density LoD) are the
   recommended setup, and LShunter's S thresholds default to -4 / -2. In `lsstats_rasters`,
   `@cendif`, `@exceed`, `@dif`, `@dif0`, `@outlier`, `@lod`, `@gradient` and `@aligned` stand for
   this run's Align outputs -- prefer `@cendif` and `@exceed` over older difference rasters.
3. **Check it** without running anything:

   ```
   Rscript R/run_site.R <site file> --configs-only
   ```

   It stops on an unknown keyword, a missing input or executable, or an unknown `@` token, and
   notes when the RIL chain uses a different DEM from the reference DTM. Then show the user the
   site file and the generated `configs\align_pipeline.txt` `[align]` section and
   `configs\lshunter_params.txt`.

## 4. Run and monitor

Start it in the background (Align alone takes about 25 min on a 370 M-cell site; the
ground-return counts and mosaics can take hours):

```
Rscript R/run_site.R <site file>
```

Every output goes to the run folder `<site_dir>\<site_run_name>\` (the report to its `report\`
unless `report_dir` is set). Logs go to `<report_dir>\logs\` (`align_pipeline.log`,
`ril_pipeline.log`, `lshunter_lsstats_pipeline.log`); Align's own console output is in
`<run folder>\align_run\align.log`. Don't poll in a tight loop -- wait for the completion
notification, checking the logs at long intervals. Give the user one line at each Align milestone:
the global loop converged (iterations, RMSE), the tiles chosen (`Selected: tiles nx x ny, prior
scaled by s`, and whether it says the choice is the finest grid tried), the tile loop converged.

On failure, read the end of the failing log and report the cause. Rerunning resumes from what
exists. To change a setting, edit the site file (never the generated configs: `run_site.R`
rewrites them on every run). **Ask before deleting any output** -- after a parameter change, the
outputs it affects must be deleted for the change to take effect (each pipeline's header says
which; for Align, deleting the outlier raster `.flt` reruns it).

## 5. Write the interpretation and render

When the run finishes, read these (all in the run folder), then write `<report_dir>\narrative.md`
(Markdown, no title -- the report supplies the heading):

- `align_run\align_run.csv` (tiles, prior scale, `finest_chosen`, iterations and convergence of
  each loop, starting/ending RMSE), `align_iterations.csv` (per iteration: RMSE, change of the
  aligned DEM P50/P99/max overall and steep, P99 ratio, predicted remaining change, nodes held),
  `align_select.csv` (held-out gain ± SE per grid and prior scale; `best`, `chosen`),
  `align_nodes.csv` (node shifts in mm relative to the global solution, with x/y);
- `<outbins>_fitFinal.csv` and `_finalCoef.csv` (the fence), `<outbins>_densityLod.csv` (density
  LoD bins and fits) and `_densityLodCorr.csv` (noise correlogram, SDz);
- LShunter's patch table `patch_<tag>_table.csv` and `_null.csv`, and its log in
  `<report_dir>\logs\` (candidate counts, sign-flip and simulated false-discovery counts);
- with deposits, `patch_<tag>_table_deposits.csv` (each deposit's scar, volume ± SD, travel
  distance beyond the last coherent loss, drop, reach angle, channel and valley-floor shares),
  `_unlinked.csv` (gain patches with no scar above them: many of these means the gain thresholds
  are low, or gains are vegetation or fill), and `_balance.csv` (per scar: its volume, its
  deposits' volume, the net with SD from the correlated noise, the deposit/scar ratio, and the
  scar its runout enters);
- the LSstats csv files `lsstats_<tag>_<i>.csv` (with each polygon's Volume, SD_Volume and
  Z_Volume) and `_landforms_<i>.csv`;
- the noise reference in `strip_noise\` (unless `noise_run: FALSE`): each survey's
  `<label>_noise_summary.txt` (pass method, and why others were rejected), `_classes.csv`,
  `_density_classes.csv` and `_pairs.csv` (strip offsets), and the report's pass-noise floor
  table.

Cover, briefly:

- the tiles and prior Align chose and how clear the choice was (the chosen gain against the best,
  in SE); if `finest_chosen`, say that finer tiles weren't tried and that `min_tile_size` set the
  limit; whether each loop converged (one that hit `iterations` deserves a sentence);
- RMSE before and after (it is measured inside each run's own fence, so it can't be compared
  with other runs), and the level of detection: how the loss limit varies with slope, canopy and
  ground-return density;
- the spatial pattern of node shifts, and the nodes held by their neighbours (too few points);
- the noise reference: which pass method each survey used (say plainly when it is none -- the LoD
  then rests on stable ground alone), strip offsets (above about 0.1 m is a systematic strip
  error), and the pass-noise floor's share of the observed variance by canopy and density: the
  floor is each survey's own noise, a lower bound and not a prediction of the DoD noise, so don't
  read agreement between them as validation. Where its share is small, the level of detection is
  set by differences between the surveys (classification, gridding, alignment) and real change,
  which more ground returns wouldn't fix;
- how many candidate landslides LShunter found, their sizes, depths and volumes ± SD, how many
  gain patches the sign-flip control found (false-alarm level), and which landforms hold them;
- caveats a reader needs (e.g. LSstats used a landform raster other than RIL's output).

Every number must come from those files -- quote, don't estimate. Then rerun
`Rscript R/run_site.R <site file>`: it skips every pipeline and re-renders the report, which
includes the narrative.

Give the user the report path (`<report_dir>\<study_name>_report.html`). It is self-contained
and can be emailed as is; offer to publish it as a private claude.ai page if they want a link.
Also give the patch polygons (`<run folder>\patch_<tag>_polygons.shp`, fields described in its
`_fields.csv`) for use in a GIS.

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
  the `<out_bins>_densityLod.csv` fits); noise reference from flight-pass overlaps (below);
  candidate landslides (`patch_raster` over a hillshade of
  `ref_dtm`, the `lshunter_table` summary, and a line naming `patch_polygons` with its feature
  count and `_fields.csv`); landforms (`landform_raster`, the LSstats csvs);
  Interpretation (include `narrative_md` if it exists, else a note that it hasn't been written);
  run details (`created`, `config`, `elapsed_min`, `log_dir`).
- Noise reference section, from `strip_noise_dir` (`NOFILE` when `noise_run: FALSE`: one note),
  with files named `<label>_noise_*` for `noise_ref_label` and `noise_align_label`. Say first that
  each survey's noise is measured between its own flight passes (no real change), so it is a floor
  under the level of detection, not a prediction of it. Per survey: the `_summary.txt` lines as
  bullets (pass method, why others were rejected, pass checks); a strip-offset table from
  `_pairs.csv` (pairs, cells, median / 90th percentile / max of |`median_m`|), noting a
  systematic strip error when the max exceeds 0.1 m; then one table merging both surveys'
  `_classes.csv` (cells, `iqr_m`, `gt10`, `gt20` by canopy and slope, in the csv's row order).
  Then the pass-noise floor table, which needs both `_density_classes.csv` files, `out_cendif`,
  `<out_lod>_iqr`, both DTMs and DSMs and `ref_ground_count`/`align_ground_count`:
  - sample `out_cendif` (about 1e6 regular cells) and keep cells whose canopy class (DSM - DTM,
    breaks read back from the class labels) is the same in both surveys;
  - each survey's ground returns per m2 at the cell: its count raster summed over 5 x 5 cells,
    divided by their area;
  - per survey and canopy x slope class, fit `iqr_m^2 = c0 + c1 * inv_density_sum` over
    `_density_classes.csv` rows with at least 5000 cells (weights = cells); one pass's noise is
    `sqrt(max(c0/2, 0) + max(c1, 0) / d)`, and the floor is both surveys' values combined in
    quadrature;
  - table by canopy class and the sparser survey's density (< 0.25, 0.25-0.75, 0.75-2, 2-4,
    >= 4 returns/m2; groups of at least 300 cells): cells, observed IQR, median of Align's
    `_iqr`, median floor, and the floor's share of the observed variance (floor^2 / IQR^2,
    at most 1);
  - a note that the remainder is what differs between the surveys (classification, gridding,
    sensor, season, alignment) plus real change, and grows as ground returns thin; where the
    floor's share is small, more points wouldn't lower the level of detection.
  Keep the comparisons' filters free of NA (test both canopy classes for NA before comparing
  them), or the grouping shifts silently.
- Each optional input (`NOFILE` in the manifest, or a file that doesn't exist) gets a one-line
  note instead of its figure or table, so the report renders for every configuration.
- Downsample rasters for display (`terra::aggregate`) so the html stays emailable.
