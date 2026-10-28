---
name: landslide-pipelines
description: Run the full landslide-mapping method for one study area -- DTM coregistration with Align (automatic tile/iteration search), RIL landforms, LShunter candidate landslides and LSstats summaries -- then write the interpretation and render the shareable html report. Use when someone wants to process a new site or pair of lidar dates, rerun a site, or regenerate its report.
---

# Landslide pipelines

This repo chains three pipelines and a report:

| Step | Script | Produces |
|---|---|---|
| 1. Align (and its inputs) | `R/run_align_pipeline.R` -> `R/align_supervisor.R` | aligned DTM, DoD, outlier raster, `align_attempts/*.csv` |
| 2. RIL (independent of 1) | `R/ril_pipeline.R` | landform raster (codes 1-10 as LSstats names them: Floodplain, Bedrock_Hollow, Inner_Gorge, Steep_/Convergent, Divergent, Planar, Other) |
| 3. LShunter + LSstats | `R/lshunter_lsstats_pipeline.R` | patch raster, `<out_csv>_<i>.csv`, `<out_csv>_landforms_<i>.csv` |
| Report | `PipelineReport.qmd` | one self-contained html file |

`R/run_landslide_pipelines.R` runs 1 and 2 side by side, then 3, then renders the report. **Run
that script; don't reimplement its steps.** Every pipeline skips outputs that already exist, so
rerunning it is how you resume after a failure and how you re-render after editing the
narrative.

The number of tiles and iterations Align uses is decided by `R/align_supervisor.R`, not by you:
Align stops iterating when the RMSE changes by less than `convergence` (default 1 mm, at most
`iterations`, default 10), and the supervisor tries up to `max_attempts` (default 3) tilings,
keeping a finer one only if it lowers the ending RMSE by at least `min_rmse_reduction` (default
2%). If the user wants different behavior, change those keywords in the `[align]` section --
never edit the attempt outputs by hand. Read the header of `R/align_supervisor.R` once before
your first run.

## 1. Collect the inputs

Ask for anything not already given, in one message:

- **Study name**, and an **output folder** for this site.
- **Study-area polygon** (any CRS) and the **target CRS** (default EPSG:2856).
- For each of the reference (older) and align (newer) lidar: the **DTM tile folder** (or one DTM
  file), the **LAZ tile folder**, and optionally the **DSM tile folder** (or file).
- The **complete DTM for RIL** (ESRI GridFloat `.flt`, named `<prefix>_<ID>.flt` -- the programs
  derive their data ID from the text after the first `_`).
- **RIL-chain input files** to use as templates: the bldgrds (enforce), ValleyFloor, and RIL
  Fortran input files from a previous site. `ril_pipeline.R` runs in mode 1 only, so these carry
  the method's parameters; you will copy them and change only the paths.
- **LShunter and LSstats parameters**: their `*_params.txt` files from a previous site, or the
  values (thresholds, gradient/accumulation rasters, rasters to summarize, landform raster).
- The **executable folders** (`Projects\<Name>\x64\release\` under the sibling repos), if they
  differ from the templates' defaults.

## 2. Check prerequisites

Run these and fix anything that fails before writing configs:

```
Rscript -e "stopifnot('convergence' %in% names(formals(TWutils::align_input)), 'minPoints' %in% names(formals(TWutils::align_input))); for (p in c('terra','sf','processx','ggplot2','knitr')) stopifnot(requireNamespace(p, quietly = TRUE))"
quarto --version
```

- If TWutils lacks `convergence`/`minPoints`, reinstall it from `..\TWutils` (`R CMD INSTALL`).
- Align must be built from `alignSlope2.f90` version 1.0.6 or later (the `CONVERGENCE` keyword).
  An older `align.exe` is caught a few seconds into the first attempt -- the supervisor stops with
  "ignored CONVERGENCE"; the fix is to rebuild `Projects\align` in Visual Studio. Tell the user;
  don't try to work around it.

## 3. Write the configs

Put everything for the site in its output folder. Start each file from the repo template and
change only what the user's inputs require; keep the templates' comments.

1. **Align** -- copy `R/run_align_pipeline_params_template.txt`. Fill `[study_area]`, `[ref_set]`/
   `[align_set]` (LAZ `input_dir`, `output_dir`), `[ref_dtm]`/`[align_dtm]` (`tile_dir` +
   `output_dtm`, or `input_dtm` + `output_dtm`), `[ref_dsm]`/`[align_dsm]` if DSMs were given,
   and `[align]` outputs. Keep `tileNx: auto`, `tileNy: auto`, `iterations: 10`,
   `convergence: 0.001` unless the user says otherwise; set `tile_size_m` only if they have a
   view on the scale of misalignment (default 2000 m).
2. **RIL** -- copy the three template input files, and in each change the `DEM:` (bldgrds:
   `DEM FILE:` too) line to the complete DTM, `SCRATCH DIRECTORY:` to the site's scratch folder,
   and every output path (e.g. `OUTPUT RIL RASTER:`) to the output folder. All three must name
   the same DEM. Then copy `R/ril_pipeline_params_template.txt` and point it at them.
3. **LShunter/LSstats** -- start from the user's params files (or `R/LShunter_params_template.txt`
   / `R/LSstats_params_template.txt`). `outlier` **must** equal `[align] outOutlier` (the driver
   stops otherwise); LSstats `polygon_raster` must equal LShunter `out_patch`. Then copy
   `R/lshunter_lsstats_pipeline_params_template.txt`.
4. **Driver** -- copy `R/run_landslide_pipelines_params_template.txt` and set `study_name` and
   the three config paths.

Show the user the `[align]` section and the paths you changed before running.

## 4. Run and monitor

Start the driver in the background (it can run for hours on large sites):

```
Rscript R/run_landslide_pipelines.R <site folder>\landslide_pipelines.txt
```

Logs go to `<report_dir>\logs\` (`align_pipeline.log`, `ril_pipeline.log`,
`lshunter_lsstats_pipeline.log`); each Align attempt also writes
`<outOutlier folder>\align_attempts\attempt_<n>\align.log`. Don't poll in a tight loop -- wait
for the completion notification, checking the logs at long intervals. When an attempt finishes,
tell the user in one line what it found (tiles, iterations used, ending RMSE, kept or not), from
`align_attempts.csv` once written, or the attempt log before that.

On failure, read the end of the failing log and report the cause. Rerunning resumes from what
exists. **Ask before deleting any output** -- after a parameter change, the outputs it affects
must be deleted for the change to take effect (each pipeline's header says which).

## 5. Write the interpretation and render

When the driver finishes, read `align_attempts.csv`, `align_iterations.csv`, `align_tiles.csv`,
`<outbins>_fitFinal.csv`, and the LSstats csv files, then write `<report_dir>\narrative.md`
(Markdown, no title -- the report supplies the heading). Cover, briefly:

- which tiling was used and why the search stopped there; whether every loop converged (a loop
  that hit `iterations` without converging deserves a sentence);
- RMSE before and after, and the level of detection (the fences in `_fitFinal.csv`) -- how it
  varies with slope;
- the spatial pattern of tile shifts, any tile filled from neighbors, standard errors that are
  large relative to the shifts;
- how many candidate landslides LShunter found, their sizes, and which landforms hold them;
- caveats a reader needs (e.g. LSstats used a landform raster other than RIL's output).

Every number must come from those files -- quote, don't estimate. Then rerun the driver: it skips
every pipeline and re-renders the report, which includes the narrative.

Give the user the report path (`<report_dir>\<study_name>_report.html`). It is self-contained
and can be emailed as is; offer to publish it as a private claude.ai page if they want a link.
