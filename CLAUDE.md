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
   GeoPackage export, count/density rasterization, mosaic, reprojection to EPSG:2856). The
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
- **run_align_pipeline.R** — driver: runs `ground_returns.R` once each for a "reference" and an
  "align" set of LAZ tiles to build the two ground-point-density rasters `align_dtms.R`'s Align call
  needs, then calls `TWutils::align()` directly with those two rasters as `refGrnd`/`alignGrnd`.
- **ground_returns.R** — assembled from `SusNRun_PointCloudWorkflow_blank_082626.qmd` (see the
  branch caveat above); the batched lidar point-cloud pre-processing workflow described in the
  document map (ground filtering → optional GeoPackage → per-batch count raster → mosaic →
  reprojection to EPSG:2856).
- **run_ground_returns_batches.R** — driver: runs `ground_returns.R` sequentially, once per named
  set in a batch config file, each as its own `Rscript` subprocess (isolating each set's
  `future::plan` and global-environment state rather than `source()`-ing repeatedly in one session).
- **RIL.R** — runs `TWutils::RIL()` (Fortran program RIL, see `..\ChannelUtilities\RIL\RIL.f90`) to
  trace the channel network from a DEM and classify valley floor / hollow / inner-gorge landforms.
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
