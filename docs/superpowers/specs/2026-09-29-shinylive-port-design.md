# PACNet web app: serverless port (Shinylive + GitHub Pages)

Date: 2026-09-29
Status: approved to proceed (user waived section-by-section review)

## Goal

Replace the AWS-hosted Shiny server deployment of the PACNet web app with a
static site on GitHub Pages that runs the existing R code in the visitor's
browser (webR, via Shinylive). The same `webApp/` directory must keep working
with `shiny::runApp()` / `shiny::runGitHub()` for anyone with R installed.

### Success criteria

1. `https://pcahan1.github.io/PACNet/` serves the app; no server-side compute.
2. Feature parity with the current app:
   - "CellNet analysis": upload metadata + expression CSVs, pick tissue,
     classification heatmap, engineered-reference comparison heatmap,
     automatic retraining when the query lacks classifier genes (with the
     "New Training Validation" tab), NIS plots, CSV downloads.
   - "Browse reference data": engineered-reference heatmap per tissue.
3. Results are identical to the current app for classification of query and
   reference-panel samples and for NIS, for every tissue (verified against
   golden outputs produced by the original code on the original data).
4. The retraining path uses `expTrain` rounded to 4 significant digits.
   Retrained results must agree with the original-data retrain on every
   query top call; score differences stay within the random-forest
   seed-to-seed noise measured in the spike (max |d| ~0.14).
5. `shiny::runApp("webApp")` works in native R with only CRAN packages; the
   CellNet and cancerCellNet packages are no longer required.
6. The AWS deploy script is removed and the README points to the new site.

### Non-goals

- UI redesign. The existing `www/index.html` template and layout stay.
- GRN status analysis (already commented out in the current app).
- Mouse support, new tissues, new features.
- Changing any numerical method.

## Background (spike results, 2026-09-28/29)

- Every CRAN dependency has a webR (wasm) binary. CellNet/cancerCellNet are
  pure R; the app uses 32 functions from cancerCellNet only.
- In webR under Node, the no-retrain analysis (liver) runs in ~25 s, peak R
  heap 590 MB. Retraining failed on memory as written; with two patches
  (free `expTrain` before training; slice the quick-pairs scoring) it runs
  in ~6.5 min, peak R heap 815 MB.
- `parallel::detectCores()` returns `NA` in wasm, breaking three functions.
- Reference data shrinks from 747 MB to ~97 MB (see Data format).
- GitHub release assets do not send CORS headers; S3 bucket does not either.
  Data is therefore served from the Pages site itself.

## Architecture

```
PACNet/
├── webApp/                      Shiny app dir; exported by shinylive (code only)
│   ├── ui.R                     unchanged except removal of dead code paths
│   ├── server.R                 Shiny wiring; calls pipeline.R, no raw data I/O
│   ├── pipeline.R               NEW: pure analysis functions (no Shiny)
│   ├── cnlib.R                  NEW: vendored cancerCellNet subset + patches
│   ├── data_access.R            NEW: pacnet_load(): local file or URL fetch
│   ├── pacnet_utils.R           pacnet_nis() adapted to slim NIS object
│   ├── plotting.R               unchanged
│   └── www/                     unchanged template/assets
├── data/                        NEW: slim .rds reference files (committed)
├── tests/
│   ├── testthat/                NEW: regression tests (native R)
│   └── golden/                  NEW: outputs of the ORIGINAL code on ORIGINAL data
├── tools/
│   └── make_slim_data.R         NEW: original S3 .rda -> data/*.rds + tests/golden
└── .github/workflows/
    ├── pages.yml                NEW: shinylive export + data -> GitHub Pages
    └── tests.yml                NEW: run tests/testthat in native R
```

Data lives outside `webApp/` because `shinylive::export()` embeds every file in
the app directory into `app.json`, which is downloaded at startup.

The published site is assembled by the Pages workflow:

```
site/                 <- shinylive::export("webApp", "site")
site/data/*.rds       <- copied from data/
```

## Units

### `webApp/data_access.R`

`pacnet_load(name)` is the only way app code reads reference data. It
returns the object stored in `name` (an `.rds` file name).

- Native R: `readRDS(file.path(getOption("pacnet.data_dir", "../data"), name))`.
  `runApp("webApp")` and `runGitHub(..., subdir = "webApp")` both run with
  `webApp/` as the working directory, so `../data` resolves.
- webR (`R.version$os == "emscripten"`): download
  `paste0(pacnet_data_url(), name)` to `tempdir()` once per session, then
  `readRDS`. `pacnet_data_url()` returns `getOption("pacnet.data_url")` if set;
  otherwise it derives the site root from the webR worker's own URL
  (`webr::eval_js("self.location.href")`, stripping `/shinylive/...`) and
  appends `data/`. This makes the same export work at
  `pcahan1.github.io/PACNet/`, on a fork, and on a local static server.
- Download failures raise an R error with the URL in the message; the server
  shows it via the existing `shinyalert` error path.

### `webApp/cnlib.R`

Single file holding the 32 cancerCellNet functions in the app's call closure
(`zscore` excluded; `pacnet_utils.R` defines it). Header notes the source
(pcahan1/cancerCellNet, MIT) and the commit it was copied from. Changes from
upstream, each marked with a `# PACNet-web:` comment:

1. Parallel code paths removed (`detectCores`, `snow::` cluster branches);
   everything runs serially. Upstream already defaults to serial
   (`coreProportion = 0`).
2. `ptGetTop()` quick-pairs branch scores gene pairs in slices of
   `sliceSize` and row-binds the per-slice statistics before
   `findBestPairs()`. `findBestPairs()` orders by `cval` only, so results
   are unchanged; peak memory drops from ~235 MB per class to ~26 MB.

### `webApp/pipeline.R`

Pure functions, no Shiny calls, used by both `server.R` and the tests:

| Function | Purpose |
|---|---|
| `normalize_query(expDat, sampTab)` | reorder columns to `sampTab$sample_name`, `downSampleW(1e5)`, `log1p` |
| `load_classifier(tissue)` | classifier + iGenes; `none_of_the_above` handled |
| `needs_retrain(broadReturn, queryExp)` | `!all(cgenes %in% rownames(queryExp))` |
| `retrain_classifier(queryExp, iGenes)` | loads `Hs_expTrain_sig4`/`Hs_stTrain`, splits (seed 39), frees the full matrix before training, trains with the current parameters, returns `list(broadReturn, classMatrix, grps, assessment)` for the validation tab |
| `classify(broadReturn, expDat, nrand)` | `broadClass_predict` wrapper |
| `load_engineered_ref(tissue)` | rank matrix + sample table (column-aligned) + top-TF order |
| `run_nis(queryExp, sampTab, iGenes, tissue)` | `logRank`, load `<tissue>_nis.rds`, subset to `iGenes`, `pacnet_nis()` |

`server.R` keeps all UI logic, progress updates, alerts, rendering and
downloads, and replaces inline data loading/subsetting with these calls.
The ~60-line GRN/trainNormParam subsetting block is removed; subsetting the
slim NIS object is a few lines inside `run_nis()`.

### `webApp/pacnet_utils.R`

`pacnet_nis()` takes the slim object instead of `grnAll` + `trainNormParam`:
`pacnet_nis(expDat, stQuery, iGenes, nis, colname_sid, relaWeight = 0)`.
Arithmetic is unchanged. `subsetGRNall()` / `subsetTrainNormParam()` are
deleted.

## Data format (`data/`)

Produced by `tools/make_slim_data.R` from the original S3 objects. All files
are `saveRDS(..., compress = "xz")`.

| File | Source | Transformation | Exact? |
|---|---|---|---|
| `<T>_broadClassifier100.rds` | same name `.rda` | none | yes |
| `<T>_studies_iGenes.rds` | same | none | yes |
| `<T>_engineeredRef_sampTab_all.rds` | same | none | yes |
| `<T>_top_TF_display_order.rds` | same | none | yes |
| `<T>_engineeredRef_ranks.rds` | `<T>_engineeredRef_normalized_expDat_all.rda` | columns aligned to sample table; integer rank within each sample (`ties.method = "min"`) | yes: classification uses only within-sample gene comparisons |
| `<T>_nis.rds` | `<T>_grnAll.rda` + `<T>_trainNormParam.rda` | `list(tfTargets = ctGRNs$tfTargets[[T]], netGenes = names(ctGRNs$geneLists[[T]]), tVals = list(mean, sd) for netGenes)` | yes: the only fields NIS reads |
| `none_of_the_above_broadClassifier100.rds`, `..._studies_iGenes.rds` | same | none | yes |
| `Hs_stTrain.rds` | `Hs_stTrain_Jun-20-2017.rda` | none | yes |
| `Hs_expTrain_sig4.rds` | `Hs_expTrain_Jun-20-2017.rda` | `signif(x, 4)` | no (user-approved) |

`T` ∈ {heart, hspc, intestine_colon, liver, lung, neuron, skeletal_muscle}.

Dropped: `*_GRN_statusQuery` (loaded in browse mode but never displayed),
`none_of_the_above_grnAll` and `none_of_the_above_trainNormParam` (NIS is
skipped for that choice).

Expected sizes: ~48 MB total excluding `expTrain`, plus 49 MB for
`Hs_expTrain_sig4.rds`. Largest per-tissue file ~14 MB.

## Golden outputs (`tests/golden/`)

`tools/make_slim_data.R` also runs the ORIGINAL pipeline (upstream
cancerCellNet functions, original `pacnet_nis`, full `grnAll`, original
double-precision reference data) on `example_data/` and saves:

- per tissue (7 + `none_of_the_above`): query classification matrix (real
  sample columns only; random-profile columns are stochastic);
- per tissue (7): reference-panel classification matrix (real columns);
- per tissue (7): NIS matrix, plus one extra run with 25 network genes
  removed from `iGenes` to exercise the subsetting path;
- one liver retrain from the original `expTrain` (query + reference
  classification, real columns).

## Testing

Native-R testthat suite (`tests/testthat/`), run in CI on every push:

- `test-data.R`: every expected file loads through `pacnet_load()`; basic
  shape/type checks (ranks are integer; NIS gene sets consistent).
- `test-classify.R`: for each tissue, `classify()` on slim data equals
  golden (`expect_equal`, real columns).
- `test-nis.R`: for each tissue, `run_nis()` equals golden, including the
  subsetting case.
- `test-retrain.R` (only when `PACNET_SLOW=1`; CI sets it):
  `retrain_classifier()` on the example query; query top calls equal the
  golden original-data retrain; max |score diff| < 0.15.

Browser verification (run during implementation; documented in the README
dev notes): serve `site/` locally, and with Playwright drive
(a) CellNet analysis on the example CSVs with tissue = liver, no retrain;
(b) a forced-retrain case (example CSV with classifier genes removed);
(c) Browse reference data for one tissue. Check that heatmaps, NIS plots
and downloads appear and that no console errors are thrown.

## Deployment

`.github/workflows/pages.yml` (on push to `main` touching `webApp/**` or
`data/**`, and manual dispatch):

1. `r-lib/actions/setup-r`; install `shinylive` from CRAN.
2. `shinylive::export("webApp", "site")` (wasm packages bundled).
3. `cp -r data site/data`; `touch site/.nojekyll`.
4. `actions/upload-pages-artifact` + `actions/deploy-pages`.

One-time manual step for the repo owner: Settings → Pages → Source:
"GitHub Actions".

`.github/workflows/tests.yml`: setup-r, install CRAN deps, run
`testthat::test_dir("tests/testthat")` with `PACNET_SLOW=1`.

Removed: `webApp/deploy_agnosticCellNet_web.sh`.

README: web-app link changes to the Pages URL; add "run locally" section
(`install.packages(...)` line + `shiny::runGitHub("PACNet", "pcahan1",
subdir = "webApp")`).

## Risks and mitigations

1. **`htmlTemplate` + bundled jQuery under Shinylive.** The template loads
   its own jQuery before Shiny's head content. If this breaks Shinylive,
   drop the template's duplicate jQuery include; if still broken, fall back
   to `fluidPage()` with the same panels. Verified early in the plan.
2. **`download.file()` inside the Shinylive worker.** Verified first, with a
   minimal Shinylive app, before touching the real app. Fallback: read via
   `url()` connection or `webr::eval_js` + `fetch`.
3. **Retraining memory for other tissues.** Liver was measured.
   `none_of_the_above` has the largest gene set; the plan includes a
   webR retrain run for it.
4. **First load cost.** webR plus packages is tens of MB on first visit
   (cached afterwards). The existing progress/spinner UI remains; README
   notes the first-load delay.
5. **Provenance.** `tools/make_slim_data.R` downloads originals from
   `cnobjects.s3.amazonaws.com`. If the bucket is retired, the originals
   should be archived first (e.g. Zenodo); the README Data table links
   still point at S3 and are unchanged by this work.
