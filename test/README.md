# Consensus-clustering test suite

## Purpose

The suite validates the canonical workflow from raw cognitive-data preparation
through base clustering, five-method consensus, cross-validation, held-out
prediction, resampling inference, holdout assignment calibration, and
known-truth simulation.
All fixtures are deterministic and synthetic. Tests write only to temporary
directories unless explicitly generating `test/test_report.md`.

## Running the tests

Fast suite from the repository root:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' test\run_tests.R
```

Execute every test file individually, including the slow known-truth pipeline,
and regenerate the validation report:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' test\generate_test_report.R
```

The second command sets `RUN_SLOW_TESTS=true` for the report run and returns a
nonzero exit code after writing the report if any validation fails.

Dependencies are restored from `renv.lock`. The validated toolchain is R 4.6.1
with `renv` activation through the repository `.Rprofile`.

## Test infrastructure

- `run_tests.R` runs all `test-*.R` files through `testthat::test_dir()`; the
  slow pipeline is skipped unless `RUN_SLOW_TESTS=true`.
- `generate_test_report.R` runs every test file separately with the slow flag
  enabled and writes per-file and per-case results to `test_report.md`.
- `helper-fixtures.R` generates complete raw cognitive data and temporary Excel
  workbooks for loader and preprocessing tests.
- `helper-simulations.R` generates independent clustered feature-space data,
  known labels, and reusable module-loading helpers.
- `baseline.md` records the original pre-refactor baseline and environment
  history; it is not an executable test.

## Test files and validation intent

### `test-base-clustering.R`

- Verifies the deterministic structured output contract, participant order,
  assignments, fitted metadata, selected cluster counts, and diagnostics for
  all requested base methods.
- Rejects duplicated/misordered IDs and infeasible candidate cluster counts.
- Confirms HDBSCAN label `0` is represented as unavailable (`NA`) instead of a
  regular cluster.
- Confirms SPEC base HCA uses average linkage/UPGMA while compatibility settings
  remain explicit.

### `test-consensus-clustering.R`

- Checks pairwise-available co-association numerators and denominators when
  method assignments contain missing values.
- Validates explicit handling and diagnostics for participant pairs with no
  method evidence.
- Confirms compatibility mode retains its historical three-method default.
- Confirms SPEC mode includes HCA, spherical k-means, GMM, affinity propagation,
  and HDBSCAN, with Ward D2 applied only to consensus distance.

### `test-cv-stability.R`

- Verifies stable labels require at least three matching assignments among the
  four training folds; 2/4 agreement and ties remain non-clustered.
- Checks deterministic modal-cluster-count reference-fold selection.
- Validates method-level fold assignments, overlap ARI, predictions, accuracy,
  coverage, stable labels, and consensus metrics returned by SPEC CV.
- Checks the per-participant, per-method fold-label table contains every
  aligned training-fold label plus final cross-fold stability diagnostics.

### `test-distance-matrix.R`

- Checks hand-calculated two-dimensional cosine task similarity.
- Validates the standardized Pearson-product fallback when only performance or
  reaction time is structurally available, plus zero for no shared dimension.
- Confirms zero task vectors retain self-similarity one and similarity zero to
  nonzero vectors.

### `test-folds-and-label-alignment.R`

- Checks nearly equal, reproducible folds balanced over sex, antidepressant
  use, and age bins.
- Validates Hungarian label alignment under a pure label permutation.
- Confirms a CV fold retains training and held-out assignments for every base
  method and consensus.

### `test-full-synthetic-pipeline.R`

- Runs the five-fold SPEC workflow with all five algorithms on ten paired task
  measures transformed into 45 task-similarity features.
- Tests final-development, stable-CV, and held-out recovery against independent
  known labels using ARI and coverage thresholds.
- Evaluates every 2-to-5-method consensus subset and verifies the oracle best
  subset is at least as accurate as the all-five subset.
- Validates the generated synthetic Markdown result, including method ranking
  and participant-level recovered labels.
- This test is opt-in for the fast runner and always runs in the generated test
  report.

### `test-gam-imputation.R`

- Checks fixed-seed chained GAM imputation is deterministic, bounded by observed
  target ranges, complete, and leaves observed values unchanged.
- Verifies the audited observed-value resampling fallback when too few target
  observations are available for a GAM fit.

### `test-generalization.R`

- Validates cohort-aligned reference and holdout proportions, simultaneous
  intervals, assigned coverage, exclusions, exact multinomial calibration, and
  reproducible Monte Carlo fallback diagnostics.

### `test-legacy-surfaces.R`

- Ensures superseded global scripts, preparation notebooks, and report artifacts
  remain deleted so obsolete analytical paths cannot silently return.

### `test-load-and-prepdata.R`

- Runs the public loader on a generated Excel workbook.
- Checks that cohort holdout sampling occurs after the missingness screen and
  that the realized holdout count equals the requested fraction of screened
  eligible participants.
- Checks participant exclusion, development/holdout disjointness, finite
  distance output, exactly 45 task-pair features, audit metadata, and temporary
  TSV/JSON artifacts.

### `test-prediction.R`

- Uses holdout sentinels to prove fitted feature preprocessing learns only from
  training participants.
- Checks the loader exposes disjoint raw development and untouched holdout rows.
- Verifies raw imputation values are training-derived and unaffected by holdout
  extremes.
- Exercises HCA, spherical k-means, GMM, AP, and HDBSCAN prediction adapters in
  the fitted 45-feature space, followed by fixed-partition consensus prediction.
- Checks pairwise-available held-out co-association, prediction scores, margins,
  contributing-method counts, rejection diagnostics, GMM posterior information,
  and HDBSCAN noise handling.

### `test-resampling.R`

- Confirms each two-variable cognitive task block uses one shared donor
  permutation and therefore preserves its within-task performance/RT pairing.
- Checks participant bootstrap shape, replacement behavior, source IDs, and
  unique execution IDs.
- Validates plus-one empirical p-values, the documented cluster-count tail, and
  2.5%/97.5% percentile confidence-interval schemas.
- Checks deterministic, non-overlapping per-iteration seeds.
- Runs two-worker PSOCK resampling and verifies identical permutation,
  bootstrap, and inference tables against sequential execution.
- Runs a reduced parallel resampling iteration through the complete
  cross-validation, base-clustering, prediction, and consensus call chain.

### `test-workflow-integration.R`

- Confirms the default output directory remains the repository-root `results/`
  directory even when the caller changes the current working directory.
- Runs the canonical compatibility entry point and checks structured RDS/CSV
  artifacts, including the `cv_final_labels.csv` schema without spurious
  serialization columns, and relational output keys.
- Runs SPEC mode with all five methods and a fixed holdout, verifying active
  methods, assigned/excluded accounting, simultaneous proportion intervals,
  and the renamed `generalization_calibration.csv` artifact.

## Reading `test_report.md`

The report records the tested Git commit, R/testthat versions, overall result,
per-file cases and assertion counts, elapsed time, and every named validation
case. A `PASS` means no failed expectation or error occurred. The synthetic
pipeline is not marked as skipped in this report because the generator enables
the slow-test flag.
## Custom-data adapter checks

`test-custom-data.R` covers arbitrary categorical balancing, optional age
balancing, unchanged legacy fold assignments, configurable age bands,
training-only transformations, and the custom loader's inference dispatch.
The inference-dispatch check substitutes a single observed evaluation; it does
not execute permutation or bootstrap resampling.

From the repository root, a separate synthetic usage example is available:

```r
source("src/custom_data_example.R")
example_result <- run_custom_data_example()
```

This example bypasses resampling and writes no outputs. See
`src/custom_data.md` for adapter input/output contracts and distance compatibility.

