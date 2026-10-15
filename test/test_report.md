# Consensus-clustering test report

## Execution metadata

- Generated: 2026-08-05 20:27:16 UTC
- Git base commit: `3df0d4a`
- Working-tree changes included in this run: `TRUE`
- R: `R version 4.6.1 (2026-06-24 ucrt)`
- testthat: `3.3.2`
- Command: `Rscript test/generate_test_report.R`
- `RUN_SLOW_TESTS=true`: every test file, including the known-truth full pipeline, was executed.
- Inputs: deterministic generated fixtures and temporary files; no production workbook was read or modified.

## Overall result

| result | files | cases | assertions | passed_expectations | failed_expectations_or_cases | skipped_cases | warnings | elapsed_seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PASS | 13 | 36 | 169 | 169 | 0 | 0 | 0 | 80.18 |

## Result by test file

| file | cases | assertions | passed | failed | skipped | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| test-base-clustering.R | 4 | 14 | 14 | 0 | 0 | 0 | 9.56 | PASS |
| test-consensus-clustering.R | 4 | 20 | 20 | 0 | 0 | 0 | 8.14 | PASS |
| test-cv-stability.R | 3 | 11 | 11 | 0 | 0 | 0 | 11.97 | PASS |
| test-distance-matrix.R | 3 | 13 | 13 | 0 | 0 | 0 | 0.57 | PASS |
| test-folds-and-label-alignment.R | 3 | 11 | 11 | 0 | 0 | 0 | 2.21 | PASS |
| test-full-synthetic-pipeline.R | 1 | 14 | 14 | 0 | 0 | 0 | 23.13 | PASS |
| test-gam-imputation.R | 2 | 11 | 11 | 0 | 0 | 0 | 0.97 | PASS |
| test-generalization.R | 3 | 15 | 15 | 0 | 0 | 0 | 0.38 | PASS |
| test-legacy-surfaces.R | 1 | 1 | 1 | 0 | 0 | 0 | 0.17 | PASS |
| test-load-and-prepdata.R | 1 | 11 | 11 | 0 | 0 | 0 | 3.24 | PASS |
| test-prediction.R | 5 | 23 | 23 | 0 | 0 | 0 | 4.98 | PASS |
| test-resampling.R | 3 | 11 | 11 | 0 | 0 | 0 | 0.56 | PASS |
| test-workflow-integration.R | 3 | 14 | 14 | 0 | 0 | 0 | 14.3 | PASS |

## Validation cases

### `test-base-clustering.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| base clustering exposes a deterministic structured contract | 8 | 8 | 0 | 0 | 7.52 | PASS |
| base clustering validates participant ordering and candidate k | 2 | 2 | 0 | 0 | 0.11 | PASS |
| HDBSCAN noise is unavailable rather than a regular cluster | 1 | 1 | 0 | 0 | 0.06 | PASS |
| SPEC mode uses UPGMA for the base HCA | 3 | 3 | 0 | 0 | 0.03 | PASS |

### `test-consensus-clustering.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| co-association uses pairwise available denominators | 7 | 7 | 0 | 0 | 0.16 | PASS |
| no-evidence policy is explicit and diagnosed | 3 | 3 | 0 | 0 | 0.05 | PASS |
| compatibility consensus defaults to the three legacy methods | 6 | 6 | 0 | 0 | 4.01 | PASS |
| SPEC consensus includes all five available methods | 4 | 4 | 0 | 0 | 3.69 | PASS |

### `test-cv-stability.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| stable labels require at least three of four assignments | 3 | 3 | 0 | 0 | 0.08 | PASS |
| reference fold uses the modal cluster count deterministically | 2 | 2 | 0 | 0 | 0.04 | PASS |
| SPEC CV returns method-level assignments and metrics | 6 | 6 | 0 | 0 | 11.66 | PASS |

### `test-distance-matrix.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| two-dimensional task pairs use cosine similarity | 3 | 3 | 0 | 0 | 0.11 | PASS |
| structural one-measure pairs use Pearson-product features | 7 | 7 | 0 | 0 | 0.28 | PASS |
| zero task vectors preserve self-similarity and zero other-similarity | 3 | 3 | 0 | 0 | 0.03 | PASS |

### `test-folds-and-label-alignment.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| balanced folds are nearly equal and reproducible | 4 | 4 | 0 | 0 | 0.06 | PASS |
| legacy alignment recovers a pure label permutation | 2 | 2 | 0 | 0 | 0.07 | PASS |
| one CV fold retains every base assignment | 5 | 5 | 0 | 0 | 1.87 | PASS |

### `test-full-synthetic-pipeline.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| full SPEC pipeline recovers known synthetic structure | 14 | 14 | 0 | 0 | 22.99 | PASS |

### `test-gam-imputation.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| GAM imputation is deterministic, bounded, and leaves observed values unchanged | 8 | 8 | 0 | 0 | 0.74 | PASS |
| too few observed values use the audited observed-value fallback | 3 | 3 | 0 | 0 | 0.08 | PASS |

### `test-generalization.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| generalization reports calibration, coverage, and intervals | 9 | 9 | 0 | 0 | 0.14 | PASS |
| Monte Carlo calibration is explicit and reproducible | 4 | 4 | 0 | 0 | 0.05 | PASS |
| auto calibration uses exact enumeration when tractable | 2 | 2 | 0 | 0 | 0.05 | PASS |

### `test-legacy-surfaces.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| superseded analysis, preparation, and report artifacts are retired | 1 | 1 | 0 | 0 | 0.03 | PASS |

### `test-load-and-prepdata.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| loader creates a finite distance matrix and documented output artifacts | 11 | 11 | 0 | 0 | 3.08 | PASS |

### `test-prediction.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| fitted feature preprocessing is unaffected by holdout sentinels | 7 | 7 | 0 | 0 | 0.21 | PASS |
| loader exposes raw development and untouched holdout rows | 4 | 4 | 0 | 0 | 0.51 | PASS |
| raw cognitive imputation is fitted on training rows only | 4 | 4 | 0 | 0 | 0.19 | PASS |
| held-out prediction stays in the fitted feature and consensus spaces | 6 | 6 | 0 | 0 | 3.86 | PASS |
| consensus prediction uses pairwise available agreement | 2 | 2 | 0 | 0 | 0.05 | PASS |

### `test-resampling.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| blocked permutation preserves within-task pairs | 4 | 4 | 0 | 0 | 0.1 | PASS |
| participant bootstrap preserves shape and unique execution IDs | 3 | 3 | 0 | 0 | 0.06 | PASS |
| empirical inference uses plus-one p-values and percentile CIs | 4 | 4 | 0 | 0 | 0.25 | PASS |

### `test-workflow-integration.R`

| validation | assertions | passed | failed | warnings | seconds | status |
| --- | --- | --- | --- | --- | --- | --- |
| default analysis output is anchored to the repository root | 2 | 2 | 0 | 0 | 0.03 | PASS |
| canonical entry point writes structured compatibility artifacts | 8 | 8 | 0 | 0 | 5.72 | PASS |
| SPEC entry point uses all methods and evaluates a fixed holdout | 4 | 4 | 0 | 0 | 8.23 | PASS |

## Interpretation

Every validation case passed, including the all-five-method synthetic pipeline and held-out consensus prediction.

Detailed test intent and coverage are documented in `test/README.md`.
