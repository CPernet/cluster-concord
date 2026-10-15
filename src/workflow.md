# Consensus clustering workflow

## Canonical entry point

The production entry point is `src/R/clustering_analysis.R`.

The completed entry point will:

1. load an Excel workbook through `load_and_prepdata.R`;
2. create the development/generalization split;
3. run five-fold base and consensus clustering;
4. derive stable labels and held-out predictions;
5. optionally run permutation and bootstrap inference; and
6. write structured tables and serialized results to the configured results directory.

## Inputs and outputs

- Raw input workbooks belong in `data/` and are never modified.
- Source code belongs in `src/R/`.
- Automated tests and deterministic fixtures belong in `test/`.
- Generated analysis artifacts belong in `results/` or a caller-provided output directory.

The public workflow accepts explicit study and participant filters, clustering methods, candidate cluster counts, fold count, resampling count, seed, worker count, and output directory. All five base algorithms are the SPEC-mode
default.

## Reproducibility policy

- Every public stochastic function accepts an explicit seed.
- Sequential execution is the Windows-safe default.
- Parallel execution is opt-in. PSOCK workers inherit the project library
  paths, source the canonical workflow modules by normalized absolute path, and
  pass a function/package preflight before resampling begins.
- Permutation data, bootstrap data, and their two analyses use distinct seeds
  derived from the master seed and iteration number, so load-balanced scheduling
  does not change results.
- Fast tests run through `Rscript test/run_tests.R`.
- Slow synthetic integration and Monte Carlo tests are enabled with `RUN_SLOW_TESTS=true`.
- The intended Windows toolchain uses `rig` to select R and `renv` to restore the committed dependency lock.
