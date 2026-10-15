# Cluster Concord

Multi-method consensus clustering in R with flexible data adapters, balanced cross-validation, stability assessment, and held-out prediction and calibration.

**Authors:** Malthe Andersen & Cyril Pernet — Neurobiology Research Unit, Denmark.

M.A. and C.P. designed the algorithm and method. Coding and testing were performed using OpenAI Codex.

The included behavioural application implements the consensus-clustering workflow defined in
[`SPEC.md`](SPEC.md). It identifies cognitive profiles from ten paired task
measures, evaluates whether those profiles are stable out of sample, and tests
their uncertainty and generalization.

The canonical entry point is
[`src/R/clustering_analysis.R`](src/R/clustering_analysis.R). This is a
script-based analysis repository, not an R package, so it is installed by
selecting the required R version and restoring the project library from
`renv.lock`.

## What the analysis does

The default `mode = "spec"` workflow:

1. Loads and screens the raw Excel data, uses the cleaned input cohort, and then
   creates a seeded 10% holdout from the screened configured cohort.
2. Fits imputation and feature transformations using development/training data
   only.
3. Represents each participant in the specified 45-feature task-similarity
   space derived from ten cognitive tasks.
4. Fits five base clustering methods:
   average-linkage hierarchical clustering (UPGMA), spherical k-means, Gaussian
   mixture models, affinity propagation, and HDBSCAN.
5. Builds a pairwise-available co-association matrix across the methods and
   obtains consensus clusters with Ward D2 clustering of `1 - co-association`.
   HDBSCAN noise and rejected assignments remain missing (`NA`), rather than
   being treated as an extra cluster.
6. Runs five balanced cross-validation folds, aligns fold labels, and defines a
   stable label only when at least three of the other four training fits agree.
7. Predicts held-out participants with each fitted base method in the
   training-derived 45-feature space, then assigns them to the fixed training
   consensus partition from their cross-method co-association evidence.
8. Fits the final model on the full development population and evaluates the
   untouched holdout's marginal assignment calibration against a cohort-aligned
   development reference. This is not treated as subject-level accuracy or an
   independent replication of the cluster structure.
9. By default, runs 1,000 blocked permutations and 1,000 participant
   bootstraps for inference.

The historical `mode = "compatibility"` path is retained only for reproduction
of the earlier analysis. New analyses should use the default SPEC mode.

## Repository layout

| Path                            | Purpose                                                                                                 |
| ------------------------------- | ------------------------------------------------------------------------------------------------------- |
| `SPEC.md`                     | Statistical and analytical requirements.                                                                |
| `src/R/clustering_analysis.R` | Public end-to-end analysis entry point.                                                                 |
| `src/R/`                      | Data preparation, clustering, prediction, resampling, generalization, and synthetic-validation modules. |
| `data/`                       | Input workbooks and synthetic input data.                                                               |
| `results/`                    | Generated analysis and validation artifacts.                                                            |
| `test/`                       | Deterministic test suite, test documentation, and generated test report.                                |
| `.Rprofile`                   | Automatically activates this repository's`renv` environment.                                          |
| `renv.lock`                   | Exact R and R-package dependency record.                                                                |

## Reproducible installation

### 1. Obtain the repository

Clone it with Git, or download and extract an archive. All following commands
must be run from the repository root.

```sh
git clone <repository-url>
cd Consensus_clustering
```

### 2. Install `rig`

[`rig`](https://github.com/r-lib/rig) manages parallel R installations and
selects the default version on Windows, macOS, and supported Linux
distributions.

Windows 10 or later, in PowerShell:

```powershell
winget install posit.rig
```

Restart the terminal if `rig` is not immediately found on `PATH`.

macOS with Homebrew:

```sh
brew tap r-lib/rig
brew trust r-lib/rig
brew install --cask rig
```

Ubuntu or Debian:

```sh
sudo curl -L https://rig.r-pkg.org/deb/rig.gpg \
  -o /etc/apt/trusted.gpg.d/rig.gpg
echo "deb http://rig.r-pkg.org/deb rig main" | \
  sudo tee /etc/apt/sources.list.d/rig.list
sudo apt update
sudo apt install r-rig
```

Fedora/RHEL, OpenSUSE/SLES, and other supported Linux distributions should use
the corresponding package from the
[`rig` installation instructions](https://github.com/r-lib/rig#-installing-rig).

### 3. Install and select R 4.6.1

These commands are the same on all three operating systems. Administrator or
`sudo` authorization may be requested by `rig`.

```sh
rig add 4.6.1
rig default 4.6.1
rig list
Rscript --version
```

`Rscript --version` must report R 4.6.1 before dependencies are restored. On
Windows, if the terminal still resolves an older R installation, either restart
it or call the validated executable directly:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' --version
```



### Alternative: install R 4.6.1 via conda (no admin rights required)

If you do not have administrator or `sudo` access to install `rig` system-wide, you can obtain R 4.6.1 through a user-owned conda environment instead. This was verified on OpenSUSE Leap 15.6 but applies to any Linux distribution where conda is available.

**3a. Create a dedicated conda environment with R 4.6.1**

Use `--override-channels` to restrict the solver to conda-forge only, which avoids conflicts with other channels (pytorch, pyg, etc.) that may be present in an existing environment:

```sh
conda create --prefix ~/conda-envs/r461 \
  -c conda-forge r-base=4.6.1 \
  --override-channels -y
```

**3b. Activate the environment and verify R**

```sh
conda activate ~/conda-envs/r461
Rscript --version   # must report R version 4.6.1
```

Keep this environment active for all subsequent steps (including
`renv::restore()` and running the analysis). The R binary used is
`~/conda-envs/r461/bin/Rscript`.

**3c. (Optional) add the environment to your shell configuration**

To activate it automatically in every new session add the following to your
`~/.bashrc` or `~/.bash_profile`:

```sh
conda activate ~/conda-envs/r461
```


### 4. Restore the locked R packages

From the repository root:

```sh
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript -e 'renv::status()'
```

The committed `.Rprofile` activates `renv` automatically. The first command
installs the package versions recorded in `renv.lock` into a project-specific
library; it does not replace packages used by unrelated R projects. The first
restore requires internet access and may require operating-system development
libraries on Linux. `renv::status()` should report that the project is
synchronized.

Do not run `renv::snapshot()` merely to resolve a version mismatch: that would
change the recorded environment. Select R 4.6.1 and run `renv::restore()`
instead. See the official
[`renv::restore()` documentation](https://rstudio.github.io/renv/reference/restore.html)
for details.

## Input data

For other datasets, see [the custom-data adapter guide](src/custom_data.md) and
the runnable [numeric-data example](src/custom_data_example.R). Users supply
data preparation and a compatible feature/similarity/distance representation;
the repository supplies consensus clustering and its evaluation workflow.
Custom distances must remain compatible with the selected methods and held-out
prediction; an arbitrary distance matrix alone is not sufficient.

The shared input is `data/data.xlsx`, a 423-participant cleaned extract containing fully randomly generated replacement participant identifiers and privacy-recoded ages. Approval to share this dataset has been obtained. The workbook is described as pseudonymized, not as confirmed anonymous.
No participant IDs are excluded by default (`exclude_ids = character(0)`).
The optional exclusion argument remains available only when explicitly supplied.
Missing observations are displayed as the literal text `NA` in Excel. The loader
uses `readxl::read_excel(..., na = c("", "NA"))`, so these import as R missing
values. Missing observations are not replaced with zero or imputed in the input.
See `data/data.json` for variable descriptions and details of the identifier and age recoding.


The normal entry point accepts a raw `.xlsx` workbook through `input_path`. At
minimum, the workbook must satisfy the schema checked by
`load_and_prepdata.R`, including:

- `CIMBI ID` as the participant key;
- `Project ID (primary)`, `Age at neuropsych test`, and `Overall comments`;
- the raw cognitive columns required to derive the ten task measures;
- the configured fold-balancing columns, which default to `Gender`,
  `BDDin Q9 - Antidepressants`, and `Age at neuropsych test`.

Participant IDs must be non-missing and unique after the loader's documented
deduplication. The source workbook is read but never modified. Generated files
go to `results/` unless another `output_dir` is supplied.

## Random replacement IDs and privacy recoding of ages

Participant IDs are fully randomly generated replacement identifiers, as confirmed by the data provider. They are not derived from the original participant codes. The dataset is described as pseudonymized because random IDs alone do not establish anonymity of the participant-level data.

Age is stored as numeric representatives of 3-year categories. The workbook includes one deliberate adjacent-category privacy recoding to remove an age/gender singleton. Representatives are released privacy categories, not exact ages or guaranteed true-age intervals. Private source-matching information is not included here.

## Run the analysis

### Full production analysis

Start R or RStudio in the repository root, then run:

```r
source("src/R/clustering_analysis.R")

result <- clustering_analysis(
  input_path = "data/data.xlsx"
)
```

Replace the example path with the intended workbook. With defaults, this uses
all five clustering algorithms, five-fold cross-validation, a 10% screened
eligible-cohort holdout, 1,000 blocked permutations, and 1,000 bootstraps. The
holdout calibration p-value is computed by exact multinomial enumeration when
tractable; otherwise it uses 9,999 Monte Carlo replicates. The run can therefore
take substantial time.

The same call can be made without an interactive R session on Windows, macOS,
or Linux:

```sh
Rscript -e 'source("src/R/clustering_analysis.R"); result <- clustering_analysis(input_path = "data/data.xlsx", workers = 20L)'
```

### Faster development run

Use this form to exercise the complete clustering and prediction path without
the expensive permutation/bootstrap stage and without writing output files:

```r
source("src/R/clustering_analysis.R")

result <- clustering_analysis(
  input_path = "data/data.xlsx",
  resamples = 0L,
  write_outputs = FALSE
)
```

Skipping resampling is appropriate for development and synthetic
known-structure validation, not for the final inferential analysis.

### Common optional arguments

| Argument                        | Default                     | Effect                                                                                                 |
| ------------------------------- | --------------------------- | ------------------------------------------------------------------------------------------------------ |
| `methods`                     | all five base methods       | Selects the methods contributing to consensus.                                                         |
| `k_candidates`                | `2:6`                     | Candidate cluster counts for bounded methods and consensus.                                            |
| `folds`                       | `5L`                      | Number of balanced outer cross-validation folds.                                                       |
| `resamples`                   | `1000L`                   | Number of blocked permutations and participant bootstraps;`0L` skips inference.                      |
| `seed`                        | `2025L`                   | Master seed used to derive deterministic stage-specific seeds.                                         |
| `workers`                     | `1L`                      | Sequential execution by default; larger values opt into parallel resampling.                           |
| `output_dir`                  | repository-root`results/` | Destination for generated RDS, CSV, and audit artifacts, independent of the current working directory. |
| `study_filter`                | `NULL`                    | Optional allowed values of`Project ID (primary)`.                                                    |
| `exclude_ids`                 | `character(0)`           | Participant IDs removed by the raw-data loader.                                                        |
| `holdout_fraction`            | `0.10`                    | Fraction sampled after screening from the eligible holdout population.                                 |
| `holdout_cohort`              | `"KOHORTE 2"`             | Eligible`Overall comments` value; `NULL` uses all screened participants.                           |
| `generalization_method`       | `"auto"`                  | Uses exact multinomial calibration when tractable, otherwise Monte Carlo.                              |
| `generalization_permutations` | `9999L`                   | Monte Carlo replicates when simulation is selected or required.                                        |
| `generalization_confidence`   | `0.95`                    | Family-wise confidence for simultaneous holdout-proportion intervals.                                  |
| `write_outputs`               | `TRUE`                    | Controls whether result artifacts are written.                                                         |

Advanced callers may provide an in-memory `data` frame and a separate
`generalization_data` frame, rename the balancing columns, or supply complete
base, consensus, distance, and feature-pipeline configurations. The header of
`src/R/clustering_analysis.R` documents those extension points.

### Parallel resampling on a server

Values of `workers > 1` use a PSOCK cluster for the blocked-permutation and
participant-bootstrap iterations. Each worker is initialized with the master's
project-library paths and sources the canonical files in `src/R` from normalized
absolute paths before resampling begins. A preflight then checks the workflow
functions and package namespaces required by the selected clustering methods.
Missing dependencies therefore fail before the first expensive resampled CV fit
and report the affected worker and its `.libPaths()`.

All workers must run on a host that can read the checked-out repository and its
restored `renv` project library at the same absolute paths. The default integer
form, such as `workers = 20L`, starts workers on the same host as the main R
process. It is not a scheduler submission mechanism for separate compute nodes.

Every resampling iteration derives four fixed seeds from the master seed and
iteration number: permutation data, bootstrap data, permutation analysis, and
bootstrap analysis. Results are therefore independent of load-balanced worker
scheduling and should match `workers = 1L` for the same runtime and package
environment.

## Outputs

Unless `write_outputs = FALSE`, `clustering_analysis()` writes the following
files to `output_dir` (by default, the repository `results/` directory).

| File                                   | Availability                       | Contents                                                                                                                                                                                                                                    |
| -------------------------------------- | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `clustering_analysis_result.rds`     | Every written run                  | Complete serializable result, including CV fits and summaries, final development fit, holdout calibration, optional inference, metadata, and file paths.                                                                                    |
| `cv_cluster_counts.csv`              | Every written run                  | One row per fold and method with the selected cluster count.                                                                                                                                                                                |
| `cv_predictions.csv`                 | Every written run                  | One row per held-out participant, fold, and method with the aligned prediction, stable-label target, and agreement diagnostics.                                                                                                             |
| `cv_final_labels.csv`                | Every written run                  | One row per participant and method.`fold_1` through `fold_F` contain aligned labels from training fits; the held-out fold is missing. `final_label`, `stability`, `stable`, and `n_assignments` summarize support across folds. |
| `cv_overlap_ari.csv`                 | Every written run                  | Pairwise adjusted Rand indices for each method on participants shared by two training folds.                                                                                                                                                |
| `permutation_p_values.csv`           | `resamples > 0`                  | Empirical blocked-permutation p-values for the CV metrics.                                                                                                                                                                                  |
| `bootstrap_confidence_intervals.csv` | `resamples > 0`                  | Participant-bootstrap percentile confidence intervals for CV metrics.                                                                                                                                                                       |
| `generalization_calibration.csv`     | SPEC mode with holdout calibration | Reference and holdout counts/proportions, simultaneous intervals, coverage, and calibration diagnostics.                                                                                                                                    |

The two assignment CSVs serve different purposes: `cv_predictions.csv` contains
the prediction made for each participant's held-out fold, while
`cv_final_labels.csv` contains the labels from the other folds where that
participant was in training and the derived cross-fold final label. They should
not be combined as duplicate observations.

`load_and_prepdata()` can also be called directly with `write_outputs = TRUE`.
That standalone preparation workflow writes `subject_distance_matrix.tsv` (the
participant distance matrix, keyed by `CIMBI`) and
`subject_distance_metadata.json` (input, seed, package, task-availability, and
feature-construction audit metadata) to its `output_dir`.

`run_synthetic_validation()` writes its human-readable validation summary to
`results/synthetic_validation_results.md` and places the same analysis
artifacts listed above in `results/synthetic_validation/`.

The returned object contains the same information in structured form:

```r
names(result)
result$cv$summary$metrics
result$cv$summary$final_labels
result$cv$summary$fold_labels
result$final_fit$consensus
result$generalization
result$inference$inference$p_values
result$output_paths
```

Entries that were intentionally skipped are `NULL`; for example,
`result$inference` is `NULL` when `resamples = 0L`.

## Synthetic known-truth validation

The synthetic module generates three known cognitive profiles, runs the full
five-method clustering and holdout-prediction pipeline, measures label recovery,
and evaluates every 2-to-5-method consensus subset against the known truth. It
deliberately skips bootstrap and permutation inference because those components
are tested separately.

```r
source("src/R/synthetic_validation.R")
validation <- run_synthetic_validation()
```

The human-readable report is written to
`results/synthetic_validation_results.md`; detailed artifacts are written under
`results/synthetic_validation/`.

## Tests

Run the fast deterministic suite from the repository root:

```sh
Rscript test/run_tests.R
```

Run every test, including the slower full synthetic pipeline, and regenerate
`test/test_report.md`:

```sh
Rscript test/generate_test_report.R
```

[`test/README.md`](test/README.md) lists every test file and explains its
validation purpose. A nonzero command exit status means at least one validation
failed.

## Cross-platform reproducibility

The validated development environment is `rig` 0.8.1, R 4.6.1, and `renv`
1.2.4. `renv.lock` pins the R-package versions, while fixed seeds make the
stochastic workflow reproducible within the limits of the implementations.

Windows, macOS, and Linux may still produce small numerical differences because
R packages can use different compilers, BLAS/LAPACK implementations, math
libraries, CPU instructions, or parallel scheduling. Those differences can
occasionally affect a borderline cluster assignment or tie. Reproducibility
checks should therefore compare documented metrics, assignments, and tolerances
rather than require byte-identical result files across operating systems.
## Reproducibility after ID replacement and row shuffling

The consensus_clustering copy orders records internally from non-ID, non-age
contents before sampling, learned preprocessing and clustering. These transient
keys are reversible encodings for ordering, not privacy protection, and are not
written as workbook columns. Identical non-ID/non-age records require an
independent persistent analysis token; the loader reports ambiguity without
deleting rows. Unique, nonmissing participant IDs remain necessary to address
outputs. Automatic deduplication and participant-ID exclusions are removed.

Holdout selection remains uniform sampling without replacement from screened
KOHORTE 2 participants. Separate seed offsets govern holdout and fold sampling.
Cross-validation randomly allocates records within sex, medication and broad age
groups, with random tie-breaking and nearly equal total fold sizes.

Age groups have configurable width 12 years and origin 18:
[18,30), [30,42), [42,54), [54,66), continuing automatically above or below
this range. Width and origin can be changed with age_bin_width and
age_bin_origin in clustering_analysis(). Missing or nonfinite ages form a
Missing group. This does not change cognitive measurements or perform age
anonymization. cv_stratum_counts.csv flags joint groups smaller than the number
of folds; cv_fold_balance.csv records their actual allocation.

Aligned 2- or 3-year age intervals nest inside the default 12-year bands when
they use origin 18 and an in-interval numeric representative. Arbitrary rounding
can cross a boundary and may change folds. Categorized age strings are not a
supported numeric representation and must be encoded explicitly before use.

This establishes a new randomization baseline. Original and pseudonymized
workbooks must agree under this implementation; agreement with the historical
ID-sorted split is a separate question. SPEC.md and renv.lock are preserved.

The metamorphic checks in test/test-order-invariance.R exercise holdout,
multiple row shuffles, ID replacement, actual model fitting and nested age
recoding, in addition to the existing regression suite.


## Input identifiers and cleaned source metadata

The reader accepts exactly one identifier column named `ID`, `CIMBI ID`, or `CIMBI`, and uses `CIMBI` internally. IDs remain unique and non-missing.
For already-cleaned workbooks, `Project ID (primary)` is optional. When absent, no source-project filter is applied; this is recorded in the preparation audit. When present, the existing source-project exclusion rule is retained.
`Overall comments` remains required because it identifies the configured KOHORTE 2 holdout population and calibration reference. The shared workbook retains this column.
To reproduce the no-resampling analysis, run `Rscript src/run_no_resampling.R` from the project root after `renv::restore()`.

## results

Using the pseudonymized data with privacy-recoded ages, the final clustering is unchanged compared to the published results, but cross-validation estimates differ materially. All 328 development cluster assignments and 24 holdout predictions are identical. The six consensus clusters still contain 80, 36, 63, 77, 24, and 48 participants. The complete holdout evaluation is unchanged, including p = 0.1870246. The stable-participant total hides a different turnover of subjects: 58 lost stability and 65 gained it. This is explained by our cross validation using age as a controlling factor ensuring folds are not biased, but binning make this sampling less stable.
