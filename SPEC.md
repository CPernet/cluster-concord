# Consensus clustering specification

This repository identifies cognitive profiles from ten paired task measures,
evaluates the stability of assignments within the development population, and
checks the marginal calibration of assignments in an untouched holdout.

## Consensus clustering

Five base algorithms are fitted: average-linkage hierarchical agglomerative
clustering (UPGMA), spherical k-means, Gaussian mixture modelling, affinity
propagation, and HDBSCAN. Each participant is represented by the 45 lower-triangle
entries of the ten-task similarity matrix.

Consensus is based on pairwise agreement among available base assignments. For
each participant pair, the co-association value is the proportion of methods
that assign the pair to the same cluster among methods that assign both
participants. HDBSCAN noise and rejected predictions remain missing and are not
treated as an additional cluster. Ward D2 clustering is applied to
`1 - co-association` to obtain the consensus partition.

## Statistical framework

### Cross-validation and assignment stability

Five outer folds are constructed to be nearly equal in size while balancing
sex, antidepressant use, and age distribution. Every learned preprocessing
operation is fitted on the current training data and applied unchanged to its
validation data. Cluster labels are aligned across fits because their numeric
values are arbitrary.

For every fold, base method, and consensus solution, the workflow records the
selected number of clusters and held-out predictions. Fold stability is
quantified with adjusted Rand index on overlapping training participants.
Prediction agreement compares each held-out assignment with a subject-level
label derived from the other four training fits. A stable label requires
agreement in at least three of those four fits; otherwise it remains missing.

These cross-validation results assess reproducibility and assignment agreement
inside the observed population. They are distinct from the holdout calibration
test below.

### Permutation inference and bootstrap confidence intervals

When `resamples > 0`, two development-population resampling procedures rerun the
five-fold pipeline:

- Blocked permutation independently shuffles each of the ten paired
  performance/reaction-time task blocks across participants. This preserves the
  within-task pair but disrupts between-task cognitive profiles. Empirical
  upper- or lower-tail p-values compare observed CV metrics with their
  permutation null distributions.
- Participant bootstrap samples development participants with replacement and
  uses percentile intervals to quantify uncertainty in the CV metrics.

These resamples concern development-set clustering evidence. They are not the
same as the multinomial calibration calculation below.

## Holdout assignment calibration

### Holdout construction

Protocol exclusions, deduplication, project filtering, and the configured
missingness screen are applied before holdout sampling. With the workbook
defaults, `holdout_fraction = 0.10` is sampled without replacement from screened
participants whose `Overall comments` value is `KOHORTE 2`. Consequently, 10%
refers to the screened eligible cohort, not necessarily 10% of every screened
participant across all cohorts. Setting `holdout_cohort = NULL` makes all
screened participants eligible.

The holdout is removed before imputation, transformation, cluster fitting, or
model selection. Preprocessing and all five clustering methods are fitted on
development data only. Each holdout participant is transformed with those
fixed preprocessing parameters, predicted by the available base methods, and
assigned to the fixed development consensus partition using cross-method
co-association evidence. Rejected assignments remain missing and are reported
through coverage.

When sampling is restricted to a configured cohort, the reference cluster
proportions are calculated from development participants in that same cohort.
When all participants are eligible, all development participants form the
reference. This aligns the calibration null population with the population
from which the holdout was sampled.

### Estimand and null hypothesis

The calibration diagnostic compares *marginal predicted cluster proportions*.
Let:

- `n` be the number of non-rejected holdout assignments;
- `K` be the number of development consensus clusters;
- `p_k` be the proportion of reference-development participants in cluster
  `k`;
- `O_k` be the number of assigned holdout participants predicted into cluster
  `k`; and
- `E_k = n p_k` be its expected count under the null.

The null hypothesis is

`H0: (O_1, ..., O_K) follows Multinomial(n, p_1, ..., p_K)`.

The discrepancy statistic is Pearson's multinomial statistic:

`X2 = sum_k ((O_k - E_k)^2 / E_k)`.

This is a calibration or distribution-shift diagnostic. It does **not** test
whether individual assignments are correct, whether the cluster structure was
independently recovered, or whether the clustering “generalizes” in every
sense. A small p-value means the marginal holdout assignment distribution is
inconsistent with the fixed reference proportions under the multinomial model;
it does not identify whether the cause is sampling variation, cohort shift,
model instability, or a genuine population difference.

### Exact and Monte Carlo p-values

`generalization_method = "auto"` is the default. The number of possible count
vectors is

`choose(n + K - 1, K - 1)`.

If that number does not exceed `generalization_exact_max_tables` (default
1,000,000), every possible multinomial count vector is enumerated. Its exact
multinomial probability is calculated, and the exact upper-tail p-value is the
sum of probabilities for count vectors whose statistic is at least the observed
statistic. This avoids simulation error for small holdouts.

If exact enumeration would exceed the limit, `auto` falls back to Monte Carlo.
Monte Carlo can also be requested explicitly with
`generalization_method = "monte_carlo"`. For each of `B` replicates, default
`B = 9,999`, the workflow:

1. generates one count vector from `Multinomial(n, p_1, ..., p_K)`;
2. computes its Pearson statistic with the same fixed expected counts; and
3. records whether it is at least as large as the observed statistic.

If `b` of the `B` simulated statistics are at least as extreme as observed, the
reported p-value is

`p = (b + 1) / (B + 1)`.

The `+1` correction prevents zero p-values and treats the observed configuration
consistently with the simulated null draws. The minimum attainable p-value is
`1 / (B + 1)`. The result also reports `b`, the approximate Monte Carlo standard
error `sqrt(p(1-p)/(B+1))`, and the p-value resolution. Increasing `B` reduces
only Monte Carlo error; it does not increase the holdout sample size, power, or
assignment accuracy.

### Proportion uncertainty and reporting

For every holdout cluster proportion, the workflow reports a conservative
simultaneous confidence interval. It uses Wilson score limits with a Bonferroni
adjustment across the `K` cluster proportions, with family-wise confidence set
by `generalization_confidence` (default 0.95). The output also reports:

- reference and holdout counts and proportions;
- assigned and rejected participant counts;
- assignment coverage;
- the exact or Monte Carlo method actually used;
- the calibration statistic and p-value;
- the reference population description; and
- Monte Carlo diagnostics when simulation is used.

Because a small holdout can produce wide intervals and results that change with
one assignment, interpretation must emphasize counts, intervals, coverage, and
effect size rather than a binary p-value threshold. Assignment stability is
assessed primarily through the leakage-free outer cross-validation. Strong
claims of external generalization require a genuinely independent replication
cohort; the internal holdout cannot substitute for one.

## Workflow

### Input

The Excel loader applies project and participant exclusions, deduplicates IDs,
screens cognitive-task missingness, creates the seeded holdout, and returns raw
development and holdout frames. Development-only preprocessing then produces
the 45-feature representation and cosine subject distance used for clustering.

### Core function

`src/R/clustering_analysis.R` is the public entry point. It runs five-fold
cross-validation, fits the fixed full-development clustering solution, predicts
the untouched holdout, evaluates holdout assignment calibration, and optionally
runs blocked permutations and participant bootstraps.

### Outputs

Unless `write_outputs = FALSE`, `clustering_analysis()` writes the following
artifacts to `output_dir` (default `results/`). The RDS is the complete
serializable result; the CSV files are relational views intended for inspection
and reporting.

| File | Availability | Contents |
| --- | --- | --- |
| `clustering_analysis_result.rds` | Every written run | Complete result object, including CV fits and summaries, final development fit, generalization result, inference result, metadata, and output paths. |
| `cv_cluster_counts.csv` | Every written run | One row per fold and method with the selected number of clusters. |
| `cv_predictions.csv` | Every written run | One row per held-out participant, fold, and method with the aligned prediction, stable-label target, and prediction-agreement diagnostics. |
| `cv_final_labels.csv` | Every written run | One row per participant and method. `fold_1` through `fold_F` are aligned labels from training fits; the participant's held-out fold is missing. `final_label` is the cross-fold stable label, with `stability`, `stable`, and `n_assignments` recording its support. |
| `cv_overlap_ari.csv` | Every written run | Pairwise adjusted Rand indices for overlapping training participants in each method's fold fits. |
| `permutation_p_values.csv` | `resamples > 0` | Empirical blocked-permutation p-values for CV metrics. |
| `bootstrap_confidence_intervals.csv` | `resamples > 0` | Participant-bootstrap percentile confidence intervals for CV metrics. |
| `generalization_calibration.csv` | When holdout calibration is run | Reference and holdout counts and proportions, simultaneous intervals, and calibration details for the fixed development partition. |

`cv_final_labels.csv` describes assignments fitted without the participant's
held-out fold, whereas `cv_predictions.csv` records that held-out fold's
prediction. These are complementary outputs and should not be treated as
duplicate measurements.
