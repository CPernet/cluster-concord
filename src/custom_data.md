# Using your own data

The consensus, cross-validation, agreement, stability, prediction and calibration
workflow can be reused with a dataset-specific adapter. The bundled
`load_and_prepdata.R` is the behavioural application, not a universal importer.
Keep it intact and supply your own functions rather than replacing the module:
it also contains shared ordering utilities used elsewhere.

## Recommended route: load data yourself

Source `src/custom_data_example.R`, then run `run_custom_data_example()` for a
complete synthetic example with no resampling and no written results. Its
`make_numeric_pipeline()` is an illustrative numeric/cosine adapter, not a
recommended representation for every scientific question.

```r
source("src/custom_data_example.R")
result <- clustering_analysis(
  data = development,
  generalization_data = holdout, # optional; NULL means no holdout evaluation
  predictor_variables = c("measure_a", "measure_b", "measure_c"),
  feature_pipeline = make_numeric_pipeline(c("measure_a", "measure_b", "measure_c")),
  balance_cols = c("site", "group"),
  age_col = "age", age_bin_width = 10, age_bin_origin = 20,
  holdout_cohort = NULL,
  resamples = 0L, write_outputs = FALSE
)
```

Each input row is one observation. Map your ID column to `CIMBI`: unique,
nonempty character IDs are the internal interface, including in output tables.
Keep development and holdout participants disjoint. Use the same feature columns
in both tables. Do not include IDs or balancing variables as predictors unless
scientifically intended. Records otherwise identical after excluding ID and age
require an independent persistent analysis token for deterministic ordering;
do not include that token among features.

## Optional custom loader

Pass `data_loader = your_loader` with `input_path` and `feature_pipeline`.
The loader receives named arguments `input_path`, `holdout_fraction`,
`holdout_cohort`, `seed`, `output_dir`, `write_outputs`, and `return`.
Accept the applicable arguments explicitly and use `...` for unused ones.
In `mode = "spec"` it must return a list with:

- `screened_development_raw_data`: a data frame of eligible development rows.
- `screened_holdout_data`: a data frame of untouched holdout rows, or a zero-row
  data frame with the same schema when there is no holdout.

Additional audit fields may be returned and are retained in `result$loaded`.
The loader owns screening and split selection. It must not fit imputation,
scaling, or other learned transformations on all rows before cross-validation.
The built-in loader additionally returns prepared data, distances and audit
objects; these extra fields are not required for a custom SPEC loader.
The historical compatibility mode instead uses `prepared_data` and is not the
recommended route for new datasets.

For custom loaders the holdout calibration reference is all development rows;
the behavioural `Overall comments` cohort rule is not applied automatically.
Do any intended cohort selection in your adapter.

## Feature pipeline contract

`feature_pipeline(train, test)` returns `list(train = ..., test = ...,
preprocessor = ...)`. Fit every learned transformation on `train` only, then
apply those fitted parameters to both tables. `preprocessor` stores the fitted
state for inspection. Each transformed table contains:

| Field | Contract |
|---|---|
| `features` | Finite numeric matrix, rows named by observation IDs, consistently named and ordered feature columns |
| `participant_ids` | Character IDs in exactly the feature-row order |
| `similarity` | Finite square within-table similarity matrix in the same observation order |
| `distance` | R `dist` object in that same order, consistent with the representation |

This is not an arbitrary-distance-only interface. The existing algorithms use
features, similarities and distances in different ways; held-out prediction
also uses features, including cosine comparisons. Changing distance alone can
make fitting and prediction inconsistent. A different geometry may require
adapting `prediction.R` and selecting compatible base methods. In SPEC mode,
the feature pipeline provides the clustering inputs; `distance_builder` alone
does not override it. There is no requirement for exactly 45 features in the
shared matrix interface; 45 is the behavioural application's representation.

## Fold balancing

- `balance_cols = NULL` preserves the historical sex/medication column defaults
  (or the names supplied through `sex_col` and `medication_col`).
- A character vector selects any number of categorical balancing columns.
- `balance_cols = character(0)` disables categorical balancing.
- `age_col = NULL` disables age balancing. Together with no categorical columns,
  this yields seeded folds with nearly equal sizes and no stratification.
- Age bands have configurable `age_bin_width = 12` and `age_bin_origin = 18`.
  Defaults mean `[18,30)`, `[30,42)`, etc., extended across the observed range.
  These are balancing bands, not age eligibility limits or privacy recoding.
- Other continuous covariates should be binned by the adapter before inclusion
  in `balance_cols`; listed columns are treated as categorical.

Missing categories form their own group. Joint strata guide fold assignment;
small strata cannot populate every fold and are reported in diagnostics.
Avoid diagnostic field names `fold`, `n`, `age_bin`, `all_rows`, and
`can_populate_every_fold` as custom balancing column names.

## Optional resampling

Use `resamples = 0L` to bypass permutation/bootstrap inference while retaining
CV, final fitting and holdout evaluation. If enabling inference, pass meaningful
`task_blocks`: a list of column-name vectors whose variables must move together
under the permutation null. Do not inherit cognitive task blocks for unrelated
data. A supplied custom feature pipeline is preserved during resampling;
the default behavioural workflow retains its original preprocessing path.
Custom callbacks used with parallel workers must carry their dependencies
(for example in their closure); begin with `workers = 1L` when validating an adapter.
