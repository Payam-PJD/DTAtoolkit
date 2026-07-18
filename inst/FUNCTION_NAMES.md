# Standardized and legacy function names

Use the `dta_*` names for new analyses. Every historical function listed below
remains exported so existing scripts continue to run.

| Preferred function | Replaces or consolidates | Purpose |
|---|---|---|
| `dta_forest()` | `forest.diag()`, `forest.diag.combined()`, `forest.diag.subgroup()`, `forest.diag.subgroup.combined()` | One pooled/subgroup forest function; use `type = "sens"`, `"spec"`, or `"combined"` |
| `dta_sroc()` | `multiple.srocs()` | One pooled/subgroup SROC function |
| `dta_lr_matrix()` | `multiple.LRmats()` | One pooled/subgroup likelihood-ratio matrix |
| `dta_nomogram()` | `nomogrammer_plus()`, `nomogrammer_subgroups()` | Data-based pooled/subgroup nomograms |
| `dta_fagan_nomogram()` | `nomogrammer()` | Scalar Fagan nomogram from sensitivity/specificity or likelihood ratios |
| `dta_publication_bias()` | `pubbias.diag()` | Publication-bias diagnostics and MVPBT bootstrap |
| `dta_influentials()` | `dta.outliers()`, `dta.outliers.multi()` | Influential-study analysis, without a figure by default |
| `remove_inf()` | `forest.diag.no()`, `forest.diag.subgroup.no()`, `multiple.srocs.no()` | Pipe-friendly influential-study removal and regeneration |
| `dta_auc_bootstrap()` | `AUC_bootstrap()`, misspelled `AUC_boot_paralell()` | Bootstrap AUC interval |
| `dta_auc_ci_cache_info()` | new | Inspect the latest in-memory SROC AUC-CI cache entry |
| `dta_clear_auc_ci_cache()` | new | Clear the in-memory SROC AUC-CI cache |
| `dta_mvpbt_bootstrap()` | `MVPBT_bootstrap()` | MVPBT bootstrap |
| `dta_mvpbt()` | `MVPBT_boot()` | MVPBT test |
| `dta_parametric_bootstrap_sample()` | `PBS3()` | Low-level parametric bootstrap sample generator |
| `dta_subgroup_comparisons()` | `reitsma.subgroup.comparisons()` | Omnibus and pairwise bivariate subgroup comparisons |
| `dta_heterogeneity_text()` | `het.string()` | Forest heterogeneity annotation |
| `dta_exclude_studies()` | `exclude_by_names()` | Exclude named studies |
| `dta_rename_counts()` | `rename()` | Copy suffixed count columns to TP/TN/FP/FN |
| `dta_rename_prefixed_counts()` | `rename.first()` | Copy prefixed count columns to TP/TN/FP/FN |
| `dta_find_repeated_studies()` | `find_repeated_studies()` | Find repeated study IDs |
| `dta_find_repeated_studies_by_subgroup()` | `find_repeated_studies_by_subgroup()` | Find repeated study IDs within subgroups |
| `dta_combine_cohorts()` | `df.mixer()` | Combine training and validation cohorts |
| `dta_palette()` | unchanged | Resolve the common plot palette |

The beta's earlier snake-case aliases are also retained with their existing
interfaces: `exclude_studies_by_name()`, `rename_count_columns()`,
`rename_prefixed_count_columns()`, `heterogeneity_string()`,
`reitsma_subgroup_comparisons()`, `forest_diag()`,
`forest_diag_combined()`, `forest_diag_subgroup()`,
`forest_diag_subgroup_combined()`, `forest_diag_no()`,
`forest_diag_subgroup_no()`, `multiple_srocs()`, `multiple_srocs_no()`,
`dta_outliers()`, `dta_outliers_multi()`, `auc_bootstrap()`, `pbs3()`,
`mvpbt_bootstrap()`, `mvpbt_boot()`, `pubbias_diag()`,
`multiple_lr_mats()`, and `df_mixer()`. The new `dta_*` interface is the only
interface that consolidates pooled/subgroup routing and supports
`remove_inf()` regeneration recipes.

## Pipe examples

```r
dat |>
  dta_forest(subgroup = Sequence, type = "combined",
             study = study, tp = tp, tn = tn, fp = fp, fn = fn) |>
  remove_inf()

dat |>
  dta_sroc(subgroup = Sequence, auc_ci = TRUE,
           study = study, tp = tp, tn = tn, fp = fp, fn = fn) |>
  remove_inf()
```

For the second example, `remove_inf()` does not start a second AUC bootstrap by
default. The reduced SROC is drawn with point AUC estimates. Use
`recompute_bootstrap = TRUE` when a new interval for the reduced dataset is
required.

Preferred column arguments use data masking. Bare column names are recommended;
character names, complete vectors, expressions, and `data$column` references
are also accepted. The preferred count arguments are `tp`, `tn`, `fp`, and
`fn`; beta 0.1's longer argument names remain supported.
