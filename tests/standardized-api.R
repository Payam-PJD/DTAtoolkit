library(DTAtoolkit)

data(dta_example)

dat <- dta_example
dat$group <- dat$subgroup
dat$row_id <- dat$unique_row_id
dat$author_id <- dat$study_id

inf_dat <- data.frame(
  study = paste0("S", seq_len(12)),
  tp = c(45, 44, 43, 42, 41, 40, 35, 34, 33, 32, 31, 1),
  fn = c(5, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 49),
  tn = c(45, 46, 44, 43, 42, 41, 35, 36, 34, 33, 32, 1),
  fp = c(5, 4, 6, 7, 8, 9, 5, 4, 6, 7, 8, 49),
  group = rep(c("A", "B"), each = 6),
  row_id = seq_len(12),
  author_id = paste0("author-", seq_len(12)),
  stringsAsFactors = FALSE
)

columns <- list(
  study_name = "study",
  true_positive = "tp",
  true_negative = "tn",
  false_positive = "fp",
  false_negative = "fn"
)

pooled <- do.call(
  dta_forest,
  c(list(data = dat, type = "sens", draw = FALSE), columns)
)
stopifnot(
  inherits(pooled, "dta_analysis"),
  identical(pooled$analysis, "forest"),
  identical(pooled$result$metaprop.sens$text.random, "Bivariate model")
)

subgroup <- do.call(
  dta_forest,
  c(list(
    data = dat, subgroup = "group", type = "spec", draw = FALSE,
    pairwise = "never"
  ), columns)
)
stopifnot(
  inherits(subgroup, "dta_analysis"),
  length(subgroup$result$reitsmas$subgroups) == 3L
)

adjusted_forest <- do.call(
  dta_forest,
  c(list(
    data = dat, subgroup = "group", type = "combined", draw = FALSE,
    left_panel_width_adjustment = 1,
    font_size = 9, pairwise = "never"
  ), columns)
)
stopifnot(
  !identical(adjusted_forest$result$layout$measurement, "fallback"),
  all(is.finite(adjusted_forest$result$layout$widths.cm)),
  !isTRUE(all.equal(
    unname(adjusted_forest$result$layout$base.ratio),
    c(0.7, 0.3),
    check.attributes = FALSE
  )),
  isTRUE(all.equal(
    adjusted_forest$result$layout$ratio,
    adjusted_forest$result$layout$base.ratio + c(0.01, -0.01),
    check.attributes = FALSE
  )),
  identical(
    adjusted_forest$result$layout$ratio.adjustment, 1
  ),
  identical(adjusted_forest$result$layout$font.size, 9)
)

# The former argument name remains a compatibility alias.
resolve_adjustment <- getFromNamespace(
  ".dta_resolve_left_panel_width_adjustment", "DTAtoolkit"
)
stopifnot(identical(
  resolve_adjustment(
    left_panel_width_adjustment = 0,
    sensitivity_width_adjustment = 1,
    preferred_supplied = FALSE,
    legacy_supplied = TRUE
  ),
  1
))

# Preferred functions use data masking: bare names, character names, vectors,
# and data$column expressions are all valid column selectors.
masked_forest <- dat |>
  dta_forest(
    subgroup = group,
    type = "sens",
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    unique_row_id = row_id,
    study_id = author_id,
    left_column_1 = group,
    pairwise = "never",
    draw = FALSE
  )
stopifnot(
  identical(masked_forest$recipe$arguments$subgroup, ".dta_subgroup"),
  identical(masked_forest$data$uniquer.row.id, dat$row_id),
  identical(masked_forest$data$study.id, dat$author_id)
)

influence_messages <- capture.output(
  masked_influence <- inf_dat |>
    dta_influentials(
      subgroup = group,
      study = study,
      tp = tp, tn = tn, fp = fp, fn = fn,
      unique_row_id = row_id,
      study_id = author_id
    )
)
stopifnot(
  identical(masked_influence$scope, "within_subgroup"),
  identical(names(masked_influence$fits), c("A", "B")),
  all(masked_influence$diagnostics$subgroup %in% c("A", "B")),
  all(c("study_name", "unique_row_id") %in%
        names(masked_influence$influential_studies)),
  any(grepl("author/study name: S12", influence_messages, fixed = TRUE)),
  any(grepl("unique row ID: 12", influence_messages, fixed = TRUE))
)

masked_reduced <- masked_forest |>
  remove_inf(draw = FALSE, verbose = FALSE)
stopifnot(
  identical(masked_reduced$influential_analysis$scope, "within_subgroup"),
  identical(
    names(masked_reduced$influential_analysis$fits),
    c("Platform A", "Platform B", "Platform C")
  ),
  all(c("study_name", "unique_row_id", "subgroup") %in%
        names(masked_reduced$removed_studies)),
  is.data.frame(masked_reduced$removed_studies)
)

influence <- do.call(
  dta_influentials,
  c(list(data = inf_dat, verbose = FALSE), columns)
)
reduced_data <- influence |> remove_inf()
stopifnot(
  inherits(influence, "dta_influentials"),
  is.data.frame(reduced_data),
  nrow(reduced_data) < nrow(dat)
)

sroc_result <- inf_dat |>
  dta_sroc(
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    draw = FALSE
  )
# Pretend that the originating call requested an AUC interval. This checks the
# no-second-bootstrap branch without spending time on a bootstrap in tests.
sroc_result$recipe$arguments$auc_ci <- TRUE
reduced_sroc <- sroc_result |>
  remove_inf(draw = FALSE, verbose = FALSE)
stopifnot(
  inherits(reduced_sroc, "dta_analysis"),
  nrow(reduced_sroc$data) < nrow(sroc_result$data),
  identical(reduced_sroc$recipe$arguments$auc_ci, FALSE),
  grepl("skipped", reduced_sroc$bootstrap_after_removal)
)

lr_result <- dat |>
  dta_lr_matrix(
    subgroup = group,
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    draw = FALSE
  )
nomogram_result <- dat |>
  dta_nomogram(
    prevalence = 0.2,
    study = study,
    tp = tp, tn = tn, fp = fp, fn = fn,
    return_list = TRUE,
    draw = FALSE
  )
stopifnot(
  identical(lr_result$recipe$function_name, "dta_lr_matrix"),
  identical(nomogram_result$recipe$function_name, "dta_nomogram"),
  identical(
    formals(dta_publication_bias)$n_boots,
    2000
  )
)
