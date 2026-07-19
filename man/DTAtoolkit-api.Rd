\name{DTAtoolkit-api}
\alias{DTAtoolkit-api}
\alias{exclude_by_names}
\alias{exclude_studies_by_name}
\alias{rename}
\alias{rename.first}
\alias{rename_count_columns}
\alias{rename_prefixed_count_columns}
\alias{het.string}
\alias{heterogeneity_string}
\alias{reitsma.subgroup.comparisons}
\alias{reitsma_subgroup_comparisons}
\alias{forest.diag}
\alias{forest_diag}
\alias{forest.diag.combined}
\alias{forest_diag_combined}
\alias{forest.diag.subgroup}
\alias{forest_diag_subgroup}
\alias{forest.diag.subgroup.combined}
\alias{forest_diag_subgroup_combined}
\alias{forest.diag.no}
\alias{forest_diag_no}
\alias{forest.diag.subgroup.no}
\alias{forest_diag_subgroup_no}
\alias{multiple.srocs}
\alias{multiple_srocs}
\alias{multiple.srocs.no}
\alias{multiple_srocs_no}
\alias{dta.outliers}
\alias{dta_outliers}
\alias{dta.outliers.multi}
\alias{dta_outliers_multi}
\alias{dta_palette}
\alias{dta_auc_bootstrap}
\alias{dta_auc_ci_cache_info}
\alias{dta_clear_auc_ci_cache}
\alias{dta_mvpbt_bootstrap}
\alias{dta_mvpbt}
\alias{dta_subgroup_comparisons}
\alias{dta_heterogeneity_text}
\alias{dta_exclude_studies}
\alias{dta_fagan_nomogram}
\alias{dta_rename_counts}
\alias{dta_rename_prefixed_counts}
\alias{dta_find_repeated_studies}
\alias{dta_find_repeated_studies_by_subgroup}
\alias{dta_combine_cohorts}
\alias{dta_parametric_bootstrap_sample}
\alias{AUC_bootstrap}
\alias{auc_bootstrap}
\alias{AUC_boot_paralell}
\alias{PBS3}
\alias{pbs3}
\alias{MVPBT_bootstrap}
\alias{mvpbt_bootstrap}
\alias{MVPBT_boot}
\alias{mvpbt_boot}
\alias{pubbias.diag}
\alias{pubbias_diag}
\alias{find_repeated_studies}
\alias{find_repeated_studies_by_subgroup}
\alias{multiple.LRmats}
\alias{multiple_lr_mats}
\alias{nomogrammer}
\alias{nomogrammer_plus}
\alias{nomogrammer_subgroups}
\alias{df.mixer}
\alias{df_mixer}
\title{DTAtoolkit Analysis Functions}
\description{
The public analysis API for data preparation, forest plots, SROC plots,
subgroup testing, influential-study analysis, bootstrap inference, likelihood-ratio
matrices, nomograms, repeated-study checks, and data mixing.
}
\details{
Use the consolidated \code{dta_*} names for new analyses. Legacy names and the
beta's earlier snake-case aliases remain exported unchanged. The complete
mapping is provided in \file{FUNCTION_NAMES.md} and the main workflows are
introduced in the package README.
}
\section{Data columns}{
Preferred functions evaluate column selectors inside the supplied data. Bare
names such as \code{tp = tp}, character names, vectors, expressions, and
\code{data$column} references are accepted. The preferred names are
\code{study}, \code{tp}, \code{tn}, \code{fp}, and \code{fn}; the beta 0.1
long-form arguments remain supported. The historical columns
\code{names}, \code{TP}, \code{TN}, \code{FP}, and \code{FN} remain defaults.
When omitted, \code{uniquer.row.id} is the row number and \code{study.id}
defaults to the study-name values.
}
\section{Compatibility}{
Prefer the \code{dta_*} names in new code. Existing scripts may continue to use
the legacy names, including the misspelled \code{AUC_boot_paralell()} alias.
\code{dta_forest()} consolidates all pooled, subgroup, single-outcome, and
combined forest functions.
}
\section{SROC AUC-CI cache}{
SROC AUC confidence intervals use a single-entry cache in session memory by
default. Its identity includes an MD5 digest of retained TP/TN/FP/FN values and
subgroup labels, the seed, and the number of bootstrap replicates. Use
\code{cache_auc_ci = FALSE} in \code{dta_sroc()} or
\code{AUC.CI.cache = FALSE} in \code{multiple.srocs()} to bypass it.
\code{dta_auc_ci_cache_info()} inspects the current entry and
\code{dta_clear_auc_ci_cache()} clears it. Mismatched input clears the previous
entry before computation, and results are stored only after all subgroups
finish successfully.
}
\section{Forest model label}{
The pooled Reitsma result in every forest display is labelled
\sQuote{Bivariate model}. The generic random-effects label is not used.
}
\section{Forest panel ratio and text size}{
Combined plots accept \code{left_panel_width_adjustment}, a percentage-point
change to the dynamically calculated left-panel share. The right-panel share
changes by the opposite amount. The deprecated
\code{sensitivity_width_adjustment} alias remains accepted, and legacy
combined functions expose \code{sensitivity.width.adjustment}.
\code{font_size} in the preferred API and \code{fontsize} in legacy forest
functions control forest text size.
}
\section{Likelihood-ratio plot layout}{
\code{multiple_lr_mats()} and \code{multiple.LRmats()} use responsive legend
placement by default. Wide devices use a separate right-hand legend panel;
narrow devices use a bottom panel and reflow its content. The placement can be
overridden with \code{legend.position}, and \code{legend.cex} can override the
automatically reduced legend size.
}
\section{Plot colours}{
SROC and likelihood-ratio matrix plots share the six-colour \code{"mada6"}
default. The \code{color.palette} argument also accepts a palette name returned
by \code{palette.pals()} or \code{hcl.pals()}, or a manually supplied character
vector of R colours. Use \code{dta_palette()} to resolve or preview colours.
Explicit legacy colour arguments take precedence for backward compatibility.
}
\section{Study marker sizes}{
Individual-study marker area in SROC, likelihood-ratio matrix, and forest plots
is proportional to \code{TP + TN + FP + FN}. A single maximum population is
used across every retained subgroup in a figure, making marker sizes directly
comparable between subgroups.
}
\seealso{\code{help(package = "DTAtoolkit")}}
