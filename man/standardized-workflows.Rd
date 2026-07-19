\name{standardized-workflows}
\alias{standardized-workflows}
\alias{dta_forest}
\alias{dta_sroc}
\alias{dta_lr_matrix}
\alias{dta_nomogram}
\alias{dta_publication_bias}
\alias{dta_influentials}
\alias{remove_inf}
\title{Standardized Diagnostic-Accuracy Workflows}
\description{
Unified, pipe-oriented functions for diagnostic-accuracy forest plots, SROC
plots, likelihood-ratio matrices, nomograms, publication-bias diagnostics, and
within-subgroup influential-study analysis.
}
\details{
The preferred functions use \code{data} as their first argument and therefore
work naturally with R's native pipe. Data-column selectors are evaluated inside
\code{data}. For example, \code{tp = tp} selects the \code{tp} column without
requiring \code{data$tp}. Character column names, vectors, expressions, and
explicit \code{data$column} expressions remain valid.

The short selectors \code{study}, \code{tp}, \code{tn}, \code{fp}, and
\code{fn} are preferred. For compatibility, \code{study_name},
\code{true_positive}, \code{true_negative}, \code{false_positive}, and
\code{false_negative} remain accepted. Supplying both names for the same
column is an error.
}
\section{Unified forest function}{
\code{dta_forest()} automatically uses a subgroup model when \code{subgroup}
is supplied and a pooled model otherwise. The \code{type} argument accepts
\code{"sens"}, \code{"spec"}, or \code{"combined"}. Forest pooled rows are
labelled \sQuote{Bivariate model}.

For combined plots, \code{left_panel_width_adjustment} changes the dynamically
calculated left-panel share in percentage points, while the right-panel share
changes by the opposite amount. For example, if automatic layout calculates
57/43, a value of \code{1} changes it to 58/42 and \code{-1} changes it to
56/44. Zero preserves the automatic ratio. The former
\code{sensitivity_width_adjustment} name remains accepted as a deprecated
compatibility alias. \code{font_size} sets the forest text size for pooled,
subgroup, and combined plots.
}
\section{Influential analysis and subgroup inheritance}{
\code{dta_influentials()} uses the bivariate mean-shift method by default and
does not draw a figure unless \code{plot = TRUE}. When \code{subgroup} is
supplied, the analysis is fitted independently inside every eligible subgroup;
there is no pooled influential test. The bivariate method requires at least
four complete studies per subgroup and the DOR method requires at least three.
Subgroups below the applicable minimum are reported in
\code{skipped_subgroups} and are not tested or modified.

Results returned by \code{dta_forest()}, \code{dta_sroc()},
\code{dta_lr_matrix()}, and \code{dta_nomogram()} retain their resolved
subgroup column. Piping one of these results to \code{remove_inf()} therefore
reuses exactly that subgroup definition. Influential messages and the
\code{influential_studies} and \code{removed_studies} tables pair author/study
name with unique row ID, study ID, and subgroup.
}
\section{Bootstrap execution order}{
R evaluates the left side of a native pipe before calling
\code{remove_inf()}. Thus \code{dta_sroc(auc_ci = TRUE) |> remove_inf()} first
runs an AUC bootstrap on the original data. By default a second AUC bootstrap
is skipped after removal; \code{recompute_bootstrap = TRUE} requests a new
interval for the reduced data.

Similarly, \code{dta_publication_bias() |> remove_inf()} runs MVPBT on the
original data and reruns it on the reduced data when rows are removed. To run
only a post-removal bootstrap, first pipe the data through
\code{dta_influentials() |> remove_inf()}, then pass the reduced data to the
bootstrap analysis.

AUC and MVPBT bootstraps use the legacy single-pass \code{foreach} execution
model. By default the package uses two fewer workers than the logical cores
reported by R; use \code{n_cores} (or legacy \code{n.cores}) to override this,
or \code{parallel = FALSE} for sequential execution. One backend is reused
across bootstrap calls and all subgroup AUC calculations, while the active
backend returns to sequential mode after each analysis. Progress is updated as
completed results return to the main R session through a \code{doSNOW}
callback, without dividing the bootstrap into batches. With
\code{progress = FALSE}, execution uses the original \code{doParallel}
backend. The running calculation can be stopped with RStudio's Stop command or
Escape.

The controls are independent: \code{progress = FALSE} hides the bar but keeps
parallel execution, while \code{parallel = FALSE} runs sequentially. If worker
initialization fails, the package discards its managed pool, creates fresh
workers, and retries once.

On Escape or another user interrupt, the package closes its managed worker
cluster and socket connections and restores sequential execution before
propagating the interrupt. Connections unrelated to DTAtoolkit are not closed.
}
\section{AUC confidence-interval cache}{
The latest complete SROC \code{AUC_CIs} object is cached in session memory.
Cache identity combines an MD5 digest of retained diagnostic counts and
subgroup labels with the seed and replicate count. Identical preferred and
legacy calls reuse it. Any mismatch clears the old entry before recomputation;
the replacement is stored only after all subgroup bootstraps succeed.
}
\value{
The plotting and publication-bias functions return a \code{dta_analysis}
object containing the fitted result, standardized data, and a regeneration
recipe. \code{dta_influentials()} returns diagnostics, subgroup-specific fits,
and an \code{influential_studies} table. \code{remove_inf()} returns filtered
data when given an influential-analysis result, or a regenerated
\code{dta_analysis} object with a \code{removed_studies} table when given an
analysis result.
}
\examples{
\dontrun{
forest_result <- dat |>
  dta_forest(
    subgroup = Sequence,
    study = author,
    tp = tp, tn = tn, fp = fp, fn = fn,
    type = "combined"
  )

reduced_forest <- forest_result |> remove_inf()

reduced_data <- dat |>
  dta_influentials(
    subgroup = Sequence,
    study = author,
    tp = tp, tn = tn, fp = fp, fn = fn
  ) |>
  remove_inf()
}
}
