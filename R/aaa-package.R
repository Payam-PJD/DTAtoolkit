# Package-level declarations.

utils::globalVariables(c(
  ".data", ".mada_subgroup", "fpr", "line", "lo_y", "sens", "tfpr",
  "tsens", "x"
))

.dta_bootstrap_state <- new.env(parent = emptyenv())
.dta_bootstrap_state$cluster <- NULL
.dta_bootstrap_state$n.cores <- NULL

.dta_auc_ci_cache_state <- new.env(parent = emptyenv())
.dta_auc_ci_cache_state$entry <- NULL

.onUnload <- function(libpath) {
  if (!is.null(.dta_bootstrap_state$cluster)) {
    backend <- list(
      managed = TRUE,
      cluster = .dta_bootstrap_state$cluster,
      n.cores = .dta_bootstrap_state$n.cores
    )
    try(.dta_stop_bootstrap_backend(backend), silent = TRUE)
  }
  invisible(NULL)
}
