# Standardized public API ---------------------------------------------------
#
# New code should use the dta_* functions in this file. The original public
# functions remain available in legacy-api.R and are not removed or renamed.

.dta_standard_data <- function(data, study_name = NULL,
                               true_positive = NULL, true_negative = NULL,
                               false_positive = NULL, false_negative = NULL,
                               unique_row_id = NULL, study_id = NULL,
                               complete_only = FALSE) {
  .dta_prepare_data(
    data,
    study.names = study_name,
    TP = true_positive,
    TN = true_negative,
    FP = false_positive,
    FN = false_negative,
    uniquer.row.id = unique_row_id,
    study.id = study_id,
    complete.only = complete_only
  )
}

# Resolve a column argument inside the data mask used by the standardized API.
# This permits bare column names (tp = tp), character names (tp = "tp"),
# explicit vectors, expressions, and the older data$tp form.
.dta_column_from_call <- function(data, call, argument, environment,
                                  required = FALSE) {
  supplied <- argument %in% names(call)
  if (!supplied) {
    if (required) {
      stop("'", argument, "' must be supplied.", call. = FALSE)
    }
    return(NULL)
  }

  expression <- call[[argument]]
  if (is.null(expression) || identical(expression, quote(NULL))) {
    if (required) {
      stop("'", argument, "' cannot be NULL.", call. = FALSE)
    }
    return(NULL)
  }
  data_mask <- if (is.data.frame(data)) data else as.data.frame(data)
  value <- tryCatch(
    eval(expression, envir = data_mask, enclos = environment),
    error = function(error) {
      stop(
        "Could not resolve column argument '", argument,
        "' inside 'data': ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (is.null(value)) {
    if (required) {
      stop("'", argument, "' cannot resolve to NULL.", call. = FALSE)
    }
    return(NULL)
  }
  if (is.character(value) && length(value) == 1L &&
      value %in% names(data_mask)) {
    value <- data_mask[[value]]
  }
  if (length(value) != nrow(data_mask)) {
    stop(
      "'", argument,
      "' must be a bare column name, a character column name, or contain ",
      "one value per row of 'data'.",
      call. = FALSE
    )
  }
  value
}

.dta_column_alias_from_call <- function(data, call, preferred, compatible,
                                        environment) {
  supplied <- c(preferred, compatible)[
    c(preferred, compatible) %in% names(call)
  ]
  if (length(supplied) > 1L) {
    stop(
      "Supply only one of '", preferred, "' and '", compatible, "'.",
      call. = FALSE
    )
  }
  if (!length(supplied)) return(NULL)
  .dta_column_from_call(data, call, supplied, environment)
}

.dta_core_columns_from_call <- function(data, call, environment) {
  list(
    study = .dta_column_alias_from_call(
      data, call, "study", "study_name", environment
    ),
    tp = .dta_column_alias_from_call(
      data, call, "tp", "true_positive", environment
    ),
    tn = .dta_column_alias_from_call(
      data, call, "tn", "true_negative", environment
    ),
    fp = .dta_column_alias_from_call(
      data, call, "fp", "false_positive", environment
    ),
    fn = .dta_column_alias_from_call(
      data, call, "fn", "false_negative", environment
    ),
    unique_row_id = .dta_column_from_call(
      data, call, "unique_row_id", environment
    ),
    study_id = .dta_column_from_call(
      data, call, "study_id", environment
    )
  )
}

.dta_internal_columns <- function() {
  list(
    study.names = "names",
    TP = "TP",
    TN = "TN",
    FP = "FP",
    FN = "FN",
    uniquer.row.id = "uniquer.row.id",
    study.id = "study.id"
  )
}

.dta_merge_args <- function(defaults, supplied, locked = character()) {
  if (length(supplied)) {
    supplied <- supplied[!names(supplied) %in% locked]
    defaults[names(supplied)] <- supplied
  }
  defaults
}

.dta_run_graphics <- function(draw, fun) {
  if (isTRUE(draw)) return(fun())
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  fun()
}

.dta_analysis_result <- function(result, data, kind, function_name,
                                 arguments, draw, call) {
  structure(
    list(
      result = result,
      data = data,
      analysis = kind,
      call = call,
      recipe = list(
        function_name = function_name,
        arguments = arguments,
        draw = draw
      ),
      influential_analysis = NULL,
      removed_data = NULL,
      removed_unique_row_ids = NULL,
      removed_studies = NULL,
      bootstrap_after_removal = NULL
    ),
    class = c("dta_analysis", "list")
  )
}

.dta_study_records <- function(data, rows = rep(TRUE, nrow(data))) {
  subgroup <- if (".dta_subgroup" %in% names(data)) {
    as.character(data$.dta_subgroup)
  } else {
    rep(NA_character_, nrow(data))
  }
  data.frame(
    unique_row_id = data$uniquer.row.id[rows],
    study_name = data$names[rows],
    study_id = data$study.id[rows],
    subgroup = subgroup[rows],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' @export
print.dta_analysis <- function(x, ...) {
  cat("DTAtoolkit", x$analysis, "analysis\n")
  cat("Studies:", nrow(x$data), "\n")
  if (inherits(x$result, c("gg", "ggplot", "patchwork"))) {
    print(x$result)
  }
  invisible(x)
}

.dta_match_forest_type <- function(type) {
  key <- tolower(as.character(type)[1])
  aliases <- c(
    combined = "combined", both = "combined",
    sensitivity = "sensitivity", sens = "sensitivity", se = "sensitivity",
    specificity = "specificity", spec = "specificity", sp = "specificity"
  )
  if (is.na(key) || !key %in% names(aliases)) {
    stop("'type' must be 'combined', 'sensitivity' (or 'sens'), or ",
         "'specificity' (or 'spec').", call. = FALSE)
  }
  unname(aliases[[key]])
}

#' Unified forest plot for diagnostic accuracy meta-analysis
#' @export
dta_forest <- function(data,
                       subgroup = NULL,
                       type = c("combined", "sensitivity", "specificity"),
                       study_name = NULL,
                       true_positive = NULL,
                       true_negative = NULL,
                       false_positive = NULL,
                       false_negative = NULL,
                       unique_row_id = NULL,
                       study_id = NULL,
                       sort_by = NULL,
                       subgroup_label = "Subgroup",
                       left_column_1 = NULL,
                       left_label_1 = NULL,
                       left_column_2 = NULL,
                       left_label_2 = NULL,
                       xlim = c(50, 100),
                       plot_heterogeneity = TRUE,
                       plot_subgroup_heterogeneity = TRUE,
                       plot_overall = TRUE,
                       minimum_subgroup_studies = 3L,
                       omnibus_alpha = 0.05,
                       pairwise = c("ask", "always", "never"),
                       p_adjust_method = "holm",
                       draw = TRUE,
                       ...,
                       study = NULL,
                       tp = NULL,
                       tn = NULL,
                       fp = NULL,
                       fn = NULL,
                       sensitivity_width_adjustment = 0,
                       font_size = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  type <- .dta_match_forest_type(type)
  pairwise <- match.arg(pairwise)
  sensitivity_width_adjustment <- .dta_validate_forest_ratio_adjustment(
    sensitivity_width_adjustment
  )
  resolved.font_size <- .dta_resolve_forest_fontsize(font_size)
  if (!is.numeric(minimum_subgroup_studies) ||
      length(minimum_subgroup_studies) != 1L ||
      !is.finite(minimum_subgroup_studies) || minimum_subgroup_studies < 1) {
    stop("'minimum_subgroup_studies' must be one positive integer.",
         call. = FALSE)
  }
  minimum_subgroup_studies <- as.integer(minimum_subgroup_studies)

  columns <- .dta_core_columns_from_call(
    data, captured_call, calling_environment
  )
  subgroup_values <- .dta_column_from_call(
    data, captured_call, "subgroup", calling_environment
  )
  sort_values <- .dta_column_from_call(
    data, captured_call, "sort_by", calling_environment
  )
  left_1 <- .dta_column_from_call(
    data, captured_call, "left_column_1", calling_environment
  )
  left_2 <- .dta_column_from_call(
    data, captured_call, "left_column_2", calling_environment
  )
  prepared <- .dta_standard_data(
    data, columns$study, columns$tp, columns$tn, columns$fp,
    columns$fn, columns$unique_row_id, columns$study_id,
    complete_only = FALSE
  )
  if (!is.null(subgroup_values)) prepared$.dta_subgroup <- subgroup_values
  if (!is.null(sort_values)) prepared$.dta_sort <- sort_values
  if (!is.null(left_1)) prepared$.dta_left_1 <- left_1
  if (!is.null(left_2)) prepared$.dta_left_2 <- left_2

  complete <- stats::complete.cases(
    prepared[, c("TP", "TN", "FP", "FN"), drop = FALSE]
  )
  use_subgroups <- !is.null(subgroup_values)
  if (use_subgroups) {
    counts <- table(prepared$.dta_subgroup[complete])
    retained <- names(counts[counts >= minimum_subgroup_studies])
    prepared <- prepared[
      as.character(prepared$.dta_subgroup) %in% retained,
      , drop = FALSE
    ]
    if (!length(retained)) {
      stop("No subgroup meets 'minimum_subgroup_studies'.", call. = FALSE)
    }
    # A one-level subgroup cannot support an interaction model. It is handled
    # as a pooled analysis while keeping the same studies.
    use_subgroups <- length(retained) > 1L
  }

  column_args <- .dta_internal_columns()
  dots <- list(...)
  if (use_subgroups) {
    common <- c(
      list(
        dat = prepared,
        subgrouping.variable = prepared$.dta_subgroup,
        sortvar = if (is.null(sort_values)) NULL else prepared$.dta_sort,
        sglabel = subgroup_label,
        lcols1 = if (is.null(left_1)) NULL else prepared$.dta_left_1,
        llab1 = left_label_1,
        lcols2 = if (is.null(left_2)) NULL else prepared$.dta_left_2,
        llab2 = left_label_2,
        object.return = TRUE,
        plot.het.overall = plot_heterogeneity,
        plot.het.subgroup = plot_subgroup_heterogeneity,
        plot.overall = plot_overall,
        only.subgroups.bigger.than.3 = FALSE,
        forest.xlim = xlim,
        omnibus.alpha = omnibus_alpha,
        pairwise = pairwise,
        p.adjust.method = p_adjust_method
      ),
      column_args
    )
    if (type == "combined") {
      common$sensitivity.width.adjustment <- sensitivity_width_adjustment
      common$fontsize <- resolved.font_size
      target <- forest.diag.subgroup.combined
    } else {
      common$sens.forest <- type == "sensitivity"
      common$spec.forest <- type == "specificity"
      common$fontsize <- resolved.font_size
      target <- forest.diag.subgroup
    }
  } else {
    common <- c(
      list(
        dat = prepared,
        lcols1 = if (is.null(left_1)) NULL else prepared$.dta_left_1,
        llab1 = left_label_1,
        lcols2 = if (is.null(left_2)) NULL else prepared$.dta_left_2,
        llab2 = left_label_2,
        object.return = TRUE,
        plot.het = plot_heterogeneity,
        xlim = xlim
      ),
      column_args
    )
    if (type == "combined") {
      common$sensitivity.width.adjustment <- sensitivity_width_adjustment
      common$fontsize <- resolved.font_size
      target <- forest.diag.combined
    } else {
      common$sens.forest <- type == "sensitivity"
      common$spec.forest <- type == "specificity"
      common$fontsize <- resolved.font_size
      target <- forest.diag
    }
  }
  call_args <- .dta_merge_args(
    common, dots,
    locked = c("dat", "object.return", "sens.forest", "spec.forest")
  )
  result <- .dta_run_graphics(draw, function() do.call(target, call_args))

  recipe_args <- c(
    list(
      subgroup = if (is.null(subgroup_values)) NULL else ".dta_subgroup",
      type = type,
      study_name = "names",
      true_positive = "TP", true_negative = "TN",
      false_positive = "FP", false_negative = "FN",
      unique_row_id = "uniquer.row.id", study_id = "study.id",
      sort_by = if (is.null(sort_values)) NULL else ".dta_sort",
      subgroup_label = subgroup_label,
      left_column_1 = if (is.null(left_1)) NULL else ".dta_left_1",
      left_label_1 = left_label_1,
      left_column_2 = if (is.null(left_2)) NULL else ".dta_left_2",
      left_label_2 = left_label_2,
      xlim = xlim,
      plot_heterogeneity = plot_heterogeneity,
      plot_subgroup_heterogeneity = plot_subgroup_heterogeneity,
      plot_overall = plot_overall,
      minimum_subgroup_studies = minimum_subgroup_studies,
      omnibus_alpha = omnibus_alpha,
      pairwise = pairwise,
      p_adjust_method = p_adjust_method,
      sensitivity_width_adjustment = sensitivity_width_adjustment,
      font_size = font_size,
      draw = draw
    ),
    dots
  )
  .dta_analysis_result(
    result, prepared, "forest", "dta_forest", recipe_args, draw,
    match.call()
  )
}

#' Inspect the session AUC-CI cache
#' @export
dta_auc_ci_cache_info <- function(include_object = FALSE) {
  if (!is.logical(include_object) || length(include_object) != 1L ||
      is.na(include_object)) {
    stop("'include_object' must be TRUE or FALSE.", call. = FALSE)
  }
  entry <- .dta_auc_ci_cache_state$entry
  if (is.null(entry)) return(NULL)
  info <- entry$metadata
  info$available <- TRUE
  if (isTRUE(include_object)) info$AUC_CIs <- entry$value
  info
}

#' Clear the session AUC-CI cache
#' @export
dta_clear_auc_ci_cache <- function() {
  .dta_auc_ci_cache_clear()
}

#' Standardized SROC analysis and plot
#' @export
dta_sroc <- function(data, subgroup = NULL,
                     study_name = NULL,
                     true_positive = NULL, true_negative = NULL,
                     false_positive = NULL, false_negative = NULL,
                     unique_row_id = NULL, study_id = NULL,
                     color_palette = "mada6",
                     main_title = "SROC curves",
                     auc_ci = FALSE, n_boots = 2000,
                     parallel = TRUE, n_cores = NULL,
                     progress = TRUE, seed = NULL,
                     cache_auc_ci = TRUE,
                     omnibus_alpha = 0.05,
                     pairwise = c("ask", "always", "never"),
                     p_adjust_method = "holm",
                     draw = TRUE, ...,
                     study = NULL,
                     tp = NULL,
                     tn = NULL,
                     fp = NULL,
                     fn = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  pairwise <- match.arg(pairwise)
  columns <- .dta_core_columns_from_call(
    data, captured_call, calling_environment
  )
  subgroup_values <- .dta_column_from_call(
    data, captured_call, "subgroup", calling_environment
  )
  prepared <- .dta_standard_data(
    data, columns$study, columns$tp, columns$tn, columns$fp,
    columns$fn, columns$unique_row_id, columns$study_id,
    complete_only = FALSE
  )
  if (!is.null(subgroup_values)) prepared$.dta_subgroup <- subgroup_values
  dots <- list(...)
  base <- c(
    list(
      dat = prepared,
      subgrouping.variable = if (is.null(subgroup_values)) NULL else prepared$.dta_subgroup,
      main.title = main_title,
      object.return = TRUE,
      AUC.CI = auc_ci,
      n.boots = n_boots,
      AUC.CI.parallel = parallel,
      AUC.CI.n.cores = n_cores,
      AUC.CI.progress = progress,
      AUC.CI.seed = seed,
      AUC.CI.cache = cache_auc_ci,
      omnibus.alpha = omnibus_alpha,
      pairwise = pairwise,
      p.adjust.method = p_adjust_method,
      color.palette = color_palette
    ),
    .dta_internal_columns()
  )
  call_args <- .dta_merge_args(base, dots, locked = c("dat", "object.return"))
  result <- .dta_run_graphics(
    draw, function() do.call(multiple.srocs, call_args)
  )
  recipe_args <- c(
    list(
      subgroup = if (is.null(subgroup_values)) NULL else ".dta_subgroup",
      study_name = "names",
      true_positive = "TP", true_negative = "TN",
      false_positive = "FP", false_negative = "FN",
      unique_row_id = "uniquer.row.id", study_id = "study.id",
      color_palette = color_palette,
      main_title = main_title,
      auc_ci = auc_ci, n_boots = n_boots,
      parallel = parallel, n_cores = n_cores,
      progress = progress, seed = seed,
      cache_auc_ci = cache_auc_ci,
      omnibus_alpha = omnibus_alpha, pairwise = pairwise,
      p_adjust_method = p_adjust_method, draw = draw
    ), dots
  )
  .dta_analysis_result(
    result, prepared, "SROC", "dta_sroc", recipe_args, draw, match.call()
  )
}

.dta_influentials_one <- function(data, method, alpha, p_adjust_method,
                                  plot) {
  madauni_fit <- metafor_fit <- influence_fit <- reitsma_fit <- mean_shift <- NULL
  if (method == "DOR") {
    if (!requireNamespace("metafor", quietly = TRUE)) {
      stop("The 'metafor' package is required for method = 'DOR'.",
           call. = FALSE)
    }
    madauni_fit <- mada::madauni(data, type = "DOR", method = "DSL")
    metafor_fit <- metafor::rma(
      yi = log(madauni_fit$descr$DOR$DOR),
      sei = madauni_fit$descr$DOR$se.lnDOR,
      weights = madauni_fit$weights,
      method = "DL"
    )
    influence_fit <- stats::influence(metafor_fit)
    if (isTRUE(plot)) graphics::plot(influence_fit)
    score <- as.numeric(influence_fit$inf$rstudent)
    diagnostics <- data.frame(
      row_index = data$.dta_source_row,
      unique_row_id = data$uniquer.row.id,
      study_id = data$study.id,
      study_name = data$names,
      statistic = score,
      df = NA_real_, p_value = NA_real_, p_adjusted = NA_real_,
      influential = abs(score) > 2,
      fit_error = NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    mean_shift <- .dta_bivariate_mean_shift(
      data, alpha = alpha, p.adjust.method = p_adjust_method
    )
    reitsma_fit <- mean_shift$null.fit
    source <- mean_shift$diagnostics
    diagnostics <- data.frame(
      row_index = data$.dta_source_row,
      unique_row_id = source$unique_row_id,
      study_id = source$study_id,
      study_name = source$study_name,
      statistic = source$chi_squared,
      df = source$df,
      p_value = source$p_value,
      p_adjusted = source$p_adjusted,
      influential = source$p_adjusted < alpha & !is.na(source$p_adjusted),
      fit_error = source$fit_error,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  list(
    diagnostics = diagnostics,
    madauni = madauni_fit,
    metafor = metafor_fit,
    influence = influence_fit,
    reitsma = reitsma_fit,
    mean_shift = mean_shift
  )
}

#' Influential-study analysis for a bivariate DTA model
#' @export
dta_influentials <- function(data, subgroup = NULL,
                             method = c("Bivariate", "DOR"),
                             alpha = 0.05,
                             p_adjust_method = "holm",
                             study_name = NULL,
                             true_positive = NULL, true_negative = NULL,
                             false_positive = NULL, false_negative = NULL,
                             unique_row_id = NULL, study_id = NULL,
                             plot = FALSE, verbose = TRUE,
                             study = NULL,
                             tp = NULL,
                             tn = NULL,
                             fp = NULL,
                             fn = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  method <- match.arg(method)
  columns <- .dta_core_columns_from_call(
    data, captured_call, calling_environment
  )
  subgroup_values <- .dta_column_from_call(
    data, captured_call, "subgroup", calling_environment
  )
  prepared <- .dta_standard_data(
    data, columns$study, columns$tp, columns$tn, columns$fp,
    columns$fn, columns$unique_row_id, columns$study_id,
    complete_only = FALSE
  )
  prepared$.dta_source_row <- seq_len(nrow(prepared))
  prepared$.original_row <- prepared$.dta_source_row
  if (!is.null(subgroup_values)) prepared$.dta_subgroup <- subgroup_values
  complete <- stats::complete.cases(
    prepared[, c("TP", "TN", "FP", "FN"), drop = FALSE]
  )
  analysed <- prepared[complete, , drop = FALSE]
  minimum <- if (method == "Bivariate") 4L else 3L
  groups <- if (is.null(subgroup_values)) {
    rep("Pooled", nrow(analysed))
  } else {
    as.character(analysed$.dta_subgroup)
  }
  group_counts <- table(groups)
  valid_groups <- sort(names(group_counts[group_counts >= minimum]))
  skipped_groups <- sort(setdiff(names(group_counts), valid_groups))
  if (!length(valid_groups)) {
    stop(method, " influential analysis requires at least ", minimum,
         " complete studies", if (!is.null(subgroup_values)) " per subgroup" else "",
         ".", call. = FALSE)
  }
  if (length(skipped_groups)) {
    warning(
      "Influential analysis was not run for subgroup(s) ",
      paste(skipped_groups, collapse = ", "), " because method = '",
      method, "' requires at least ", minimum,
      " complete studies per subgroup.",
      call. = FALSE
    )
  }

  fits <- diagnostics <- list()
  for (group in valid_groups) {
    group_data <- analysed[which(groups == group), , drop = FALSE]
    fit <- .dta_influentials_one(
      group_data, method, alpha, p_adjust_method, plot
    )
    fit$diagnostics$subgroup <- if (is.null(subgroup_values)) NA_character_ else group
    fits[[group]] <- fit
    diagnostics[[group]] <- fit$diagnostics
  }
  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL
  flagged <- diagnostics$influential %in% TRUE
  influential_studies <- diagnostics[
    flagged,
    c("unique_row_id", "study_name", "study_id", "subgroup",
      "statistic", "p_adjusted"),
    drop = FALSE
  ]

  if (isTRUE(verbose)) {
    if (any(flagged)) {
      for (i in which(flagged)) {
        subgroup_text <- if (is.na(diagnostics$subgroup[i])) "" else {
          paste0(", subgroup: ", diagnostics$subgroup[i])
        }
        p_text <- if (method == "Bivariate") {
          paste0(" (adjusted p = ",
                 format.pval(diagnostics$p_adjusted[i], digits = 3), ")")
        } else ""
        cat(paste0(
          "Potentially influential study - author/study name: ",
          diagnostics$study_name[i],
          ", unique row ID: ", diagnostics$unique_row_id[i],
          subgroup_text, p_text, "\n"
        ))
      }
    }
    cat("Influential analysis method:", method, "\n")
    cat("Number of potentially influential rows:", sum(flagged), "\n")
  }

  structure(
    list(
      data = prepared,
      diagnostics = diagnostics,
      fits = fits,
      influential_indices = diagnostics$row_index[flagged],
      influential_unique_row_ids = diagnostics$unique_row_id[flagged],
      influential_study_ids = diagnostics$study_id[flagged],
      influential_studies = influential_studies,
      method = method,
      alpha = alpha,
      p_adjust_method = p_adjust_method,
      subgroup = if (is.null(subgroup_values)) NULL else ".dta_subgroup",
      subgroup_variable = if (is.null(subgroup_values)) NULL else {
        paste(deparse(captured_call$subgroup), collapse = "")
      },
      scope = if (is.null(subgroup_values)) "pooled" else "within_subgroup",
      analyzed_subgroups = valid_groups,
      skipped_subgroups = skipped_groups
    ),
    class = c("dta_influentials", "list")
  )
}

#' Remove influential studies and regenerate an analysis
#' @export
remove_inf <- function(x,
                       method = c("Bivariate", "DOR"),
                       alpha = 0.05,
                       p_adjust_method = "holm",
                       subgroups = NULL,
                       remove_study_across_subgroups = FALSE,
                       recompute_bootstrap = FALSE,
                       draw = NULL,
                       verbose = TRUE,
                       ...) {
  method <- match.arg(method)
  if (inherits(x, "dta_influentials")) {
    influence <- x
    data <- x$data
    recipe <- NULL
  } else if (inherits(x, "dta_analysis")) {
    data <- x$data
    recipe <- x$recipe
    inherited_subgroup <- recipe$arguments$subgroup
    influence <- dta_influentials(
      data,
      subgroup = inherited_subgroup,
      method = method, alpha = alpha,
      p_adjust_method = p_adjust_method,
      study_name = "names", true_positive = "TP", true_negative = "TN",
      false_positive = "FP", false_negative = "FN",
      unique_row_id = "uniquer.row.id", study_id = "study.id",
      plot = FALSE, verbose = verbose
    )
  } else if (is.data.frame(x)) {
    influence <- dta_influentials(
      x, method = method, alpha = alpha,
      p_adjust_method = p_adjust_method, verbose = verbose, ...
    )
    data <- influence$data
    recipe <- NULL
  } else {
    stop("'x' must be a result from a dta_* analysis, a ",
         "dta_influentials() result, or a diagnostic data frame.",
         call. = FALSE)
  }

  diagnostics <- influence$diagnostics
  selected <- diagnostics$influential %in% TRUE
  if (!is.null(subgroups)) {
    available <- unique(stats::na.omit(diagnostics$subgroup))
    if (is.numeric(subgroups)) {
      if (any(subgroups < 1L | subgroups > length(available))) {
        stop("'subgroups' contains an invalid subgroup index.", call. = FALSE)
      }
      subgroups <- available[subgroups]
    }
    unknown <- setdiff(as.character(subgroups), available)
    if (length(unknown)) {
      stop("Unknown subgroup(s): ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    selected <- selected & diagnostics$subgroup %in% subgroups
  }
  selected_ids <- unique(diagnostics$unique_row_id[selected])
  selected_studies <- unique(diagnostics$study_id[selected])
  remove <- if (isTRUE(remove_study_across_subgroups)) {
    data$study.id %in% selected_studies
  } else {
    data$uniquer.row.id %in% selected_ids
  }
  filtered <- data[!remove, , drop = FALSE]
  removed_studies <- .dta_study_records(data, remove)
  if (isTRUE(verbose) && nrow(removed_studies)) {
    for (i in seq_len(nrow(removed_studies))) {
      subgroup_text <- if (is.na(removed_studies$subgroup[i])) "" else {
        paste0(", subgroup: ", removed_studies$subgroup[i])
      }
      cat(paste0(
        "Removing influential study - author/study name: ",
        removed_studies$study_name[i],
        ", unique row ID: ", removed_studies$unique_row_id[i],
        subgroup_text, "\n"
      ))
    }
  }

  if (is.null(recipe)) return(filtered)
  if (!any(remove)) {
    x$influential_analysis <- influence
    x$removed_data <- data[FALSE, , drop = FALSE]
    x$removed_unique_row_ids <- data$uniquer.row.id[FALSE]
    x$removed_studies <- .dta_study_records(data, rep(FALSE, nrow(data)))
    x$bootstrap_after_removal <- "not needed; no influential rows were removed"
    return(x)
  }

  rerun_args <- recipe$arguments
  rerun_args$data <- filtered
  if (!is.null(draw)) rerun_args$draw <- isTRUE(draw)
  overrides <- list(...)
  if (length(overrides)) rerun_args[names(overrides)] <- overrides

  bootstrap_status <- "not applicable"
  if (identical(x$analysis, "SROC") && isTRUE(rerun_args$auc_ci)) {
    if (isTRUE(recompute_bootstrap)) {
      bootstrap_status <- "recomputed for the reduced dataset"
    } else {
      rerun_args$auc_ci <- FALSE
      bootstrap_status <- paste(
        "skipped for the reduced dataset; set",
        "recompute_bootstrap = TRUE to calculate a new AUC interval"
      )
      if (isTRUE(verbose)) cat("AUC bootstrap", bootstrap_status, "\n")
    }
  }

  target <- get(recipe$function_name, envir = asNamespace("DTAtoolkit"))
  regenerated <- do.call(target, rerun_args)
  regenerated$influential_analysis <- influence
  regenerated$removed_data <- data[remove, , drop = FALSE]
  regenerated$removed_unique_row_ids <- data$uniquer.row.id[remove]
  regenerated$removed_studies <- removed_studies
  regenerated$bootstrap_after_removal <- bootstrap_status
  regenerated
}

#' Standardized likelihood-ratio matrix
#' @export
dta_lr_matrix <- function(data, subgroup = NULL,
                          study_name = NULL,
                          true_positive = NULL, true_negative = NULL,
                          false_positive = NULL, false_negative = NULL,
                          unique_row_id = NULL, study_id = NULL,
                          color_palette = "mada6",
                          legend_position = c("auto", "right", "bottom", "inside", "none"),
                          draw = TRUE, ...,
                          study = NULL,
                          tp = NULL,
                          tn = NULL,
                          fp = NULL,
                          fn = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  legend_position <- match.arg(legend_position)
  columns <- .dta_core_columns_from_call(
    data, captured_call, calling_environment
  )
  subgroup_values <- .dta_column_from_call(
    data, captured_call, "subgroup", calling_environment
  )
  prepared <- .dta_standard_data(
    data, columns$study, columns$tp, columns$tn, columns$fp,
    columns$fn, columns$unique_row_id, columns$study_id,
    complete_only = FALSE
  )
  if (!is.null(subgroup_values)) prepared$.dta_subgroup <- subgroup_values
  dots <- list(...)
  base <- c(
    list(
      dat = prepared,
      subgrouping.variable = if (is.null(subgroup_values)) NULL else prepared$.dta_subgroup,
      color.palette = color_palette,
      legend.position = legend_position
    ),
    .dta_internal_columns()
  )
  call_args <- .dta_merge_args(base, dots, locked = "dat")
  result <- .dta_run_graphics(
    draw, function() do.call(multiple.LRmats, call_args)
  )
  recipe_args <- c(
    list(
      subgroup = if (is.null(subgroup_values)) NULL else ".dta_subgroup",
      study_name = "names",
      true_positive = "TP", true_negative = "TN",
      false_positive = "FP", false_negative = "FN",
      unique_row_id = "uniquer.row.id", study_id = "study.id",
      color_palette = color_palette,
      legend_position = legend_position,
      draw = draw
    ), dots
  )
  .dta_analysis_result(
    result, prepared, "likelihood-ratio matrix", "dta_lr_matrix",
    recipe_args, draw, match.call()
  )
}

#' Standardized data-based Fagan nomogram
#' @export
dta_nomogram <- function(data, prevalence, subgroup = NULL,
                         study_name = NULL,
                         true_positive = NULL, true_negative = NULL,
                         false_positive = NULL, false_negative = NULL,
                         unique_row_id = NULL, study_id = NULL,
                         return_list = FALSE, draw = TRUE, ...,
                         study = NULL,
                         tp = NULL,
                         tn = NULL,
                         fp = NULL,
                         fn = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  columns <- .dta_core_columns_from_call(
    data, captured_call, calling_environment
  )
  subgroup_values <- .dta_column_from_call(
    data, captured_call, "subgroup", calling_environment
  )
  prepared <- .dta_standard_data(
    data, columns$study, columns$tp, columns$tn, columns$fp,
    columns$fn, columns$unique_row_id, columns$study_id,
    complete_only = FALSE
  )
  if (!is.null(subgroup_values)) prepared$.dta_subgroup <- subgroup_values
  dots <- list(...)
  base <- c(
    list(
      dat = prepared,
      prevalence = prevalence,
      return.list = return_list
    ),
    .dta_internal_columns()
  )
  target <- if (is.null(subgroup_values)) {
    nomogrammer_plus
  } else {
    base$subgrouping.variable <- prepared$.dta_subgroup
    nomogrammer_subgroups
  }
  call_args <- .dta_merge_args(base, dots, locked = "dat")
  result <- do.call(target, call_args)
  if (isTRUE(draw)) print(result)
  recipe_args <- c(
    list(
      prevalence = prevalence,
      subgroup = if (is.null(subgroup_values)) NULL else ".dta_subgroup",
      study_name = "names",
      true_positive = "TP", true_negative = "TN",
      false_positive = "FP", false_negative = "FN",
      unique_row_id = "uniquer.row.id", study_id = "study.id",
      return_list = return_list, draw = draw
    ), dots
  )
  .dta_analysis_result(
    result, prepared, "nomogram", "dta_nomogram", recipe_args, draw,
    match.call()
  )
}

#' Standardized publication-bias analysis
#' @export
dta_publication_bias <- function(data, n_boots = 2000,
                                 parallel = TRUE, n_cores = NULL,
                                 progress = TRUE, seed = NULL,
                                 study_name = NULL,
                                 true_positive = NULL, true_negative = NULL,
                                 false_positive = NULL, false_negative = NULL,
                                 unique_row_id = NULL, study_id = NULL,
                                 draw = TRUE,
                                 study = NULL,
                                 tp = NULL,
                                 tn = NULL,
                                 fp = NULL,
                                 fn = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  columns <- .dta_core_columns_from_call(
    data, captured_call, calling_environment
  )
  prepared <- .dta_standard_data(
    data, columns$study, columns$tp, columns$tn, columns$fp,
    columns$fn, columns$unique_row_id, columns$study_id,
    complete_only = FALSE
  )
  call_args <- c(
    list(
      dat = prepared, n.boots = n_boots, parallel = parallel,
      n.cores = n_cores, progress = progress, seed = seed
    ),
    .dta_internal_columns()
  )
  result <- .dta_run_graphics(
    draw, function() do.call(pubbias.diag, call_args)
  )
  recipe_args <- list(
    n_boots = n_boots, parallel = parallel, n_cores = n_cores,
    progress = progress, seed = seed,
    study_name = "names",
    true_positive = "TP", true_negative = "TN",
    false_positive = "FP", false_negative = "FN",
    unique_row_id = "uniquer.row.id", study_id = "study.id",
    draw = draw
  )
  .dta_analysis_result(
    result, prepared, "publication bias", "dta_publication_bias",
    recipe_args, draw, match.call()
  )
}

#' Standardized scalar Fagan nomogram
#' @export
dta_fagan_nomogram <- function(prevalence,
                               sensitivity = NULL,
                               specificity = NULL,
                               positive_lr = NULL,
                               negative_lr = NULL,
                               detail = TRUE,
                               null_line = TRUE,
                               label_size = 12 / 5,
                               verbose = FALSE,
                               annotation_x = 0.75,
                               annotation_y = 2,
                               draw = TRUE) {
  arguments <- list(
    Prevalence = prevalence,
    Detail = detail,
    NullLine = null_line,
    LabelSize = label_size,
    Verbose = verbose,
    x_var = annotation_x,
    y_var = annotation_y
  )
  if (!is.null(sensitivity)) arguments$Sens <- sensitivity
  if (!is.null(specificity)) arguments$Spec <- specificity
  if (!is.null(positive_lr)) arguments$Plr <- positive_lr
  if (!is.null(negative_lr)) arguments$Nlr <- negative_lr
  result <- do.call(nomogrammer, arguments)
  if (isTRUE(draw)) print(result)
  invisible(result)
}

# Data-masked standardized utility wrappers.
dta_exclude_studies <- function(data, studies, study = NULL,
                                study_name = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  study_values <- .dta_column_alias_from_call(
    data, captured_call, "study", "study_name", parent.frame()
  )
  if (is.null(study_values)) {
    study_values <- .dta_column_values(
      data, NULL, "names", "study", required = TRUE
    )
  }
  data[!study_values %in% studies, , drop = FALSE]
}

dta_subgroup_comparisons <- function(
    data, subgroup,
    omnibus_alpha = 0.05,
    pairwise = c("ask", "always", "never"),
    p_adjust_method = "holm",
    digits = 4,
    study_name = NULL,
    true_positive = NULL, true_negative = NULL,
    false_positive = NULL, false_negative = NULL,
    unique_row_id = NULL, study_id = NULL,
    study = NULL, tp = NULL, tn = NULL, fp = NULL, fn = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  columns <- .dta_core_columns_from_call(
    data, captured_call, calling_environment
  )
  subgroup_values <- .dta_column_from_call(
    data, captured_call, "subgroup", calling_environment, required = TRUE
  )
  prepared <- .dta_standard_data(
    data, columns$study, columns$tp, columns$tn, columns$fp,
    columns$fn, columns$unique_row_id, columns$study_id,
    complete_only = FALSE
  )
  reitsma.subgroup.comparisons(
    dat = prepared,
    subgrouping.variable = subgroup_values,
    omnibus.alpha = omnibus_alpha,
    pairwise = match.arg(pairwise),
    p.adjust.method = p_adjust_method,
    digits = digits,
    study.names = "names", TP = "TP", TN = "TN", FP = "FP", FN = "FN",
    uniquer.row.id = "uniquer.row.id", study.id = "study.id"
  )
}

dta_find_repeated_studies <- function(data, study = NULL,
                                      unique_row_id = NULL,
                                      study_id = NULL,
                                      study_name = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  study_values <- .dta_column_alias_from_call(
    data, captured_call, "study", "study_name", calling_environment
  )
  unique_values <- .dta_column_from_call(
    data, captured_call, "unique_row_id", calling_environment
  )
  study_ids <- .dta_column_from_call(
    data, captured_call, "study_id", calling_environment
  )
  find_repeated_studies(
    data,
    uniquer.row.id = unique_values,
    study.id = study_ids,
    study.names = study_values
  )
}

dta_find_repeated_studies_by_subgroup <- function(
    data, subgroup, study = NULL, unique_row_id = NULL,
    study_id = NULL, study_name = NULL) {
  captured_call <- match.call(expand.dots = FALSE)
  calling_environment <- parent.frame()
  subgroup_values <- .dta_column_from_call(
    data, captured_call, "subgroup", calling_environment, required = TRUE
  )
  study_values <- .dta_column_alias_from_call(
    data, captured_call, "study", "study_name", calling_environment
  )
  unique_values <- .dta_column_from_call(
    data, captured_call, "unique_row_id", calling_environment
  )
  study_ids <- .dta_column_from_call(
    data, captured_call, "study_id", calling_environment
  )
  find_repeated_studies_by_subgroup(
    data,
    subgroup_var = subgroup_values,
    uniquer.row.id = unique_values,
    study.id = study_ids,
    study.names = study_values
  )
}

# Remaining standardized utility names are direct aliases because they do not
# accept data-column selectors.
dta_rename_counts <- rename
dta_rename_prefixed_counts <- rename.first
dta_heterogeneity_text <- het.string
dta_auc_bootstrap <- AUC_bootstrap
dta_mvpbt_bootstrap <- MVPBT_bootstrap
dta_mvpbt <- MVPBT_boot
dta_parametric_bootstrap_sample <- PBS3
dta_combine_cohorts <- df.mixer
