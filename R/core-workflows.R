# Legacy-compatible API migrated from My_Mada_functions.R.
#
# Package dependencies are imported through NAMESPACE or referenced with
# package-qualified calls. Loading DTAtoolkit therefore does not attach its
# dependencies to the user's search path.

# Diagnostic-table columns can be supplied either as column names or as vectors.
# The historical defaults (names, TP, TN, FP, FN) remain available. Internally,
# data are standardised to those names because mada/meta expect them.
.dta_column_values <- function(dat, value = NULL, default.column = NULL,
                               argument, required = TRUE) {
  if (is.null(value)) {
    if (!is.null(default.column) && default.column %in% names(dat)) {
      value <- dat[[default.column]]
    } else if (required) {
      stop(
        "'", argument, "' was not supplied and the default column '",
        default.column, "' is not present in 'dat'."
      )
    } else {
      return(NULL)
    }
  } else if (is.character(value) && length(value) == 1L && value %in% names(dat)) {
    value <- dat[[value]]
  }

  if (length(value) != nrow(dat)) {
    stop("'", argument, "' must be a column name or contain one value per row of 'dat'.")
  }
  value
}

.dta_prepare_data <- function(dat, study.names = NULL, TP = NULL, TN = NULL,
                              FP = NULL, FN = NULL, uniquer.row.id = NULL,
                              study.id = NULL, complete.only = TRUE) {
  if (!is.data.frame(dat)) {
    dat <- as.data.frame(dat)
  }

  resolved <- list(
    names = .dta_column_values(dat, study.names, "names", "study.names"),
    TP = .dta_column_values(dat, TP, "TP", "TP"),
    TN = .dta_column_values(dat, TN, "TN", "TN"),
    FP = .dta_column_values(dat, FP, "FP", "FP"),
    FN = .dta_column_values(dat, FN, "FN", "FN")
  )
  resolved$uniquer.row.id <- if (is.null(uniquer.row.id)) {
    seq_len(nrow(dat))
  } else {
    .dta_column_values(dat, uniquer.row.id, NULL, "uniquer.row.id")
  }
  resolved$study.id <- if (is.null(study.id)) {
    resolved$names
  } else {
    .dta_column_values(dat, study.id, NULL, "study.id")
  }

  if (anyNA(resolved$uniquer.row.id) || anyDuplicated(resolved$uniquer.row.id)) {
    stop("'uniquer.row.id' must contain unique, non-missing values.")
  }

  dat[["names"]] <- resolved$names
  dat[["TP"]] <- resolved$TP
  dat[["TN"]] <- resolved$TN
  dat[["FP"]] <- resolved$FP
  dat[["FN"]] <- resolved$FN
  dat[["uniquer.row.id"]] <- resolved$uniquer.row.id
  dat[["study.id"]] <- resolved$study.id

  if (isTRUE(complete.only)) {
    dat <- dat[stats::complete.cases(dat[c("TP", "TN", "FP", "FN")]), , drop = FALSE]
  }
  dat
}

exclude_by_names <- function(dat, study_names, study.name.column = NULL){
  study.values <- .dta_column_values(
    dat, study.name.column, "names", "study.name.column"
  )
  dat_no <- dat[!study.values %in% study_names, , drop = FALSE]
  return(dat_no)
}

rename <- function(dat, option){
  
  for (metric in c('TP', 'TN', 'FP', 'FN')){
    oldvar <- paste(metric, option)
    dat[[metric]] <- round(dat[[oldvar]])
  }
  return(dat)
}

rename.first <- function(dat, option){
  
  for (metric in c('TP', 'TN', 'FP', 'FN')){
    oldvar <- paste(option, metric)
    dat[[metric]] <- round(dat[[oldvar]])
  }
  return(dat)
}


# het.string <- function(reitsma.fit){
#   summary.reitsma <- summary(reitsma.fit)
#   Zhou <- percent(summary.reitsma$i2$Zhou, accuracy = 0.1)
#   Holling.ua <- paste("[",
#                       percent(summary.reitsma$i2$HollingUnadjusted3, accuracy = 0.1),
#                       "-", 
#                       percent(summary.reitsma$i2$HollingUnadjusted2, accuracy = 0.1),
#                       "]",
#                       sep = "")
#   Holling.a <- paste("[", 
#                      percent(summary.reitsma$i2$HollingAdjusted3, accuracy = 0.1),
#                      "-", 
#                      percent(summary.reitsma$i2$HollingAdjusted2, accuracy = 0.1),
#                      "]",
#                      sep = "")
#   combined.string <- paste(
#     "I^2: ",
#     "ZD: ",
#     Zhou,
#     ", Holling(ua): ",
#     Holling.ua,
#     ", Holling(a): ",
#     Holling.a,
#     sep = ""
#   )
#   return(combined.string)
# }

het.string <- function(reitsma.fit){
  summary.reitsma <- summary(reitsma.fit)
  Holling.ua <- paste("(",
                      percent(summary.reitsma$i2$HollingUnadjusted3, accuracy = 0.1),
                      "-", 
                      percent(summary.reitsma$i2$HollingUnadjusted2, accuracy = 0.1),
                      ")",
                      sep = "")
  combined.string <- paste(
    "I^2 range (Holling's): ",
    Holling.ua,
    sep = ""
  )
  return(combined.string)
}


# Test whether sensitivity and false-positive rate differ between subgroups in
# a bivariate Reitsma meta-regression. If the joint omnibus test is significant,
# optional pairwise Wald contrasts are calculated from the full model. A
# specificity contrast is the negative of the corresponding false-positive-
# rate contrast, so both have the same two-sided p-value.
reitsma.subgroup.comparisons <- function(
    dat,
    subgrouping.variable,
    omnibus.alpha = 0.05,
    pairwise = c("ask", "always", "never"),
    p.adjust.method = "holm",
    digits = 4,
    study.names = NULL, TP = NULL, TN = NULL, FP = NULL, FN = NULL,
    uniquer.row.id = NULL, study.id = NULL) {
  pairwise <- match.arg(pairwise)

  if (!is.numeric(omnibus.alpha) || length(omnibus.alpha) != 1L ||
      is.na(omnibus.alpha) || omnibus.alpha <= 0 || omnibus.alpha >= 1) {
    stop("'omnibus.alpha' must be one number strictly between 0 and 1.")
  }
  if (!p.adjust.method %in% p.adjust.methods) {
    stop("Unknown p-value adjustment method: ", p.adjust.method)
  }
  if (length(subgrouping.variable) != nrow(dat)) {
    stop("'subgrouping.variable' must have one value per row of 'dat'.")
  }

  dat[[".dta_subgroup"]] <- subgrouping.variable
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  subgrouping.variable <- dat[[".dta_subgroup"]]

  required <- c("TP", "FN", "FP", "TN")
  missing.columns <- setdiff(required, names(dat))
  if (length(missing.columns) > 0L) {
    stop("Missing required column(s): ", paste(missing.columns, collapse = ", "))
  }

  keep <- complete.cases(dat[, required, drop = FALSE]) &
    !is.na(subgrouping.variable)
  analysis.dat <- dat[keep, , drop = FALSE]
  analysis.dat$.mada_subgroup <- droplevels(factor(subgrouping.variable[keep]))
  subgroup.names <- levels(analysis.dat$.mada_subgroup)

  if (length(subgroup.names) < 2L) {
    stop("At least two non-empty subgroups are required for an interaction test.")
  }

  full.fit <- mada::reitsma(
    data = analysis.dat,
    formula = cbind(tsens, tfpr) ~ .mada_subgroup,
    method = "ml"
  )
  intercept.fit <- mada::reitsma(
    data = analysis.dat,
    formula = cbind(tsens, tfpr) ~ 1,
    method = "ml"
  )
  omnibus <- stats::anova(full.fit, intercept.fit)
  omnibus.p <- unname(omnibus$statistic["Pr(>Chi-squared)"])
  omnibus.significant <- is.finite(omnibus.p) && omnibus.p < omnibus.alpha

  cat("\nOmnibus subgroup interaction (bivariate Reitsma model)\n")
  print(data.frame(
    Chi.squared = unname(omnibus$statistic["Chi-squared"]),
    df = unname(omnibus$statistic["Df"]),
    p.value = omnibus.p,
    check.names = FALSE
  ), row.names = FALSE, digits = digits)

  run.pairwise <- FALSE
  if (omnibus.significant) {
    if (pairwise == "always") {
      run.pairwise <- TRUE
    } else if (pairwise == "ask") {
      if (interactive()) {
        answer <- readline(
          "The omnibus interaction is significant. Run pairwise comparisons? [y/N]: "
        )
        run.pairwise <- tolower(trimws(answer)) %in% c("y", "yes")
      } else {
        message(
          "The omnibus interaction is significant. Pairwise comparisons were not run ",
          "because this is a non-interactive session; use pairwise = \"always\" to run them."
        )
      }
    }
  } else {
    cat("No pairwise comparisons requested because the omnibus interaction was not significant.\n")
  }

  pairwise.results <- NULL
  if (run.pairwise) {
    coefficient.matrix <- stats::coef(full.fit)
    covariance.matrix <- stats::vcov(full.fit)
    model.terms <- rownames(coefficient.matrix)

    group.design <- stats::model.matrix(
      ~ .mada_subgroup,
      data = data.frame(
        .mada_subgroup = factor(subgroup.names, levels = subgroup.names)
      )
    )
    group.design <- group.design[, model.terms, drop = FALSE]
    rownames(group.design) <- subgroup.names

    subgroup.pairs <- utils::combn(subgroup.names, 2L, simplify = FALSE)
    contrast.rows <- lapply(subgroup.pairs, function(pair) {
      contrast <- group.design[pair[1L], ] - group.design[pair[2L], ]

      outcome.row <- function(model.outcome, reported.outcome, sign = 1) {
        parameter.names <- paste0(model.outcome, ".", model.terms)
        if (!all(parameter.names %in% rownames(covariance.matrix))) {
          stop("The Reitsma coefficient and covariance layouts are incompatible.")
        }

        contrast.full <- setNames(
          numeric(nrow(covariance.matrix)),
          rownames(covariance.matrix)
        )
        contrast.full[parameter.names] <- sign * contrast
        estimate <- sum(sign * contrast * coefficient.matrix[, model.outcome])
        variance <- drop(
          t(contrast.full) %*% covariance.matrix %*% contrast.full
        )
        standard.error <- sqrt(variance)
        z.value <- estimate / standard.error

        data.frame(
          comparison = paste(pair[1L], "vs", pair[2L]),
          outcome = reported.outcome,
          logit.difference = estimate,
          odds.ratio = exp(estimate),
          std.error = standard.error,
          z.value = z.value,
          p.value = 2 * stats::pnorm(abs(z.value), lower.tail = FALSE),
          stringsAsFactors = FALSE
        )
      }

      rbind(
        outcome.row("tsens", "Sensitivity"),
        outcome.row("tfpr", "Specificity", sign = -1)
      )
    })

    pairwise.results <- do.call(rbind, contrast.rows)
    pairwise.results$p.adjusted <- ave(
      pairwise.results$p.value,
      pairwise.results$outcome,
      FUN = function(x) stats::p.adjust(x, method = p.adjust.method)
    )
    rownames(pairwise.results) <- NULL

    cat(
      "\nPairwise subgroup comparisons (full bivariate model; ",
      p.adjust.method,
      " adjustment within each outcome)\n",
      sep = ""
    )
    print(pairwise.results, row.names = FALSE, digits = digits)
  }

  invisible(list(
    full.fit = full.fit,
    intercept.fit = intercept.fit,
    omnibus = omnibus,
    omnibus.p.value = omnibus.p,
    omnibus.significant = omnibus.significant,
    pairwise = pairwise.results
  ))
}








forest.diag <- function(dat,
                        lcols1 =NULL,
                        llab1 = NULL,
                        lcols2 =NULL,
                        llab2 = NULL,
                        object.return = F,
                        sens.forest = T,
                        spec.forest = F,
                        plot.het = T,
                        calcwidth.addline.opt = T,
                        xlim = c(50,100),
                        study.names = NULL, TP = NULL, TN = NULL,
                        FP = NULL, FN = NULL, uniquer.row.id = NULL,
                        study.id = NULL,
                        fontsize = NULL)
  {
  fontsize <- .dta_resolve_forest_fontsize(fontsize)
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  reitsma.fit <- reitsma(dat, method ='ml')
  summary.reitsma <- summary(reitsma.fit)
  metaprop.sens <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TP"]],
                            n = dat[["TP"]] + dat[["FN"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "sensitivity",
                            method.ci = "WS",
                            common = FALSE)
  metaprop.spec <- metaprop(data = dat,
                            studlab = dat[["names"]],
                            event = dat[["TN"]],
                            n = dat[["TN"]] + dat[["FP"]], method.tau =  "ml",
                            sm = "PRAW",
                            outclab = "specificity",
                            method.ci = "WS",
                            common = FALSE)
  # These pooled estimates come from the bivariate Reitsma model. Override
  # meta's generic random-effects label so the forest plot describes the model
  # that was actually fitted.
  metaprop.sens$text.random <- "Bivariate model"
  metaprop.spec$text.random <- "Bivariate model"
  
  
  madad.dat <- madad(dat)
  metaprop.sens$upper <- madad.dat$sens$sens.ci[,2]
  metaprop.sens$lower <- madad.dat$sens$sens.ci[,1]
  metaprop.sens$TE.random <- summary.reitsma$coefficients["sensitivity", "Estimate"]
  metaprop.sens$lower.random <- summary.reitsma$coefficients["sensitivity", "95%ci.lb"]
  metaprop.sens$upper.random <- summary.reitsma$coefficients["sensitivity", "95%ci.ub"]
  # metaprop.sens$lcols <- dat[["lcols"]]
  metaprop.spec$upper <- madad.dat$spec$spec.ci[,2]
  metaprop.spec$lower <- madad.dat$spec$spec.ci[,1]
  metaprop.spec$TE.random <- 1- summary.reitsma$coefficients[4,1]
  metaprop.spec$upper.random <- 1- summary.reitsma$coefficients[4,5]
  metaprop.spec$lower.random <- 1- summary.reitsma$coefficients[4,6]
  study.population <- .dta_population(dat)
  metaprop.sens <- .dta_use_population_weights(
    metaprop.sens, study.population
  )
  metaprop.spec <- .dta_use_population_weights(
    metaprop.spec, study.population
  )
  # metaprop.spec$lcols <- dat[["lcols"]]
  combined.set <- list()
  combined.set$reitsma <- reitsma.fit
  combined.set$summary.reitsma <- summary.reitsma
  combined.set$metaprop.sens <- metaprop.sens
  combined.set$metaprop.spec <- metaprop.spec
  if (plot.het){
    plot.string.het.overall <- het.string(reitsma.fit)
  }else{
    plot.string.het.overall <- " "
  }
  
  
  if (is.null(lcols1)){
    lcols1.string <- NULL
  }else{
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)){
    lcols2.string <- NULL
  }else{
    lcols2.string <- "lcols2"
  }
  if (sens.forest){
    meta::forest(metaprop.sens,
           xlim = xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Sensitivity %", "95% C.I.     "),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study", llab1, llab2),
           xlab = "Sensitivity", smlab = "",
           weight.study = "common", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           fontsize = fontsize,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           digits = 1,
           hetstat = FALSE,
           text.addline1 = plot.string.het.overall,
           ref = 100 * combined.set$metaprop.sens$TE.random,
           calcwidth.addline = calcwidth.addline.opt,
           just = "center"
    )
  }
  if (spec.forest){
    meta::forest(metaprop.spec,
           xlim = xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Specificty %", "95% C.I.     "),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study", llab1, llab2),
           xlab = "Specifcity", smlab = "",
           weight.study = "common", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           fontsize = fontsize,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           digits = 1,
           hetstat = FALSE,
           text.addline1 = plot.string.het.overall,
           ref = 100 * combined.set$metaprop.spec$TE.random,
           calcwidth.addline = calcwidth.addline.opt,
           just = "center"
    )
  }
    # Calculate sensitivity estimates and confidence intervals
  sens_estimate <- 100 * combined.set$metaprop.sens$TE.random
  sens_lower <- 100 * combined.set$metaprop.sens$lower.random
  sens_upper <- 100 * combined.set$metaprop.sens$upper.random
  
  # Format to 2 decimal places
  sens_estimate_formatted <- sprintf("%.2f", sens_estimate)
  sens_lower_formatted <- sprintf("%.2f", sens_lower)
  sens_upper_formatted <- sprintf("%.2f", sens_upper)
  
  # Print sensitivity
  cat("Sensitivity:", sens_estimate_formatted, "% (95% CI:", sens_lower_formatted, "-", sens_upper_formatted, "%)\n")
  
  # Calculate specificity estimates and confidence intervals
  spec_estimate <- 100 * combined.set$metaprop.spec$TE.random
  spec_lower <- 100 * combined.set$metaprop.spec$lower.random
  spec_upper <- 100 * combined.set$metaprop.spec$upper.random
  
  # Format to 2 decimal places
  spec_estimate_formatted <- sprintf("%.2f", spec_estimate)
  spec_lower_formatted <- sprintf("%.2f", spec_lower)
  spec_upper_formatted <- sprintf("%.2f", spec_upper)
  
  # Print specificity
  cat("Specificity:", spec_estimate_formatted, "% (95% CI:", spec_lower_formatted, "-", spec_upper_formatted, "%)\n")
  if (object.return){
    return(combined.set)
  }
  }


.dta_forest_grob_width <- function(grob, units = "cm") {
  if (is.null(grob)) return(NA_real_)

  convert_width <- function(x) {
    value <- tryCatch(
      grid::convertWidth(x, units, valueOnly = TRUE),
      error = function(e) NA_real_
    )
    value <- as.numeric(value)[1]
    if (is.finite(value) && value > 0) value else NA_real_
  }

  # grid.grabExpr() returns a gTree whose own grobWidth is zero. The natural
  # width of a meta::forest() plot is held in the layout of its captured root
  # viewport, including text-derived grobwidth units. Traverse the viewport
  # tree so automatic layout also works with meta versions predating
  # meta::forest_dims().
  viewport_widths <- numeric()
  visit_viewports <- function(x) {
    if (is.null(x)) return(invisible(NULL))
    if (inherits(x, "vpTree")) {
      visit_viewports(x$parent)
      visit_viewports(x$children)
      return(invisible(NULL))
    }
    if (inherits(x, "vpList")) {
      for (i in seq_along(x)) visit_viewports(x[[i]])
      return(invisible(NULL))
    }
    if (inherits(x, "viewport")) {
      layout_widths <- x$layout$widths
      if (!is.null(layout_widths)) {
        value <- convert_width(sum(layout_widths))
        if (is.finite(value)) viewport_widths <<- c(viewport_widths, value)
      }
    }
    invisible(NULL)
  }
  visit_viewports(grob$childrenvp)

  candidates <- viewport_widths
  if (!is.null(grob$widths)) {
    candidates <- c(candidates, convert_width(sum(grob$widths)))
  }
  candidates <- c(candidates, convert_width(grid::grobWidth(grob)))
  candidates <- candidates[is.finite(candidates) & candidates > 0]
  if (length(candidates)) max(candidates) else NA_real_
}

.dta_forest_width <- function(meta.object, forest.args, grob = NULL,
                              units = "cm") {
  if ("forest_dims" %in% getNamespaceExports("meta")) {
    forest.dims <- getExportedValue("meta", "forest_dims")
    dims <- tryCatch(
      do.call(
        forest.dims,
        c(list(x = meta.object, units = units), forest.args)
      ),
      error = function(e) NULL
    )
    if (!is.null(dims) && is.finite(dims$width) && dims$width > 0) {
      return(list(width = as.numeric(dims$width), method = "meta::forest_dims"))
    }
  }

  width <- .dta_forest_grob_width(grob, units = units)
  if (is.finite(width) && width > 0) {
    return(list(width = width, method = "captured grob"))
  }

  list(width = NA_real_, method = "unavailable")
}

.dta_forest_dims <- function(meta.object, forest.args, units = "in") {
  if ("forest_dims" %in% getNamespaceExports("meta")) {
    forest.dims <- getExportedValue("meta", "forest_dims")
    dims <- tryCatch(
      do.call(
        forest.dims,
        c(list(x = meta.object, units = units), forest.args)
      ),
      error = function(e) NULL
    )
    if (!is.null(dims) && all(is.finite(c(dims$width, dims$height))) &&
        dims$width > 0 && dims$height > 0) {
      return(dims)
    }
  }
  NULL
}

.dta_resolve_forest_fontsize <- function(fontsize = NULL) {
  if (is.null(fontsize)) {
    fontsize <- tryCatch(as.numeric(meta::gs("fontsize")),
                         error = function(e) 12)
  }
  if (!is.numeric(fontsize) || length(fontsize) != 1L ||
      !is.finite(fontsize) || fontsize <= 0) {
    stop("'fontsize' must be NULL or one positive finite number.",
         call. = FALSE)
  }
  as.numeric(fontsize)
}

.dta_validate_forest_ratio_adjustment <- function(
    sensitivity.width.adjustment = 0) {
  if (!is.numeric(sensitivity.width.adjustment) ||
      length(sensitivity.width.adjustment) != 1L ||
      is.na(sensitivity.width.adjustment) ||
      !is.finite(sensitivity.width.adjustment)) {
    stop(
      "'sensitivity.width.adjustment' must be one finite number.",
      call. = FALSE
    )
  }
  as.numeric(sensitivity.width.adjustment)
}

.dta_adjust_forest_ratio <- function(ratio,
                                     sensitivity.width.adjustment = 0) {
  sensitivity.width.adjustment <- .dta_validate_forest_ratio_adjustment(
    sensitivity.width.adjustment
  )
  ratio <- ratio / sum(ratio)
  change <- sensitivity.width.adjustment / 100
  adjusted <- c(
    sensitivity = unname(ratio[1]) + change,
    specificity = unname(ratio[2]) - change
  )
  if (any(adjusted <= 0) || any(adjusted >= 1)) {
    lower <- -100 * unname(ratio[1])
    upper <- 100 * unname(ratio[2])
    stop(
      "'sensitivity.width.adjustment' must be greater than ",
      format(lower), " and less than ", format(upper),
      " percentage points for the calculated ratio.",
      call. = FALSE
    )
  }
  adjusted
}

.dta_capture_forest <- function(plot.function, dims = NULL,
                                height = NULL) {
  muffle_display_warning <- function(w) {
    if (grepl(
      "unsupported operation on the graphics display list",
      conditionMessage(w),
      fixed = TRUE
    )) {
      invokeRestart("muffleWarning")
    }
  }
  if (!is.null(dims)) {
    if (is.null(height)) height <- dims$height
    return(withCallingHandlers(
      grid::grid.grabExpr(
        plot.function(),
        warn = 0,
        width = dims$width,
        height = height
      ),
      warning = muffle_display_warning
    ))
  }
  withCallingHandlers(
    ggplotify::as.grob(plot.function),
    warning = muffle_display_warning
  )
}

.dta_fit_forests_to_device <- function(meta.sens, sens.args,
                                       meta.spec, spec.args,
                                       grob.sens = NULL,
                                       grob.spec = NULL,
                                       enabled = TRUE,
                                       device.padding.in = 0.15,
                                       minimum.scale = 0.6) {
  sens.dims <- .dta_forest_dims(meta.sens, sens.args)
  spec.dims <- .dta_forest_dims(meta.spec, spec.args)
  device.width <- tryCatch(
    as.numeric(grDevices::dev.size("in")[1]),
    error = function(e) NA_real_
  )
  available.width <- device.width - 2 * device.padding.in
  natural.widths <- c(
    sensitivity = if (!is.null(sens.dims)) sens.dims$width else
      .dta_forest_grob_width(grob.sens, "in"),
    specificity = if (!is.null(spec.dims)) spec.dims$width else
      .dta_forest_grob_width(grob.spec, "in")
  )
  natural.width <- if (all(is.finite(natural.widths))) {
    sum(natural.widths)
  } else {
    NA_real_
  }
  scale <- 1

  if (isTRUE(enabled) && is.finite(available.width) &&
      is.finite(natural.width) && natural.width > available.width) {
    scale <- max(minimum.scale, available.width / natural.width)
    sens.fontsize <- .dta_resolve_forest_fontsize(sens.args$fontsize)
    spec.fontsize <- .dta_resolve_forest_fontsize(spec.args$fontsize)
    sens.args$fontsize <- sens.fontsize * scale
    spec.args$fontsize <- spec.fontsize * scale
    scale_plotwidth <- function(plotwidth) {
      width.cm <- if (is.null(plotwidth)) {
        6
      } else {
        tryCatch(
          grid::convertWidth(plotwidth, "cm", valueOnly = TRUE),
          error = function(e) 6
        )
      }
      grid::unit(as.numeric(width.cm)[1] * scale, "cm")
    }
    sens.args$plotwidth <- scale_plotwidth(sens.args$plotwidth)
    spec.args$plotwidth <- scale_plotwidth(spec.args$plotwidth)
    if (!is.null(sens.args$fs.hetstat)) {
      sens.args$fs.hetstat <- sens.args$fs.hetstat * scale
    }
    if (!is.null(spec.args$fs.hetstat)) {
      spec.args$fs.hetstat <- spec.args$fs.hetstat * scale
    }
    sens.dims <- .dta_forest_dims(meta.sens, sens.args)
    spec.dims <- .dta_forest_dims(meta.spec, spec.args)
  }

  list(
    sens.args = sens.args,
    spec.args = spec.args,
    sens.dims = sens.dims,
    spec.dims = spec.dims,
    scale = scale,
    natural.widths.in = natural.widths,
    device.width.in = device.width,
    available.width.in = available.width
  )
}

.dta_capture_and_fit_forests <- function(meta.sens, sens.args,
                                         meta.spec, spec.args,
                                         enabled = TRUE,
                                         minimum.scale = 0.6,
                                         maximum.iterations = 3L) {
  capture_pair <- function() {
    sens.dims <- .dta_forest_dims(meta.sens, sens.args)
    spec.dims <- .dta_forest_dims(meta.spec, spec.args)
    common.height <- if (!is.null(sens.dims) && !is.null(spec.dims)) {
      max(sens.dims$height, spec.dims$height)
    } else {
      NULL
    }
    list(
      sens.dims = sens.dims,
      spec.dims = spec.dims,
      grob.sens = .dta_capture_forest(
        function() do.call(meta::forest, c(list(x = meta.sens), sens.args)),
        sens.dims, common.height
      ),
      grob.spec = .dta_capture_forest(
        function() do.call(meta::forest, c(list(x = meta.spec), spec.args)),
        spec.dims, common.height
      )
    )
  }

  captured <- capture_pair()
  cumulative.scale <- 1
  device.fit <- NULL
  for (iteration in seq_len(maximum.iterations)) {
    # The first pass respects the public readability floor. Later passes only
    # correct the small residual caused by text and fixed gaps not scaling
    # perfectly with font size.
    iteration.minimum <- if (iteration == 1L) minimum.scale else 0.1
    device.fit <- .dta_fit_forests_to_device(
      meta.sens = meta.sens,
      sens.args = sens.args,
      meta.spec = meta.spec,
      spec.args = spec.args,
      grob.sens = captured$grob.sens,
      grob.spec = captured$grob.spec,
      enabled = enabled,
      minimum.scale = iteration.minimum
    )
    if (!isTRUE(enabled) || device.fit$scale >= 0.995) break

    sens.args <- device.fit$sens.args
    spec.args <- device.fit$spec.args
    cumulative.scale <- cumulative.scale * device.fit$scale
    captured <- capture_pair()
  }

  if (is.null(device.fit)) {
    stop("Internal forest device fitting failed.", call. = FALSE)
  }
  synced <- .dta_sync_forest_heights(
    captured$grob.sens, captured$grob.spec
  )

  list(
    sens.args = sens.args,
    spec.args = spec.args,
    sens.dims = captured$sens.dims,
    spec.dims = captured$spec.dims,
    grob.sens = synced$grob1,
    grob.spec = synced$grob2,
    scale = cumulative.scale,
    device.width.in = device.fit$device.width.in,
    available.width.in = device.fit$available.width.in
  )
}

.dta_draw_combined_forests <- function(grob.sens, grob.spec,
                                       forest.layout,
                                       inner.trim.cm = 0) {
  auto.widths <- forest.layout$mode == "auto" &&
    all(is.finite(forest.layout$widths.cm))

  if (auto.widths) {
    panel.widths <- pmax(
      forest.layout$widths.cm - inner.trim.cm,
      1
    )
    if (!is.null(forest.layout$ratio.adjustment) &&
        forest.layout$ratio.adjustment != 0) {
      panel.widths <- sum(panel.widths) * forest.layout$ratio
    }
    arranged <- gridExtra::arrangeGrob(
      grob.sens, grob.spec,
      ncol = 2,
      widths = grid::unit(panel.widths, "cm")
    )
    total.width <- sum(panel.widths)
    device.width <- tryCatch(
      2.54 * as.numeric(grDevices::dev.size("in")[1]),
      error = function(e) NA_real_
    )
    viewport.width <- if (is.finite(device.width) && total.width < device.width) {
      grid::unit(total.width, "cm")
    } else {
      grid::unit(1, "npc")
    }
    combined <- grid::grobTree(
      arranged,
      vp = grid::viewport(width = viewport.width, height = grid::unit(1, "npc"))
    )
  } else {
    arranged <- gridExtra::arrangeGrob(
      grob.sens, grob.spec,
      ncol = 2,
      widths = forest.layout$ratio
    )
    combined <- arranged
  }

  grid::grid.newpage()
  grid::grid.draw(combined)
  combined
}

.dta_sync_forest_heights <- function(grob1, grob2) {
  if (!is.null(grob1$heights) && !is.null(grob2$heights) &&
      length(grob1$heights) == length(grob2$heights)) {
    common.heights <- grid::unit.pmax(grob1$heights, grob2$heights)
    grob1$heights <- common.heights
    grob2$heights <- common.heights
  }
  list(grob1 = grob1, grob2 = grob2)
}

.dta_resolve_forest_layout <- function(meta.sens, sens.args, grob.sens,
                                       meta.spec, spec.args, grob.spec,
                                       layout.mode = c("auto", "manual"),
                                       ratio = NULL) {
  layout.mode <- match.arg(layout.mode)

  # Supplying the legacy ratio remains an explicit manual override.
  if (!is.null(ratio)) {
    if (!is.numeric(ratio) || length(ratio) != 2L ||
        any(!is.finite(ratio)) || any(ratio <= 0)) {
      stop("'ratio' must contain two positive finite numbers.", call. = FALSE)
    }
    layout.mode <- "manual"
  }

  if (layout.mode == "manual") {
    if (is.null(ratio)) ratio <- c(0.7, 0.3)
    ratio <- ratio / sum(ratio)
    return(list(
      mode = "manual",
      ratio = ratio,
      widths.cm = c(sensitivity = NA_real_, specificity = NA_real_),
      measurement = "manual"
    ))
  }

  sens.width <- .dta_forest_width(meta.sens, sens.args, grob.sens)
  spec.width <- .dta_forest_width(meta.spec, spec.args, grob.spec)
  widths <- c(sensitivity = sens.width$width, specificity = spec.width$width)

  if (any(!is.finite(widths)) || any(widths <= 0)) {
    warning(
      "Automatic forest-width measurement failed; using c(0.7, 0.3).",
      call. = FALSE
    )
    widths <- c(sensitivity = NA_real_, specificity = NA_real_)
    ratio <- c(0.7, 0.3)
    method <- "fallback"
  } else {
    ratio <- widths / sum(widths)
    method <- paste(unique(c(sens.width$method, spec.width$method)), collapse = "; ")
  }

  list(
    mode = "auto",
    ratio = ratio,
    widths.cm = widths,
    measurement = method
  )
}


forest.diag.combined <- function(dat,
                                 lcols1 = NULL,
                                 llab1 = NULL,
                                 lcols2 = NULL,
                                 llab2 = NULL,
                                 object.return = FALSE,
                                 plot.het = TRUE,
                                 calcwidth.addline.opt = TRUE,
                                 xlim = c(50, 100),
                                 leftspace = "    ",
                                 rightspace = "        ",
                                 ratio = NULL,
                                 layout.mode = c("auto", "manual"),
                                 auto.fit.device = TRUE,
                                 auto.minimum.scale = 0.6,
                                 auto.inner.trim.cm = 0,
                                 study.names = NULL, TP = NULL, TN = NULL,
                                 FP = NULL, FN = NULL, uniquer.row.id = NULL,
                                 study.id = NULL,
                                 sensitivity.width.adjustment = 0,
                                 fontsize = NULL,
                                 ...
                                 ) {
  layout.mode <- match.arg(layout.mode)
  fontsize <- .dta_resolve_forest_fontsize(fontsize)
  # Prepare data
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  
  # Fit Reitsma model
  reitsma.fit <- reitsma(dat, method = 'ml')
  summary.reitsma <- summary(reitsma.fit)
  
  # Meta-analysis for sensitivity
  metaprop.sens <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TP"]],
    n = dat[["TP"]] + dat[["FN"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Sensitivity",
    method.ci = "WS",
    common = FALSE
  )
  
  # Meta-analysis for specificity
  metaprop.spec <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TN"]],
    n = dat[["TN"]] + dat[["FP"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Specificity",
    method.ci = "WS",
    common = FALSE
  )
  
  # Update metaprop objects with confidence intervals from mada
  madad.dat <- madad(dat)
  
  # Sensitivity
  metaprop.sens$upper <- madad.dat$sens$sens.ci[, 2]
  metaprop.sens$lower <- madad.dat$sens$sens.ci[, 1]
  metaprop.sens$TE.random <- summary.reitsma$coefficients["sensitivity", "Estimate"]
  metaprop.sens$lower.random <- summary.reitsma$coefficients["sensitivity", "95%ci.lb"]
  metaprop.sens$upper.random <- summary.reitsma$coefficients["sensitivity", "95%ci.ub"]
  
  # Specificity (using numeric indices as per your instructions)
  metaprop.spec$upper <- madad.dat$spec$spec.ci[, 2]
  metaprop.spec$lower <- madad.dat$spec$spec.ci[, 1]
  metaprop.spec$TE.random <- 1 - summary.reitsma$coefficients[4, 1]
  metaprop.spec$lower.random <- 1 - summary.reitsma$coefficients[4, 6]
  metaprop.spec$upper.random <- 1 - summary.reitsma$coefficients[4, 5]
  study.population <- .dta_population(dat)
  metaprop.sens <- .dta_use_population_weights(
    metaprop.sens, study.population
  )
  metaprop.spec <- .dta_use_population_weights(
    metaprop.spec, study.population
  )
  
  # Prepare left column labels
  if (is.null(lcols1)) {
    lcols1.string <- NULL
  } else {
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)) {
    lcols2.string <- NULL
  } else {
    lcols2.string <- "lcols2"
  }
  
  leftcols_sens <- c("studlab", lcols1.string, lcols2.string)
  leftlabs_sens <- c("Study", llab1, llab2)
  
  # Remove NULL values from leftcols and leftlabs
  leftcols_sens <- leftcols_sens[!sapply(leftcols_sens, is.null)]
  leftlabs_sens <- leftlabs_sens[!sapply(leftlabs_sens, is.null)]
  
  # Generate the sensitivity forest plot and capture it as a grob
  metaprop.sens$text.random <- "Bivariate model"
  sens.forest.args <- list(
    xlim = xlim,
    pscale = 100,
    just.addcols.right = "center",
    rightcols = c("effect", "ci"),
    rightlabs = c("Sensitivity %", "95% C.I."),
    leftcols = leftcols_sens,
    leftlabs = leftlabs_sens,
    xlab = "Sensitivity",
    smlab = "",
    weight.study = "common",
    squaresize = 0.7,
    col.square = "navy",
    col.square.lines = "navy",
    col.diamond = "maroon",
    col.diamond.lines = "maroon",
    pooled.totals = FALSE,
    comb.fixed = FALSE,
    fs.hetstat = 10,
    fontsize = fontsize,
    print.tau2 = FALSE,
    print.Q = FALSE,
    print.pval.Q = FALSE,
    print.I2 = FALSE,
    digits = 1,
    hetstat = FALSE,
    text.addline1 = if (plot.het) het.string(reitsma.fit) else " ",
    ref = 100 * metaprop.sens$TE.random,
    calcwidth.addline = calcwidth.addline.opt,
    just = "center"
  )
  # Generate the specificity forest plot and capture it as a grob
  
  manual.layout <- layout.mode == "manual" || !is.null(ratio)
  if (manual.layout) {
    metaprop.spec$leftspace <- rep(leftspace, length(metaprop.spec$studlab))
    metaprop.spec$rightspace <- rep(rightspace, length(metaprop.spec$studlab))
  }
  metaprop.spec$text.random <- ""
  spec.forest.args <- list(
    xlim = xlim,
    pscale = 100,
    just.addcols.right = "center",
    rightcols = if (manual.layout) c("effect", "ci", "rightspace") else c("effect", "ci"),
    rightlabs = if (manual.layout) c("Specificity %", "95% C.I.", "") else c("Specificity %", "95% C.I."),
    leftcols = if (manual.layout) "leftspace" else FALSE,
    leftlabs = NULL,
    studlab = FALSE,
    xlab = "Specificity",
    smlab = "",
    weight.study = "common",
    squaresize = 0.7,
    col.square = "navy",
    col.square.lines = "navy",
    col.diamond = "maroon",
    col.diamond.lines = "maroon",
    pooled.totals = FALSE,
    comb.fixed = FALSE,
    fs.hetstat = 10,
    fontsize = fontsize,
    print.tau2 = FALSE,
    print.Q = FALSE,
    print.pval.Q = FALSE,
    print.I2 = FALSE,
    digits = 1,
    hetstat = FALSE,
    text.addline1 = "",
    ref = 100 * metaprop.spec$TE.random,
    calcwidth.addline = calcwidth.addline.opt,
    just = "center"
  )
  device.fit <- .dta_capture_and_fit_forests(
    meta.sens = metaprop.sens,
    sens.args = sens.forest.args,
    meta.spec = metaprop.spec,
    spec.args = spec.forest.args,
    enabled = auto.fit.device && layout.mode == "auto" && is.null(ratio),
    minimum.scale = auto.minimum.scale
  )
  sens.forest.args <- device.fit$sens.args
  spec.forest.args <- device.fit$spec.args
  grob.sens <- device.fit$grob.sens
  grob.spec <- device.fit$grob.spec

  forest.layout <- .dta_resolve_forest_layout(
    meta.sens = metaprop.sens,
    sens.args = sens.forest.args,
    grob.sens = grob.sens,
    meta.spec = metaprop.spec,
    spec.args = spec.forest.args,
    grob.spec = grob.spec,
    layout.mode = layout.mode,
    ratio = ratio
  )
  forest.layout$device.scale <- device.fit$scale
  forest.layout$device.width.in <- device.fit$device.width.in
  forest.layout$available.width.in <- device.fit$available.width.in
  forest.layout$base.ratio <- forest.layout$ratio
  forest.layout$ratio <- .dta_adjust_forest_ratio(
    forest.layout$ratio, sensitivity.width.adjustment
  )
  forest.layout$ratio.adjustment <- sensitivity.width.adjustment
  forest.layout$font.size <- fontsize

  if (forest.layout$mode == "auto") {
    layout.note <- forest.layout$measurement
    if (is.finite(forest.layout$device.scale) &&
        forest.layout$device.scale < 0.999) {
      layout.note <- paste0(
        layout.note,
        "; device scale ",
        sprintf("%.1f%%", 100 * forest.layout$device.scale)
      )
    }
    cat(
      "Automatic forest layout:",
      sprintf("sensitivity %.1f%%, specificity %.1f%%",
              100 * forest.layout$ratio[1], 100 * forest.layout$ratio[2]),
      sprintf("(%s)\n", layout.note)
    )
  }

  # Combine the two grobs in a fixed-width centred viewport. Extra device
  # width is therefore placed at the outer margins, not between the panels.
  combined <- .dta_draw_combined_forests(
    grob.sens, grob.spec,
    forest.layout = forest.layout,
    inner.trim.cm = auto.inner.trim.cm
  )
  
  # Print overall estimates
  sens_estimate <- 100 * metaprop.sens$TE.random
  sens_lower <- 100 * metaprop.sens$lower.random
  sens_upper <- 100 * metaprop.sens$upper.random
  
  spec_estimate <- 100 * metaprop.spec$TE.random
  spec_lower <- 100 * metaprop.spec$lower.random
  spec_upper <- 100 * metaprop.spec$upper.random
  
  cat("Sensitivity:", sprintf("%.2f", sens_estimate), "% (95% CI:", sprintf("%.2f", sens_lower), "-", sprintf("%.2f", sens_upper), "%)\n")
  cat("Specificity:", sprintf("%.2f", spec_estimate), "% (95% CI:", sprintf("%.2f", spec_lower), "-", sprintf("%.2f", spec_upper), "%)\n")
  
  # Return objects if requested
  if (object.return) {
    combined.set <- list(
      reitsma = reitsma.fit,
      summary.reitsma = summary.reitsma,
      metaprop.sens = metaprop.sens,
      metaprop.spec = metaprop.spec,
      layout = forest.layout,
      plot = combined
    )
    return(combined.set)
  }
}



forest.diag.subgroup <- function (dat, # Diagnostic-accuracy dataframe; columns are configurable below.
                                  subgrouping.variable, #dat$subgrouping.variable
                                  sortvar = NULL ,
                                  sglabel = "subgroup" , #string for subgroup labe;
                                  lcols1 = NULL, #Variable to be shown in left side dat$lcols
                                  llab1 = NULL, #a string for the label of lcols
                                  lcols2 = NULL, #Variable to be shown in left side dat$lcols
                                  llab2 = NULL, #a string for the label of lcols
                                  object.return = T, # IF True it returns an object any object used in the function
                                  sens.forest = T, #If true, it draws  sensitivity forest plot
                                  spec.forest = F,# If true it draws specificity forrest plot
                                  only.subgroups.bigger.than.3 = T, # If true it excludes subgroups smaller than 3 subjets.
                                  plot.het.overall =T,
                                  plot.het.subgroup =T,
                                  plot.overall = T,
                                  calcwidth.shet.opt =F,
                                  forest.xlim = c(50,100),
                                  omnibus.alpha = 0.05,
                                  pairwise = c("ask", "always", "never"),
                                  p.adjust.method = "holm",
                                  study.names = NULL, TP = NULL, TN = NULL,
                                  FP = NULL, FN = NULL, uniquer.row.id = NULL,
                                  study.id = NULL,
                                  fontsize = NULL
                                  )
{
  pairwise <- match.arg(pairwise)
  fontsize <- .dta_resolve_forest_fontsize(fontsize)
  dat[["subgrouping.variable"]] <- subgrouping.variable
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  if (!is.null(sortvar)) dat[[".dta_sortvar"]] <- sortvar
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  if (is.null(sortvar)) {
    dat <- dat[order(dat[["subgrouping.variable"]]), ]
  } else {
    dat <- dat[order(dat[["subgrouping.variable"]], dat[[".dta_sortvar"]]), ]
  }
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  if (only.subgroups.bigger.than.3){
    subgroup.list <- names(counts[counts >= 3])
  }
  subgroup.list <- sort(subgroup.list)
  valid.subgroup.list <- make.names(subgroup.list, unique = T)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  reitsmas <- list()
  reitsmas$reitsma.overall <- reitsma(data = dat, method = "ml")
  reitsmas$subgroup.tests <- reitsma.subgroup.comparisons(
    dat = dat,
    subgrouping.variable = dat[["subgrouping.variable"]],
    omnibus.alpha = omnibus.alpha,
    pairwise = pairwise,
    p.adjust.method = p.adjust.method
  )
  reitsmas$reitsma.reg.fit <- reitsmas$subgroup.tests$full.fit
  reitsmas$reitsma.intercept <- reitsmas$subgroup.tests$intercept.fit
  reitsmas$anova.reitsma <- reitsmas$subgroup.tests$omnibus
  reitsmas$subgroups <- list()
  summaries <- list()
  summaries$reitsma.overalll <- summary(reitsmas$reitsma.overall)
  summaries$reitsma.reg.fit <- summary(reitsmas$reitsma.reg.fit)
  summaries$subgroups <- list() 
  madads <- list()
  madads$overall <- madad(dat)
  madads$subgroups <- list()
  props <- list()
  props$sens <- list()
  props$spec <- list()
  props$sens$overall <-  metaprop(data = dat,
                                  studlab = dat[["names"]],
                                  event = dat[["TP"]],
                                  n = dat[["TP"]] + dat[["FN"]], method.tau =  "ml",
                                  sm = "PRAW",
                                  outclab = "sensitivity",
                                  method.ci = "WS",
                                  common = FALSE,
                                  subgroup = dat[["subgrouping.variable"]],
                                  subgroup.name = sglabel)
  props$sens$overall$upper <- madads$overall$sens$sens.ci[,2]
  props$sens$overall$lower <- madads$overall$sens$sens.ci[,1]
  props$sens$overall$TE.random <- summaries$reitsma.overall$coefficients[3,1]
  props$sens$overall$lower.random <- summaries$reitsma.overall$coefficients[3,5]
  props$sens$overall$upper.random <- summaries$reitsma.overall$coefficients[3,6]
  props$spec$overall <- metaprop(data = dat,
                                 studlab = dat[["names"]],
                                 event = dat[["TN"]],
                                 n = dat[["FP"]] + dat[["TN"]], method.tau =  "ml",
                                 sm = "PRAW",
                                 outclab = "specificity",
                                 method.ci = "WS",
                                 common = FALSE,
                                 subgroup = dat[["subgrouping.variable"]],
                                 subgroup.name = sglabel,
  )
  props$spec$overall$upper <- madads$overall$spec$spec.ci[,2]
  props$spec$overall$lower <- madads$overall$spec$spec.ci[,1]
  props$spec$overall$TE.random <- 1 - summaries$reitsma.overall$coefficients[4,1]
  props$spec$overall$upper.random <- 1 - summaries$reitsma.overall$coefficients[4,5]
  props$spec$overall$lower.random <- 1 - summaries$reitsma.overall$coefficients[4,6]
  study.population <- .dta_population(dat)
  props$sens$overall <- .dta_use_population_weights(
    props$sens$overall, study.population
  )
  props$spec$overall <- .dta_use_population_weights(
    props$spec$overall, study.population
  )
  hets.list <- list()
  for (sg in 1:length(subgroup.list)){
    reitsma.sg <- reitsma(data = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ], method = 'ml')
    madad.sg <- madad(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ])
    name <- valid.subgroup.list[sg]
    reitsmas$subgroups[[name]] <- reitsma.sg 
    madads$subgroups[[name]] <- madad.sg
    summary.sg <- summary(reitsma.sg)
    summaries$subgroups[[name]] <- summary.sg
    props$sens$overall$TE.random.w[subgroup.list[sg]] <- summary.sg$coefficients[3,1]
    props$sens$overall$upper.random.w[subgroup.list[sg]] <- summary.sg$coefficients[3,6]
    props$sens$overall$lower.random.w[subgroup.list[sg]] <- summary.sg$coefficients[3,5]
    props$spec$overall$TE.random.w[subgroup.list[sg]] <- 1 - summary.sg$coefficients[4,1]
    props$spec$overall$upper.random.w[subgroup.list[sg]] <- 1 - summary.sg$coefficients[4,5]
    props$spec$overall$lower.random.w[subgroup.list[sg]] <- 1 - summary.sg$coefficients[4,6]
    if (plot.het.subgroup) {
      hets.list[sg] <- het.string(reitsma.sg)
    }
  }
  if (plot.het.overall){
    overall.hetstring <- het.string(reitsmas$reitsma.overall)
  } else {
    overall.hetstring <- " "
  }
  if (!plot.het.subgroup){
    hets.list <- "Subgroup pooled effect"
  }
  if (is.null(lcols1)){
    lcols1.string <- NULL
  }else{
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)){
    lcols2.string <- NULL
  }else{
    lcols2.string <- "lcols2"
  }
  if (sens.forest){
    meta::forest(props$sens$overall,
           xlim = forest.xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Sensitivity %", "95% C.I."),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study                                                   ", llab1 , llab2),
           xlab = "Sensitivity", smlab = "",
           weight.study = "common", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           fontsize = fontsize,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           # overall = T,
           lty.random = 0,
           digits = 1,
           subgroup.hetstat = F,
           hetstat = FALSE,
           test.subgroup = FALSE,
           text.addline2 =  if (plot.overall) het.string(reitsmas$reitsma.overall) else " ",
           text.addline1 = paste("Between-group difference (p): ", round(reitsmas$anova.reitsma$statistic[3], digits = 3)),
           text.random.w = hets.list,
           text.random = "Bivariate model",
           ref = if (plot.overall) 100 * props$sens$overall$TE.random else NA,
           calcwidth.random = calcwidth.shet.opt,
           just = "center",
           overall=plot.overall
    
    )
  }
  
  if (spec.forest){
    meta::forest(props$spec$overall,
           xlim = forest.xlim,
           pscale = 100,
           just.addcols.right = "center",
           rightcols = c("effect", "ci"),
           rightlabs = c("Specificity %", "95% C.I."),
           leftcols = c("studlab", lcols1.string, lcols2.string),
           leftlabs = c("Study                                                   ", llab1, llab2),
           xlab = "Specifcity", smlab = "",
           weight.study = "common", squaresize = 0.7, col.square = "navy",
           col.square.lines = "navy",
           col.diamond = "maroon", col.diamond.lines = "maroon",
           pooled.totals = FALSE,
           comb.fixed = FALSE,
           fs.hetstat = 10,
           fontsize = fontsize,
           print.tau2 = FALSE,
           print.Q = FALSE,
           print.pval.Q = FALSE,
           print.I2 = FALSE,
           # overall = T,
           lty.random = 0,
           digits = 1,
           subgroup.hetstat = F,
           hetstat = FALSE,
           test.subgroup = FALSE,
           text.addline2 =  if (plot.overall) het.string(reitsmas$reitsma.overall) else " ",
           text.addline1 = paste("Between-group difference: p = ", round(reitsmas$anova.reitsma$statistic[3], digits = 3)),
           text.random.w = hets.list,
           text.random = "Bivariate model",
           ref = if (plot.overall) 100 * props$spec$overall$TE.random else NA,
           calcwidth.random = calcwidth.shet.opt,
           just = "center",
           overall = plot.overall
    )
  }



  # Loop over subgroups
for (sg in 1:length(subgroup.list)) {
    subgroup_name <- subgroup.list[sg]
    valid_name <- valid.subgroup.list[sg]
    
    # Extract sensitivity estimates and confidence intervals
    sens_estimate <- 100 * props$sens$overall$TE.random.w[subgroup_name]
    sens_lower <- 100 * props$sens$overall$lower.random.w[subgroup_name]
    sens_upper <- 100 * props$sens$overall$upper.random.w[subgroup_name]
    
    # Format to 2 decimal places
    sens_estimate_formatted <- sprintf("%.2f", sens_estimate)
    sens_lower_formatted <- sprintf("%.2f", sens_lower)
    sens_upper_formatted <- sprintf("%.2f", sens_upper)
    
    # Extract specificity estimates and confidence intervals
    spec_estimate <- 100 * props$spec$overall$TE.random.w[subgroup_name]
    spec_lower <- 100 * props$spec$overall$lower.random.w[subgroup_name]
    spec_upper <- 100 * props$spec$overall$upper.random.w[subgroup_name]
    
    # Format to 2 decimal places
    spec_estimate_formatted <- sprintf("%.2f", spec_estimate)
    spec_lower_formatted <- sprintf("%.2f", spec_lower)
    spec_upper_formatted <- sprintf("%.2f", spec_upper)
    
    # Print subgroup name
    cat("\nSubgroup:", subgroup_name, "\n")
    # Print sensitivity
    cat("  Sensitivity:", sens_estimate_formatted, "% (95% CI:", sens_lower_formatted, "-", sens_upper_formatted, "%)\n")
    # Print specificity
    cat("  Specificity:", spec_estimate_formatted, "% (95% CI:", spec_lower_formatted, "-", spec_upper_formatted, "%)\n")
}
  
  # The fitted model is returned in the result object. Avoid printing mada's
  # generic model-family heading here because the forest display uses the more
  # precise "Bivariate model" label.
  invisible(summary(reitsmas$reitsma.reg.fit))
          
  if (object.return){
    returned.object <- list()
    returned.object$reitsmas <- reitsmas
    returned.object$madads <- madads
    returned.object$metaprops <- props
    returned.object$summaries <- summaries
    returned.object$valid.subgroup.names <- valid.subgroup.list
    returned.object$subgroup.names <- subgroup.list
    return(returned.object)
  }
}




forest.diag.subgroup.combined <- function(dat, # Diagnostic-accuracy dataframe; columns are configurable below.
                                          subgrouping.variable, # dat$subgrouping.variable
                                          sortvar = NULL,
                                          sglabel = "Subgroup", # String for subgroup label
                                          lcols1 = NULL,        # Variable to be shown in left side dat$lcols1
                                          llab1 = NULL,         # Label for lcols1
                                          lcols2 = NULL,        # Variable to be shown in left side dat$lcols2
                                          llab2 = NULL,         # Label for lcols2
                                          object.return = TRUE, # If TRUE, returns an object used in the function
                                          plot.het.overall = TRUE,
                                          plot.het.subgroup = TRUE,
                                          plot.overall = TRUE,
                                          calcwidth.shet.opt = FALSE,
                                           forest.xlim = c(50, 100),
                                           leftspace = "    ",
                                           rightspace = "        ",
                                           ratio = NULL,
                                           layout.mode = c("auto", "manual"),
                                           auto.fit.device = TRUE,
                                           auto.minimum.scale = 0.6,
                                           auto.inner.trim.cm = 0,
                                           only.subgroups.bigger.than.3 = T,
                                           omnibus.alpha = 0.05,
                                           pairwise = c("ask", "always", "never"),
                                           p.adjust.method = "holm",
                                           study.names = NULL, TP = NULL,
                                           TN = NULL, FP = NULL, FN = NULL,
                                           uniquer.row.id = NULL,
                                           study.id = NULL,
                                           sensitivity.width.adjustment = 0,
                                           fontsize = NULL,
                                            ...) {
  pairwise <- match.arg(pairwise)
  layout.mode <- match.arg(layout.mode)
  fontsize <- .dta_resolve_forest_fontsize(fontsize)
# Prepare data
  dat[["subgrouping.variable"]] <- subgrouping.variable
  dat[["lcols1"]] <- lcols1
  dat[["lcols2"]] <- lcols2
  if (!is.null(sortvar)) dat[[".dta_sortvar"]] <- sortvar
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  if (is.null(sortvar)) {
    dat <- dat[order(dat[["subgrouping.variable"]]), ]
  } else {
    dat <- dat[order(dat[["subgrouping.variable"]], dat[[".dta_sortvar"]]), ]
  }
  
  # Get list of subgroups
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  if (only.subgroups.bigger.than.3) {
  subgroup.list <- names(counts[counts >= 3])
    }
  subgroup.list <- sort(subgroup.list)
  valid.subgroup.list <- make.names(subgroup.list, unique = TRUE)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  
  # Initialize lists to store results
  reitsmas <- list()
  reitsmas$subgroups <- list()
  summaries <- list()
  summaries$subgroups <- list()
  madads <- list()
  madads$subgroups <- list()
  props <- list()
  props$sens <- list()
  props$spec <- list()
  
  # Overall Reitsma model
  reitsmas$reitsma.overall <- reitsma(data = dat, method = "ml")
  reitsmas$subgroup.tests <- reitsma.subgroup.comparisons(
    dat = dat,
    subgrouping.variable = dat[["subgrouping.variable"]],
    omnibus.alpha = omnibus.alpha,
    pairwise = pairwise,
    p.adjust.method = p.adjust.method
  )
  reitsmas$reitsma.reg.fit <- reitsmas$subgroup.tests$full.fit
  reitsmas$reitsma.intercept <- reitsmas$subgroup.tests$intercept.fit
  reitsmas$anova.reitsma <- reitsmas$subgroup.tests$omnibus
  summaries$reitsma.overall <- summary(reitsmas$reitsma.overall)
  madads$overall <- madad(dat)
  
  # Meta-analysis for sensitivity and specificity (overall)
  props$sens$overall <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TP"]],
    n = dat[["TP"]] + dat[["FN"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Sensitivity",
    method.ci = "WS",
    common = FALSE,
    subgroup = dat[["subgrouping.variable"]],
    subgroup.name = sglabel
  )

  props$spec$overall <- metaprop(
    data = dat,
    studlab = dat[["names"]],
    event = dat[["TN"]],
    n = dat[["TN"]] + dat[["FP"]],
    method.tau = "ML",
    sm = "PRAW",
    outclab = "Specificity",
    method.ci = "WS",
    common = FALSE,
    subgroup = dat[["subgrouping.variable"]],
    subgroup.name = sglabel
  )
  
  # Update confidence intervals from mada
  props$sens$overall$upper <- madads$overall$sens$sens.ci[, 2]
  props$sens$overall$lower <- madads$overall$sens$sens.ci[, 1]
  props$sens$overall$TE.random <- summaries$reitsma.overall$coefficients[3, 1]
  props$sens$overall$lower.random <- summaries$reitsma.overall$coefficients[3, 5]
  props$sens$overall$upper.random <- summaries$reitsma.overall$coefficients[3, 6]

  props$spec$overall$upper <- madads$overall$spec$spec.ci[, 2]
  props$spec$overall$lower <- madads$overall$spec$spec.ci[, 1]
  props$spec$overall$TE.random <- 1 - summaries$reitsma.overall$coefficients[4, 1]
  props$spec$overall$lower.random <- 1 - summaries$reitsma.overall$coefficients[4, 6]
  props$spec$overall$upper.random <- 1 - summaries$reitsma.overall$coefficients[4, 5]
  study.population <- .dta_population(dat)
  props$sens$overall <- .dta_use_population_weights(
    props$sens$overall, study.population
  )
  props$spec$overall <- .dta_use_population_weights(
    props$spec$overall, study.population
  )

  
  hets.list <- list()
  
  # Loop over subgroups
  for (sg in seq_along(subgroup.list)) {
    subgroup_name <- subgroup.list[sg]
    valid_name <- valid.subgroup.list[sg]
    subgroup_data <- dat[dat[["subgrouping.variable"]] == subgroup_name, ]
    
    # Reitsma model for subgroup
    reitsma.sg <- reitsma(data = subgroup_data, method = 'ml')
    reitsmas$subgroups[[valid_name]] <- reitsma.sg
    summaries$subgroups[[valid_name]] <- summary(reitsma.sg)
    madads$subgroups[[valid_name]] <- madad(subgroup_data)
    
    # Sensitivity estimates for subgroup
    props$sens$overall$TE.random.w[[subgroup_name]] <- summaries$subgroups[[valid_name]]$coefficients[3, 1]
    props$sens$overall$lower.random.w[[subgroup_name]] <- summaries$subgroups[[valid_name]]$coefficients[3, 5]
    props$sens$overall$upper.random.w[[subgroup_name]] <- summaries$subgroups[[valid_name]]$coefficients[3, 6]
    
    # Specificity estimates for subgroup
    props$spec$overall$TE.random.w[[subgroup_name]] <- 1 - summaries$subgroups[[valid_name]]$coefficients[4, 1]
    props$spec$overall$lower.random.w[[subgroup_name]] <- 1 - summaries$subgroups[[valid_name]]$coefficients[4, 6]
    props$spec$overall$upper.random.w[[subgroup_name]] <- 1 - summaries$subgroups[[valid_name]]$coefficients[4, 5]
    
    if (plot.het.subgroup) {
      hets.list[[subgroup_name]] <- het.string(reitsma.sg)
    }
  }

  # Prepare left column labels
  if (is.null(lcols1)) {
    lcols1.string <- NULL
  } else {
    lcols1.string <- "lcols1"
  }
  if (is.null(lcols2)) {
    lcols2.string <- NULL
  } else {
    lcols2.string <- "lcols2"
  }
  
  leftcols <- c("studlab", lcols1.string, lcols2.string)
  leftlabs <- c("Study", llab1, llab2)
  
  # Remove NULL values from leftcols and leftlabs
  leftcols <- leftcols[!sapply(leftcols, is.null)]
  leftlabs <- leftlabs[!sapply(leftlabs, is.null)]

  sens.subgroup.text <- if (plot.het.subgroup) {
    unname(unlist(hets.list, use.names = FALSE))
  } else {
    rep("", length(subgroup.list))
  }
  
  # Generate the sensitivity forest plot and capture it as a grob
  sens.forest.args <- list(
    xlim = forest.xlim,
    pscale = 100,
    just.addcols.right = "center",
    rightcols = c("effect", "ci"),
    rightlabs = c("Sensitivity %", "95% C.I."),
    leftcols = leftcols,
    leftlabs = leftlabs,
    xlab = "Sensitivity",
    smlab = "",
    weight.study = "common",
    squaresize = 0.7,
    col.square = "navy",
    col.square.lines = "navy",
    col.diamond = "maroon",
    col.diamond.lines = "maroon",
    pooled.totals = FALSE,
    comb.fixed = FALSE,
    fs.hetstat = 10,
    print.tau2 = FALSE,
    print.Q = FALSE,
    print.pval.Q = FALSE,
    print.I2 = FALSE,
    digits = 1,
    subgroup.hetstat = FALSE,
    hetstat = FALSE,
    test.subgroup = FALSE,
    text.random = if (plot.overall) "Bivariate model" else "",
    text.random.w = sens.subgroup.text,
    text.addline2 = if (plot.overall && plot.het.overall) {
      het.string(reitsmas$reitsma.overall)
    } else {
      ""
    },
    text.addline1 = paste("Between-group difference (p): ", round(reitsmas$anova.reitsma$statistic[3], digits = 3)),
    ref = if (plot.overall) 100 * props$sens$overall$TE.random else NA,
    calcwidth.random = calcwidth.shet.opt,
    fontsize = fontsize,
    just = "center",
    overall = plot.overall
  )
  # Generate the specificity forest plot and capture it as a grob
  manual.layout <- layout.mode == "manual" || !is.null(ratio)
  if (manual.layout) {
    props$spec$overall$leftspace <- rep(leftspace, length(props$spec$overall$studlab))
    props$spec$overall$rightspace <- rep(rightspace, length(props$spec$overall$studlab))
  }
  # Retain the original unique subgroup levels. Their text is hidden below
  # with col.subgroup, while the rows remain available for panel alignment.
  empty_list <- rep("", length(subgroup.list))
  spec.forest.args <- list(
    xlim = forest.xlim,
    pscale = 100,
    just.addcols.right = "center",
    rightcols = if (manual.layout) c("effect", "ci", "rightspace") else c("effect", "ci"),
    rightlabs = if (manual.layout) c("Specificity %", "95% C.I.", "") else c("Specificity %", "95% C.I."),
    leftcols = if (manual.layout) "leftspace" else FALSE,
    leftlabs = NULL,
    studlab = FALSE,
    xlab = "Specificity",
    smlab = "",
    weight.study = "common",
    squaresize = 0.7,
    col.square = "navy",
    col.square.lines = "navy",
    col.diamond = "maroon",
    col.diamond.lines = "maroon",
    pooled.totals = FALSE,
    comb.fixed = FALSE,
    fs.hetstat = 10,
    print.tau2 = FALSE,
    print.Q = FALSE,
    print.pval.Q = FALSE,
    print.I2 = FALSE,
    digits = 1,
    subgroup.hetstat = FALSE,
    hetstat = FALSE,
    test.subgroup = FALSE,
    text.random.w = empty_list,
    text.random = "",
    text.addline2 = "",
    text.addline1 = "",
    ref = if (plot.overall) 100 * props$spec$overall$TE.random else NA,
    calcwidth.random = calcwidth.shet.opt,
    fontsize = fontsize,
    just = "center",
    overall = plot.overall,
    col.subgroup = "white"
  )
  device.fit <- .dta_capture_and_fit_forests(
    meta.sens = props$sens$overall,
    sens.args = sens.forest.args,
    meta.spec = props$spec$overall,
    spec.args = spec.forest.args,
    enabled = auto.fit.device && layout.mode == "auto" && is.null(ratio),
    minimum.scale = auto.minimum.scale
  )
  sens.forest.args <- device.fit$sens.args
  spec.forest.args <- device.fit$spec.args
  grob.sens <- device.fit$grob.sens
  grob.spec <- device.fit$grob.spec

  forest.layout <- .dta_resolve_forest_layout(
    meta.sens = props$sens$overall,
    sens.args = sens.forest.args,
    grob.sens = grob.sens,
    meta.spec = props$spec$overall,
    spec.args = spec.forest.args,
    grob.spec = grob.spec,
    layout.mode = layout.mode,
    ratio = ratio
  )
  forest.layout$device.scale <- device.fit$scale
  forest.layout$device.width.in <- device.fit$device.width.in
  forest.layout$available.width.in <- device.fit$available.width.in
  forest.layout$base.ratio <- forest.layout$ratio
  forest.layout$ratio <- .dta_adjust_forest_ratio(
    forest.layout$ratio, sensitivity.width.adjustment
  )
  forest.layout$ratio.adjustment <- sensitivity.width.adjustment
  forest.layout$font.size <- fontsize

  if (forest.layout$mode == "auto") {
    layout.note <- forest.layout$measurement
    if (is.finite(forest.layout$device.scale) &&
        forest.layout$device.scale < 0.999) {
      layout.note <- paste0(
        layout.note,
        "; device scale ",
        sprintf("%.1f%%", 100 * forest.layout$device.scale)
      )
    }
    cat(
      "Automatic forest layout:",
      sprintf("sensitivity %.1f%%, specificity %.1f%%",
              100 * forest.layout$ratio[1], 100 * forest.layout$ratio[2]),
      sprintf("(%s)\n", layout.note)
    )
  }
  
  # Combine the two grobs in a fixed-width centred viewport. Extra device
  # width is therefore placed at the outer margins, not between the panels.
  combined <- .dta_draw_combined_forests(
    grob.sens, grob.spec,
    forest.layout = forest.layout,
    inner.trim.cm = auto.inner.trim.cm
  )
  
  # Print overall estimates
  sens_estimate <- 100 * props$sens$overall$TE.random
  sens_lower <- 100 * props$sens$overall$lower.random
  sens_upper <- 100 * props$sens$overall$upper.random
  
  spec_estimate <- 100 * props$spec$overall$TE.random
  spec_lower <- 100 * props$spec$overall$lower.random
  spec_upper <- 100 * props$spec$overall$upper.random
  
  cat("Overall Sensitivity:", sprintf("%.2f", sens_estimate), "% (95% CI:", sprintf("%.2f", sens_lower), "-", sprintf("%.2f", sens_upper), "%)\n")
  cat("Overall Specificity:", sprintf("%.2f", spec_estimate), "% (95% CI:", sprintf("%.2f", spec_lower), "-", sprintf("%.2f", spec_upper), "%)\n")
  
  # Loop over subgroups to print estimates
  for (sg in seq_along(subgroup.list)) {
    subgroup_name <- subgroup.list[sg]
    valid_name <- valid.subgroup.list[sg]
    
    # Extract sensitivity estimates and confidence intervals
    sens_estimate <- 100 * props$sens$overall$TE.random.w[[subgroup_name]]
    sens_lower <- 100 * props$sens$overall$lower.random.w[[subgroup_name]]
    sens_upper <- 100 * props$sens$overall$upper.random.w[[subgroup_name]]
    
    # Format to 2 decimal places
    sens_estimate_formatted <- sprintf("%.2f", sens_estimate)
    sens_lower_formatted <- sprintf("%.2f", sens_lower)
    sens_upper_formatted <- sprintf("%.2f", sens_upper)
    
    # Extract specificity estimates and confidence intervals
    spec_estimate <- 100 * props$spec$overall$TE.random.w[[subgroup_name]]
    spec_lower <- 100 * props$spec$overall$lower.random.w[[subgroup_name]]
    spec_upper <- 100 * props$spec$overall$upper.random.w[[subgroup_name]]
    
    # Format to 2 decimal places
    spec_estimate_formatted <- sprintf("%.2f", spec_estimate)
    spec_lower_formatted <- sprintf("%.2f", spec_lower)
    spec_upper_formatted <- sprintf("%.2f", spec_upper)
    
    # Print subgroup name
    cat("\nSubgroup:", subgroup_name, "\n")
    # Print sensitivity
    cat("  Sensitivity:", sens_estimate_formatted, "% (95% CI:", sens_lower_formatted, "-", sens_upper_formatted, "%)\n")
    # Print specificity
    cat("  Specificity:", spec_estimate_formatted, "% (95% CI:", spec_lower_formatted, "-", spec_upper_formatted, "%)\n")
  }
  
  # Print summary of the Reitsma regression
  invisible(summary(reitsmas$reitsma.overall))
  
  # Return objects if requested
  if (object.return) {
    returned.object <- list()
    returned.object$reitsmas <- reitsmas
    returned.object$madads <- madads
    returned.object$metaprops <- props
    returned.object$summaries <- summaries
    returned.object$valid.subgroup.names <- valid.subgroup.list
    returned.object$subgroup.names <- subgroup.list
    returned.object$layout <- forest.layout
    returned.object$plot <- combined
    return(returned.object)
  }
}






.dta_bootstrap_cores <- function(parallel = TRUE, n.cores = NULL) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("'parallel' must be TRUE or FALSE.")
  }

  detected <- parallel::detectCores(logical = TRUE)
  if (is.na(detected) || detected < 1L) {
    detected <- 1L
  }

  if (is.null(n.cores)) {
    n.cores <- max(1L, detected - 2L)
  }
  if (!is.numeric(n.cores) || length(n.cores) != 1L || is.na(n.cores) ||
      !is.finite(n.cores) || n.cores < 1L || n.cores != as.integer(n.cores)) {
    stop("'n.cores' must be a single positive integer or NULL.")
  }

  n.cores <- min(as.integer(n.cores), detected)
  if (!parallel) 1L else n.cores
}

.dta_progress_reporter <- function(total, label = "Bootstrap", enabled = TRUE) {
  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
    stop("'progress' must be TRUE or FALSE.")
  }
  if (!enabled) {
    return(list(update = function(value) invisible(value),
                close = function() invisible(NULL)))
  }

  cat(label, " progress:\n", sep = "")
  bar <- utils::txtProgressBar(min = 0, max = total, initial = 0, style = 3)
  list(
    update = function(value) utils::setTxtProgressBar(bar, value),
    close = function() close(bar)
  )
}

.dta_start_bootstrap_backend <- function(n.cores, live.progress = FALSE) {
  if (n.cores <= 1L) {
    return(list(
      created = FALSE, managed = FALSE, deactivate = FALSE,
      cluster = NULL, n.cores = 1L
    ))
  }

  # Reuse an already registered foreach backend. This mirrors the original
  # implementation and, importantly, avoids starting a new group of R workers
  # for every subgroup bootstrap.
  if (foreach::getDoParRegistered() && foreach::getDoParWorkers() > 1L) {
    backend.name <- foreach::getDoParName()
    return(list(
      created = FALSE,
      managed = FALSE,
      deactivate = FALSE,
      cluster = NULL,
      n.cores = as.integer(foreach::getDoParWorkers()),
      live.progress = grepl("doSNOW", backend.name, fixed = TRUE)
    ))
  }

  cluster <- .dta_bootstrap_state$cluster
  reuse <- !is.null(cluster) &&
    identical(.dta_bootstrap_state$n.cores, n.cores)
  if (reuse) {
    reuse <- isTRUE(tryCatch(
      all(unlist(parallel::clusterCall(cluster, function() TRUE))),
      error = function(e) FALSE
    ))
  }
  if (!reuse && !is.null(cluster)) {
    try(parallel::stopCluster(cluster), silent = TRUE)
    .dta_bootstrap_state$cluster <- NULL
    .dta_bootstrap_state$n.cores <- NULL
  }
  if (!reuse) {
    cluster <- parallel::makePSOCKcluster(n.cores)
    library.paths <- .libPaths()
    invisible(parallel::clusterCall(
      cluster, function(paths) .libPaths(paths), paths = library.paths
    ))
    .dta_bootstrap_state$cluster <- cluster
    .dta_bootstrap_state$n.cores <- n.cores
  }
  if (isTRUE(live.progress)) {
    doSNOW::registerDoSNOW(cluster)
  } else {
    doParallel::registerDoParallel(cluster)
  }
  list(
    created = !reuse,
    managed = TRUE,
    deactivate = TRUE,
    cluster = cluster,
    n.cores = n.cores,
    live.progress = isTRUE(live.progress)
  )
}

.dta_deactivate_bootstrap_backend <- function(backend) {
  if (isTRUE(backend$deactivate)) {
    foreach::registerDoSEQ()
  }
  invisible(NULL)
}

.dta_stop_bootstrap_backend <- function(backend) {
  if (isTRUE(backend$managed) && !is.null(backend$cluster)) {
    try(parallel::stopCluster(backend$cluster), silent = TRUE)
    .dta_bootstrap_state$cluster <- NULL
    .dta_bootstrap_state$n.cores <- NULL
    foreach::registerDoSEQ()
  }
  invisible(NULL)
}

.dta_reset_bootstrap_backend <- function() {
  base::suspendInterrupts({
    cluster <- .dta_bootstrap_state$cluster
    .dta_bootstrap_state$cluster <- NULL
    .dta_bootstrap_state$n.cores <- NULL
    foreach::registerDoSEQ()
    if (!is.null(cluster)) {
      try(parallel::stopCluster(cluster), silent = TRUE)
    }
  })
  invisible(NULL)
}

.dta_abort_bootstrap <- function(condition) {
  had.managed.workers <- !is.null(.dta_bootstrap_state$cluster)
  .dta_reset_bootstrap_backend()
  if (had.managed.workers) {
    message(
      "Bootstrap interrupted; DTAtoolkit parallel workers and connections ",
      "were closed."
    )
  } else {
    message(
      "Bootstrap interrupted; the active foreach backend was detached. ",
      "No DTAtoolkit-managed worker pool was active; user-managed workers ",
      "were left untouched."
    )
  }
  stop(condition)
}

.dta_compact_bootstrap_function <- function(FUN) {
  globals <- codetools::findGlobals(FUN, merge = FALSE)
  globals <- unique(c(globals$variables, globals$functions))
  source.environment <- environment(FUN)
  local.globals <- globals[vapply(
    globals,
    exists,
    logical(1L),
    envir = source.environment,
    inherits = TRUE
  )]
  values <- lapply(
    local.globals,
    get,
    envir = source.environment,
    inherits = TRUE
  )
  compact.environment <- list2env(
    stats::setNames(values, local.globals),
    parent = baseenv()
  )
  compact.function <- FUN
  environment(compact.function) <- compact.environment
  compact.function
}

.dta_bootstrap_map <- function(B, FUN, parallel = TRUE, n.cores = NULL,
                               progress = TRUE, label = "Bootstrap",
                               seed = NULL) {
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || !is.finite(B) ||
      B < 1L || B > .Machine$integer.max || B != as.integer(B)) {
    stop("'B' must be a single positive integer.")
  }
  B <- as.integer(B)
  n.cores <- .dta_bootstrap_cores(parallel, n.cores)
  requested.cores <- n.cores

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed)) {
      stop("'seed' must be a single finite number or NULL.")
    }
    set.seed(seed)
  }
  replicate.seeds <- sample.int(.Machine$integer.max, B, replace = TRUE)

  reporter <- .dta_progress_reporter(B, label, progress)
  on.exit(reporter$close(), add = TRUE)

  if (n.cores <= 1L) {
    results <- vector("list", B)
    for (i in seq_len(B)) {
      set.seed(replicate.seeds[[i]])
      results[[i]] <- tryCatch(FUN(i), error = function(e) NA_real_)
      reporter$update(i)
    }
  } else {
    backend <- .dta_start_bootstrap_backend(
      n.cores, live.progress = isTRUE(progress)
    )
    if (isTRUE(backend$deactivate)) {
      on.exit(.dta_deactivate_bootstrap_backend(backend), add = TRUE)
    }
    n.cores <- backend$n.cores

    bootstrap.fun <- .dta_compact_bootstrap_function(FUN)
    run.foreach <- function(active.backend) {
      seeds <- replicate.seeds
      bootstrap.worker <- bootstrap.fun
      live.progress <- isTRUE(progress) &&
        isTRUE(active.backend$live.progress)
      progress.callback <- function(completed) reporter$update(completed)
      snow.options <- if (live.progress) {
        list(progress = progress.callback, preschedule = FALSE)
      } else {
        list()
      }

      # This is the original foreach execution model: submit the complete set
      # of bootstrap replicates once. With progress enabled, doSNOW calls the
      # callback as each result arrives, without scheduling separate batches.
      values <- foreach::foreach(
        i = seq_len(B),
        .combine = c,
        .init = list(),
        .inorder = TRUE,
        .multicombine = TRUE,
        .maxcombine = 100L,
        .errorhandling = "pass",
        .options.snow = snow.options
      ) %dopar% {
        list(tryCatch({
          set.seed(seeds[[i]])
          bootstrap.worker(i)
        }, error = function(e) NA_real_))
      }
      if (isTRUE(progress) && !live.progress) reporter$update(B)
      values
    }

    results <- tryCatch(
      run.foreach(backend),
      interrupt = function(e) .dta_abort_bootstrap(e),
      error = identity
    )
    interrupted.error <- inherits(results, "error") && grepl(
      "user interrupt|interrupted|execution halted by user",
      conditionMessage(results),
      ignore.case = TRUE
    )
    if (interrupted.error) {
      .dta_abort_bootstrap(results)
    }
    recoverable <- inherits(results, "error") && grepl(
      paste0(
        "worker initialization failed|error reading from connection|",
        "invalid connection|unserialize.*connection"
      ),
      conditionMessage(results),
      ignore.case = TRUE
    )
    if (recoverable) {
      message(
        "Parallel workers failed to initialize; rebuilding the worker pool ",
        "and retrying once."
      )
      .dta_reset_bootstrap_backend()
      backend <- .dta_start_bootstrap_backend(
        requested.cores, live.progress = isTRUE(progress)
      )
      if (isTRUE(backend$deactivate)) {
        on.exit(.dta_deactivate_bootstrap_backend(backend), add = TRUE)
      }
      n.cores <- backend$n.cores
      results <- tryCatch(
        run.foreach(backend),
        interrupt = function(e) .dta_abort_bootstrap(e)
      )
    } else if (inherits(results, "error")) {
      stop(results)
    }
  }

  list(
    values = results,
    parallel = n.cores > 1L,
    n.cores = n.cores,
    B = B
  )
}

AUC_bootstrap <- function(TP, FP, FN, TN, B = 2000, alpha = 0.95,
                          parallel = TRUE, n.cores = NULL, progress = TRUE,
                          seed = NULL) {
  counts <- data.frame(TP = TP, FP = FP, FN = FN, TN = TN)
  numeric.columns <- vapply(counts, is.numeric, logical(1L))
  if (nrow(counts) < 1L || !all(numeric.columns) || anyNA(counts) ||
      any(!is.finite(as.matrix(counts))) || any(counts < 0) ||
      any(as.matrix(counts) != floor(as.matrix(counts)))) {
    stop(paste(
      "TP, FP, FN, and TN must be equally sized, non-missing,",
      "non-negative integer-count vectors."
    ))
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be between 0 and 1.")
  }

  N <- nrow(counts)
  n1 <- counts$TP + counts$FN
  n2 <- counts$TN + counts$FP
  fit0 <- mada::reitsma(counts)
  auc <- summary(fit0)$AUC$AUC
  mu1 <- as.numeric(fit0$coefficients)
  G1 <- fit0$Psi

  one.bootstrap <- function(b) {
    tryCatch({
      latent <- MASS::mvrnorm(n = N, mu = mu1, Sigma = G1)
      if (N == 1L) {
        latent <- matrix(latent, nrow = 1L)
      }
      probabilities <- stats::plogis(latent)
      TPb <- stats::rbinom(N, size = n1, prob = probabilities[, 1L])
      FPb <- stats::rbinom(N, size = n2, prob = probabilities[, 2L])
      boot.data <- data.frame(
        TP = TPb,
        FP = FPb,
        FN = n1 - TPb,
        TN = n2 - FPb
      )
      fit.pb <- suppressWarnings(mada::reitsma(boot.data))
      as.numeric(summary(fit.pb)$AUC$AUC)
    }, error = function(e) NA_real_)
  }

  bootstrap <- .dta_bootstrap_map(
    B = B,
    FUN = one.bootstrap,
    parallel = parallel,
    n.cores = n.cores,
    progress = progress,
    label = "AUC bootstrap",
    seed = seed
  )
  auc.pb <- as.numeric(unlist(bootstrap$values, use.names = FALSE))
  successful <- sum(is.finite(auc.pb))
  if (successful == 0L) {
    stop("All AUC bootstrap model fits failed.")
  }
  if (successful < B) {
    warning(B - successful, " of ", B,
            " AUC bootstrap fits failed and were omitted from the confidence interval.")
  }

  probs <- c(0.5 * (1 - alpha), 1 - 0.5 * (1 - alpha))
  CI <- stats::quantile(auc.pb, probs, na.rm = TRUE, names = TRUE)
  list(
    AUC = auc,
    CI = CI,
    bootstrap.AUC = auc.pb,
    n.boots = B,
    n.successful = successful,
    parallel = bootstrap$parallel,
    n.cores = bootstrap$n.cores
  )
}

# Backward-compatible alias for the previous function name.
AUC_boot_paralell <- function(TP, FP, FN, TN, B = 2000, alpha = 0.95,
                              parallel = TRUE, n.cores = NULL,
                              progress = TRUE, seed = NULL) {
  AUC_bootstrap(
    TP = TP, FP = FP, FN = FN, TN = TN, B = B, alpha = alpha,
    parallel = parallel, n.cores = n.cores, progress = progress, seed = seed
  )
}

.dta_auc_ci_cache_descriptor <- function(dat, n.boots, seed) {
  if (!is.numeric(n.boots) || length(n.boots) != 1L || is.na(n.boots) ||
      !is.finite(n.boots) || n.boots < 1L ||
      n.boots > .Machine$integer.max || n.boots != as.integer(n.boots)) {
    stop("'n.boots' must be a single positive integer.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
       !is.finite(seed))) {
    stop("'AUC.CI.seed' must be a single finite number or NULL.",
         call. = FALSE)
  }
  data.signature <- list(
    counts = cbind(
      TP = as.numeric(dat[["TP"]]),
      TN = as.numeric(dat[["TN"]]),
      FP = as.numeric(dat[["FP"]]),
      FN = as.numeric(dat[["FN"]])
    ),
    subgroup = enc2utf8(as.character(dat[["subgrouping.variable"]]))
  )
  settings <- list(
    n.boots = as.integer(n.boots),
    seed = if (is.null(seed)) NULL else as.numeric(seed),
    confidence.level = 0.95,
    bootstrap.model = "reitsma-parametric-v1"
  )
  data.md5 <- digest::digest(data.signature, algo = "md5", serialize = TRUE)
  key <- digest::digest(
    list(data.md5 = data.md5, settings = settings),
    algo = "md5",
    serialize = TRUE
  )
  list(
    key = key,
    data.md5 = data.md5,
    data.signature = data.signature,
    settings = settings,
    subgroups = sort(unique(data.signature$subgroup)),
    n.rows = nrow(data.signature$counts)
  )
}

.dta_auc_ci_cache_lookup <- function(descriptor) {
  entry <- .dta_auc_ci_cache_state$entry
  had.entry <- !is.null(entry)
  hit <- had.entry &&
    identical(entry$key, descriptor$key) &&
    identical(entry$data.signature, descriptor$data.signature) &&
    identical(entry$settings, descriptor$settings)
  if (!hit) .dta_auc_ci_cache_state$entry <- NULL
  list(
    hit = hit,
    had.entry = had.entry,
    value = if (hit) entry$value else NULL,
    metadata = if (hit) entry$metadata else NULL
  )
}

.dta_auc_ci_cache_store <- function(descriptor, value) {
  metadata <- list(
    cache.key = descriptor$key,
    data.md5 = descriptor$data.md5,
    n.boots = descriptor$settings$n.boots,
    seed = descriptor$settings$seed,
    subgroups = descriptor$subgroups,
    n.rows = descriptor$n.rows,
    created.at = Sys.time()
  )
  .dta_auc_ci_cache_state$entry <- list(
    key = descriptor$key,
    data.signature = descriptor$data.signature,
    settings = descriptor$settings,
    value = value,
    metadata = metadata
  )
  metadata
}

.dta_auc_ci_cache_clear <- function() {
  had.entry <- !is.null(.dta_auc_ci_cache_state$entry)
  .dta_auc_ci_cache_state$entry <- NULL
  invisible(had.entry)
}





multiple.srocs <- function(dat, # Diagnostic-accuracy dataframe; columns are configurable below.
                           subgrouping.variable = NULL, #dat$subgrouping variable
                           sroc.colors = c("blue", "maroon", "black", "skyblue", "#20cb20", "red"), #colors for the sroc and summary estimates and ellipse 
                           points.colors = c("#0000FF20", "#A52A2A20","#1d1c1c20" ,"#87ceeb30", "#20cb2020", "#ff000350"),#colors for the point estimates
                           pch.list = c(16, 15, 17, 18, 15, 13), #shape for point estimate default is : c(1, 0, 2, 5, 7, 13)
                           summary.pch.list = c(16, 15, 17, 18, 15, 13), #shape for summary estimate
                           plot.ellipse = T,
                           plot.points = T,
                           plot.legend =T,
                           main.title = "SROC curves for subgroups",
                           object.return = T, #if true, it returns an object with all used objects
                           legend.AUC = T,
                           AUC.CI = F, # If true it also returns confidence interval for AUC.
                           n.boots = 2000, # number of bootstraps for AUC CI
                           AUC.CI.object = NULL, # if you have already done auc ci calculation and have the results  as an object, place it here for faster implementation
                           AUC.CI.parallel = TRUE,
                           AUC.CI.n.cores = NULL,
                           AUC.CI.progress = TRUE,
                           AUC.CI.seed = NULL,
                           magnify = 2, # argument to magnify weight size of point estimates in the sroc
                           extrapolate = TRUE,
                           omnibus.alpha = 0.05,
                           pairwise = c("ask", "always", "never"),
                           p.adjust.method = "holm",
                           study.names = NULL, TP = NULL, TN = NULL,
                           FP = NULL, FN = NULL, uniquer.row.id = NULL,
                           study.id = NULL,
                           color.palette = "mada6",
                           point.alpha = 0.18,
                           point.size.range = NULL,
                           AUC.CI.cache = TRUE
){
  legacy.sroc.colors <- if (missing(sroc.colors)) NULL else sroc.colors
  legacy.points.colors <- if (missing(points.colors)) NULL else points.colors
  pairwise <- match.arg(pairwise)
  ellipse.lty <- 0
  if (plot.ellipse){
    ellipse.lty <- 2
  }
  if (is.null(subgrouping.variable)){
    dat[["subgrouping.variable"]] <- "Pooled"
  }else{
    dat[["subgrouping.variable"]] <- subgrouping.variable  
  }
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  subgroup.list <- unique(dat[["subgrouping.variable"]])
  counts <- table(dat[["subgrouping.variable"]])
  subgroup.list <- names(counts[counts >= 3])
  subgroup.list <- sort(subgroup.list)
  if (length(subgroup.list) == 0L) {
    stop("At least one subgroup with three or more complete studies is required.",
         call. = FALSE)
  }
  valid.subgroup.list <- make.names(subgroup.list, unique = T)
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  plot.colors <- .dta_resolve_plot_colours(
    n = length(subgroup.list),
    color.palette = color.palette,
    legacy.primary = legacy.sroc.colors,
    legacy.points = legacy.points.colors,
    point.alpha = point.alpha
  )
  sroc.colors <- plot.colors$primary
  points.colors <- plot.colors$points
  pch.list <- rep(pch.list, length.out = length(subgroup.list))
  summary.pch.list <- rep(
    summary.pch.list, length.out = length(subgroup.list)
  )
  if (!is.numeric(magnify) || length(magnify) != 1L ||
      !is.finite(magnify) || magnify <= 0) {
    stop("'magnify' must be one positive number.", call. = FALSE)
  }
  if (is.null(point.size.range)) {
    point.size.range <- c(max(0.35, 0.35 * magnify), magnify)
  }
  dat[["population"]] <- .dta_population(dat)
  dat[["point.size"]] <- .dta_population_point_sizes(
    dat[["population"]], point.size.range
  )
  reitsmas <- list()
  reitsmas$reitsma.overall <- reitsma(data = dat, method = "ml")
  reitsmas$subgroups <- list()
  summaries <- list()
  if (length(subgroup.list)>1){
    summaries$reitsma.overall <- summary(reitsmas$reitsma.overall)
    reitsmas$subgroup.tests <- reitsma.subgroup.comparisons(
      dat = dat,
      subgrouping.variable = dat[["subgrouping.variable"]],
      omnibus.alpha = omnibus.alpha,
      pairwise = pairwise,
      p.adjust.method = p.adjust.method
    )
    reitsmas$reitsma.reg.fit <- reitsmas$subgroup.tests$full.fit
    reitsmas$reitsma.intercept <- reitsmas$subgroup.tests$intercept.fit
    reitsmas$anova.reitsma <- reitsmas$subgroup.tests$omnibus
  }
  summaries$subgroups <- list()
  madauni.main <- madauni(dat, type = "DOR", method = "DSL")
  datasets <- list()
  datasets$main <- dat
  datasets$subgroups <- list()
  for (sg in 1:length(subgroup.list)){
    reitsma.sg <- reitsma(data = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ], method = 'ml')
    name <- valid.subgroup.list[sg]
    reitsmas$subgroups[[name]] <- reitsma.sg 
    summaries$subgroups[[name]] <- summary(reitsma.sg)
    datasets$subgroups[[name]] <- dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]
  }
  if (length(subgroup.list)>1) {
    plot(reitsmas$subgroups[[valid.subgroup.list[1]]],
         extrapolate = FALSE,
         sroclwd = 2,
         predict = FALSE,
         pch = summary.pch.list[1],
         plotsumm = FALSE,
         main = main.title,
         sub = if (reitsmas$anova.reitsma$statistic[3]<0.01) "Difference between subgroups (bivariate model): p < 0.01" else paste("Difference between subgroups (bivariate model): p =" , round(reitsmas$anova.reitsma$statistic[3], digits = 3))
    )
  }else{
    plot(reitsmas$subgroups[[valid.subgroup.list[1]]],
         extrapolate = FALSE,
         sroclwd = 2,
         predict = FALSE,
         pch = summary.pch.list[1],
         plotsumm = FALSE,
         main = main.title,
    )
  }
  
  
  ROCellipse(reitsmas$subgroups[[valid.subgroup.list[1]]],
             lty = ellipse.lty,
             pch = summary.pch.list[1],
             add = TRUE,
             col = sroc.colors[1]
  )
  # Modify the lines(sroc(...)) call for the first subgroup
  if (extrapolate == FALSE) {
    fpr_values <- fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[1]), ])
    min_fpr <- max(min(fpr_values), 0.01) # Ensure min_fpr is at least 0.01
    max_fpr <- min(max(fpr_values), 0.99) # Ensure max_fpr is at most 0.99
    fpr_seq <- seq(min_fpr, max_fpr, length.out = 99)
    lines(sroc(reitsmas$subgroups[[valid.subgroup.list[1]]],
               fpr = fpr_seq),
          lty = 1,
          col = sroc.colors[1],
          lwd = 2
    )
  } else {
    lines(sroc(reitsmas$subgroups[[valid.subgroup.list[1]]]),
          lty = 1,
          col = sroc.colors[1],
          lwd = 2
    )
  }
  if (plot.points){
    points(fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[1]), ]),
           sens(dat[which(dat[["subgrouping.variable"]] == subgroup.list[1]), ]),
           pch = pch.list[1],
           cex = dat[which(dat[["subgrouping.variable"]] == subgroup.list[1]), ][["point.size"]],
           col = points.colors[1]
    )
  }
  
  if (length(subgroup.list)>1){
    for (sg in 2:length(subgroup.list)){
      # Modify the lines(sroc(...)) call for other subgroups
      if (extrapolate == FALSE) {
        fpr_values <- fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ])
        min_fpr <- max(min(fpr_values), 0.01) # Ensure min_fpr is at least 0.01
        max_fpr <- min(max(fpr_values), 0.99) # Ensure max_fpr is at most 0.99
        fpr_seq <- seq(min_fpr, max_fpr, length.out = 99)
        lines(sroc(reitsmas$subgroups[[valid.subgroup.list[sg]]],
                   fpr = fpr_seq),
              lty = 1,
              col = sroc.colors[sg],
              lwd = 2
        )
      } else {
        lines(sroc(reitsmas$subgroups[[valid.subgroup.list[sg]]]),
              lty = 1,
              col = sroc.colors[sg],
              lwd = 2
        )
      }
      ROCellipse(reitsmas$subgroups[[valid.subgroup.list[sg]]],
                 lty = ellipse.lty,
                 pch = summary.pch.list[sg],
                 add = TRUE,
                 col = sroc.colors[sg],
                 cex = 1.5
      )
      if(plot.points){
        points(fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
               sens(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
               pch = pch.list[sg],
               cex = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ][["point.size"]],
               col = points.colors[sg]
        )
      }
      
    }
  }
  
  no.subgroups <- length(subgroup.list)
  legend.pch <- summary.pch.list[1:no.subgroups]
  legend.col <- sroc.colors[1:no.subgroups]
  subgroups.auc.list <- list()
  AUC_CIs <- list()
  if (!is.logical(AUC.CI.cache) || length(AUC.CI.cache) != 1L ||
      is.na(AUC.CI.cache)) {
    stop("'AUC.CI.cache' must be TRUE or FALSE.", call. = FALSE)
  }
  auc.cache.descriptor <- NULL
  auc.cache.metadata <- list(
    enabled = isTRUE(AUC.CI.cache),
    hit = FALSE,
    source = if (isTRUE(AUC.CI)) "computed" else "not requested"
  )
  bootstrap.AUC.CI <- isTRUE(AUC.CI) && is.null(AUC.CI.object)
  if (isTRUE(AUC.CI) && !is.null(AUC.CI.object)) {
    AUC_CIs <- AUC.CI.object
    bootstrap.AUC.CI <- FALSE
    auc.cache.metadata$source <- "supplied AUC.CI.object"
  } else if (bootstrap.AUC.CI && isTRUE(AUC.CI.cache)) {
    auc.cache.descriptor <- .dta_auc_ci_cache_descriptor(
      dat = dat, n.boots = n.boots, seed = AUC.CI.seed
    )
    cache.lookup <- .dta_auc_ci_cache_lookup(auc.cache.descriptor)
    if (cache.lookup$hit) {
      AUC_CIs <- cache.lookup$value
      bootstrap.AUC.CI <- FALSE
      auc.cache.metadata <- utils::modifyList(
        cache.lookup$metadata,
        list(enabled = TRUE, hit = TRUE, source = "session memory")
      )
      message(
        "AUC CI cache hit (data MD5 ",
        auc.cache.descriptor$data.md5,
        "): reusing the previous ", n.boots, "-replicate result."
      )
    } else {
      auc.cache.metadata <- list(
        enabled = TRUE,
        hit = FALSE,
        source = "computed",
        cache.key = auc.cache.descriptor$key,
        data.md5 = auc.cache.descriptor$data.md5,
        n.boots = as.integer(n.boots),
        seed = if (is.null(AUC.CI.seed)) NULL else as.numeric(AUC.CI.seed),
        previous.entry.cleared = isTRUE(cache.lookup$had.entry)
      )
      if (isTRUE(cache.lookup$had.entry)) {
        message(
          "AUC CI cache mismatch: the previous entry was cleared; ",
          "running a fresh bootstrap."
        )
      }
    }
  } else if (bootstrap.AUC.CI) {
    auc.cache.metadata$source <- "computed without cache"
  }
  for (sg in 1:length(subgroup.list)){
    summary.sg <- summaries$subgroups[[valid.subgroup.list[sg]]]
    if(!AUC.CI){
      subgroups.auc.list[sg] <-  if (legend.AUC) paste(subgroup.list[sg],
                                                       " (AUC: ",
                                                       round(summary.sg$AUC$AUC, digits = 2),
                                                       ")",
                                                       sep = ""
      ) else subgroup.list[sg]
      
    } else {
      if (bootstrap.AUC.CI) {
        sg.dat <- dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]
        subgroup.seed <- if (is.null(AUC.CI.seed)) NULL else AUC.CI.seed + sg - 1L
        AUC_CIs[[valid.subgroup.list[sg]]] <- AUC_bootstrap(
          TP = sg.dat[["TP"]],
          FP = sg.dat[["FP"]],
          FN = sg.dat[["FN"]],
          TN = sg.dat[["TN"]],
          B = n.boots,
          parallel = AUC.CI.parallel,
          n.cores = AUC.CI.n.cores,
          progress = AUC.CI.progress,
          seed = subgroup.seed
        )
      }
      AUC_CIs_sg <- AUC_CIs[[valid.subgroup.list[sg]]]
      subgroups.auc.list[sg] <- paste(subgroup.list[sg],
                                      " - AUC: ",
                                      round(AUC_CIs_sg$AUC, digits = 2),
                                      " [",
                                      round(AUC_CIs_sg$CI[1], digits = 2),
                                      " - ",
                                      round(AUC_CIs_sg$CI[2], digits = 2),
                                      "]",
                                      sep = ""
      )
    }
    
  }
  if (isTRUE(AUC.CI) && bootstrap.AUC.CI && isTRUE(AUC.CI.cache)) {
    previous.entry.cleared <- isTRUE(
      auc.cache.metadata$previous.entry.cleared
    )
    stored.metadata <- .dta_auc_ci_cache_store(
      auc.cache.descriptor, AUC_CIs
    )
    auc.cache.metadata <- utils::modifyList(
      stored.metadata,
      list(
        enabled = TRUE,
        hit = FALSE,
        source = "computed and cached",
        previous.entry.cleared = previous.entry.cleared
      )
    )
  }
  if (plot.legend){
    legend("bottomright",
           legend = subgroups.auc.list,
           pch = legend.pch,
           lty = 1,
           col = legend.col
           , lwd =1.5
    )
  }
  
  # Check if there are multiple subgroups
  if (length(subgroup.list) > 1) {
    # Loop over each subgroup
    for (sg in 1:length(subgroup.list)) {
      subgroup_name <- subgroup.list[sg]
      valid_name <- valid.subgroup.list[sg]
      summary.sg <- summaries$subgroups[[valid_name]]
      
      # Extract AUC estimate
      auc_estimate <- summary.sg$AUC$AUC
      
      if (AUC.CI) {
        # Confidence intervals are stored in AUC_CIs
        auc_ci <- AUC_CIs[[valid_name]]$CI
        auc_lower <- auc_ci[1]
        auc_upper <- auc_ci[2]
        
        # Format to 2 decimal places
        auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
        auc_lower_formatted <- sprintf("%.2f", auc_lower)
        auc_upper_formatted <- sprintf("%.2f", auc_upper)
        
        # Print AUC with confidence intervals
        cat("\nSubgroup:", subgroup_name, "\n")
        cat("  AUC:", auc_estimate_formatted, "(95% CI:", auc_lower_formatted, "-", auc_upper_formatted, ")\n")
      } else {
        # Format AUC estimate to 2 decimal places
        auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
        
        # Print AUC without confidence intervals
        cat("\nSubgroup:", subgroup_name, "\n")
        cat("  AUC:", auc_estimate_formatted, "\n")
      }
    }
  } else {
    # Only one subgroup or overall data
    summary.overall <- summaries$subgroups[[valid.subgroup.list[1]]]
    auc_estimate <- summary.overall$AUC$AUC
    
    if (AUC.CI) {
      # Confidence intervals are stored in AUC_CIs
      auc_ci <- AUC_CIs[[valid.subgroup.list[1]]]$CI
      auc_lower <- auc_ci[1]
      auc_upper <- auc_ci[2]
      
      # Format to 2 decimal places
      auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
      auc_lower_formatted <- sprintf("%.2f", auc_lower)
      auc_upper_formatted <- sprintf("%.2f", auc_upper)
      
      # Print AUC with confidence intervals
      cat("AUC:", auc_estimate_formatted, "(95% CI:", auc_lower_formatted, "-", auc_upper_formatted, ")\n")
    } else {
      # Format AUC estimate to 2 decimal places
      auc_estimate_formatted <- sprintf("%.2f", auc_estimate)
      
      # Print AUC without confidence intervals
      cat("AUC:", auc_estimate_formatted, "\n")
    }
  }       
  if (object.return){
    returned.object <- list()
    returned.object$reitsmas <- reitsmas
    returned.object$summaries <- summaries
    returned.object$madauni <- madauni.main
    if(length(subgroup.list)>1){
      returned.object$anova <- reitsmas$anova.reitsma
    }
    returned.object$datasets <- datasets
    returned.object$AUC_CIs <-AUC_CIs
    returned.object$AUC_CI.cache <- auc.cache.metadata
    returned.object$plot.colors <- plot.colors
    returned.object$point.size.range <- point.size.range
    return(returned.object)
  }
  for (sg  in 1:length(subgroup.list)){
    print(paste(subgroup.list[sg], het.string(reitsmas$subgroups[[valid.subgroup.list[sg]]]), sep = " : "))
    
    if(plot.points){
      points(fpr(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
             sens(dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]),
             pch = pch.list[sg],
             cex = dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ][["point.size"]],
             col = points.colors[sg]
      )
    }
  }
}

                                 

# Test each study for a two-dimensional mean shift in transformed sensitivity
# and false-positive rate. Both the null and alternative Reitsma models use ML,
# as required for likelihood-ratio comparisons with different fixed effects.
.dta_bivariate_mean_shift <- function(dat,
                                      alpha = 0.05,
                                      p.adjust.method = "holm") {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number between 0 and 1.")
  }
  if (!(p.adjust.method %in% stats::p.adjust.methods)) {
    stop("Unknown p-value adjustment method: ", p.adjust.method)
  }
  if (nrow(dat) < 4L) {
    stop("Bivariate mean-shift testing requires at least four complete studies.")
  }

  null.fit <- mada::reitsma(data = dat, method = "ml")
  alternative.fits <- vector("list", nrow(dat))
  tests <- vector("list", nrow(dat))
  statistic <- df <- p.value <- rep(NA_real_, nrow(dat))
  fit.error <- rep(NA_character_, nrow(dat))

  for (i in seq_len(nrow(dat))) {
    shifted.dat <- dat
    shifted.dat$.mada_mean_shift <- as.numeric(seq_len(nrow(dat)) == i)

    result <- tryCatch({
      alternative.fit <- mada::reitsma(
        data = shifted.dat,
        formula = cbind(tsens, tfpr) ~ .mada_mean_shift,
        method = "ml"
      )
      test <- stats::anova(null.fit, alternative.fit)
      list(fit = alternative.fit, test = test)
    }, error = function(e) e)

    if (inherits(result, "error")) {
      fit.error[i] <- conditionMessage(result)
    } else {
      alternative.fits[[i]] <- result$fit
      tests[[i]] <- result$test
      statistic[i] <- as.numeric(result$test$statistic[1])
      df[i] <- as.numeric(result$test$statistic[2])
      p.value[i] <- as.numeric(result$test$statistic[3])
    }
  }

  p.adjusted <- rep(NA_real_, nrow(dat))
  valid.p <- which(is.finite(p.value))
  if (length(valid.p) > 0L) {
    p.adjusted[valid.p] <- stats::p.adjust(
      p.value[valid.p],
      method = p.adjust.method
    )
  }

  diagnostics <- data.frame(
    row_index = dat$.original_row,
    unique_row_id = dat[["uniquer.row.id"]],
    study_id = dat[["study.id"]],
    study_name = dat[["names"]],
    chi_squared = statistic,
    df = df,
    p_value = p.value,
    p_adjusted = p.adjusted,
    outlier = !is.na(p.adjusted) & p.adjusted < alpha,
    fit_error = fit.error,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  failed <- sum(!is.na(fit.error))
  if (failed > 0L) {
    warning(
      failed,
      " study-specific mean-shift model(s) could not be fitted; ",
      "their diagnostic results are NA.",
      call. = FALSE
    )
  }

  list(
    null.fit = null.fit,
    alternative.fits = alternative.fits,
    tests = tests,
    diagnostics = diagnostics,
    alpha = alpha,
    p.adjust.method = p.adjust.method
  )
}


dta.outliers <- function(dat,
                         object.return = FALSE,
                         method = c("Bivariate", "DOR"),
                         alpha = 0.05,
                         p.adjust.method = "holm",
                         study.names = NULL, TP = NULL, TN = NULL,
                         FP = NULL, FN = NULL, uniquer.row.id = NULL,
                         study.id = NULL) {
  method <- match.arg(method)
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id,
    complete.only = FALSE
  )
  dat$.original_row <- seq_len(nrow(dat))
  complete <- stats::complete.cases(dat[, c("TP", "TN", "FP", "FN"), drop = FALSE])
  dat <- dat[complete, , drop = FALSE]

  madauni.dat <- NULL
  metafor.object <- NULL
  inf.object <- NULL
  reitsma.fit <- NULL
  mean.shift <- NULL

  if (method == "DOR") {
    # Legacy implementation retained unchanged for backward compatibility.
    madauni.dat <- mada::madauni(x = dat, type = "DOR", method = "DSL")
    metafor.object <- metafor::rma(
      yi = log(madauni.dat$descr$DOR$DOR),
      sei = madauni.dat$descr$DOR$se.lnDOR,
      weights = madauni.dat$weights,
      method = "DL"
    )
    inf.object <- stats::influence(metafor.object)
    graphics::plot(inf.object)
    local.outlier.indices <- which(abs(inf.object$inf$rstudent) > 2)
  } else {
    mean.shift <- .dta_bivariate_mean_shift(
      dat = dat,
      alpha = alpha,
      p.adjust.method = p.adjust.method
    )
    reitsma.fit <- mean.shift$null.fit
    inf.object <- mean.shift$diagnostics
    local.outlier.indices <- which(mean.shift$diagnostics$outlier)
  }

  outlier_indices <- dat$.original_row[local.outlier.indices]
  outlier_ids <- dat[["uniquer.row.id"]][local.outlier.indices]
  outlier_study_ids <- dat[["study.id"]][local.outlier.indices]

  for (std in local.outlier.indices) {
    cat(
      paste0(
        "Unique row ID ",
        dat[["uniquer.row.id"]][std],
        ".   ",
        dat[["names"]][std],
        " is a potential outlier",
        if (method == "Bivariate") {
          paste0(" (adjusted p = ",
                 format.pval(mean.shift$diagnostics$p_adjusted[std], digits = 3),
                 ")")
        } else ""
      ),
      "\n"
    )
  }
  cat("Outlier detection method:", method, "\n")
  cat("Number of potential outliers:", length(local.outlier.indices), "\n")

  if (object.return) {
    # Existing fields are preserved. For the Bivariate method, madauni is NULL
    # and inf contains the study-level mean-shift diagnostic table.
    returned.object <- list(
      madauni = madauni.dat,
      inf = inf.object,
      outlier_indices = outlier_indices,
      outlier_ids = outlier_ids,
      outlier_unique_row_ids = outlier_ids,
      outlier_study_ids = outlier_study_ids,
      method = method,
      metafor = metafor.object,
      reitsma = reitsma.fit,
      mean_shift = mean.shift
    )
    return(returned.object)
  }
}

forest.diag.no <- function(dat,
                           combined = T,
                           ...,
                           outlier.method = c("Bivariate", "DOR"),
                           outlier.alpha = 0.05,
                           outlier.p.adjust.method = "holm",
                           study.names = NULL, TP = NULL, TN = NULL,
                           FP = NULL, FN = NULL, uniquer.row.id = NULL,
                           study.id = NULL) {
  outlier.method <- match.arg(outlier.method)
  # Capture additional arguments
  args_list <- list(...)
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id,
    complete.only = FALSE
  )
  
  # First, run dta.outliers and capture the result
  outliers_result <- dta.outliers(
    dat,
    object.return = TRUE,
    method = outlier.method,
    alpha = outlier.alpha,
    p.adjust.method = outlier.p.adjust.method,
    uniquer.row.id = "uniquer.row.id",
    study.id = "study.id"
  )
  
  # The dta.outliers function already prints outlier information
  # Now, extract the indices of outlier studies
  outlier_indices <- outliers_result$outlier_indices
  
  # Remove outlier studies from the dataframe
  if (length(outlier_indices) > 0) {
    dat_no_outliers <- dat[-outlier_indices, ]
  } else {
    dat_no_outliers <- dat
  }
  
  # Adjust any data-dependent arguments to match the modified dataframe
  adjust_args <- function(arg_value) {
    if (length(outlier_indices) > 0L &&
        is.vector(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[-outlier_indices])
    } else {
      return(arg_value)
    }
  }
  
  # Apply the adjustment to all arguments
  args_list <- lapply(args_list, adjust_args)
  args_list$uniquer.row.id <- "uniquer.row.id"
  args_list$study.id <- "study.id"
  
  # Now, run forest.diag with the modified dataframe and adjusted arguments
  if (combined){
do.call(forest.diag.combined, c(list(dat = dat_no_outliers), args_list))
    } else {
  do.call(forest.diag, c(list(dat = dat_no_outliers), args_list))
    }
}

         

dta.outliers.multi <- function(dat,
                               subgrouping.variable,
                               object.return = FALSE,
                               method = c("Bivariate", "DOR"),
                               alpha = 0.05,
                               p.adjust.method = "holm",
                               study.names = NULL, TP = NULL, TN = NULL,
                               FP = NULL, FN = NULL, uniquer.row.id = NULL,
                               study.id = NULL) {
  method <- match.arg(method)
  dat[["subgrouping.variable"]] <- subgrouping.variable
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  counts <- table(dat[["subgrouping.variable"]])
  minimum.studies <- if (method == "Bivariate") 4L else 3L
  subgroup.list <- names(counts[counts >= minimum.studies])
  subgroup.list <- sort(subgroup.list)
  # Use original subgroup names for consistency
  valid.subgroup.list <- subgroup.list
  dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
  
  metafors <- list()
  madaunis <- list()
  infs <- list()
  reitsmas <- list()
  mean.shifts <- list()
  outlier_ids <- list()  # Store unique identifiers of outliers
  outlier_study_ids <- list()
  
  for (sg in seq_along(subgroup.list)) {
    datsg <- dat[which(dat[["subgrouping.variable"]] == subgroup.list[sg]), ]
    name.sg <- valid.subgroup.list[sg]
    cat(paste0("*** In subgroup [", subgroup.list[sg], "] :\n"))

    subgroup.result <- dta.outliers(
      dat = datsg,
      object.return = TRUE,
      method = method,
      alpha = alpha,
      p.adjust.method = p.adjust.method,
      uniquer.row.id = "uniquer.row.id",
      study.id = "study.id"
    )

    # Single-bracket assignment preserves subgroup names even when a stored
    # value is NULL or no outliers were detected.
    metafors[name.sg] <- list(subgroup.result$metafor)
    madaunis[name.sg] <- list(subgroup.result$madauni)
    infs[name.sg] <- list(subgroup.result$inf)
    reitsmas[name.sg] <- list(subgroup.result$reitsma)
    mean.shifts[name.sg] <- list(subgroup.result$mean_shift)
    outlier_ids[name.sg] <- list(subgroup.result$outlier_ids)
    outlier_study_ids[name.sg] <- list(subgroup.result$outlier_study_ids)
  }
  
  if (object.return) {
    returned.object <- list()
    returned.object$metafors <- metafors
    returned.object$madaunis <- madaunis
    returned.object$infs <- infs
    returned.object$reitsmas <- reitsmas
    returned.object$mean_shifts <- mean.shifts
    returned.object$outlier_ids <- outlier_ids
    returned.object$outlier_unique_row_ids <- outlier_ids
    returned.object$outlier_study_ids <- outlier_study_ids
    returned.object$method <- method
    # Include subgroup names for reference
    returned.object$subgroup.list <- subgroup.list
    return(returned.object)
  }
}

# Modified forest.diag.subgroup.no function
forest.diag.subgroup.no <- function(dat, 
                                    subgrouping.variable,
                                    combined = TRUE,
                                    ..., 
                                    only.subgroups.bigger.than.3 = TRUE,
                                    exclude.outliers.in.subgroups = NULL,
                                    exclude_from_all_subgroups = FALSE,
                                    outlier.method = c("Bivariate", "DOR"),
                                    outlier.alpha = 0.05,
                                    outlier.p.adjust.method = "holm",
                                    study.names = NULL, TP = NULL, TN = NULL,
                                    FP = NULL, FN = NULL,
                                    uniquer.row.id = NULL, study.id = NULL) {
  outlier.method <- match.arg(outlier.method)
  # Add subgrouping.variable to dat
  dat[["subgrouping.variable"]] <- subgrouping.variable
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id,
    complete.only = FALSE
  )
  
  # Capture additional arguments
  args_list <- list(...)
  
  # Run dta.outliers.multi and capture the result
  outliers_result <- dta.outliers.multi(dat, 
                                        subgrouping.variable = subgrouping.variable, 
                                        object.return = TRUE,
                                        method = outlier.method,
                                        alpha = outlier.alpha,
                                        p.adjust.method = outlier.p.adjust.method,
                                        uniquer.row.id = "uniquer.row.id",
                                        study.id = "study.id")
  
  # Extract outlier identifiers and subgroup names
  outlier_ids_list <- outliers_result$outlier_ids  # Named list with subgroup names
  outlier_study_ids_list <- outliers_result$outlier_study_ids
  subgroup_names <- names(outlier_ids_list)
  
  # Determine which subgroups to exclude outliers from
  if (is.null(exclude.outliers.in.subgroups)) {
    # Exclude outliers from all subgroups (default behavior)
    outlier_ids <- unlist(outlier_ids_list, use.names = FALSE)
    outlier_study_ids <- unique(unlist(outlier_study_ids_list, use.names = FALSE))
  } else {
    # Exclude outliers only from specified subgroups
    if (is.numeric(exclude.outliers.in.subgroups)) {
      # Subgroup indices provided
      subgroup_indices <- exclude.outliers.in.subgroups
      # Validate indices
      if (any(subgroup_indices < 1 | subgroup_indices > length(subgroup_names))) {
        stop("exclude.outliers.in.subgroups contains invalid subgroup indices.")
      }
      subgroup_names_to_exclude <- subgroup_names[subgroup_indices]
    } else if (is.character(exclude.outliers.in.subgroups)) {
      # Subgroup names provided
      subgroup_names_to_exclude <- exclude.outliers.in.subgroups
      # Validate names
      invalid_names <- setdiff(subgroup_names_to_exclude, subgroup_names)
      if (length(invalid_names) > 0) {
        stop(paste("The following subgroups are invalid:", paste(invalid_names, collapse = ", ")))
      }
    } else {
      stop("exclude.outliers.in.subgroups must be numeric indices or character names.")
    }
    
    # Subset outlier_ids_list to only include specified subgroups
    outlier_ids_list_subset <- outlier_ids_list[subgroup_names %in% subgroup_names_to_exclude]
    outlier_study_ids_list_subset <- outlier_study_ids_list[subgroup_names %in% subgroup_names_to_exclude]
    outlier_ids <- unlist(outlier_ids_list_subset, use.names = FALSE)
    outlier_study_ids <- unique(unlist(outlier_study_ids_list_subset, use.names = FALSE))
  }
  
  # Remove outlier studies from the dataframe
  if (length(outlier_ids) > 0) {
    if(exclude_from_all_subgroups){
      outlier_ids <- dat[["uniquer.row.id"]][dat[["study.id"]] %in% outlier_study_ids]
    }
    keep.rows <- !dat[["uniquer.row.id"]] %in% outlier_ids
    dat_no_outliers <- dat[keep.rows, , drop = FALSE]
    } else {
    keep.rows <- rep(TRUE, nrow(dat))
    dat_no_outliers <- dat
  }
  
  # Extract the adjusted subgrouping.variable
  subgrouping.variable_no_outliers <- dat_no_outliers[["subgrouping.variable"]]
  
  # Adjust any data-dependent arguments to match the modified dataframe
  adjust_args <- function(arg_value) {
    if (is.vector(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[keep.rows])
    } else if (is.list(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[keep.rows])
    } else {
      return(arg_value)
    }
  }
  
  # Apply the adjustment to all arguments
  args_list <- lapply(args_list, adjust_args)
  args_list$uniquer.row.id <- "uniquer.row.id"
  args_list$study.id <- "study.id"
  
  # Now, run forest.diag.subgroup with the modified dataframe and adjusted arguments
  if (combined) {
    do.call(forest.diag.subgroup.combined, c(list(dat = dat_no_outliers, 
                                                  subgrouping.variable = subgrouping.variable_no_outliers, 
                                                  only.subgroups.bigger.than.3 = only.subgroups.bigger.than.3), 
                                             args_list))
  } else {
    do.call(forest.diag.subgroup, c(list(dat = dat_no_outliers, 
                                         subgrouping.variable = subgrouping.variable_no_outliers, 
                                         only.subgroups.bigger.than.3 = only.subgroups.bigger.than.3), 
                                    args_list))
  }
}









multiple.srocs.no <- function(dat, 
                              subgrouping.variable = NULL, 
                              ..., 
                              object.return = TRUE,
                              AUC.CI.object = NULL,
                              exclude.outliers.in.subgroups = NULL,
                              exclude_from_all_subgroups = F,
                              outlier.method = c("Bivariate", "DOR"),
                              outlier.alpha = 0.05,
                              outlier.p.adjust.method = "holm",
                              study.names = NULL, TP = NULL, TN = NULL,
                              FP = NULL, FN = NULL,
                              uniquer.row.id = NULL, study.id = NULL) {
  outlier.method <- match.arg(outlier.method)
  # Capture additional arguments
  args_list <- list(...)
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id,
    complete.only = FALSE
  )
  
  # Add default values for data-dependent arguments if not provided.
  # Colours are deliberately left unresolved so multiple.srocs() can apply
  # its shared palette standard (including named base-R palettes).
  if (is.null(args_list$pch.list)) {
    args_list$pch.list <- c(16, 15, 17, 18, 15, 13)
  }
  if (is.null(args_list$summary.pch.list)) {
    args_list$summary.pch.list <- c(16, 15, 17, 18, 15, 13)
  }
  
  # Determine the number of subgroups
  if (is.null(subgrouping.variable)) {
    subgrouping.variable <- rep("Pooled", nrow(dat))
  }
  dat[["subgrouping.variable"]] <- subgrouping.variable
  subgroup.list <- unique(subgrouping.variable)
  
  # Run appropriate outlier detection function
  if (length(subgroup.list) <= 1) {
    # Only one subgroup, use dta.outliers
    outliers_result <- dta.outliers(
      dat,
      object.return = TRUE,
      method = outlier.method,
      alpha = outlier.alpha,
      p.adjust.method = outlier.p.adjust.method,
      uniquer.row.id = "uniquer.row.id",
      study.id = "study.id"
    )
    outlier_ids <- outliers_result$outlier_ids  # Use unique identifiers
  } else {
    # Multiple subgroups, use dta.outliers.multi
    outliers_result <- dta.outliers.multi(dat, 
                                          subgrouping.variable = subgrouping.variable, 
                                          object.return = TRUE,
                                          method = outlier.method,
                                          alpha = outlier.alpha,
                                          p.adjust.method = outlier.p.adjust.method,
                                          uniquer.row.id = "uniquer.row.id",
                                          study.id = "study.id")
    
    # Extract outlier identifiers and subgroup names
    outlier_ids_list <- outliers_result$outlier_ids  # Named list with subgroup names
    outlier_study_ids_list <- outliers_result$outlier_study_ids
    subgroup_names <- names(outlier_ids_list)
    
    # Determine which subgroups to exclude outliers from
    if (is.null(exclude.outliers.in.subgroups)) {
      # Exclude outliers from all subgroups (default behavior)
      outlier_ids <- unlist(outlier_ids_list, use.names = FALSE)
      outlier_study_ids <- unique(unlist(outlier_study_ids_list, use.names = FALSE))
    } else {
      # Exclude outliers only from specified subgroups
      if (is.numeric(exclude.outliers.in.subgroups)) {
        # Subgroup indices provided
        subgroup_indices <- exclude.outliers.in.subgroups
        # Validate indices
        if (any(subgroup_indices < 1 | subgroup_indices > length(subgroup_names))) {
          stop("exclude.outliers.in.subgroups contains invalid subgroup indices.")
        }
        subgroup_names_to_exclude <- subgroup_names[subgroup_indices]
      } else if (is.character(exclude.outliers.in.subgroups)) {
        # Subgroup names provided
        subgroup_names_to_exclude <- exclude.outliers.in.subgroups
        # Validate names
        invalid_names <- setdiff(subgroup_names_to_exclude, subgroup_names)
        if (length(invalid_names) > 0) {
          stop(paste("The following subgroups are invalid:", paste(invalid_names, collapse = ", ")))
        }
      } else {
        stop("exclude.outliers.in.subgroups must be numeric indices or character names.")
      }
      
      # Subset outlier_ids_list to only include specified subgroups
      outlier_ids_list_subset <- outlier_ids_list[subgroup_names %in% subgroup_names_to_exclude]
      outlier_study_ids_list_subset <- outlier_study_ids_list[subgroup_names %in% subgroup_names_to_exclude]
      outlier_ids <- unlist(outlier_ids_list_subset, use.names = FALSE)
      outlier_study_ids <- unique(unlist(outlier_study_ids_list_subset, use.names = FALSE))
    }
    
  }
  
  # Remove outlier studies from the dataframe
  if (length(outlier_ids) > 0) {
    if(exclude_from_all_subgroups && length(subgroup.list) > 1){
      outlier_ids <- dat[["uniquer.row.id"]][dat[["study.id"]] %in% outlier_study_ids]
    }
    keep.rows <- !dat[["uniquer.row.id"]] %in% outlier_ids
    dat_no_outliers <- dat[keep.rows, , drop = FALSE]
  } else {
    keep.rows <- rep(TRUE, nrow(dat))
    dat_no_outliers <- dat
  }
  
  # Extract the adjusted subgrouping.variable
  subgrouping.variable_no_outliers <- dat_no_outliers[["subgrouping.variable"]]
  
  # Adjust any data-dependent arguments to match the modified dataframe
  adjust_args <- function(arg_value) {
    if (is.vector(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[keep.rows])
    } else if (is.list(arg_value) && length(arg_value) == nrow(dat)) {
      return(arg_value[keep.rows])
    } else {
      return(arg_value)
    }
  }
  
  # Apply the adjustment to all arguments except subgrouping.variable
  args_list <- lapply(args_list, adjust_args)
  args_list$uniquer.row.id <- "uniquer.row.id"
  args_list$study.id <- "study.id"
  
  # Now, run multiple.srocs with the modified dataframe and adjusted arguments
  do.call(multiple.srocs, c(list(dat = dat_no_outliers, 
                                 subgrouping.variable = subgrouping.variable_no_outliers, 
                                 object.return = object.return, 
                                 AUC.CI.object = AUC.CI.object), 
                            args_list))
}

         
PBS3 <- function(y,S,b0,V0){
  
  N <- dim(y)[1]
  p <- dim(y)[2]
  
  y.pb <- matrix(numeric(N*p),N)
  
  for(i in 1:N){
    
    yi <- y[i,]
    Si <- matrix(c(S[i,1],S[i,2],S[i,2],S[i,3]),p)
    
    Vi <- Si + V0
    
    Xi <- diag( sqrt(diag(Si) + diag(V0))^-1 )
    Psii <- Xi %*% Vi %*% t(Xi)
    
    mui <- Xi %*% b0
    
    y.pb[i,] <- MASS::ginv(Xi) %*% MASS::mvrnorm(1, mui, Psii)
    
  }
  
  return(y.pb)
  
}


MVPBT_bootstrap <- function(y, S, B = 2000, parallel = TRUE,
                            n.cores = NULL, progress = TRUE, seed = NULL) {
  y <- as.matrix(y)
  S <- as.matrix(S)
  if (nrow(y) < 1L || ncol(y) != 2L) {
    stop("'y' must be a two-column matrix with at least one row.")
  }
  if (nrow(S) != nrow(y) || ncol(S) != 3L) {
    stop("'S' must have one row per study and three covariance columns.")
  }

  V0 <- mvmeta::mvmeta(y, S)$Psi
  Q0 <- MVPBT::MVPBT2(y, S)
  pbs3 <- PBS3

  one.bootstrap <- function(b) {
    tryCatch({
      y.pb <- pbs3(y, S, Q0$b0, V0)
      Q.b <- MVPBT::MVPBT2(y.pb, S)
      as.numeric(Q.b$T)
    }, error = function(e) NA_real_)
  }

  bootstrap <- .dta_bootstrap_map(
    B = B,
    FUN = one.bootstrap,
    parallel = parallel,
    n.cores = n.cores,
    progress = progress,
    label = "MVPBT bootstrap",
    seed = seed
  )
  T.b <- as.numeric(unlist(bootstrap$values, use.names = FALSE))
  successful <- sum(is.finite(T.b))
  if (successful == 0L) {
    stop("All MVPBT bootstrap model fits failed.")
  }
  if (successful < B) {
    warning(B - successful, " of ", B,
            " MVPBT bootstrap fits failed and were omitted from the p-value.")
  }
  T.b.valid <- T.b[is.finite(T.b)]

  # Equivalent to the original upper-tail empirical calculation when there are no ties.
  P <- sum(T.b.valid > as.numeric(Q0$T)) / (length(T.b.valid) + 1L)
  list(
    T.b = T.b,
    T = Q0$T,
    P = P,
    n.boots = B,
    n.successful = successful,
    parallel = bootstrap$parallel,
    n.cores = bootstrap$n.cores
  )
}

# Backward-compatible name retained for existing analyses.
MVPBT_boot <- function(y, S, B = 2000, parallel = TRUE,
                       n.cores = NULL, progress = TRUE, seed = NULL) {
  MVPBT_bootstrap(
    y = y, S = S, B = B, parallel = parallel, n.cores = n.cores,
    progress = progress, seed = seed
  )
}


pubbias.diag <- function(dat, n.boots = 2000, parallel = TRUE,
                         n.cores = NULL, progress = TRUE, seed = NULL,
                         study.names = NULL, TP = NULL, TN = NULL,
                         FP = NULL, FN = NULL, uniquer.row.id = NULL,
                         study.id = NULL){
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  dta.dat <- MVPBT::edta(
    TP = dat[["TP"]], FN = dat[["FN"]],
    TN = dat[["TN"]], FP = dat[["FP"]]
  )
  y <- dta.dat[["y"]]
  S <- dta.dat[["S"]]

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  res1 <- metafor::rma(y[, 1], S[, 1])
  metafor::funnel(res1, main = "(b) Funnel plot for logit(Se)")
  res2 <- metafor::rma(y[, 2], S[, 3])
  metafor::funnel(res2, main = "(c) Funnel plot for logit(FPR)")

  MVPBT3.dat <- MVPBT_bootstrap(
    y = y, S = S, B = n.boots, parallel = parallel, n.cores = n.cores,
    progress = progress, seed = seed
  )
  print(MVPBT3.dat[["P"]])
  return(MVPBT3.dat)
}



find_repeated_studies <- function(dat, uniquer.row.id = NULL,
                                  study.id = NULL, study.names = NULL) {
  dat[["uniquer.row.id"]] <- if (is.null(uniquer.row.id)) {
    seq_len(nrow(dat))
  } else {
    .dta_column_values(dat, uniquer.row.id, NULL, "uniquer.row.id")
  }
  dat[["study.id"]] <- if (is.null(study.id)) {
    study.names <- .dta_column_values(dat, study.names, "names", "study.names")
    study.names
  } else {
    .dta_column_values(dat, study.id, NULL, "study.id")
  }
  if (anyNA(dat[["uniquer.row.id"]]) || anyDuplicated(dat[["uniquer.row.id"]])) {
    stop("'uniquer.row.id' must contain unique, non-missing values.")
  }

  repeated_studies <- dat %>%
    group_by(.data[["study.id"]]) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(.data[["study.id"]])
  
  # For each repeated study ID, report its unique row IDs.
  if (length(repeated_studies) > 0) {
    for(study in repeated_studies) {
      associated_rows <- dat %>%
        filter(.data[["study.id"]] == study) %>%
        pull(.data[["uniquer.row.id"]])
      
      cat("Study ID", study, "has the following unique row IDs:\n")
      print(unique(associated_rows))
      cat("\n")
    }
  } else {
    print("No repeated studies found.")
  }
  
  # Return the repeated study IDs.
  return(repeated_studies)
}

find_repeated_studies_by_subgroup <- function(dat, subgroup_var,
                                              uniquer.row.id = NULL,
                                              study.id = NULL,
                                              study.names = NULL) {
  dat[["uniquer.row.id"]] <- if (is.null(uniquer.row.id)) {
    seq_len(nrow(dat))
  } else {
    .dta_column_values(dat, uniquer.row.id, NULL, "uniquer.row.id")
  }
  dat[["study.id"]] <- if (is.null(study.id)) {
    study.names <- .dta_column_values(dat, study.names, "names", "study.names")
    study.names
  } else {
    .dta_column_values(dat, study.id, NULL, "study.id")
  }
  if (anyNA(dat[["uniquer.row.id"]]) || anyDuplicated(dat[["uniquer.row.id"]])) {
    stop("'uniquer.row.id' must contain unique, non-missing values.")
  }
  subgroup.values <- .dta_column_values(dat, subgroup_var, NULL, "subgroup_var")
  dat[[".dta_subgroup"]] <- subgroup.values

  # Unique subgroups
  subgroups <- unique(dat[[".dta_subgroup"]])
  
  for (subgroup in subgroups) {
    cat("For subgroup:", subgroup, "\n")
    cat("-------------------------------\n")
    
    subgroup_data <- dat %>% filter(.data[[".dta_subgroup"]] == subgroup)
    
    # Identify study IDs that have been repeated within the subgroup.
    repeated_studies <- subgroup_data %>%
      group_by(.data[["study.id"]]) %>%
      tally() %>%
      filter(n > 1) %>%
      pull(.data[["study.id"]])
    
    # For each repeated study ID, report its unique row IDs.
    if (length(repeated_studies) > 0) {
      for(study in repeated_studies) {
        associated_rows <- subgroup_data %>%
          filter(.data[["study.id"]] == study) %>%
          pull(.data[["uniquer.row.id"]])
        
        cat("Study ID", study, "has the following unique row IDs:\n")
        print(unique(associated_rows))
        cat("\n")
      }
    } else {
      print(paste("No repeated studies found in subgroup", subgroup))
      cat("\n")
    }
    cat("-------------------------------\n\n")
  }
}

.dta_lr_legend_content <- function(subgroups, sum.colors) {
  decision.labels <- c(
    "LUQ  Exclusion + confirmation  (LR+ > 10, LR- < 0.1)",
    "RUQ  Confirmation only  (LR+ > 10, LR- > 0.1)",
    "LLQ  Exclusion only  (LR+ < 10, LR- < 0.1)",
    "RLQ  Neither  (LR+ < 10, LR- > 0.1)"
  )
  subgroup.labels <- c(
    subgroups,
    "",
    "Summary estimate with 95% CI",
    "Study estimate"
  )
  list(
    decision = list(
      legend = decision.labels,
      pch = rep(NA_integer_, length(decision.labels)),
      col = rep("black", length(decision.labels)),
      pt.cex = rep(1, length(decision.labels)),
      title = "Decision regions"
    ),
    subgroup = list(
      legend = subgroup.labels,
      pch = c(rep(18, length(subgroups)), NA, 18, 16),
      col = c(sum.colors, NA, "#555555", "#999999"),
      pt.cex = c(rep(1.5, length(subgroups)), 1, 1.5, 1),
      title = "Subgroups and symbols"
    )
  )
}

.dta_draw_lr_legend_panel <- function(content, arrangement = c("vertical", "columns"),
                                      cex = 0.78) {
  arrangement <- match.arg(arrangement)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))

  legend.args <- function(spec, cex, plot = FALSE) {
    c(
      list(
        x = 0, y = 1, xjust = 0, yjust = 1, cex = cex,
        bty = "n", y.intersp = 1.1, x.intersp = 0.8,
        title.adj = 0, plot = plot
      ),
      spec
    )
  }

  measured <- NULL
  for (i in seq_len(20L)) {
    decision <- do.call(graphics::legend, legend.args(content$decision, cex))
    subgroup <- do.call(graphics::legend, legend.args(content$subgroup, cex))
    if (arrangement == "columns") {
      fits <- max(decision$rect$w, subgroup$rect$w) <= 0.46 &&
        max(decision$rect$h, subgroup$rect$h) <= 0.94
    } else {
      fits <- max(decision$rect$w, subgroup$rect$w) <= 0.94 &&
        decision$rect$h + subgroup$rect$h + 0.03 <= 0.96
    }
    measured <- list(decision = decision, subgroup = subgroup)
    if (fits || cex <= 0.48) {
      break
    }
    cex <- max(0.48, cex * 0.92)
  }

  if (arrangement == "columns") {
    decision.args <- legend.args(content$decision, cex, plot = TRUE)
    decision.args$x <- 0.02
    decision.args$y <- 0.98
    subgroup.args <- legend.args(content$subgroup, cex, plot = TRUE)
    subgroup.args$x <- 0.52
    subgroup.args$y <- 0.98
  } else {
    decision.args <- legend.args(content$decision, cex, plot = TRUE)
    decision.args$x <- 0.03
    decision.args$y <- 0.98
    subgroup.args <- legend.args(content$subgroup, cex, plot = TRUE)
    subgroup.args$x <- 0.03
    subgroup.args$y <- 0.98 - measured$decision$rect$h - 0.03
  }
  do.call(graphics::legend, decision.args)
  do.call(graphics::legend, subgroup.args)
  invisible(cex)
}

multiple.LRmats <- function(dat,
                            subgrouping.variable = NULL,
                            sum.colors = c("blue", "maroon", "black", "skyblue", "#20cb20", "red"),
                            points.colors = c("#0000FF20", "#A52A2A20", "#1d1c1c20", "#87ceeb30", "#20cb2020", "#ff000050"),
                            inset_var = -.4,
                            study.names = NULL, TP = NULL, TN = NULL,
                            FP = NULL, FN = NULL, uniquer.row.id = NULL,
                            study.id = NULL,
                            legend.position = c("auto", "right", "bottom", "inside", "none"),
                            legend.cex = NULL,
                            min.plot.width = 4.2,
                            xlim = c(0.01, 1),
                            ylim = c(1, 100),
                            color.palette = "mada6",
                            point.alpha = 0.18,
                            point.size.range = c(0.7, 1.8)) {
  legacy.sum.colors <- if (missing(sum.colors)) NULL else sum.colors
  legacy.points.colors <- if (missing(points.colors)) NULL else points.colors
  legend.position <- match.arg(legend.position)
  if (!is.null(subgrouping.variable) && length(subgrouping.variable) != nrow(dat)) {
    stop("'subgrouping.variable' must contain one value per row of 'dat'.")
  }
  if (length(xlim) != 2L || any(!is.finite(xlim)) || any(xlim <= 0) ||
      xlim[1] >= xlim[2]) {
    stop("'xlim' must contain two increasing positive values for the log scale.")
  }
  if (length(ylim) != 2L || any(!is.finite(ylim)) || any(ylim <= 0) ||
      ylim[1] >= ylim[2]) {
    stop("'ylim' must contain two increasing positive values for the log scale.")
  }
  if (!is.numeric(min.plot.width) || length(min.plot.width) != 1L ||
      !is.finite(min.plot.width) || min.plot.width <= 0) {
    stop("'min.plot.width' must be one positive number in inches.")
  }

  dat[["subgrouping.variable"]] <- if (is.null(subgrouping.variable)) {
    "Pooled"
  } else {
    subgrouping.variable
  }
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )

  counts <- table(dat[["subgrouping.variable"]])
  subgroup.list <- sort(names(counts[counts >= 3L]))
  if (length(subgroup.list) == 0L) {
    stop("At least one subgroup with three or more complete studies is required.")
  }
  valid.subgroup.list <- make.names(subgroup.list, unique = TRUE)
  dat <- dat[dat[["subgrouping.variable"]] %in% subgroup.list, , drop = FALSE]
  plot.colors <- .dta_resolve_plot_colours(
    n = length(subgroup.list),
    color.palette = color.palette,
    legacy.primary = legacy.sum.colors,
    legacy.points = legacy.points.colors,
    point.alpha = point.alpha
  )
  sum.colors <- plot.colors$primary
  points.colors <- plot.colors$points
  dat[["population"]] <- .dta_population(dat)
  dat[["point.size"]] <- .dta_population_point_sizes(
    dat[["population"]], point.size.range
  )

  reitsmas <- list(subgroups = list())
  pLRs <- list(subgroups = list())
  nLRs <- list(subgroups = list())
  summaries <- list(subgroups = list())
  datasets <- list(main = dat, subgroups = list())

  for (sg in seq_along(subgroup.list)) {
    subgroup.data <- dat[
      dat[["subgrouping.variable"]] == subgroup.list[sg], , drop = FALSE
    ]
    name <- valid.subgroup.list[sg]
    reitsma.sg <- reitsma(data = subgroup.data, method = "ml")
    reitsmas$subgroups[[name]] <- reitsma.sg
    summaries$subgroups[[name]] <- summary(SummaryPts(reitsma.sg))
    subgroup.madad <- madad(subgroup.data)
    pLRs$subgroups[[name]] <- subgroup.madad$posLR$posLR
    nLRs$subgroups[[name]] <- subgroup.madad$negLR$negLR
    datasets$subgroups[[name]] <- subgroup.data
  }

  device.size <- grDevices::dev.size("in")
  if (is.null(legend.cex)) {
    legend.cex <- max(0.62, min(0.82, 0.88 - 0.03 * length(subgroup.list)))
  }
  if (!is.numeric(legend.cex) || length(legend.cex) != 1L ||
      !is.finite(legend.cex) || legend.cex <= 0) {
    stop("'legend.cex' must be NULL or one positive number.")
  }
  legend.content <- .dta_lr_legend_content(subgroup.list, sum.colors)
  longest.label <- max(nchar(c(
    legend.content$decision$legend,
    legend.content$subgroup$legend
  )))
  estimated.legend.width <- longest.label * graphics::par("cin")[1] *
    legend.cex * 0.72 + 0.65

  resolved.position <- legend.position
  if (resolved.position == "auto") {
    resolved.position <- if (
      device.size[1] >= min.plot.width + estimated.legend.width + 0.5
    ) "right" else "bottom"
  }
  bottom.arrangement <- if (device.size[1] >= 6.4) "columns" else "vertical"

  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::layout(matrix(1))
    suppressWarnings(graphics::par(old.par))
  }, add = TRUE)

  if (resolved.position == "right") {
    legend.width <- min(
      max(estimated.legend.width, 2.6),
      max(2.2, device.size[1] - min.plot.width)
    )
    plot.width <- max(2.2, device.size[1] - legend.width)
    graphics::layout(matrix(c(1, 2), nrow = 1L), widths = c(plot.width, legend.width))
    graphics::par(mar = c(5.0, 4.6, 1.1, 0.5) + 0.1, xpd = FALSE)
  } else if (resolved.position == "bottom") {
    legend.rows <- if (bottom.arrangement == "columns") {
      max(5L, length(subgroup.list) + 4L)
    } else {
      5L + length(subgroup.list) + 4L
    }
    legend.height <- min(
      max(1.45, legend.rows * graphics::par("cin")[2] * legend.cex * 1.15),
      max(1.45, device.size[2] * 0.42)
    )
    plot.height <- max(2.5, device.size[2] - legend.height)
    graphics::layout(matrix(c(1, 2), ncol = 1L), heights = c(plot.height, legend.height))
    graphics::par(mar = c(4.5, 4.6, 1.0, 1.0) + 0.1, xpd = FALSE)
  } else {
    graphics::layout(matrix(1))
    graphics::par(mar = c(5.0, 4.6, 1.1, 1.1) + 0.1, xpd = FALSE)
  }

  graphics::plot(
    nLRs$subgroups[[valid.subgroup.list[1]]],
    pLRs$subgroups[[valid.subgroup.list[1]]],
    log = "xy", type = "n", xlim = xlim, ylim = ylim,
    xlab = "Negative Likelihood Ratio",
    ylab = "Positive Likelihood Ratio",
    cex.lab = 1.1, cex.axis = 0.75, las = 1
  )
  usr <- graphics::par("usr")
  graphics::segments(
    x0 = 10^usr[1], x1 = 10^usr[2], y0 = 10, y1 = 10,
    col = "darkgray", lty = 2, lwd = 3
  )
  graphics::segments(
    x0 = 0.1, x1 = 0.1, y0 = 10^usr[3], y1 = 10^usr[4],
    col = "darkgray", lty = 2, lwd = 3
  )

  for (sg in seq_along(subgroup.list)) {
    name <- valid.subgroup.list[sg]
    subgroup.summary <- summaries$subgroups[[name]]
    graphics::points(
      nLRs$subgroups[[name]], pLRs$subgroups[[name]],
      pch = 19, col = points.colors[sg],
      cex = datasets$subgroups[[name]][["point.size"]]
    )
    graphics::points(
      subgroup.summary["negLR", 1], subgroup.summary["posLR", 1],
      pch = 18, col = sum.colors[sg], cex = 2.6
    )
    ci.color <- grDevices::adjustcolor(sum.colors[sg], alpha.f = 0.6)
    graphics::arrows(
      subgroup.summary["negLR", 1], subgroup.summary["posLR", 3],
      subgroup.summary["negLR", 1], subgroup.summary["posLR", 4],
      length = 0.05, col = ci.color, angle = 90, code = 3, lwd = 2
    )
    graphics::arrows(
      subgroup.summary["negLR", 3], subgroup.summary["posLR", 1],
      subgroup.summary["negLR", 4], subgroup.summary["posLR", 1],
      length = 0.05, col = ci.color, angle = 90, code = 3, lwd = 2
    )
  }

  effective.legend.cex <- legend.cex
  if (resolved.position %in% c("right", "bottom")) {
    graphics::par(mar = rep(0.1, 4), xpd = NA)
    effective.legend.cex <- .dta_draw_lr_legend_panel(
      legend.content,
      arrangement = if (resolved.position == "right") "vertical" else bottom.arrangement,
      cex = legend.cex
    )
  } else if (resolved.position == "inside") {
    compact.legend <- c(
      legend.content$decision$legend,
      "",
      legend.content$subgroup$legend
    )
    compact.pch <- c(
      legend.content$decision$pch,
      NA,
      legend.content$subgroup$pch
    )
    compact.col <- c(
      legend.content$decision$col,
      NA,
      legend.content$subgroup$col
    )
    compact.pt.cex <- c(
      legend.content$decision$pt.cex,
      1,
      legend.content$subgroup$pt.cex
    )
    inside.inset <- if (identical(inset_var, -0.4)) 0.02 else inset_var
    graphics::legend(
      "topright", inset = c(inside.inset, 0.02), legend = compact.legend,
      pch = compact.pch, col = compact.col, pt.cex = compact.pt.cex,
      cex = legend.cex, bty = "o", box.lwd = 1, box.col = "black",
      y.intersp = 1.1, xpd = NA
    )
  }

  list(
    reitsmas = reitsmas,
    pLRs = pLRs,
    nLRs = nLRs,
    summaries = summaries,
    datasets = datasets,
    plot.colors = plot.colors,
    point.size.range = point.size.range,
    plot.layout = list(
      requested = legend.position,
      resolved = resolved.position,
      bottom.arrangement = if (resolved.position == "bottom") bottom.arrangement else NULL,
      device.inches = stats::setNames(device.size, c("width", "height")),
      legend.cex = effective.legend.cex
    )
  )
}













nomogrammer <- function(Prevalence,
                        Sens = NULL,
                        Spec = NULL,
                        Plr = NULL,
                        Nlr = NULL,
                        Detail = T,
                        NullLine = T,
                        LabelSize = (12/5),
                        Verbose = FALSE,
                        x_var = .75,
                        y_var = 2){
  
  ## Function inputs:
  # Prevalence (prior probability) as a number between 0 and 1
  # Either
  # Sens & Spec
  # model sensitivity and specificity as a number between 0 and 1
  # Or
  # Likelihood ratios
  # Positive and Negative LRs (numeric)
  
  ## Function options:
  # Detail: If true, will overlay key statistics onto the plot
  # NullLine: If true, will add a line from prior prob through LR = 1
  # LabelSize: Tweak this number to change the label sizes
  # Verbose: Print out relevant metrics in the conso
  
  
  
  
  
  ## Helper functions
  ##   (defined inside nomogrammer, so remain local only & wont clutter user env)
  odds         <- function(p){
    # Function converts probability into odds
    o <- p/(1-p)
    return(o)
  }
  
  logodds      <- function(p){
    # Function returns logodds for a probability
    lo <- log10(p/(1-p))
    return(lo)
  }
  
  logodds_to_p <- function(lo){
    # Function goes from logodds back to a probability
    o <- 10^lo
    p <- o/(1+o)
    return(p)
  }
  
  p2percent <- function(p){
    # Function turns numeric probability into string percentage
    # e.g. 0.6346111 -> 63.5% 
    # scales::percent(signif(p, digits = 4))
    round(p*100, digits = 2)
    }
  
  
  ######################################
  ########## Calculations     ##########
  ######################################
  
  ## Checking inputs
  
  ## Prevalence
  # needs to exist
  if(missing(Prevalence)){
    stop("Prevalence is missing")
  }
  # needs to be numeric
  if(!is.numeric(Prevalence)){stop("Prevalence should be numeric")}
  # needs to be a prob not a percent
  if((Prevalence > 1) | (Prevalence <= 0)){stop("Prevalence should be a probability (did you give a %?)")}
  
  # Did user give sens & spec?
  if(missing(Sens) | missing(Spec)){
    sensspec <- FALSE
  } else{ sensspec <- TRUE}
  # if yes, make sure they are numbers
  if(sensspec == TRUE){
    if(!is.numeric(Sens)){stop("Sensitivity should be numeric")}
    if(!is.numeric(Spec)){stop("Specificity should be numeric")}
    # numbers that are probabilities not percentages
    if((Sens > 1) | (Sens <= 0)){stop("Sensitivity should be a probability (did you give a %?)")}
    if((Spec > 1) | (Spec <= 0)){stop("Specificity should be a probability (did you give a %?)")}
  }
  
  
  # Did user give PLR & NLR?
  if(missing(Plr) | missing(Nlr)){
    plrnlr <- FALSE
  } else{plrnlr <- TRUE}
  # if yes, make sure they are numbers
  if(plrnlr == TRUE){
    if(!is.numeric(Plr)){stop("PLR should be numeric")}
    if(!is.numeric(Nlr)){stop("NLR should be numeric")}
    # numbers that vaguely make sense
    if(Plr < 1){stop("PLR shouldn't be less than 1")}
    if(Nlr < 0){stop("NLR shouldn't be below zero")}
    if(Nlr > 1){stop("NLR shouldn't be more than 1")}
  }
  
  # Did they give a valid sensspec and plrnlr? If yes, ignore the LRs and tell them
  if((sensspec == TRUE) && (plrnlr == TRUE) ){
    warning("You provided sens/spec as well as likelihood ratios-- I ignored the LRs!")
  }
  
  
  ## If sens/spec provided, we calculate posterior probabilities & odds using sens & spec
  ##  otherwise, if plr and nlr provided, we calculate posteriors using them
  ##  if neither exist, then return an error
  if(sensspec == TRUE){
    prior_prob  <- Prevalence
    prior_odds  <- odds(prior_prob)
    sensitivity <- Sens
    specificity <- Spec
    PLR <- sensitivity/(1-specificity)
    NLR <- (1-sensitivity)/specificity
    post_odds_pos  <- prior_odds * PLR
    post_odds_neg  <- prior_odds * NLR
    post_prob_pos  <- post_odds_pos/(1+post_odds_pos)
    post_prob_neg  <- post_odds_neg/(1+post_odds_neg)
  } else if(plrnlr == TRUE){
    prior_prob  <- Prevalence
    prior_odds  <- odds(prior_prob)
    PLR <- Plr
    NLR <- Nlr
    sensitivity <- (PLR*(1-NLR))/(PLR-NLR)    ## TODO: check Adam's math! 
    specificity <- (1-PLR)/(NLR-PLR)          ## TODO: check Adam's math! 
    post_odds_pos  <- prior_odds * PLR
    post_odds_neg  <- prior_odds * NLR
    post_prob_pos  <- post_odds_pos/(1+post_odds_pos)
    post_prob_neg  <- post_odds_neg/(1+post_odds_neg)
  } else{
    stop("Couldn't find sens & spec, or positive & negative likelihood ratios")
  }
  
  
  
  ######################################
  ########## Plotting (prep)  ##########
  ######################################
  
  
  ## Set common theme preferences up front
  theme_set(theme_bw() +
              theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_text(angle = 90, face = "bold"),
                    axis.title.y.right = element_text(angle = 90, face = "bold"),
                    axis.line = element_blank(),
                    panel.grid = element_blank(),
                    legend.position = "none",
                    panel.background = element_rect(fill = "gray90"),
                    plot.title = element_text(hjust = 0.5)
              )
  )
  
  ## Setting up the points of interest along the y-axes
  
  # Select probabilities of interest (nb as percentages)
  ticks_prob    <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30,
                     40, 50, 60, 70, 80, 90, 95, 99, 99.5, 99.9)
  # Convert % to odds
  ticks_odds    <- odds(ticks_prob/100)
  # Convert % to logodds 
  ticks_logodds <- logodds(ticks_prob/100)
  
  # Select the likelihood ratios of interest (for the middle y-axis)
  ticks_lrs     <- sort(c(10^(-3:3), 2*(10^(-3:2)), 5*(10^(-3:2))))
  # Log10 them since plot is in logodds space
  ticks_log_lrs <- log10(ticks_lrs)
  
  
  
  
  ## Fixing particular x-coordinates
  left     <- 0
  right    <- 1
  middle   <- 0.5
  midright <- 0.75
  
  ## Lay out the four key plot points
  ##  (the start and finish of the positive and negative lines)
  
  # Initially these are expressed as probabilities
  df <- data.frame(x=c(left, right, left, right), 
                   y=c(prior_prob, post_prob_pos, prior_prob, post_prob_neg), 
                   line = c("pos", "pos", "neg", "neg"))
  
  adj_min      <- range(ticks_logodds)[1]
  adj_max      <- range(ticks_logodds)[2]
  adj_diff     <- adj_max - adj_min
  scale_factor <- abs(adj_min) - adj_diff/2
  #df$lo_y <- ifelse(df$x==left,(10/adj_diff)*logodds(1-df$y)-1,logodds(df$y))
  
  # Convert probabilities to logodds for plotting
  df$lo_y  <- ifelse(df$x==left,logodds(1-df$y)-scale_factor,logodds(df$y))
  # zero         <- data.frame(x = c(left,right),
  #                            y = c(0,0),
  #                            line = c('pos','pos'),
  #                            lo_y = c(-scale_factor,0))
  
  df$lo_y <- ifelse(is.na(df$lo_y)|is.infinite(df$lo_y)|is.null(df$lo_y), 0,df$lo_y)
  
  
  
  rescale   <- range(ticks_logodds) + abs(adj_min) - adj_diff/2
  rescale_x_breaks  <- ticks_logodds + abs(adj_min) - adj_diff/2  
  
  
  
  ######################################
  ########## Plot             ##########
  ######################################
  
  
  p <- ggplot(df) +
    geom_line(aes(x = x, y = lo_y, color = line), size = 1) +
    geom_vline(xintercept = middle) +
    annotate(geom = "text",
             x = rep(middle+.075, length(ticks_log_lrs)),
             y = (ticks_log_lrs-scale_factor)/2,
             label = ticks_lrs,
             size = rel(LabelSize)) +
    annotate(geom="point",
             x = rep(middle, length(ticks_log_lrs)),
             y = (ticks_log_lrs-scale_factor)/2,
             size = 1) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0),
                       limits = rescale,
                       breaks = -rescale_x_breaks,
                       labels = ticks_prob,
                       name = "Pre-test probablity (%)",
                       sec.axis = sec_axis(trans = ~.,
                                           name = "Post-test probablity (%)",
                                           labels = ticks_prob,
                                           breaks = ticks_logodds))
  
  ## Optional overlay text: prevalence, PLR/NLR, and posterior probabilities
  detailedAnnotation <- paste(
    paste("prevalence = ", p2percent(prior_prob), "%", sep= ""),
    paste("PLR =", signif(PLR, 3),", NLR =", signif(NLR, 3)),
    paste("post. pos = ", p2percent(post_prob_pos),
          "% , neg = ", p2percent(post_prob_neg), "%", sep = ""),
    sep = "\n")
  
  
  ## Optional amendments to the plot
  
  ## Do we add the null line i.e. LR = 1, illustrating an uninformative model
  if(NullLine == TRUE){
    ## If yes, first calculate the start and end points
    uninformative <- data.frame(
      x = c(left,right),
      lo_y = c( (logodds(1-prior_prob) - scale_factor) , logodds(prior_prob))
    ) 
    
    p <- p + geom_line(aes(x = x, y = lo_y), data = uninformative,
                       color = "gray20", 
                       lty = 2,
                       inherit.aes = FALSE)
  }
  
  
  ## Do we add the detailed stats to the top right?
  if(Detail == TRUE){
    p <- p + annotate(geom = "text",
                      x = x_var,
                      y = y_var,
                      label = detailedAnnotation,
                      size = rel(LabelSize))
  }
  
  if(Verbose == TRUE){
    writeLines(
      text = c(
        paste0("prevalence = ", p2percent(prior_prob)),
        paste("PLR =", signif(PLR, 3)),
        paste("NLR =", signif(NLR, 3)),
        paste("posterior probability (positive) =", p2percent(post_prob_pos)),
        paste("posterior probability (negative) =", p2percent(post_prob_neg)),
        paste("sensitivity =", p2percent(sensitivity)),
        paste("specificity =", p2percent(specificity))
        # sep = "\n"
      )
    )
  }
  
  
  return(p)
  
}



nomogrammer_plus <- function(dat, prevalence, x_var = .75, y_var = 2,
                             return.list = F, alphabet = T,
                             study.names = NULL, TP = NULL, TN = NULL,
                             FP = NULL, FN = NULL, uniquer.row.id = NULL,
                             study.id = NULL) {
  dat <- .dta_prepare_data(
    dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
  )
  posLR <- summary(SummaryPts(reitsma(dat, method = "ml")))["posLR", 1]
  negLR <- summary(SummaryPts(reitsma(dat, method = "ml")))["negLR", 1]
  if (alphabet){
      alphabet <- c("A", "B", "C", "D", "E", "F")
    } else {
      alphabet <- c("", "", "", "", "", "")
    }
  
  plots <- list() # List to store ggplot objects
  for(PR in seq_along(prevalence)) {
    plot <- nomogrammer(Prevalence = prevalence[PR], 
                        Plr = posLR,
                        Nlr = negLR,
                        x_var = if (length(x_var) == 1) x_var else x_var[PR] ,
                        y_var = if(length(y_var) == 1) y_var else y_var[PR]) + 
      labs(title = paste(alphabet[PR],". Pre-test probablity", ": ", prevalence[PR]*100, "%", sep = ""))
    plots[[PR]] <- plot
   }
  if (return.list){
    return(plots)
   } else {
  # Combine plots using patchwork
  combined_plot <- patchwork::wrap_plots(plots, ncol = length(prevalence))
  return(combined_plot)
  print(combined_plot) # Print the combined plot
   }
  }
  








nomogrammer_subgroups <- function(dat,
                                  subgrouping.variable,
                                  prevalence,
                                  x_var = .75,
                                  y_var = 2,
                                  return.list = F,
                                  study.names = NULL, TP = NULL, TN = NULL,
                                  FP = NULL, FN = NULL,
                                  uniquer.row.id = NULL, study.id = NULL
                                  ){
dat[["subgrouping.variable"]] <- subgrouping.variable  
dat <- .dta_prepare_data(
  dat, study.names, TP, TN, FP, FN, uniquer.row.id, study.id
)
subgroup.list <- unique(dat[["subgrouping.variable"]])
counts <- table(dat[["subgrouping.variable"]])
subgroup.list <- names(counts[counts >= 3])
subgroup.list <- sort(subgroup.list)
valid.subgroup.list <- make.names(subgroup.list, unique = T)
dat <- subset(dat, dat[["subgrouping.variable"]] %in% subgroup.list)
plots <- list()
for (sg in 1:length(subgroup.list)){
  sgdat <- dat[which(dat$subgrouping.variable == subgroup.list[sg]), ] 
  plr <- summary(SummaryPts(reitsma(sgdat, method = "ml")))["posLR", 1]
  nlr <- summary(SummaryPts(reitsma(sgdat, method = "ml")))["negLR", 1]
  plot <- nomogrammer(Prevalence = prevalence, 
                      Plr = plr,
                      Nlr = nlr,
                      x_var = if (length(x_var) == 1) x_var else x_var[sg] ,
                      y_var = if (length(y_var) == 1) y_var else y_var[sg]) + 
    labs(title = paste("Test: ",subgroup.list[sg], sep = ""))
  plots[[sg]] <- plot
 }
if(return.list){
 return(plots)
 } else {
  combined_plot <- patchwork::wrap_plots(plots, ncol = length(subgroup.list))
  return(combined_plot)
  print(combined_plot)
  }
}























df.mixer <- function (training, intval, extval){
  training[["Cohort"]] <- "Training"
  intval[["Cohort"]] <- "Internal Validation"
  extval[["Cohort"]] <- "External Validation"
  overall <- rbind(training,intval,extval)
  return(overall)
}
