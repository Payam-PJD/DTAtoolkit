# Shared visual standards for diagnostic-accuracy plots.

.dta_default_palette <- c(
  "blue", "maroon", "black", "skyblue", "#20cb20", "red"
)

#' Resolve a colour palette used by DTAtoolkit plots
#'
#' @param palette Either `"mada6"`, the name of a base-R palette from
#'   [grDevices::palette.pals()] or [grDevices::hcl.pals()], or a character
#'   vector of colours for manual entry.
#' @param n Number of colours to return.
#' @param alpha Optional opacity between zero and one.
#' @param reverse Reverse the resolved palette.
#'
#' @return A character vector containing `n` valid R colours.
#' @export
dta_palette <- function(palette = "mada6", n = 6L, alpha = NULL,
                        reverse = FALSE) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != as.integer(n)) {
    stop("'n' must be one positive whole number.", call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.character(palette) || length(palette) < 1L || anyNA(palette)) {
    stop(
      "'palette' must be a palette name or a non-empty character vector of colours.",
      call. = FALSE
    )
  }
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    stop("'reverse' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(alpha) &&
      (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
       alpha < 0 || alpha > 1)) {
    stop("'alpha' must be NULL or one number between zero and one.", call. = FALSE)
  }

  if (length(palette) == 1L && tolower(palette) == "mada6") {
    colours <- rep(.dta_default_palette, length.out = n)
  } else if (length(palette) == 1L) {
    palette.pals <- grDevices::palette.pals()
    palette.index <- match(tolower(palette), tolower(palette.pals))
    hcl.pals <- grDevices::hcl.pals()
    hcl.index <- match(tolower(palette), tolower(hcl.pals))

    if (!is.na(palette.index)) {
      colours <- grDevices::palette.colors(
        n, palette = palette.pals[palette.index], recycle = TRUE
      )
    } else if (!is.na(hcl.index)) {
      colours <- grDevices::hcl.colors(n, palette = hcl.pals[hcl.index])
    } else {
      colour.valid <- tryCatch({
        grDevices::col2rgb(palette)
        TRUE
      }, error = function(e) FALSE)
      if (!colour.valid) {
        stop(
          "Unknown palette or colour '", palette, "'. Use 'mada6', a name from ",
          "palette.pals()/hcl.pals(), or a vector of R colours.",
          call. = FALSE
        )
      }
      colours <- rep(palette, n)
    }
  } else {
    colours.valid <- tryCatch({
      grDevices::col2rgb(palette)
      TRUE
    }, error = function(e) FALSE)
    if (!colours.valid) {
      stop("Every manually supplied palette entry must be a valid R colour.",
           call. = FALSE)
    }
    colours <- rep(palette, length.out = n)
  }

  if (isTRUE(reverse)) colours <- rev(colours)
  if (!is.null(alpha)) {
    colours <- grDevices::adjustcolor(colours, alpha.f = alpha)
  }
  unname(colours)
}

.dta_resolve_plot_colours <- function(n, color.palette = "mada6",
                                      legacy.primary = NULL,
                                      legacy.points = NULL,
                                      point.alpha = 0.18) {
  primary <- dta_palette(
    if (is.null(legacy.primary)) color.palette else legacy.primary,
    n = n
  )
  points <- if (is.null(legacy.points)) {
    dta_palette(primary, n = n, alpha = point.alpha)
  } else {
    dta_palette(legacy.points, n = n)
  }
  list(primary = primary, points = points)
}

.dta_population <- function(dat) {
  columns <- c("TP", "TN", "FP", "FN")
  counts <- lapply(dat[columns], function(x) {
    if (!is.numeric(x)) {
      stop("TP, TN, FP, and FN must be numeric to calculate population sizes.",
           call. = FALSE)
    }
    x
  })
  population <- Reduce(`+`, counts)
  if (any(!is.finite(population)) || any(population <= 0)) {
    stop("Every included row must have a positive, finite total population.",
         call. = FALSE)
  }
  as.numeric(population)
}

.dta_population_point_sizes <- function(population,
                                        size.range = c(0.7, 2)) {
  if (!is.numeric(size.range) || length(size.range) != 2L ||
      any(!is.finite(size.range)) || any(size.range <= 0) ||
      size.range[1] > size.range[2]) {
    stop("'point.size.range' must contain two positive increasing numbers.",
         call. = FALSE)
  }
  if (!is.numeric(population) || length(population) < 1L ||
      any(!is.finite(population)) || any(population <= 0)) {
    stop("Population values must be positive and finite.", call. = FALSE)
  }

  # Base graphics cex controls marker diameter. Using sqrt(n / max(n)) makes
  # marker area proportional to population, with a legibility floor.
  sizes <- size.range[2] * sqrt(population / max(population))
  pmax(size.range[1], sizes)
}

.dta_use_population_weights <- function(meta.object, population) {
  if (length(population) != length(meta.object$TE)) {
    stop("Population weights do not align with the forest-plot studies.",
         call. = FALSE)
  }
  weights <- population / sum(population)
  meta.object$w.common <- weights
  meta.object$w.fixed <- weights
  meta.object$w.random <- weights

  subgroup <- meta.object$subgroup
  if (!is.null(subgroup) && length(subgroup) == length(weights)) {
    subgroup.weights <- tapply(weights, subgroup, sum)
    meta.object$w.common.w <- subgroup.weights
    meta.object$w.fixed.w <- subgroup.weights
    meta.object$w.random.w <- subgroup.weights
  }

  meta.object$population <- population
  meta.object$weight.study.source <- "population"
  meta.object
}
