# summary / print methods for fit_coxrtl objects ----------------------------------
#
# fit_coxrtl() (#13) fits one model per call and always returns an S3 object of
# class "fit_coxrtl": a list(beta, ci) holding the point estimate (a named
# numeric coefficient vector) in $beta and, when bootstrap intervals were
# requested, the interval summary list(se, lower, upper) in $ci (ci = NULL
# otherwise). The inference settings (bootstrap, conf.level, nboot) and the
# light-weight "fit_info" attribute (model name, data sizes, selected lambdas,
# deviance ratios) are attributes of the list. The helpers below read the
# coefficients, the (optional) intervals and fit_info; they never refit.

# --- small read-only helpers (not exported) --------------------------------------

# Does the object carry confidence intervals (a list(beta, ci) with ci present)?
.has_ci_fit <- function(x) {
  inherits(x, "fit_coxrtl") && is.list(x) && !is.null(x$ci)
}

# The coefficient vector of a fit: for a fit_coxrtl object (a list with a $beta
# component) its $beta, for a bare coefficient vector the vector itself. This
# lets every consumer (metrics, evaluation, prediction) accept a fit object or
# a plain vector interchangeably.
.coef_fit <- function(x) {
  if (is.list(x) && !is.null(x$beta)) x$beta else x
}

# The per-coefficient interval summary (list(se, lower, upper)), or NULL.
.ci_fit <- function(x) {
  if (.has_ci_fit(x)) x$ci else NULL
}

# The "fit_info" attribute attached to a fit_coxrtl object by fit_coxrtl(), or
# NULL.
.info_fit <- function(x) attr(x, "fit_info")

# Format a scalar/vector for the console with a fixed number of digits.
.fmt_val <- function(v, digits) {
  if (length(v) == 0L || all(is.na(v))) return("NA")
  paste(formatC(v, digits = digits, format = "g"), collapse = ", ")
}

# Resolve the 'head' argument of the summary and print methods: NULL (the
# default) means "show every row" and is passed through unchanged, otherwise a
# single positive whole number of rows is required. Fractional values are
# truncated downwards rather than rejected, so head = 10.7 is read as the same
# request as head = 10.
.resolve_head <- function(head) {
  if (is.null(head)) return(NULL)
  if (!is.numeric(head) || length(head) != 1L || is.na(head) ||
      !is.finite(head) || head < 1) {
    stop("'head' must be NULL (show all coefficients) or a single positive ",
         "whole number of coefficients to show.")
  }
  as.integer(head)
}

# Print the shared header lines (model, data sizes, tuning, inference) used by
# both the print and the summary methods.
.print_header <- function(model, weighted, p, n_target, n_external, n_events_t,
                          n_events_s, stats, alpha, alpha1, alpha2, bootstrap,
                          conf.level, nboot, digits) {
  # Both models use weights, so every fit reports its weighting state; they
  # differ in which weights they consume (w_t alone vs. w_t and w_s), and in
  # whether they report a single lambda or the two steps of the transfer model.
  transfer <- identical(model, "transfer")
  cat("CoxRTL fit: model '", model, "'", sep = "")
  if (!is.null(model)) {
    cat(if (isTRUE(weighted)) " (reweighted)" else " (unweighted)")
  }
  cat("\n", sep = "")
  cat("  p = ", p, " covariates", sep = "")
  if (!is.null(n_target)) {
    cat(" | target: n = ", n_target, if (!is.null(n_events_t)) {
      paste0(" (", n_events_t, " events)")
    }, sep = "")
  }
  # The nontransfer model uses the target data only, so it has no external arm
  # to report.
  if (!is.null(n_external) && n_external > 0) {
    cat(" | external: n = ", n_external, if (!is.null(n_events_s)) {
      paste0(" (", n_events_s, " events)")
    }, sep = "")
  }
  cat("\n")

  if (!transfer) {
    cat("  alpha = ", .fmt_val(alpha, digits),
        ", lambda = ", .fmt_val(stats$lambda, digits),
        "  (nzero = ", stats$nzero, ", dev.ratio = ",
        .fmt_val(stats$dev.ratio, digits), ")\n", sep = "")
  } else {
    cat("  alpha1 = ", .fmt_val(alpha1, digits), ", alpha2 = ",
        .fmt_val(alpha2, digits), "\n", sep = "")
    # One line per step, mirroring the single-step format above, so that the
    # two lambdas and their statistics stay readable side by side.
    cat("  lambda1 = ", .fmt_val(stats$lambda1, digits), "  (nzero1 = ",
        stats$nzero1, ", dev.ratio1 = ", .fmt_val(stats$dev.ratio1, digits),
        ")\n", sep = "")
    cat("  lambda2 = ", .fmt_val(stats$lambda2, digits), "  (nzero2 = ",
        stats$nzero2, ", dev.ratio2 = ", .fmt_val(stats$dev.ratio2, digits),
        ")\n", sep = "")
  }

  if (isTRUE(bootstrap)) {
    cat("  confidence intervals: bootstrap",
        ", conf.level = ", .fmt_val(conf.level, digits),
        if (!is.null(nboot)) paste0(", nboot = ", nboot) else "",
        "\n", sep = "")
  }
  invisible()
}

# --------------------------------------------------------------------------------

#' Summary Method for \code{fit_coxrtl} Objects
#'
#' Produces a compact numerical summary of a single model fitted by
#' \code{\link{fit_coxrtl}}: the data sizes, the selected lambda(s) and
#' deviance ratios and, when confidence intervals were requested, the
#' coefficient-level estimates with their standard errors and interval bounds.
#'
#' The coefficient table holds the active (non-zero) covariates ordered by
#' decreasing absolute estimate, so that the most influential covariates come
#' first and a shortened table is always the most influential part of the
#' model. Ties are kept in their original covariate order.
#'
#' @param object An object of class \code{"fit_coxrtl"} as returned by
#'   \code{\link{fit_coxrtl}()}.
#' @param head Integer or \code{NULL}. Maximum number of rows of the
#'   coefficient table to keep, so that \code{head = 10} shows only the ten
#'   coefficients with the largest absolute estimate. The default \code{NULL}
#'   keeps every active coefficient.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{"summary.fit_coxrtl"} (a list) holding the
#'   pieces needed by \code{\link[base]{print}}: the model name, the data and
#'   model-overview information, the number of active coefficients
#'   (\code{active}, which \code{head} does not affect), the number of rows kept
#'   (\code{shown}) and the coefficient table of the active (non-zero)
#'   covariates, ordered by decreasing absolute estimate, with their estimates
#'   and, when available, standard errors and confidence-interval bounds.
#'
#' @examples
#' \dontrun{
#' fit <- fit_coxrtl(X_t = X_t, X_s = X_s, y_t = y_t, y_s = y_s,
#'                   transfer = TRUE, alpha = 1, seed = 123)
#' summary(fit)
#' summary(fit, head = 10)   # only the ten most influential coefficients
#' }
#'
#' @seealso \code{\link{fit_coxrtl}}, \code{\link{print.fit_coxrtl}},
#'   \code{\link{print.summary.fit_coxrtl}}
#'
#' @method summary fit_coxrtl
#' @export
# --------------------------------------------------------------------------------

summary.fit_coxrtl <- function(object, head = NULL, ...) {
  if (!inherits(object, "fit_coxrtl")) {
    stop("'object' must be an object of class \"fit_coxrtl\", ",
         "as returned by fit_coxrtl().")
  }

  head <- .resolve_head(head)

  beta <- .coef_fit(object)
  ci   <- .ci_fit(object)
  info <- .info_fit(object)

  p <- if (!is.null(info) && !is.null(info$p)) info$p else length(beta)
  active <- sum(beta != 0, na.rm = TRUE)

  # Rank by decreasing absolute estimate so that truncating the table drops the
  # least influential covariates; which() already returns increasing positions,
  # so the index as a second key keeps ties in the original covariate order.
  keep <- which(beta != 0)
  keep <- keep[order(-abs(beta[keep]), seq_along(keep))]
  if (!is.null(head) && length(keep) > head) {
    keep <- keep[seq_len(head)]
  }

  if (length(keep) > 0L) {
    tab <- data.frame(
      variable = names(beta)[keep],
      estimate = as.numeric(beta[keep]),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!is.null(ci)) {
      tab$std.error <- ci$se[keep]
      tab$lower     <- ci$lower[keep]
      tab$upper     <- ci$upper[keep]
    }
  } else {
    tab <- NULL
  }

  structure(
    list(
      model       = if (!is.null(info)) info$model else NULL,
      weighted    = if (!is.null(info)) info$weighted else NULL,
      p           = p,
      active      = active,
      shown       = length(keep),
      n_target    = if (!is.null(info)) info$n_target else NULL,
      n_external  = if (!is.null(info)) info$n_external else NULL,
      n_events_t  = if (!is.null(info)) info$n_events_t else NULL,
      n_events_s  = if (!is.null(info)) info$n_events_s else NULL,
      stats       = if (!is.null(info)) info$stats else NULL,
      alpha       = if (!is.null(info)) info$alpha else NULL,
      alpha1      = if (!is.null(info)) info$alpha1 else NULL,
      alpha2      = if (!is.null(info)) info$alpha2 else NULL,
      bootstrap   = isTRUE(attr(object, "bootstrap")),
      conf.level  = attr(object, "conf.level"),
      nboot       = attr(object, "nboot"),
      coefficient = tab
    ),
    class = "summary.fit_coxrtl"
  )
}

#' Print Method for \code{summary.fit_coxrtl} Objects
#'
#' Prints the summary produced by \code{\link{summary.fit_coxrtl}()}: the data
#' sizes and model overview followed by the coefficient table of the active
#' covariates, ordered by decreasing absolute estimate, with standard errors and
#' confidence-interval bounds when \code{\link{fit_coxrtl}} was called with
#' \code{bootstrap = TRUE}. Whenever the table is shorter than the number of
#' active coefficients, a line noting how many of the top coefficients are shown
#' is printed between the two.
#'
#' @param x An object of class \code{"summary.fit_coxrtl"} as returned by
#'   \code{\link{summary.fit_coxrtl}()}.
#' @param digits Integer. Number of significant digits to print.
#' @param head Integer or \code{NULL}. Maximum number of rows of the
#'   coefficient table to print, applied on top of any \code{head} already
#'   passed to \code{\link{summary.fit_coxrtl}}. The default \code{NULL} prints
#'   every row the summary holds.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{fit_coxrtl}}, \code{\link{summary.fit_coxrtl}}
#'
#' @method print summary.fit_coxrtl
#' @export
# --------------------------------------------------------------------------------

print.summary.fit_coxrtl <- function(x,
                                     digits = max(3L, getOption("digits") - 3L),
                                     head = NULL,
                                     ...) {
  .print_header(x$model, x$weighted, x$p, x$n_target, x$n_external,
                x$n_events_t, x$n_events_s, x$stats, x$alpha, x$alpha1,
                x$alpha2, x$bootstrap, x$conf.level, x$nboot, digits)

  tab <- x$coefficient
  head <- .resolve_head(head)
  if (!is.null(head) && !is.null(tab) && nrow(tab) > head) {
    tab <- tab[seq_len(head), , drop = FALSE]
  }

  # The summary may have dropped rows itself, so the number actually printed is
  # read off the table rather than from the requested 'head'.
  shown <- if (is.null(tab)) 0L else nrow(tab)
  if (shown < x$active) {
    cat("  showing the top ", shown, " of ", x$active,
        " coefficients by |estimate|\n", sep = "")
  }

  if (is.null(tab)) {
    cat("  (all coefficients are zero)\n")
  } else {
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], round, digits = digits)
    print.data.frame(tab, row.names = FALSE, print.gap = 1)
  }

  invisible(x)
}

#' Print Method for \code{fit_coxrtl} Objects
#'
#' Prints a compact overview of a single model fitted by
#' \code{\link{fit_coxrtl}}: the data sizes, the selected lambda(s) and the
#' deviance ratios, including the number of non-zero coefficients selected at
#' the chosen lambda(s). Use \code{\link{summary.fit_coxrtl}} for the
#' coefficient-level table.
#'
#' @param x An object of class \code{"fit_coxrtl"} as returned by
#'   \code{\link{fit_coxrtl}()}.
#' @param digits Integer. Number of significant digits to print.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{fit_coxrtl}}, \code{\link{summary.fit_coxrtl}}
#'
#' @method print fit_coxrtl
#' @export
# --------------------------------------------------------------------------------

print.fit_coxrtl <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  s <- summary.fit_coxrtl(x)
  .print_header(s$model, s$weighted, s$p, s$n_target, s$n_external,
                s$n_events_t, s$n_events_s, s$stats, s$alpha, s$alpha1,
                s$alpha2, s$bootstrap, s$conf.level, s$nboot, digits)
  cat("Use summary() for the coefficient-level table.\n")
  invisible(x)
}
