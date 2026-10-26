# summary / print methods for calculate_weights objects --------------------------
#
# calculate_weights() (#5) returns an S3 object of class "calculate_weights":
# a list(lambda, w, theta, clip) in which $w is the n x 1 matrix of
# density-ratio weights for the denominator dataset. The weights have already
# been truncated at the percentiles given by $clip (by default the 2.5th and
# 97.5th) and normalized to mean 1 (so they sum to n), which is why the summary
# below reports both their dispersion and the effective sample size they
# correspond to. The helpers read $w, $lambda and $clip; nothing is recomputed.

# --- small read-only helpers (not exported) --------------------------------------

# The weight vector of a calculate_weights object. $w is stored as an n x 1
# matrix, so it is dropped to a plain numeric vector here.
.w_vec <- function(x) as.numeric(x$w)

# Effective sample size under the Kish (1965) definition: the size of an equally
# weighted sample carrying the same information,
#   ESS = (sum(w))^2 / sum(w^2).
# For mean-1 weights this equals n / (1 + CV^2), so 0 < ESS <= n, with equality
# exactly when all weights are equal.
.ess_w <- function(w) sum(w)^2 / sum(w^2)

# Like .fmt_val() (summary_13.R), but for a scalar embedded in a sentence:
# formatC(..., format = "g") right-pads short values, so a threshold of 0 comes
# out as "    0" and reads badly next to the numbers it is compared with. The
# padding is trimmed away here.
.fmt_scalar <- function(v, digits) trimws(.fmt_val(v, digits))

# The truncation interval of a calculate_weights object as a percentile label,
# e.g. "10th/90th" for clip = c(0.1, 0.9) or "2.5th/97.5th" for
# clip = c(0.025, 0.975); NULL when the weights were not truncated. Percentages
# are formatted with up to four significant digits, which keeps 0.1 -> "10" and
# 0.025 -> "2.5" equally tidy. Each bound is formatted on its own because
# formatC() pads the elements of a vector to a common width; the padding is
# trimmed for the same reason as in .fmt_scalar().
.fmt_clip <- function(clip) {
  if (is.null(clip)) {
    return(NULL)
  }
  pct <- vapply(100 * clip,
                function(v) trimws(formatC(v, format = "g", digits = 4)),
                character(1))
  paste0(pct[1], "th/", pct[2], "th percentiles")
}

# The weight diagnostics reported by the summary and print methods.
.w_stats <- function(w) {
  n  <- length(w)
  mn <- mean(w)
  s  <- sd(w)
  c(min = min(w),
    Q1 = unname(quantile(w, 0.25)),
    median = median(w),
    Q3 = unname(quantile(w, 0.75)),
    max = max(w),
    mean = mn,
    sd = s,
    CV = if (mn != 0) s / mn else NA_real_)
}

# --------------------------------------------------------------------------------

#' Summary Method for \code{calculate_weights} Objects
#'
#' Produces a compact summary of the density-ratio weights estimated by
#' \code{\link{calculate_weights}}: the size of the denominator dataset, the
#' regularization parameter used, the distribution of the weights (minimum,
#' quartiles, maximum, mean, standard deviation and coefficient of variation)
#' and the effective sample size the weights correspond to.
#'
#' The effective sample size follows Kish (1965),
#' \eqn{\mathrm{ESS} = (\sum_i w_i)^2 / \sum_i w_i^2}. Because the weights
#' returned by \code{\link{calculate_weights}} are normalized to mean 1, this
#' reduces to \eqn{n / (1 + \sigma^2)}, where \eqn{\sigma^2} is the
#' \emph{population} variance of the weights, \eqn{(1/n)\sum_i (w_i - 1)^2}.
#' Note that the \code{sd} and \code{CV} reported below are the sample
#' versions, so the identity is not exact in those terms. It equals \eqn{n} when
#' all weights are equal and shrinks as they become more dispersed, so it
#' quantifies how much of the denominator sample actually survives the
#' reweighting.
#'
#' @param object An object of class \code{"calculate_weights"} as returned by
#'   \code{\link{calculate_weights}()}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{"summary.calculate_weights"} (a list) with
#'   the components \code{n} (denominator sample size), \code{lambda},
#'   \code{ess}, \code{ess_pct} (\code{ess / n}), \code{weights} (a named
#'   numeric vector of weight diagnostics) and \code{clip} (the truncation
#'   probabilities applied to the weights, or \code{NULL} if there were none).
#'
#' @examples
#' \dontrun{
#' w <- calculate_weights(X_denominator = X_s, X_numerator = X_calibrate,
#'                        cv_lambda = FALSE)
#' summary(w)$ess          # effective sample size
#' print(summary(w))       # the same, formatted
#' }
#'
#' @seealso \code{\link{calculate_weights}}, \code{\link{print.calculate_weights}},
#'   \code{\link{print.summary.calculate_weights}}
#'
#' @method summary calculate_weights
#' @export
# --------------------------------------------------------------------------------

summary.calculate_weights <- function(object, ...) {
  if (!inherits(object, "calculate_weights")) {
    stop("'object' must be an object of class \"calculate_weights\", ",
         "as returned by calculate_weights().")
  }

  w   <- .w_vec(object)
  n   <- length(w)
  ess <- .ess_w(w)

  structure(
    list(
      n       = n,
      lambda  = object$lambda,
      ess     = ess,
      ess_pct = ess / n,
      weights = .w_stats(w),
      clip    = object$clip
    ),
    class = "summary.calculate_weights"
  )
}

#' Print Method for \code{summary.calculate_weights} Objects
#'
#' Prints the weight summary produced by
#' \code{\link{summary.calculate_weights}}: the denominator sample size, the
#' regularization parameter, the distribution of the weights and the effective
#' sample size.
#'
#' @param x An object of class \code{"summary.calculate_weights"} as returned by
#'   \code{\link{summary.calculate_weights}()}.
#' @param digits Integer. Number of significant digits to print.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{calculate_weights}}, \code{\link{summary.calculate_weights}}
#'
#' @method print summary.calculate_weights
#' @export
# --------------------------------------------------------------------------------

print.summary.calculate_weights <- function(x,
                                            digits = max(3L, getOption("digits") - 3L),
                                            ...) {
  cat("Density-ratio weights (calculate_weights)\n")
  cat("  n (denominator): ", x$n, "\n", sep = "")
  cat("  lambda         : ", .fmt_scalar(x$lambda, digits), "\n", sep = "")

  cat("\n  weight distribution:\n")
  tab <- as.data.frame(as.list(round(x$weights, digits)))
  print.data.frame(tab, row.names = FALSE, print.gap = 1)

  cat("\n  effective sample size (Kish): ", .fmt_scalar(x$ess, digits),
      " / ", x$n, "  (", .fmt_scalar(100 * x$ess_pct, digits), "%)\n", sep = "")
  if (is.null(x$clip)) {
    cat("  weights are not truncated, but are normalized to mean 1\n")
  } else {
    cat("  weights are truncated at the ", .fmt_clip(x$clip),
        ", then normalized to mean 1\n", sep = "")
  }

  invisible(x)
}

#' Print Method for \code{calculate_weights} Objects
#'
#' Prints a compact overview of the weights estimated by
#' \code{\link{calculate_weights}}: the denominator sample size, the
#' regularization parameter and the effective sample size. Use
#' \code{\link{summary.calculate_weights}} for the full weight distribution.
#'
#' @param x An object of class \code{"calculate_weights"} as returned by
#'   \code{\link{calculate_weights}()}.
#' @param digits Integer. Number of significant digits to print.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{calculate_weights}}, \code{\link{summary.calculate_weights}}
#'
#' @method print calculate_weights
#' @export
# --------------------------------------------------------------------------------

print.calculate_weights <- function(x,
                                    digits = max(3L, getOption("digits") - 3L),
                                    ...) {
  s <- summary.calculate_weights(x)

  cat("Density-ratio weights (calculate_weights)\n")
  cat("  n (denominator): ", s$n, "\n", sep = "")
  cat("  lambda         : ", .fmt_scalar(s$lambda, digits), "\n", sep = "")
  cat("  effective sample size (Kish): ", .fmt_scalar(s$ess, digits),
      " / ", s$n, "  (", .fmt_scalar(100 * s$ess_pct, digits), "%)\n", sep = "")
  cat("Use summary() for the weight distribution.\n")

  invisible(x)
}
