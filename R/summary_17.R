# summary / print methods for detect_transferable objects -------------------------
#
# detect_transferable() (#17) returns an S3 object of class
# "detect_transferable": a list holding the per-dataset target-domain gains
# (delta_C, delta_C_median), the per-fold gains (fold_delta), the number of
# evaluable folds (n_valid), the decision under the selected rule (transferable
# and selected) and the settings that produced it (V, threshold, rule). The
# helpers below only re-express that decision; they never refit or re-screen.
#
# The decision implemented by detect_transferable() is a strict comparison
# against 'threshold' for the gain criteria and a non-strict comparison for the
# remaining "robust" criteria:
#   rule = "mean"   : delta_C > threshold
#   rule = "median" : delta_C_median > threshold
#   rule = "robust" : the two above, plus prop_positive >= min_prop_positive
#                     and n_valid >= V - 1
# Each dataset's verdict is carried by the 'transferable' column of the
# per-dataset table, so the criterion is stated once in the header rather than
# repeated per dataset. Note that 'min_prop_positive' is not part of the
# returned object, so the criterion can be stated only symbolically.

# --- small read-only helpers (not exported) --------------------------------------

# These helpers format scalars with .fmt_scalar() from R/summary_5.R, the
# trimmed variant of .fmt_val() from R/summary_13.R.

# The condition(s) a dataset has to satisfy under the fitted rule, as text.
.dt_rule_text <- function(rule) {
  switch(rule,
         mean   = "delta_C > threshold",
         median = "delta_C_median > threshold",
         robust = paste0("delta_C > threshold and delta_C_median > threshold, ",
                         "prop_positive >= min_prop_positive, ",
                         "n_valid >= V - 1"),
         "the selected rule")
}

# --------------------------------------------------------------------------------

#' Summary Method for \code{detect_transferable} Objects
#'
#' Produces a compact summary of the transferability screen returned by
#' \code{\link{detect_transferable}}: the decision rule and threshold used and,
#' for every external dataset, the evidence it was judged on (mean and median
#' fold gain, proportion of positive folds, number of evaluable folds) together
#' with the resulting verdict.
#'
#' The verdict sits in the \code{transferable} column of the per-dataset table,
#' next to the numbers it was derived from, rather than in a separate listing.
#' It is a re-expression of the \code{transferable} component of the screened
#' object; nothing is re-fitted. Under \code{rule = "mean"} a dataset is
#' transferable when \code{delta_C > threshold}, under \code{rule = "median"}
#' when \code{delta_C_median > threshold} (the mean is ignored), and under
#' \code{rule = "robust"} when both hold and, in addition,
#' \code{prop_positive >= min_prop_positive} and \code{n_valid >= V - 1}. A
#' dataset with no evaluable fold has \code{NA} gain summaries and is never
#' transferable.
#'
#' @param object An object of class \code{"detect_transferable"} as returned by
#'   \code{\link{detect_transferable}()}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return An object of class \code{"summary.detect_transferable"} (a list) with
#'   the components \code{sources}, \code{n_datasets}, \code{transferable},
#'   \code{selected}, \code{n_transferable}, \code{V}, \code{threshold},
#'   \code{rule} and \code{table} (a one-row-per-dataset data.frame whose
#'   \code{transferable} column holds that dataset's verdict).
#'
#' @examples
#' \dontrun{
#' screen <- detect_transferable(X_t = X_t, y_t = y_t,
#'                               X_s_list = X_s_list, y_s_list = y_s_list,
#'                               rule = "robust", V = 5, seed = 123)
#' summary(screen)$selected     # positions of the transferable datasets
#' summary(screen)$table        # gains and verdicts, one row per dataset
#' print(summary(screen))       # the same, formatted
#' }
#'
#' @seealso \code{\link{detect_transferable}},
#'   \code{\link{print.detect_transferable}},
#'   \code{\link{print.summary.detect_transferable}}
#'
#' @method summary detect_transferable
#' @export
# --------------------------------------------------------------------------------

summary.detect_transferable <- function(object, ...) {
  if (!inherits(object, "detect_transferable")) {
    stop("'object' must be an object of class \"detect_transferable\", ",
         "as returned by detect_transferable().")
  }

  table <- data.frame(
    source        = object$sources,
    delta_C       = as.numeric(object$delta_C),
    delta_C_median = as.numeric(object$delta_C_median),
    prop_positive = as.numeric(object$prop_positive),
    n_valid       = as.integer(object$n_valid),
    transferable  = as.logical(object$transferable),
    check.names   = FALSE,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      sources         = object$sources,
      n_datasets      = length(object$sources),
      transferable    = object$transferable,
      selected        = object$selected,
      n_transferable  = length(which(object$transferable)),
      V               = object$V,
      threshold       = object$threshold,
      rule            = object$rule,
      table           = table
    ),
    class = "summary.detect_transferable"
  )
}

#' Print Method for \code{summary.detect_transferable} Objects
#'
#' Prints the screen summary produced by
#' \code{\link{summary.detect_transferable}}: the rule, the threshold and the
#' fold count, the per-dataset table of gain summaries with a
#' \code{transferable} column holding each dataset's verdict, and the resulting
#' set of transferable datasets.
#'
#' @param x An object of class \code{"summary.detect_transferable"} as returned
#'   by \code{\link{summary.detect_transferable}()}.
#' @param digits Integer. Number of significant digits to print.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{detect_transferable}},
#'   \code{\link{summary.detect_transferable}}
#'
#' @method print summary.detect_transferable
#' @export
# --------------------------------------------------------------------------------

print.summary.detect_transferable <- function(x,
                                              digits = max(3L, getOption("digits") - 3L),
                                              ...) {
  cat("Transferability screening (detect_transferable)\n")
  cat("  rule      : ", x$rule, "\n", sep = "")
  cat("  criterion : ", .dt_rule_text(x$rule), "\n", sep = "")
  cat("  threshold : ", .fmt_scalar(x$threshold, digits), "\n", sep = "")
  cat("  folds (V) : ", x$V, "\n", sep = "")

  cat("\n  per-dataset gain:\n")
  tab <- x$table
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], round, digits = digits)
  print.data.frame(tab, row.names = FALSE, print.gap = 1)

  sel <- x$sources[which(x$transferable)]
  cat("\n  Transferable: ")
  if (length(sel) > 0L) {
    cat(paste(sel, collapse = ", "), "  (", x$n_transferable, " of ",
        x$n_datasets, ")\n", sep = "")
  } else {
    cat("none  (0 of ", x$n_datasets, ")\n", sep = "")
  }

  invisible(x)
}

#' Print Method for \code{detect_transferable} Objects
#'
#' Prints a compact overview of the screen produced by
#' \code{\link{detect_transferable}}: the rule, threshold and fold count, and
#' the set of external datasets declared transferable. Use
#' \code{\link{summary.detect_transferable}} for the per-dataset verdicts and
#' the comparison with the threshold.
#'
#' @param x An object of class \code{"detect_transferable"} as returned by
#'   \code{\link{detect_transferable}()}.
#' @param digits Integer. Number of significant digits to print.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @seealso \code{\link{detect_transferable}},
#'   \code{\link{summary.detect_transferable}}
#'
#' @method print detect_transferable
#' @export
# --------------------------------------------------------------------------------

print.detect_transferable <- function(x,
                                      digits = max(3L, getOption("digits") - 3L),
                                      ...) {
  s <- summary.detect_transferable(x)

  cat("Transferability screening (detect_transferable)\n")
  cat("  rule      : ", s$rule, "\n", sep = "")
  cat("  threshold : ", .fmt_scalar(s$threshold, digits), "\n", sep = "")
  cat("  folds (V) : ", s$V, "\n", sep = "")

  sel <- s$sources[which(s$transferable)]
  cat("  transferable: ")
  if (length(sel) > 0L) {
    cat(paste(sel, collapse = ", "), "  (", s$n_transferable, " of ",
        s$n_datasets, ")\n", sep = "")
  } else {
    cat("none  (0 of ", s$n_datasets, ")\n", sep = "")
  }
  cat("Use summary() for the per-dataset verdicts.\n")

  invisible(x)
}
