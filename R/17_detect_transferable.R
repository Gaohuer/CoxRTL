# detect_transferable Function Documentation ------------------------------------

#' Detect Transferable External Datasets by Target-Domain Predictive Gain
#'
#' Screens each candidate external dataset by its cross-validated predictive gain
#' on the target domain and returns the subset of external datasets that are
#' worth borrowing from. This is the transferability-detection step of the
#' CoxRTL framework, used to mitigate negative transfer before any external
#' information is pooled or reweighted.
#'
#' @param X_t A numeric matrix (or data.frame of numeric columns) of predictors
#'   for the target dataset.
#' @param y_t A survival response object for the target dataset, typically
#'   created using \code{\link[survival]{Surv}()}.
#' @param X_s_list A list of numeric matrices, one per candidate external
#'   dataset, each with the same columns as \code{X_t} (in the same order).
#' @param y_s_list A list of survival response objects, one per candidate
#'   external dataset, aligned with \code{X_s_list}.
#' @param alpha Numeric. Elastic net mixing parameter shared by both steps of
#'   the screening fits (\code{alpha = 0} is ridge, \code{alpha = 1} is lasso,
#'   \code{alpha = 0.5} is the elastic net). Default is \code{0.5}.
#' @param lambda Character string or numeric value giving the \code{s} used to
#'   extract coefficients in each fit. Common choices are \code{"lambda.min"}
#'   (default) and \code{"lambda.1se"}. The penalty level itself is selected by
#'   cross-validation within the corresponding fitting data.
#' @param V Integer. Number of folds used to partition the target dataset.
#'   Default is 5.
#' @param threshold Numeric. Gain an external dataset must exceed to be declared
#'   transferable. Default is \code{0}, the criterion used in the original
#'   method. The threshold is meant to be fixed a priori: it is not estimated
#'   from the data, and tuning it on the observed gains would bias the
#'   downstream model. Setting a small positive value (for example \code{0.01})
#'   is a way to ignore differences that are negligible in practice.
#' @param rule Character. Which decision rule to apply, one of \code{"mean"}
#'   (default), \code{"median"} or \code{"robust"}; see \emph{Details}. Partial
#'   matching is allowed. \code{"mean"} reproduces the published criterion and
#'   the previous behaviour of this function.
#' @param min_prop_positive Numeric in \eqn{[0, 1]}. Minimum proportion of
#'   evaluable folds that must show a positive gain. Only used by
#'   \code{rule = "robust"}. Default is \code{0}, so by default the rule imposes
#'   no consistency requirement beyond the mean and the median; set it to a
#'   positive value (for example \code{0.8}) to demand that most folds improve.
#' @param type.measure Character string specifying the cross-validation loss
#'   used in \code{\link[glmnet]{cv.glmnet}}. For Cox models, common choices are
#'   \code{"deviance"} (default) and \code{"C"}.
#' @param seed Integer or NULL. Random seed for the fold partition of the target
#'   dataset and, through it, the whole screening run. Default is NULL.
#' @param ncore Integer. Number of cores forwarded to the underlying
#'   \code{\link{no_transfer_learning}} / \code{\link{transfer_learning}} fits to
#'   parallelize their cross-validation. \code{<= 1} (default) runs sequentially.
#' @param nfolds Integer. Number of cross-validation folds used by the
#'   underlying \code{\link[glmnet]{cv.glmnet}} call. Default is 5.
#' @param standardize Logical. Whether to standardize the predictor columns
#'   before fitting, passed to the underlying \code{cv.glmnet} calls. Default is
#'   \code{FALSE}.
#' @param cox.ties Character. Method for handling tied survival times in the Cox
#'   partial likelihood: \code{"breslow"} (default) or \code{"efron"}.
#'
#' @return An object of class \code{"detect_transferable"}, a named list:
#' \item{\code{delta_C}}{Named numeric vector of the average target-domain gain
#' \eqn{\Delta C(s)} of each external dataset. \code{NA} if no fold could be
#' evaluated for that dataset.}
#' \item{\code{delta_C_median}}{Named numeric vector of the median fold gain, a
#' location summary that is more robust than \code{delta_C} when a single fold
#' is atypical.}
#' \item{\code{prop_positive}}{Named numeric vector with the proportion of
#' evaluable folds whose gain is positive. \code{NA} when no fold could be
#' evaluated. With \code{V = 5} the possible values are 0, 0.2, ..., 1, so the
#' resolution of this criterion is set by \code{V}.}
#' \item{\code{transferable}}{Named logical vector giving the decision under the
#' selected \code{rule}.}
#' \item{\code{selected}}{Named integer vector giving the positions in
#' \code{X_s_list} of the transferable external datasets, i.e. the set
#' \eqn{\mathcal{T}}, named by \code{sources}. Use it directly to subset
#' \code{X_s_list}/\code{y_s_list} (and \code{data_s}).}
#' \item{\code{sources}}{Character vector of dataset labels used to name the
#' results (the names of \code{X_s_list} when available, otherwise
#' \code{"S1"}, \code{"S2"}, ...).}
#' \item{\code{fold_delta}}{A \code{V} by \code{K} matrix of fold-specific gains
#' \eqn{C(\hat\beta^{(\mathrm{TL},s)};\mathcal{D}^{(0)}_\upsilon) -
#' C(\hat\beta^{(0)};\mathcal{D}^{(0)}_\upsilon)}, with \code{NA} for folds that
#' could not be evaluated. Useful for diagnosing how stable the gain is.}
#' \item{\code{n_valid}}{Named integer vector with the number of folds actually
#' contributing to \code{delta_C} for each external dataset.}
#' \item{\code{V}, \code{threshold}, \code{rule}}{The fold count, decision
#' threshold and decision rule used.}
#'
#' Use \code{\link{summary.detect_transferable}} for the per-dataset verdicts and
#' how they relate to \code{threshold}.
#'
#' @details
#' For each candidate external dataset \eqn{\mathcal{D}^{(s)}}, the target
#' dataset \eqn{\mathcal{D}^{(0)}} is randomly partitioned into \code{V}
#' approximately equal folds. For fold \eqn{\upsilon}, the fitting set
#' \eqn{\mathcal{D}^{(0)}_{-\upsilon}} and the held-out validation fold
#' \eqn{\mathcal{D}^{(0)}_{\upsilon}} are formed, and two elastic-net penalized
#' Cox estimators are fitted, each with its penalty level selected by
#' cross-validation inside the corresponding fitting data:
#' \itemize{
#'   \item a \strong{target-only} estimator on \eqn{\mathcal{D}^{(0)}_{-\upsilon}},
#'   computed by \code{\link{no_transfer_learning}};
#'   \item a \strong{two-step transfer} estimator that first pools
#'   \eqn{\mathcal{D}^{(s)}} with \eqn{\mathcal{D}^{(0)}_{-\upsilon}} and then
#'   applies a target-specific debiasing correction using
#'   \eqn{\mathcal{D}^{(0)}_{-\upsilon}} as an offset, computed by
#'   \code{\link{transfer_learning}} with unit weights.
#' }
#' Both estimators are evaluated on the held-out fold \eqn{\mathcal{D}^{(0)}_{\upsilon}}
#' with Harrell's concordance index, and the per-fold gain is their difference.
#' The average gain over folds,
#' \deqn{
#' \Delta C(s) = \frac{1}{V}\sum_{\upsilon=1}^{V}\left[
#'   \mathcal{C}(\hat\beta^{(\mathrm{TL},s)}_{-\upsilon}; \mathcal{D}^{(0)}_{\upsilon}) -
#'   \mathcal{C}(\hat\beta^{(0)}_{-\upsilon}; \mathcal{D}^{(0)}_{\upsilon}) \right],
#' }
#' measures how much external dataset \eqn{s} improves target-domain prediction.
#'
#' Dataset \eqn{s} is declared transferable by comparing a summary of its
#' per-fold gains against \code{threshold}. The \code{rule} argument chooses the
#' summary and any auxiliary requirement:
#' \itemize{
#'   \item \code{"mean"} (default): the average gain exceeds \code{threshold},
#'   giving \eqn{\mathcal{T} = \{s : \Delta C(s) > 0\}} for the default
#'   \code{threshold = 0}. This reproduces the published criterion.
#'   \item \code{"median"}: the median fold gain exceeds \code{threshold}. This
#'   is robust to a single atypical fold, which matters because with a small
#'   \code{V} one unlucky fold can move the average across \code{threshold}.
#'   \item \code{"robust"}: additionally requires the median to exceed
#'   \code{threshold}, all but at most one fold (\code{n_valid >= V - 1}) to be
#'   evaluable, and, when \code{min_prop_positive} is set above \code{0}, at
#'   least that proportion of the evaluable folds to be positive. This guards
#'   against accepting a dataset whose apparent gain rests on one fold. Because
#'   the proportion is computed on at most \code{V} folds, a small \code{V} makes
#'   the requirement coarse: with \code{V = 5}, \code{min_prop_positive = 0.8}
#'   demands at least 4 of 5 positive folds, and with \code{V = 2} it demands
#'   both.
#' }
#' \code{delta_C}, \code{delta_C_median} and \code{prop_positive} are always
#' returned, so the decision can be audited, or recomputed under a different
#' rule, without refitting anything.
#'
#' Because the target-only estimator does not depend on \eqn{s}, it is fitted
#' once per fold and reused across external datasets, so the comparison for every
#' dataset is paired on the same fold partition. Screening uses \emph{unit}
#' weights: density-ratio weighting toward the deployment population is applied
#' downstream (see \code{\link{calculate_weights}} and \code{\link{fit_coxrtl}}),
#' after the transferable set has been determined.
#'
#' A fold is skipped when a fit fails or when Harrell's C-index is not estimable
#' on the held-out fold (for example, no comparable pairs); the number of usable
#' folds is reported in \code{n_valid}. A dataset with no usable fold is treated
#' as non-transferable and triggers a warning.
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' # small simulated example (see the vignette for the full-scale setup)
#' set.seed(123)
#' sim <- gen_surv_data(
#'   N = 5000, Ns = 500, Nt = 100, nt = 500,
#'   kappa_t = 1, kappa_s = 1,
#'   ratio_scale_TeS = 1, ratio_scale_TeTr = 1,
#'   p = 20, num_Xs = 3
#' )
#'
#' data_target <- sim$data_train
#' data_s      <- sim$data_s
#'
#' p <- 20
#' X_t <- as.matrix(data_target[, 1:p])
#' y_t <- Surv(data_target$time, data_target$status)
#'
#' X_s_list <- lapply(data_s, function(d) as.matrix(d[, 1:p]))
#' y_s_list <- lapply(data_s, function(d) Surv(d$time, d$status))
#'
#' # screen the three external datasets against the target domain
#' screen <- detect_transferable(
#'   X_t = X_t, y_t = y_t,
#'   X_s_list = X_s_list, y_s_list = y_s_list,
#'   alpha = 0.5, V = 5, seed = 123
#' )
#' screen$delta_C        # average predictive gain of each external dataset
#' screen$selected       # positions of the transferable external datasets
#'
#' # the per-fold evidence behind the decision
#' round(screen$fold_delta, 4)
#' screen$delta_C_median
#' screen$prop_positive
#'
#' # a stricter screen: the average and the median must both be positive and at
#' # most one fold may be unusable; adding min_prop_positive also requires most
#' # folds to improve
#' screen_strict <- detect_transferable(
#'   X_t = X_t, y_t = y_t,
#'   X_s_list = X_s_list, y_s_list = y_s_list,
#'   alpha = 0.5, V = 5, seed = 123, rule = "robust", min_prop_positive = 0.8
#' )
#' screen_strict$selected
#'
#' # keep only the transferable datasets and refit on the pooled sources
#' if (length(screen$selected) > 0) {
#'   X_s <- do.call(rbind, X_s_list[screen$selected])
#'   time_s <- unlist(lapply(data_s[screen$selected], `[[`, "time"))
#'   status_s <- unlist(lapply(data_s[screen$selected], `[[`, "status"))
#'   fit <- fit_coxrtl(
#'     X_t = X_t, X_s = X_s, y_t = y_t, y_s = Surv(time_s, status_s),
#'     transfer = TRUE, alpha = 0.5, seed = 123
#'   )
#' }
#' }
#'
#' @seealso \code{\link{fit_coxrtl}}, \code{\link{transfer_learning}},
#'   \code{\link{no_transfer_learning}}, \code{\link{calculate_weights}}
#'
#' @importFrom stats median setNames
#'
#' @export
# --------------------------------------------------------------------------------

detect_transferable <- function(
    X_t,
    y_t,
    X_s_list,
    y_s_list,
    alpha = 0.5,
    lambda = "lambda.min",
    V = 5L,
    threshold = 0,
    rule = c("mean", "median", "robust"),
    min_prop_positive = 0,
    type.measure = "deviance",
    seed = NULL,
    ncore = 1,
    nfolds = 5L,
    standardize = FALSE,
    cox.ties = "breslow"
){

  # ------------------------------------------------------------
  # Step 1. Validate the inputs
  # ------------------------------------------------------------

  X_t <- as.matrix(X_t)
  if (!is.numeric(X_t)) {
    stop("'X_t' must be numeric (a matrix or a data.frame of numeric columns).")
  }
  if (!inherits(y_t, "Surv")) {
    stop("'y_t' must be a 'Surv' object, typically created with survival::Surv().")
  }

  if (!is.list(X_s_list) || !is.list(y_s_list)) {
    stop("'X_s_list' and 'y_s_list' must be lists with one element per ",
         "candidate external dataset.")
  }
  K <- length(X_s_list)
  if (K < 1L) {
    stop("'X_s_list' must contain at least one candidate external dataset.")
  }
  if (length(y_s_list) != K) {
    stop("'X_s_list' and 'y_s_list' must have the same length.")
  }

  X_s_list <- lapply(X_s_list, as.matrix)
  if (!all(vapply(X_s_list, is.numeric, logical(1)))) {
    stop("Every element of 'X_s_list' must be numeric.")
  }
  if (!all(vapply(y_s_list, inherits, logical(1), what = "Surv"))) {
    stop("Every element of 'y_s_list' must be a 'Surv' object.")
  }

  # Every external dataset must line up with the target covariates, otherwise
  # the pooled fits below would silently use misaligned columns.
  p <- ncol(X_t)
  if (any(vapply(X_s_list, ncol, integer(1)) != p)) {
    stop("Every external dataset in 'X_s_list' must have the same number of ",
         "columns as 'X_t' (", p, ").")
  }
  cn_t <- colnames(X_t)
  for (k in seq_len(K)) {
    cn_s <- colnames(X_s_list[[k]])
    if (is.null(cn_t) != is.null(cn_s) ||
        (!is.null(cn_t) && !identical(cn_t, cn_s))) {
      stop("Covariate columns of 'X_t' and of external dataset ", k,
           " must match (same names in the same order, or unnamed in both).")
    }
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single number in [0, 1].")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("'threshold' must be a single finite number.")
  }
  rule <- match.arg(rule)
  if (!is.numeric(min_prop_positive) || length(min_prop_positive) != 1L ||
      is.na(min_prop_positive) || min_prop_positive < 0 ||
      min_prop_positive > 1) {
    stop("'min_prop_positive' must be a single number in [0, 1].")
  }

  n <- nrow(X_t)
  if (!is.numeric(V) || length(V) != 1L || is.na(V) || V < 2 || V != round(V)) {
    stop("'V' must be a single integer >= 2.")
  }
  V <- as.integer(V)
  if (V > n) {
    stop("'V' cannot exceed the number of rows in 'X_t' (", n, ").")
  }

  # Labels used to name the results
  sources <- names(X_s_list)
  if (is.null(sources) || any(!nzchar(sources))) {
    sources <- paste0("S", seq_len(K))
  }

  # ------------------------------------------------------------
  # Step 2. Partition the target dataset into V folds
  # ------------------------------------------------------------

  # The partition is drawn once, under the user's seed, and shared by every
  # external dataset so that the per-dataset comparisons are paired. The inner
  # fitters are called with seed = NULL so they do not reset the RNG stream.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  fold_id <- rep(seq_len(V), length.out = n)
  valid_idx <- split(sample.int(n), fold_id)

  # ------------------------------------------------------------
  # Step 3. Screen each external dataset by its target-domain gain
  # ------------------------------------------------------------

  fold_delta <- matrix(
    NA_real_, nrow = V, ncol = K,
    dimnames = list(paste0("fold", seq_len(V)), sources)
  )

  for (v in seq_len(V)) {

    idx_val <- valid_idx[[v]]
    idx_fit <- setdiff(seq_len(n), idx_val)
    X_fit <- X_t[idx_fit, , drop = FALSE]
    y_fit <- y_t[idx_fit]
    X_val <- X_t[idx_val, , drop = FALSE]
    y_val <- y_t[idx_val]

    # Target-only estimator, fitted once per fold and shared across datasets
    beta_0 <- tryCatch(
      no_transfer_learning(
        X = X_fit, y = y_fit,
        alpha = alpha, lambda = lambda,
        seed = NULL, type.measure = type.measure, ncore = ncore,
        nfolds = nfolds, standardize = standardize, cox.ties = cox.ties
      )$beta,
      error = function(e) NULL
    )
    if (is.null(beta_0) || anyNA(beta_0)) {
      next
    }

    c_0 <- harrell_c_index(
      time = y_val[, 1], status = y_val[, 2],
      risk = drop(X_val %*% beta_0)
    )
    if (is.na(c_0)) {
      next
    }

    for (k in seq_len(K)) {

      # Two-step transfer estimator with the candidate external dataset
      beta_tl <- tryCatch(
        transfer_learning(
          X_t = X_fit, X_s = X_s_list[[k]],
          y_t = y_fit, y_s = y_s_list[[k]],
          w_t = NULL, w_s = NULL,
          alpha1 = alpha, alpha2 = alpha,
          lambda1 = lambda, lambda2 = lambda,
          seed = NULL, type.measure = type.measure, ncore = ncore,
          nfolds = nfolds, standardize = standardize, cox.ties = cox.ties
        )$beta,
        error = function(e) NULL
      )
      if (is.null(beta_tl) || anyNA(beta_tl)) {
        next
      }

      c_tl <- harrell_c_index(
        time = y_val[, 1], status = y_val[, 2],
        risk = drop(X_val %*% beta_tl)
      )
      if (is.na(c_tl)) {
        next
      }

      fold_delta[v, k] <- c_tl - c_0
    }
  }

  # ------------------------------------------------------------
  # Step 4. Summarise the per-fold gains and apply the decision rule
  # ------------------------------------------------------------

  # All three summaries are computed for every dataset regardless of the rule,
  # so the decision can be audited and can be revisited without refitting.
  n_valid <- setNames(integer(K), sources)
  delta_C <- setNames(rep(NA_real_, K), sources)
  delta_C_median <- setNames(rep(NA_real_, K), sources)
  prop_positive <- setNames(rep(NA_real_, K), sources)
  for (k in seq_len(K)) {
    d <- fold_delta[, k]
    n_valid[k] <- sum(!is.na(d))
    if (n_valid[k] > 0L) {
      delta_C[k] <- mean(d, na.rm = TRUE)
      delta_C_median[k] <- median(d, na.rm = TRUE)
      prop_positive[k] <- sum(d > 0, na.rm = TRUE) / n_valid[k]
    }
  }

  if (any(n_valid == 0L)) {
    warning("No fold could be evaluated for external dataset(s): ",
            paste(sources[n_valid == 0L], collapse = ", "),
            ". They are treated as non-transferable; consider a larger ",
            "target dataset or fewer folds (smaller 'V').")
  }

  # Each criterion is FALSE (rather than NA) when its summary is unavailable, so
  # a dataset with no usable fold is never declared transferable.
  ok_mean <- !is.na(delta_C) & delta_C > threshold
  transferable <- switch(
    rule,
    mean = ok_mean,
    median = !is.na(delta_C_median) & delta_C_median > threshold,
    robust = ok_mean &
      !is.na(delta_C_median) & delta_C_median > threshold &
      !is.na(prop_positive) & prop_positive >= min_prop_positive &
      n_valid >= V - 1L
  )
  names(transferable) <- sources

  # ------------------------------------------------------------
  # Step 5. Return the screening result
  # ------------------------------------------------------------

  structure(
    list(
      delta_C         = delta_C,
      delta_C_median  = delta_C_median,
      prop_positive   = prop_positive,
      transferable    = transferable,
      selected        = which(transferable),
      sources         = sources,
      fold_delta      = fold_delta,
      n_valid         = n_valid,
      V               = V,
      threshold       = threshold,
      rule            = rule
    ),
    class = "detect_transferable"
  )
}

# --------------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------------

#' Harrell's Concordance Index for a Linear Predictor (Internal)
#'
#' Computes Harrell's C-index of a risk score with respect to an observed
#' survival outcome, returning \code{NA} rather than an error when the index is
#' not estimable. Used inside \code{\link{detect_transferable}} to score the
#' held-out target folds.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param status Numeric/integer vector of event indicators (\code{1} = event,
#'   \code{0} = censored).
#' @param risk Numeric vector of risk scores (higher means higher hazard), on the
#'   same length as \code{time}.
#'
#' @return A single numeric value: Harrell's C-index, or \code{NA_real_} when the
#'   index cannot be computed (fewer than two usable observations, or no
#'   comparable pair). A value above 0.5 means the risk score ranks subjects in
#'   the correct order.
#'
#' @details
#' Harrell's C-index is computed directly from its definition as the proportion
#' of comparable pairs that are concordant, with ties in the risk score counting
#' as half. A pair is comparable when the subject with the shorter observed time
#' had an event; pairs with tied times or with a censored shorter time carry no
#' information and are dropped from both the numerator and the denominator.
#'
#' The pair counting is implemented explicitly rather than through
#' \code{\link[survival]{concordance}} because the latter takes its numeric
#' right-hand side to be a predictor on the \emph{survival time} scale (larger
#' values meaning longer survival), so passing a risk score to it returns one
#' minus Harrell's C-index. Counting the pairs directly keeps the orientation
#' unambiguous and independent of that convention.
#'
#' @keywords internal
# --------------------------------------------------------------------------------

harrell_c_index <- function(time, status, risk) {

  ok <- is.finite(time) & is.finite(status) & is.finite(risk)
  time <- time[ok]
  status <- status[ok]
  risk <- risk[ok]
  n <- length(time)
  if (n < 2L) {
    return(NA_real_)
  }

  # Order by follow-up time so that every comparable partner of an event sits
  # after it; ties in time then fall outside the strict ">" comparison below.
  ord <- order(time)
  time <- time[ord]
  risk <- risk[ord]
  status <- status[ord]

  n_concordant <- 0
  n_comparable <- 0
  for (i in which(status == 1L)) {
    if (i == n) {
      next
    }
    later <- (i + 1L):n
    later <- later[time[later] > time[i]]
    if (length(later) == 0L) {
      next
    }
    risk_later <- risk[later]
    n_concordant <- n_concordant +
      sum(risk[i] > risk_later) + 0.5 * sum(risk[i] == risk_later)
    n_comparable <- n_comparable + length(later)
  }

  if (n_comparable == 0L) {
    return(NA_real_)
  }
  n_concordant / n_comparable
}
