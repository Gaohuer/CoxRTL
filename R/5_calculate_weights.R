# calculate_weights Function Documentation --------------------------------------

#' Estimate Density-Ratio Weights
#'
#' Computes density-ratio weights for covariate shift correction by estimating
#' an exponential-tilting model between a denominator dataset and a numerator
#' dataset. The regularization parameter \code{lambda} can either be selected
#' by cross-validation or fixed at 1, depending on user specification.
#'
#' @param X_denominator A numeric matrix or data.frame containing the covariates
#'   of the denominator dataset, i.e., the dataset to be reweighted.
#' @param X_numerator A numeric matrix or data.frame containing the covariates
#'   of the numerator dataset, i.e., the target distribution for reweighting.
#' @param cv_lambda Logical. If \code{TRUE}, \code{lambda} is selected by cross-validation via
#'   \code{get_theta()}. If \code{FALSE}, \code{lambda} is fixed at 1.
#' @param penalty Character string specifying the penalty type used in the
#'   density-ratio objective passed to \code{density_opt()}. One of
#'   \code{"lasso"}, \code{"ridge"} or \code{"elasticnet"} (the default).
#' @param ncore Integer. Number of cores for the lambda-grid cross-validation
#'   in \code{get_theta()} when \code{cv_lambda = TRUE}. \code{<= 1} (default)
#'   runs sequentially; \code{> 1} runs the grid in parallel.
#' @param clip Numeric or \code{NULL}. Truncation applied to the estimated
#'   weights before they are normalized, given either as the two truncation
#'   probabilities \code{c(lower, upper)} or as a single central coverage
#'   probability. The default \code{c(0.025, 0.975)} clips the weights at their 2.5th
#'   and 97.5th percentiles. A single value \code{p} in \code{(0.5, 1)} is
#'   shorthand for the symmetric interval \code{c((1 - p) / 2,
#'   1 - (1 - p) / 2)}, so \code{clip = 0.95} clips at the 2.5th and 97.5th
#'   percentiles. Use \code{NULL} to skip truncation altogether.
#'
#' @return An object of class \code{"calculate_weights"}, a list containing:
#' \item{\code{lambda}}{The selected or fixed regularization parameter.}
#' \item{\code{w}}{The normalized density-ratio weights for the denominator dataset.}
#' \item{\code{theta}}{The estimated parameter vector in the exponential-tilting model.}
#' \item{\code{clip}}{The truncation probabilities actually applied, as
#' \code{c(lower, upper)}, or \code{NULL} when truncation was skipped.}
#'
#' Use \code{\link{summary.calculate_weights}} for the weight distribution and the
#' effective sample size.
#'
#' @details
#' The function augments both datasets with an intercept term and estimates the
#' parameter vector \code{theta} by minimizing the density-ratio objective using
#' \code{nlm()}. Estimated weights are then truncated at the percentiles given
#' by \code{clip} (by default the 2.5th and 97.5th percentiles) to improve stability and
#' normalized to have mean 1.
#'
#' @examples
#' \dontrun{
#' wt_cv <- calculate_weights(X_denominator = Xs, X_numerator = Xtest, cv_lambda = TRUE)
#' wt_fixed <- calculate_weights(X_denominator = Xs, X_numerator = Xtest, cv_lambda = FALSE)
#'
#' # a wider truncation interval: the central 95% of the weights
#' wt_95 <- calculate_weights(X_denominator = Xs, X_numerator = Xtest,
#'                            cv_lambda = FALSE, clip = 0.95)
#' }
#'
#' @export
# --------------------------------------------------------------------------------

calculate_weights <- function(X_denominator, X_numerator, cv_lambda = FALSE,
                              penalty = "elasticnet", ncore = 1,
                              clip = c(0.025, 0.975)) {

  # ------------------------------------------------------------
  # Step 1. Check user input
  # ------------------------------------------------------------

  # if (missing(cv_lambda)) {
  #   stop("Argument 'cv_lambda' must be explicitly specified as TRUE or FALSE.")
  # }
  if (!is.logical(cv_lambda) || length(cv_lambda) != 1 || is.na(cv_lambda)) {
    stop("Argument 'cv_lambda' must be a single TRUE or FALSE value.")
  }

  clip <- resolve_clip(clip)

  # ------------------------------------------------------------
  # Step 2. Prepare augmented design matrices
  # ------------------------------------------------------------
  
  N_denominator <- nrow(X_denominator)
  N_numerator <- nrow(X_numerator)
  theta_ini <- rep(0, ncol(X_denominator) + 1)
  
  X_denominator <- as.matrix(cbind(rep(1, N_denominator), X_denominator))
  X_numerator <- as.matrix(cbind(rep(1, N_numerator), X_numerator))
  
  # ------------------------------------------------------------
  # Step 3. Determine the regularization parameter lambda
  # ------------------------------------------------------------
  
  if (cv_lambda) {
    theta_result <- get_theta(X_denominator, X_numerator, ncore = ncore)
    lambda_use <- theta_result$lambda.min
  } else {
    lambda_use <- 1
  }
  
  # ------------------------------------------------------------
  # Step 4. Estimate theta using the density-ratio objective
  # ------------------------------------------------------------
  
  w_opt <- nlm(
    density_opt(X_denominator, X_numerator, lambda_use, penalty = penalty),
    theta_ini,
    iterlim = 1000
  )
  theta <- w_opt$estimate
  
  # ------------------------------------------------------------
  # Step 5. Compute, truncate, and normalize weights
  # ------------------------------------------------------------
  
  w <- exp(X_denominator %*% matrix(theta, ncol = 1))
  if (!is.null(clip)) {
    w <- pmin(pmax(w, quantile(w, clip[1])), quantile(w, clip[2]))
  }
  w <- w / mean(w)

  # ------------------------------------------------------------
  # Step 6. Return results
  # ------------------------------------------------------------

  return(structure(
    list(
      lambda = lambda_use,
      w = w,
      theta = theta,
      clip = clip
    ),
    class = "calculate_weights"
  ))
}

# Internal helper: resolve the 'clip' argument to the pair of probabilities
# [lower, upper] at which the weights are truncated, or NULL when no truncation
# was requested. A single value is read as a central coverage probability, which
# is why it has to exceed 0.5: a single value below 0.5 is far more likely to be
# meant as one of the two truncation bounds, so it is rejected rather than
# silently reinterpreted.
resolve_clip <- function(clip) {
  if (is.null(clip)) {
    return(NULL)
  }
  if (!is.numeric(clip) || anyNA(clip) || !length(clip) %in% c(1L, 2L)) {
    stop("'clip' must be NULL, a single central coverage probability, or two ",
         "truncation probabilities c(lower, upper).")
  }

  if (length(clip) == 1L) {
    if (clip <= 0.5 || clip >= 1) {
      stop("A single 'clip' value is a central coverage probability and must ",
           "lie in (0.5, 1); to set the truncation probabilities directly, ",
           "pass two values such as clip = c(0.025, 0.975).")
    }
    clip <- c((1 - clip) / 2, 1 - (1 - clip) / 2)
  }

  if (clip[1] < 0 || clip[2] > 1 || clip[1] >= clip[2]) {
    stop("'clip' must satisfy 0 <= lower < upper <= 1.")
  }

  c(clip[1], clip[2])
}