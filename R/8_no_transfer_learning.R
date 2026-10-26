# no_transfer_learning Function Documentation -----------------------------------

#' Penalized Cox Regression without Transfer Learning (Internal)
#'
#' Internal workhorse behind \code{\link{fit_coxrtl}} for its
#' \code{transfer = FALSE} model. Not exported; fit this model with
#' \code{\link{fit_coxrtl}}.
#'
#' Fits a penalized Cox regression model using only the target dataset, without
#' borrowing information from any external dataset. This function serves as a
#' target-only baseline model in simulation studies and comparative analyses.
#'
#' @param X A numeric matrix of predictors for the target dataset.
#' @param y A survival response object for the target dataset, typically created
#'   using \code{Surv()}.
#' @param alpha Numeric. Elastic net mixing parameter used in the penalized
#'   Cox regression model.
#' @param lambda Character string or numeric value specifying which lambda to use
#'   when extracting coefficients from the fitted model. Common choices are
#'   \code{"lambda.min"} and \code{"lambda.1se"}.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#' @param type.measure Character string specifying the cross-validation loss
#'   used in \code{cv.glmnet()}. For Cox models, common choices are
#'   \code{"deviance"} and \code{"C"}.
#' @param ncore Integer. Number of cores for the fold-wise parallel
#'   cross-validation in \code{cv.glmnet()}. \code{<= 1} (default) runs
#'   sequentially; \code{> 1} registers a local \code{doParallel} backend for
#'   this call.
#' @param nfolds Integer. Number of cross-validation folds used in
#'   \code{cv.glmnet()}. Default is 5.
#' @param standardize Logical. Whether to standardize the predictor columns
#'   before fitting (passed to \code{cv.glmnet()}). Default is \code{FALSE}.
#' @param cox.ties Character. Method for handling tied survival times in the Cox
#'   partial likelihood: \code{"breslow"} (default) or \code{"efron"}.
#' @param w Optional numeric vector of observation weights of length
#'   \code{nrow(X)}, passed to \code{cv.glmnet()}. These are the density-ratio
#'   weights of \code{\link{calculate_weights}}, which reweight the target data
#'   toward the deployment population. Defaults to \code{NULL}, i.e. unit
#'   weights.
#'
#' @return A named list with the estimated coefficients and the model summary
#' statistics evaluated at the selected \code{lambda} only:
#' \item{\code{beta}}{Numeric vector of regression coefficients at the selected
#' \code{lambda}.}
#' \item{\code{lambda}}{The resolved numeric value of \code{lambda} actually used.}
#' \item{\code{nzero}}{Number of non-zero coefficients at \code{lambda}.}
#' \item{\code{dev.ratio}}{Fraction of deviance explained at \code{lambda}.}
#'
#' @details
#' The model is fitted using \code{cv.glmnet()} with the Cox family. Cross-validation
#' is used to construct the regularization path, and coefficients are extracted at
#' the user-specified \code{lambda}.
#'
#' @examples
#' \dontrun{
#' beta_ntl <- no_transfer_learning(
#'   X = X_target,
#'   y = y_target,
#'   alpha = 1,
#'   lambda = "lambda.1se",
#'   seed = 123
#' )
#' }
#'
#' @keywords internal
# --------------------------------------------------------------------------------

no_transfer_learning <- function(X, y, alpha, lambda, seed = NULL,
                                 type.measure = "deviance", ncore = 1,
                                 nfolds = 5L, standardize = FALSE,
                                 cox.ties = "breslow", w = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }

  if (missing(alpha) || missing(lambda)) {
    stop("Arguments 'alpha' and 'lambda' must be explicitly specified.")
  }

  # Set up a local parallel backend when ncore > 1 is requested
  backend <- start_parallel_backend(ncore)
  if (backend$owned) on.exit(stop_parallel_backend(backend), add = TRUE)

  # ------------------------------------------------------------
  # Step 1. Fit the penalized Cox model on the target dataset
  # ------------------------------------------------------------

  # Weighted and unweighted fits are two separate calls because cv.glmnet()'s
  # 'weights' formal has no default and is tested with missing(): passing
  # weights = NULL explicitly would not be treated as "no weights".
  pen_com <- if (is.null(w)) {
    cv.glmnet(
      x = X,
      y = y,
      standardize = standardize,
      type.measure = type.measure,
      alpha = alpha,
      family = "cox",
      nfolds = nfolds,
      cox.ties = cox.ties,
      parallel = backend$active
    )
  } else {
    cv.glmnet(
      x = X,
      y = y,
      weights = w,
      standardize = standardize,
      type.measure = type.measure,
      alpha = alpha,
      family = "cox",
      nfolds = nfolds,
      cox.ties = cox.ties,
      parallel = backend$active
    )
  }

  # ------------------------------------------------------------
  # Step 2. Extract coefficients at the specified lambda
  # ------------------------------------------------------------
  
  coef_com <- coef(pen_com, s = lambda)
  beta <- coef_com[, 1]
  s1 <- glmnet_lambda_stats(pen_com, lambda)

  # ------------------------------------------------------------
  # Step 3. Return results at the selected lambda
  # ------------------------------------------------------------

  list(
    beta      = beta,
    lambda    = s1$lambda,
    nzero     = s1$nzero,
    dev.ratio = s1$dev.ratio
  )
}

# Internal helper: extract, for a single requested lambda, the glmnet
# statistics evaluated at that lambda (the resolved lambda value, the number
# of non-zero coefficients and the deviance ratio). Named rules such as
# "lambda.min" are resolved to their numeric value first.
glmnet_lambda_stats <- function(fit, s) {
  if (is.character(s)) s <- fit[[s]]
  idx <- which.min(abs(fit$lambda - s))
  dev <- fit$glmnet.fit$dev.ratio[which.min(abs(fit$glmnet.fit$lambda - s))]
  list(
    lambda    = fit$lambda[idx],
    nzero     = fit$nzero[idx],
    dev.ratio = dev
  )
}
