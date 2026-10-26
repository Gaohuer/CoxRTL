# transfer_learning Function Documentation --------------------------------------

#' Two-Step Transfer Learning for Penalized Cox Regression (Internal)
#'
#' Internal workhorse behind \code{\link{fit_coxrtl}} for its
#' \code{transfer = TRUE} model. Not exported; fit this model
#' with \code{\link{fit_coxrtl}}.
#'
#' Fits a two-step transfer learning Cox model by first estimating a penalized
#' Cox model on the combined target and external datasets, and then applying
#' a target-specific debiasing step using the target dataset only.
#'
#' This function is designed for survival transfer learning settings where
#' information from an external dataset is borrowed to improve estimation in
#' the target dataset, while allowing target-specific correction through an
#' offset-based refinement step.
#'
#' @param X_t A numeric matrix of predictors for the target dataset.
#' @param X_s A numeric matrix of predictors for the external dataset.
#' @param y_t A survival response object for the target dataset, typically
#'   created using \code{Surv()}.
#' @param y_s A survival response object for the external dataset, typically
#'   created using \code{Surv()}.
#' @param w_t Numeric vector of observation weights for the target dataset.
#'   Defaults to unit weights when \code{NULL}.
#' @param w_s Numeric vector of observation weights for the external dataset.
#'   Defaults to unit weights when \code{NULL}.
#' @param alpha1 Numeric. Elastic net mixing parameter used in the first-step
#'   penalized Cox model fitted on the combined target and external datasets.
#' @param alpha2 Numeric. Elastic net mixing parameter used in the second-step
#'   target-specific debiasing model.
#' @param lambda1 Character string indicating the lambda selection rule used to
#'   extract coefficients from the first-step combined model, typically
#'   \code{"lambda.min"} or \code{"lambda.1se"}.
#' @param lambda2 Character string indicating the lambda selection rule used to
#'   extract coefficients from the second-step debiasing model, typically
#'   \code{"lambda.min"} or \code{"lambda.1se"}.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#' @param type.measure Character string specifying the cross-validation loss
#'   used in \code{cv.glmnet()}. For Cox models, common choices are
#'   \code{"deviance"} and \code{"C"}.
#' @param ncore Integer. Number of cores for the fold-wise parallel
#'   cross-validation in \code{cv.glmnet()} for both steps. \code{<= 1}
#'   (default) runs sequentially; \code{> 1} registers a local
#'   \code{doParallel} backend for this call.
#' @param nfolds Integer. Number of cross-validation folds used in
#'   \code{cv.glmnet()} for both steps. Default is 5.
#' @param standardize Logical. Whether to standardize the predictor columns
#'   before fitting (passed to \code{cv.glmnet()}). Default is \code{FALSE}.
#' @param cox.ties Character. Method for handling tied survival times in the Cox
#'   partial likelihood: \code{"breslow"} (default) or \code{"efron"}.
#'
#' @return A named list with the final transfer-learning estimate and the model
#' summary statistics of both steps, evaluated at the selected
#' \code{lambda1} (first/pooled step) and \code{lambda2} (second/debiasing step)
#' only:
#' \item{\code{beta}}{Numeric vector of final regression coefficients
#' (\code{coef_com + delta}).}
#' \item{\code{coef_com}, \code{delta}}{Numeric vectors of the first-step
#' coefficients and of the target-specific debiasing correction.}
#' \item{\code{lambda1}, \code{lambda2}}{Resolved numeric lambda values used in
#' each step.}
#' \item{\code{nzero1}, \code{dev.ratio1}}{Model summary of the first step at
#' \code{lambda1} (same meaning as in \code{\link{no_transfer_learning}}).}
#' \item{\code{nzero2}, \code{dev.ratio2}}{The corresponding quantities of the
#' second step at \code{lambda2}.}
#'
#' @details
#' The procedure consists of two steps:
#' \itemize{
#'   \item Step 1: Fit a weighted penalized Cox model on the pooled target and
#'   external datasets to obtain an initial coefficient estimate.
#'   \item Step 2: Use the linear predictor from Step 1 as an offset and fit
#'   a weighted penalized Cox model on the target dataset only to estimate a
#'   target-specific correction term.
#' }
#' The final coefficient estimate is obtained by summing the initial coefficient
#' vector and the estimated correction term.
#'
#' @examples
#' \dontrun{
#' beta_tl <- transfer_learning(
#'   X_t = X_target,
#'   X_s = X_source,
#'   y_t = y_target,
#'   y_s = y_source,
#'   w_t = w_target,
#'   w_s = w_source,
#'   alpha1 = 1,
#'   alpha2 = 1,
#'   lambda1 = "lambda.1se",
#'   lambda2 = "lambda.1se",
#'   seed = 123
#' )
#' }
#'
#' @keywords internal
# --------------------------------------------------------------------------------

transfer_learning <- function(
    X_t, X_s, y_t, y_s, w_t = NULL, w_s = NULL,
    alpha1, alpha2, lambda1 = "lambda.min", lambda2 = "lambda.min",
    seed = NULL, type.measure = "deviance", ncore = 1,
    nfolds = 5L, standardize = FALSE, cox.ties = "breslow"
){
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(w_t)) w_t<-rep(1,length(y_t))
  if(is.null(w_s)) w_s<-rep(1,length(y_s))

  # Set up a local parallel backend when ncore > 1 is requested
  backend <- start_parallel_backend(ncore)
  if (backend$owned) on.exit(stop_parallel_backend(backend), add = TRUE)

  # ------------------------------------------------------------
  # Step 1. Combine target and external survival outcomes
  # ------------------------------------------------------------
  
  y_all <- c(y_t, y_s)
  
  # ------------------------------------------------------------
  # Step 2. Fit the first-step penalized Cox model on pooled data
  # ------------------------------------------------------------
  
  pen_com <- cv.glmnet(
    x = rbind(X_t, X_s),
    y = y_all,
    standardize = standardize,
    type.measure = type.measure,
    alpha = alpha1,
    family = "cox",
    nfolds = nfolds,
    parallel = backend$active,
    cox.ties = cox.ties,
    weights = c(w_t, w_s)
  )
  
  coef_com <- coef(pen_com, s = lambda1)
  
  # ------------------------------------------------------------
  # Step 3. Construct the target-specific offset
  # ------------------------------------------------------------
  
  offset <- as.matrix(X_t %*% coef_com)
  
  # ------------------------------------------------------------
  # Step 4. Fit the second-step debiasing model on target data
  # ------------------------------------------------------------
  
  pen_d <- cv.glmnet(
    x = X_t,
    y = y_t,
    standardize = standardize,
    type.measure = type.measure,
    alpha = alpha2,
    family = "cox",
    nfolds = nfolds,
    offset = offset,
    parallel = backend$active,
    cox.ties = cox.ties,
    weights = w_t
  )
  
  delta <- as.matrix(coef(pen_d, s = lambda2))
  beta <- (coef_com + delta)[, 1]
  a1 <- glmnet_lambda_stats(pen_com, lambda1)
  a2 <- glmnet_lambda_stats(pen_d, lambda2)

  # ------------------------------------------------------------
  # Step 5. Return results at the selected lambda1 and lambda2
  # ------------------------------------------------------------

  list(
    beta       = beta,
    coef_com   = coef_com[, 1],
    delta      = delta[, 1],
    lambda1    = a1$lambda,
    nzero1     = a1$nzero,
    dev.ratio1 = a1$dev.ratio,
    lambda2    = a2$lambda,
    nzero2     = a2$nzero,
    dev.ratio2 = a2$dev.ratio
  )
}
