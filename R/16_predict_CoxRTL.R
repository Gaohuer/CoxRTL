# predict_coxrtl Function Documentation ----------------------------------------

#' Predict Risk Scores from a Fitted Cox-Type Survival Model
#'
#' Computes the risk scores of new predictor rows for a single model fitted by
#' \code{\link{fit_coxrtl}}: the linear predictor \eqn{X\hat\beta} and its
#' exponential (relative risk). No metric is computed here.
#'
#' @param coef A numeric coefficient vector whose names match covariate columns
#'   of \code{X}, or an object of class \code{"fit_coxrtl"} as returned by
#'   \code{\link{fit_coxrtl}()} (either shape: the bare coefficient vector of a
#'   point fit, or the \code{list(beta, ci)} of a \code{bootstrap = TRUE} fit).
#'   In the latter case the \code{beta} component is used.
#' @param X A numeric matrix or data.frame of predictor values for the subjects
#'   to predict (typically the deployment / evaluation covariates). Rows are
#'   subjects; the columns whose names match \code{names(coef)} are used, in the
#'   order of \code{coef}. Extra columns of \code{X} are ignored.
#'
#' @return A list with two numeric vectors of length \code{nrow(X)}:
#' \item{\code{pred}}{The linear predictor \eqn{\eta = X\hat\beta} (the risk
#' score: higher values mean higher hazard / shorter expected survival).}
#' \item{\code{risk}}{The relative risk \eqn{\exp(\eta)}.}
#' When \code{X} has row names, both vectors inherit them.
#'
#' @details
#' A penalized Cox model fitted with \code{cv.glmnet()} returns only the
#' coefficient vector: the baseline hazard is absorbed nonparametrically and is
#' not part of the fitted object. Because \code{predict_coxrtl} receives only
#' predictors and coefficients (no training \code{time}/\code{status}), it
#' cannot estimate the cumulative baseline hazard and therefore returns the
#' linear-predictor risk score rather than survival probabilities.
#'
#' @examples
#' # from a plain named coefficient vector (a small non-survival sanity check)
#' X <- matrix(rnorm(12), nrow = 4,
#'             dimnames = list(paste0("s", 1:4), paste0("x", 1:3)))
#' coef_vec <- setNames(c(0.5, -0.2, 0), colnames(X))
#' predict_coxrtl(coef = coef_vec, X = X)
#'
#' \dontrun{
#' # from a fit_coxrtl object (either a point fit or a bootstrap fit)
#' fit <- fit_coxrtl(X_t = X_target, y_t = y_target,
#'                   transfer = FALSE, alpha = 1, seed = 123)
#' predict_coxrtl(coef = fit, X = X_evaldata)
#' }
#'
#' @export
# --------------------------------------------------------------------------------

predict_coxrtl <- function(coef, X){

  # ------------------------------------------------------------
  # Step 1. Resolve the coefficient vector (vector or fit_coxrtl object)
  # ------------------------------------------------------------

  beta <- .coef_fit(coef)

  if (!is.numeric(beta) || length(beta) == 0L || anyNA(beta)) {
    stop("'coef' must be a numeric coefficient vector without missing values, ",
         "or a fit_coxrtl object as returned by fit_coxrtl().")
  }
  if (is.null(names(beta)) || any(!nzchar(names(beta)))) {
    stop("'coef' must be named, with names matching covariate columns of 'X'.")
  }

  # ------------------------------------------------------------
  # Step 2. Align the predictor matrix with the covariates of coef
  # ------------------------------------------------------------

  Xm <- as.matrix(X)
  if (!is.numeric(Xm)) {
    stop("'X' must be numeric (a matrix or a data.frame of numeric columns).")
  }
  if (is.null(colnames(Xm)) || any(!nzchar(colnames(Xm)))) {
    stop("'X' must have column names matching the covariates in 'coef'.")
  }

  missing_x <- setdiff(names(beta), colnames(Xm))
  if (length(missing_x) > 0) {
    stop("Covariate(s) in 'coef' missing from 'X': ",
         paste(missing_x, collapse = ", "))
  }

  Xm <- Xm[, names(beta), drop = FALSE]

  # ------------------------------------------------------------
  # Step 3. Compute the linear predictor and the relative risk
  # ------------------------------------------------------------

  pred <- drop(Xm %*% beta)
  if (!is.null(rownames(Xm))) names(pred) <- rownames(Xm)

  list(pred = pred, risk = exp(pred))
}
