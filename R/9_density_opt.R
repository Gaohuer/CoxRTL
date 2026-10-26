# density_opt Function Documentation --------------------------------------------

#' Construct the Penalized Density-Ratio Objective Function
#'
#' Returns an objective function for density-ratio estimation under covariate shift.
#' The objective consists of an exponential-tilting loss plus a user-specified
#' penalty term. Supported penalty types include lasso (L1), ridge (L2), and
#' elastic net penalties.
#'
#' This function is used within \code{calculate_weights()} to estimate the
#' density-ratio parameter vector \code{theta}.
#'
#' @param Xs A numeric matrix containing the augmented covariates of the
#'   denominator dataset.
#' @param Xt A numeric matrix containing the augmented covariates of the
#'   numerator dataset.
#' @param lam Numeric. Regularization parameter controlling the strength of
#'   the penalty term.
#' @param penalty Character string specifying the penalty type. Supported values
#'   are \code{"lasso"}, \code{"ridge"}, and \code{"elasticnet"}. Default is
#'   \code{"elasticnet"}.
#' @param alpha Numeric. Elastic net mixing parameter used only when
#'   \code{penalty = "elasticnet"}. Must lie in \eqn{[0, 1]}. Here,
#'   \code{alpha = 1} corresponds to lasso and \code{alpha = 0} corresponds
#'   to ridge. Default is 0.5.
#'
#' @return A function of \code{theta_0} that evaluates the penalized density-ratio
#'   objective.
#'
#' @details
#' The returned objective function is of the form
#' \deqn{
#' \frac{1}{n_s}\sum_{i=1}^{n_s}\exp(X_{s,i}^T\theta)
#' - \frac{1}{n_t}\sum_{i=1}^{n_t} X_{t,i}^T\theta
#' + \lambda \cdot \mathrm{Penalty}(\theta),
#' }
#' where the penalty is chosen according to \code{penalty}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' Xs <- cbind(1, matrix(rnorm(50), ncol = 2))   # augmented denominator covariates
#' Xt <- cbind(1, matrix(rnorm(50), ncol = 2))   # augmented numerator covariates
#'
#' fn <- density_opt(Xs, Xt, lam = 0.1, penalty = "ridge")
#' fn(rep(0, 3))
#' }
#'
#' @keywords internal
# --------------------------------------------------------------------------------

density_opt <- function(Xs, Xt, lam, penalty = "elasticnet", alpha = 0.5){
  
  # ------------------------------------------------------------
  # Step 1. Validate penalty arguments
  # ------------------------------------------------------------
  
  penalty <- match.arg(penalty, c("lasso", "ridge", "elasticnet"))
  
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("Argument 'alpha' must be a single numeric value between 0 and 1.")
  }
  
  # ------------------------------------------------------------
  # Step 2. Construct the penalized objective function
  # ------------------------------------------------------------
  
  mm <- function(theta_0){
    theta_0 <- as.matrix(theta_0, ncol = 1)
    
    base_loss <- mean(exp(Xs %*% theta_0)) - mean(Xt %*% theta_0)
    
    penalty_term <- switch(
      penalty,
      lasso = sum(abs(theta_0)),
      ridge = sum(theta_0^2),
      elasticnet = (1 - alpha) * 0.5 * sum(theta_0^2) + alpha * sum(abs(theta_0))
    )
    
    ff <- base_loss + lam * penalty_term
    return(ff)
  }
  
  return(mm)
}