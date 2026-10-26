# make_theta Function Documentation --------------------------------------------

#' Construct a Sparse Density-Ratio Parameter Vector
#'
#' Creates a parameter vector \code{theta} of length \code{p + 1} for use in
#' density-ratio or exponential-tilting models based on the augmented covariate
#' vector \eqn{Z = (1, X^T)^T}. The first element corresponds to the intercept,
#' while selected covariate positions are assigned user-specified coefficients.
#'
#' This function is primarily used to define covariate shift patterns between
#' datasets in the simulation design.
#'
#' @param p Integer. Number of covariates in \code{X}.
#' @param idx Integer vector. Indices of covariates in \code{X} to which
#'   nonzero coefficients are assigned.
#' @param coef Numeric vector. Coefficient values assigned to the covariates
#'   specified by \code{idx}. Must have the same length as \code{idx}.
#' @param intercept Numeric. Intercept term in the augmented parameter vector.
#'   Default is 0.
#'
#' @return A numeric vector of length \code{p + 1}, where:
#' \item{\code{theta[1]}}{is the intercept term;}
#' \item{\code{theta[idx + 1]}}{contains the user-specified coefficients;}
#' \item{All remaining elements}{are zero.}
#'
#' @details
#' The returned parameter vector is defined with respect to the augmented design
#' vector \eqn{Z = (1, X^T)^T}. Therefore, the coefficient for covariate
#' \code{X_j} is stored at position \code{j + 1} in \code{theta}.
#'
#' @examples
#' \dontrun{
#' theta <- make_theta(
#'   p = 20,
#'   idx = c(2, 12),
#'   coef = c(2, -0.5),
#'   intercept = 0
#' )
#' }
#'
#' @keywords internal
# --------------------------------------------------------------------------------

make_theta <- function(p, idx, coef, intercept = 0) {
  stopifnot(length(idx) == length(coef))
  theta <- numeric(p + 1)
  theta[1] <- intercept          # Intercept term for Z = (1, X^T)^T
  theta[idx + 1] <- coef         # Shift by +1 because theta includes the intercept
  theta
}