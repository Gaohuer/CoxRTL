# softmax Function Documentation ------------------------------------------------

#' Compute Row-Wise Softmax Probabilities
#'
#' Applies the softmax transformation to each row of an input matrix, converting
#' row-wise scores into normalized probabilities that sum to 1. A numerical
#' stabilization step is used by subtracting the row maximum before exponentiation.
#'
#' This function is used to convert class-specific sampling scores into multinomial
#' probabilities in the simulation-based dataset partitioning procedure.
#'
#' @param M A numeric matrix of dimension \eqn{n \times K}, where each row contains
#'   the unnormalized scores for \eqn{K} classes.
#'
#' @return A numeric matrix of the same dimension as \code{M}, where each row
#'   contains softmax-transformed probabilities summing to 1.
#'
#' @details
#' For each row \eqn{i}, the softmax transformation is defined as
#' \deqn{
#' P_{ij} = \frac{\exp(M_{ij})}{\sum_{k=1}^K \exp(M_{ik})}.
#' }
#' To improve numerical stability, the maximum value in each row is subtracted
#' before exponentiation.
#'
#' @examples
#' \dontrun{
#' M <- matrix(c(1, 2, 3,
#'               0, 1, 0), nrow = 2, byrow = TRUE)
#' softmax(M)
#' }
#'
#' @keywords internal
# --------------------------------------------------------------------------------

softmax <- function(M) {
  # Subtract the row-wise maximum for numerical stability
  M <- sweep(M, 1, apply(M, 1, max), "-")
  
  # Exponentiate the stabilized scores
  EM <- exp(M)
  
  # Normalize each row so that probabilities sum to 1
  EM / rowSums(EM)
}