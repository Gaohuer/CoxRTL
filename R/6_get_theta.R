# get_theta Function Documentation ----------------------------------------------

#' Select the Regularization Parameter for Density-Ratio Estimation
#'
#' Performs K-fold cross-validation to select the regularization parameter
#' \code{lambda} used in the density-ratio estimation objective. For each
#' candidate \code{lambda}, the function fits the exponential-tilting model
#' on the training folds and evaluates the unpenalized validation loss on
#' the held-out folds.
#'
#' This function is used within \code{calculate_weights()} to determine the
#' penalty level for density-ratio weighting under covariate shift.
#'
#' @param X_denominator A numeric matrix containing the augmented covariates
#'   of the denominator dataset. Rows correspond to subjects and columns
#'   correspond to the intercept and covariates.
#' @param X_numerator A numeric matrix containing the augmented covariates
#'   of the numerator dataset. Rows correspond to subjects and columns
#'   correspond to the intercept and covariates.
#' @param K Integer. Number of folds for cross-validation. Default is 5.
#' @param iterlim Integer. Maximum number of iterations passed to \code{nlm()}.
#'   Default is 500.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#' @param ncore Integer. Number of cores for parallelizing the candidate
#'   \code{lambda} grid across folds. \code{<= 1} (default) runs sequentially;
#'   \code{> 1} registers a local \code{doParallel} backend for this call.
#'
#' @return A list containing:
#' \item{\code{lambda.min}}{The value of \code{lambda} achieving the smallest
#' mean cross-validation loss.}
#' \item{\code{lambda.1se}}{The largest \code{lambda} whose mean cross-validation
#' loss is within one standard error of the minimum.}
#' \item{\code{cv_mean}}{Mean cross-validation loss for each candidate \code{lambda}.}
#' \item{\code{cv_se}}{Standard error of the cross-validation loss for each candidate \code{lambda}.}
#' \item{\code{cv_table}}{A matrix of fold-specific validation losses, with rows
#' corresponding to candidate \code{lambda} values and columns to folds.}
#'
#' @details
#' Separate K-fold partitions are created for the denominator and numerator datasets.
#' For each fold and each candidate \code{lambda}, the density-ratio model is fitted
#' on the training subsets using \code{nlm()}, and validation performance is assessed
#' using the unpenalized loss:
#' \deqn{
#' \frac{1}{n_s}\sum_{i=1}^{n_s}\exp(Z_{s,i}^T\theta) -
#' \frac{1}{n_t}\sum_{i=1}^{n_t} Z_{t,i}^T\theta.
#' }
#' Parallel computation is used across candidate \code{lambda} values within each fold.
#'
#' @examples
#' \dontrun{
#' theta_cv <- get_theta(
#'   X_denominator = Xs_aug,
#'   X_numerator = Xt_aug,
#'   K = 5,
#'   iterlim = 500,
#'   seed = 123
#' )
#' }
#'
#' @keywords internal
# --------------------------------------------------------------------------------

get_theta <- function(X_denominator, X_numerator, K = 5, iterlim = 500, seed = NULL,
                      ncore = 1){
  
  # ------------------------------------------------------------
  # Step 1. Set up parallel computing
  # ------------------------------------------------------------
  
  # Use ncore > 1 to enable parallelism; otherwise stay fully sequential.
  # A backend owned by this call is shut down again before the function returns.
  backend <- start_parallel_backend(ncore)
  if (backend$owned) on.exit(stop_parallel_backend(backend), add = TRUE)
  
  # ------------------------------------------------------------
  # Step 2. Initialize inputs and create fold partitions
  # ------------------------------------------------------------
  
  theta_ini <- rep(0, ncol(X_denominator))
  ns <- nrow(X_denominator)
  nt <- nrow(X_numerator)
  
  # Create approximately balanced K-fold splits separately
  # for the denominator and numerator datasets
  if(!is.null(seed)){
    set.seed(seed)
  }
  folds_s <- split(sample.int(ns), rep(1:K, length.out = ns))
  # if(!is.null(seed)){
  #   set.seed(seed)
  # }
  folds_t <- split(sample.int(nt), rep(1:K, length.out = nt))
  
  # ------------------------------------------------------------
  # Step 3. Define candidate lambda values and storage objects
  # ------------------------------------------------------------
  
  lambdas <- c(10^seq(-2, 1, length.out = 20))
  cv_mat <- matrix(NA_real_, nrow = length(lambdas), ncol = K)
  rownames(cv_mat) <- paste0("lam=", signif(lambdas, 3))
  
  # Evaluate larger lambdas first
  order_idx <- order(lambdas, decreasing = TRUE)
  thetas_warm <- vector("list", length(lambdas))
  names(thetas_warm) <- as.character(seq_along(lambdas))
  
  # ------------------------------------------------------------
  # Step 4. Perform K-fold cross-validation
  # ------------------------------------------------------------
  
  `%par%` <- if (backend$active) foreach::`%dopar%` else foreach::`%do%`

  for (k in 1:K) {
    
    idx_s_val <- folds_s[[k]]
    idx_t_val <- folds_t[[k]]
    
    Zs_tr <- X_denominator[-idx_s_val, , drop = FALSE]
    Zt_tr <- X_numerator[-idx_t_val, , drop = FALSE]
    Zs_val <- X_denominator[idx_s_val, , drop = FALSE]
    Zt_val <- X_numerator[idx_t_val, , drop = FALSE]
    
    # Fit the density-ratio model for each candidate lambda
    # and evaluate the unpenalized validation loss
    val_losses <- foreach(
      pos = order_idx,
      .combine = "c",
      .export = "density_opt"
    ) %par% {
      
      lam <- lambdas[pos]
      
      # Training: penalized objective
      fn_train <- density_opt(Zs_tr, Zt_tr, lam)
      fit <- nlm(fn_train, theta_ini, iterlim = iterlim)
      theta_hat <- as.numeric(fit$estimate)
      
      # Store fitted theta (reserved for potential warm starts)
      thetas_warm[[pos]] <- theta_hat
      
      # Validation: unpenalized loss
      val_loss <- mean(exp(drop(Zs_val %*% theta_hat))) -
        mean(drop(Zt_val %*% theta_hat))
      
      val_loss
    }
    
    cv_mat[order_idx, k] <- val_losses
  }
  
  # ------------------------------------------------------------
  # Step 5. Select lambda using minimum loss and 1-SE rule
  # ------------------------------------------------------------
  
  cv_mean <- rowMeans(cv_mat)
  cv_se <- apply(cv_mat, 1, function(x) sd(x) / sqrt(K))
  
  best_id <- which.min(cv_mean)
  lambda.min <- lambdas[best_id]
  
  # 1-SE rule: choose the largest lambda within one standard error
  # of the minimum cross-validation loss
  se_rule_id <- max(which(cv_mean <= cv_mean[best_id] + cv_se[best_id]))
  lambda.1se <- lambdas[se_rule_id]
  
  # ------------------------------------------------------------
  # Step 6. Return cross-validation results
  # ------------------------------------------------------------
  
  list(
    lambda.min = lambda.min,
    lambda.1se = lambda.1se,
    cv_mean = cv_mean,
    cv_se = cv_se,
    cv_table = cv_mat
  )
}