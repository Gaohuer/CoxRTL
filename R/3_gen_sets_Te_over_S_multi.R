# gen_sets_Te_over_S_multi Function Documentation -------------------------------

#' Partition a Covariate Pool into Target, Deployment, and Multiple External Datasets
#'
#' Splits an initial covariate pool into one target dataset, one deployment dataset,
#' and multiple external datasets under predefined covariate shift settings.
#' Sampling probabilities are constructed from exponential-tilting models defined
#' by user-specified density-ratio parameters, and dataset-specific intercepts are
#' calibrated to approximately match the desired sample sizes.
#'
#' @param X_pool A numeric matrix or data.frame containing the initial covariate pool.
#'   Each row represents one subject and each column represents one covariate.
#' @param theta_TeS A list of numeric vectors, where each element specifies the
#'   density-ratio parameter for one external dataset relative to the deployment
#'   dataset. Each vector should have length \code{p + 1}, corresponding to the
#'   augmented covariate vector \eqn{Z = (1, X^T)^T}.
#' @param theta_TeTr A numeric vector of length \code{p + 1} specifying the
#'   density-ratio parameter for the target dataset relative to the deployment dataset.
#' @param Ns An integer vector giving the desired sample sizes for the external
#'   datasets. Its length determines the number of external datasets generated.
#' @param Nt Integer. Desired sample size for the target dataset.
#' @param Nte Integer. Desired sample size for the deployment dataset.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing:
#' \item{\code{Xs1, Xs2, ..., XsK}}{External datasets sampled from \code{X_pool},
#' where \code{K = length(Ns)}.}
#' \item{\code{Xtrain}}{Target dataset sampled from \code{X_pool}.}
#' \item{\code{Xtest}}{Deployment dataset sampled from \code{X_pool}.}
#' \item{\code{aS}}{Estimated intercept terms for the external datasets, used to
#' calibrate class probabilities to the target sample sizes.}
#' \item{\code{aTr}}{Estimated intercept term for the target dataset.}
#' \item{\code{idx}}{A named list of row indices from \code{X_pool} assigned to
#' each dataset.}
#'
#' @examples
#' \dontrun{
#' theta_s1 <- make_theta(p = 20, idx = c(2, 12), coef = c(2, -0.5))
#' theta_s2 <- make_theta(p = 20, idx = c(3, 10), coef = c(1.5, -0.3))
#' theta_tr <- make_theta(p = 20, idx = c(2, 11), coef = c(1, -0.5))
#'
#' sets <- gen_sets_Te_over_S_multi(
#'   X_pool = X,
#'   theta_TeS = list(theta_s1, theta_s2),
#'   theta_TeTr = theta_tr,
#'   Ns = c(300, 300),
#'   Nt = 100,
#'   Nte = 500,
#'   seed = 123
#' )
#'
#' Xs1 <- sets$Xs1
#' Xtarget <- sets$Xtrain
#' Xdeploy <- sets$Xtest
#' }
#'
#' @keywords internal
# --------------------------------------------------------------------------------
gen_sets_Te_over_S_multi <- function(
    X_pool,
    theta_TeS,    # A list of theta_TeS vectors, allowing multiple external datasets
    theta_TeTr,
    Ns,           # A vector of sample sizes for the external datasets
    Nt, Nte,
    seed = NULL
){
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  # ------------------------------------------------------------
  # Step 1. Prepare the augmented covariate matrix and checks
  # ------------------------------------------------------------
  
  n <- nrow(X_pool)
  Z <- cbind(1, X_pool)  # Augmented covariate matrix with intercept
  
  # Determine the number of external datasets
  num_classes <- length(Ns)
  names(Ns) <- paste0("S", 1:num_classes)
  
  stopifnot(num_classes <= length(theta_TeS))  # theta_TeS must match the number of external datasets
  stopifnot(sum(Ns) + Nt + Nte <= n)           # Total requested sample size cannot exceed pool size
  
  # ------------------------------------------------------------
  # Step 2. Define target sampling proportions and calibration loss
  # ------------------------------------------------------------
  
  # S1, S2, ..., S_K correspond to external datasets,
  # Tr corresponds to the target dataset,
  # and Te corresponds to the deployment dataset
  target <- c(Ns / n, Tr = Nt / n)
  
  # Loss function for calibrating intercept terms
  loss <- function(par){
    aS <- par[1:num_classes]
    aTr <- par[num_classes + 1]
    
    scores <- c()
    for(i in 1:num_classes){
      scores <- cbind(
        scores,
        aS[i] - drop(Z %*% theta_TeS[[i]])  # Score for each external dataset
      )
    }
    scores <- cbind(
      scores,
      aTr - drop(Z %*% theta_TeTr),  # Score for target dataset
      0                              # Reference score for deployment dataset
    )
    
    P <- softmax(scores)
    sum((colMeans(P)[1:(num_classes + 1)] - target)^2)
  }
  
  # ------------------------------------------------------------
  # Step 3. Optimize intercept terms for probability calibration
  # ------------------------------------------------------------
  
  ahat <- optim(rep(0, num_classes + 1), loss)$par
  aS <- ahat[1:num_classes]
  aTr <- ahat[num_classes + 1]
  
  # ------------------------------------------------------------
  # Step 4. Compute final class probabilities
  # ------------------------------------------------------------
  
  scores <- c()
  for(i in 1:num_classes){
    scores <- cbind(
      scores,
      aS[i] - drop(Z %*% theta_TeS[[i]])  # Score for each external dataset
    )
  }
  scores <- cbind(
    scores,
    aTr - drop(Z %*% theta_TeTr),  # Score for target dataset
    0                              # Reference score for deployment dataset
  )
  
  P <- softmax(scores)  # Column order: S1, S2, ..., SK, Tr, Te
  
  # ------------------------------------------------------------
  # Step 5. Perform quota-based weighted sampling without replacement
  # ------------------------------------------------------------
  
  idx_all <- seq_len(n)
  quotas <- c(Ns, Tr = Nt, Te = Nte)
  
  # Sample smaller groups first to improve stability
  order_classes <- names(sort(quotas, decreasing = FALSE))
  
  # Store sampled row indices for each dataset
  class_idx <- vector("list", length = num_classes + 2)
  names(class_idx) <- c(paste0("S", 1:num_classes), "Tr", "Te")
  
  remaining <- idx_all
  
  pick_one_class <- function(cls, size, remaining, P) {
    if (size == 0) return(integer(0))
    col <- match(cls, c(paste0("S", 1:num_classes), "Tr", "Te"))
    sample(remaining, size = size, prob = P[remaining, col], replace = FALSE)
  }
  
  for (cls in order_classes) {
    need <- quotas[[cls]]
    if (need == 0) next
    
    got <- pick_one_class(cls, need, remaining, P)
    
    if (length(got) < need) {
      stop(sprintf(
        "Sampling failed: insufficient quota for class %s (%d/%d obtained). Consider increasing the pool size or relaxing the sampling strategy.",
        cls, length(got), need
      ))
    }
    
    class_idx[[cls]] <- got
    remaining <- setdiff(remaining, got)
  }
  
  # ------------------------------------------------------------
  # Step 6. Return the partitioned datasets and calibration results
  # ------------------------------------------------------------
  
  result <- list()
  for (i in 1:num_classes) {
    result[[paste0("Xs", i)]] <- X_pool[class_idx[[paste0("S", i)]], , drop = FALSE]
  }
  result$Xtrain <- X_pool[class_idx$Tr, , drop = FALSE]
  result$Xtest  <- X_pool[class_idx$Te, , drop = FALSE]
  result$aS <- aS
  result$aTr <- aTr
  result$idx <- class_idx
  
  return(result)
}