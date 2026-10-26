# gen_surv_data Function Documentation -------------------------------------

#' Simulated Survival Dataset Generation under Covariate Shift
#'
#' Generates synthetic survival datasets for evaluating transfer learning methods
#' under covariate shift. The function creates one target dataset, one deployment
#' dataset, and one or more external datasets, all sampled from an initial covariate
#' pool under predefined density-ratio shift settings. Survival outcomes are then
#' generated from Weibull models with potentially different regression coefficients
#' and baseline hazard parameters across target/deployment and external datasets.
#'
#' @param N Integer. Size of the initial covariate pool from which all datasets
#'   are sampled.
#' @param Ns Integer. Sample size for each external dataset.
#' @param Nt Integer. Sample size for the target dataset.
#' @param nt Integer. Sample size for the deployment dataset.
#' @param kappa_t Numeric. Baseline hazard-related parameter used for generating
#'   survival times in the target and deployment datasets.
#' @param kappa_s Numeric. Baseline hazard-related parameter used for generating
#'   survival times in the external datasets.
#' @param ratio_scale_TeS Numeric. Scaling factor controlling the magnitude of
#'   covariate shift between the deployment dataset and the external datasets.
#' @param ratio_scale_TeTr Numeric. Scaling factor controlling the magnitude of
#'   covariate shift between the deployment dataset and the target dataset.
#' @param p Integer. Number of covariates.
#' @param num_Xs Integer. Number of external datasets to generate. Default is 1.
#'
#' @return A list with the following components:
#' \item{\code{data_s}}{A list of external datasets. Each element is a data.frame
#' containing covariates (\code{X1, ..., Xp}), observed survival time
#' (\code{time}), and event indicator (\code{status}).}
#' \item{\code{data_train}}{A data.frame for the target dataset, including
#' covariates, observed survival time, and event indicator.}
#' \item{\code{data_test}}{A data.frame for the deployment dataset, including
#' covariates, observed survival time, and event indicator.}
#'
#' @examples
#' \dontrun{
#' sim_data <- gen_surv_data(
#'   N = 5000,
#'   Ns = 500,
#'   Nt = 100,
#'   nt = 500,
#'   kappa_t = 1,
#'   kappa_s = 1,
#'   ratio_scale_TeS = 1,
#'   ratio_scale_TeTr = 1,
#'   p = 20,
#'   num_Xs = 2
#' )
#'
#' target_data <- sim_data$data_train
#' deploy_data <- sim_data$data_test
#' external_data_list <- sim_data$data_s
#' }
#'
#' @export
# ------------------------------------------------------------------------------------


gen_surv_data <- function(
    N,                 # Size of the initial covariate pool
    Ns,                # Sample size for each external dataset
    Nt,                # Sample size for the target dataset
    nt,                # Sample size for the deployment dataset
    kappa_t,           # Baseline hazard parameter for target/deployment data
    kappa_s,           # Baseline hazard parameter for external data
    ratio_scale_TeS,   # Magnitude of covariate shift between deployment and external data
    ratio_scale_TeTr,  # Magnitude of covariate shift between deployment and target data
    p,                 # Number of covariates
    num_Xs = 1         # Number of external datasets
){
  
  # set.seed(1)
  
  # ------------------------------------------------------------
  # Step 1. Generate the initial covariate pool X
  # ------------------------------------------------------------
  
  # X1 and X2 follow a bivariate normal distribution
  X12 <- mvrnorm(N, c(2, 1), diag(c(2, 4)))
  
  # X3 is binary
  X3 <- rbinom(N, 1, 0.4)
  
  # Remaining covariates X4, ..., Xp follow a multivariate normal distribution
  X.Sigma <- diag(rep(1, p - 3))
  for (i in 1:(p - 3)) {
    for (j in 1:(p - 3)) {
      X.Sigma[i, j] <- 0^abs(i - j)
    }
  }
  
  X <- cbind(X12, X3, mvrnorm(N, c(rep(0, p - 3)), X.Sigma))
  colnames(X) <- paste0("X", 1:p)
  
  # ------------------------------------------------------------
  # Step 2. Define density-ratio parameters governing
  #         covariate shift across data sources
  # ------------------------------------------------------------
  
  # Covariate shift between deployment and each external dataset
  theta_TeS_list <- list()
  for (i in 1:num_Xs) {
    theta_TeS_list[[i]] <- make_theta(
      p,
      idx = c(2, 12),
      coef = c(2, -0.5) * ratio_scale_TeS
    )
  }
  
  # Covariate shift between deployment and target dataset
  theta_TeTr <- make_theta(
    p,
    idx = c(2, 11),
    coef = c(1, -0.5) * ratio_scale_TeTr
  )
  
  # ------------------------------------------------------------
  # Step 3. Split the covariate pool into:
  #         - external datasets
  #         - target dataset
  #         - deployment dataset
  # ------------------------------------------------------------
  
  sets <- gen_sets_Te_over_S_multi(
    X,
    theta_TeS_list,       # Multiple shift settings for external datasets
    theta_TeTr,           # Shift setting for target dataset
    rep(Ns, num_Xs),      # Sample sizes for external datasets
    Nt,                   # Sample size for target dataset
    nt,                   # Sample size for deployment dataset
    seed = NULL
  )
  
  # Extract partitioned datasets
  Xs_list  <- lapply(1:num_Xs, function(i) sets[[paste0("Xs", i)]])
  Xtrain <- sets$Xtrain
  Xtest <- sets$Xtest
  
  # ------------------------------------------------------------
  # Step 4. Define the Weibull scale parameter
  # ------------------------------------------------------------
  
  get_scale<-function(X,coef,kappa){
    XB<-as.matrix(cbind(X[,1],X[,2]^2,X[,3],X[,4]^2,X[,5:p])) %*% coef
    scale <- sqrt(2 / kappa) * exp(-XB/ 2)
    scale
  }
  
  coef_t <- c(c(-0.5, -0.3, 1, -0.4, 0.3,0.5, -1, 1, 0.5, 0.2),rep(0,p-10))  # train test
  
  scale_train<-get_scale(Xtrain,coef_t,kappa_t)
  scale_test<-get_scale(Xtest,coef_t,kappa_t)
  scale_s_list<-lapply(1:num_Xs, function(i){
    coef_s<-c(c(0.2, -0.3, 1, -0.4, 0.3,0.5, -1, 1, -0.2, 0.2),rep(0,p-10)) 
    get_scale(Xs_list[[i]],coef_s,kappa_s)
  })
  
  # ------------------------------------------------------------
  # Step 5. Generate survival time and censoring status
  # ------------------------------------------------------------
  
  time_status_generator <- function(X, scale, kappa) {
    T_s <- rweibull(nrow(X), shape = 2, scale = scale)
    C_s <- runif(nrow(X), 0, 3.55)
    time <- ifelse(T_s < C_s, T_s, C_s)
    status <- ifelse(T_s < C_s, 1, 0)
    as.data.frame(cbind(X, time, status))
  }
  
  Xtrain <- time_status_generator(Xtrain, scale_train, kappa_t)
  Xtest <- time_status_generator(Xtest, scale_test, kappa_t)
  
  Xs_list <- lapply(1:num_Xs, function(i) {
    time_status_generator(Xs_list[[i]], scale_s_list[[i]], kappa_s)
  })
  
  # ------------------------------------------------------------
  # Step 6. Return simulated datasets
  # ------------------------------------------------------------
  
  result<-list()
  result$data_s <- Xs_list
  result$data_train <- Xtrain
  result$data_test <- Xtest
  
  return(result)
}
