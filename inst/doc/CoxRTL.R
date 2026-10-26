## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 5)

## ----setup--------------------------------------------------------------------
library(CoxRTL)
library(survival)
library(survcomp)
library(timeROC)
library(survAUC)
set.seed(123)

## ----params-------------------------------------------------------------------
# Sample sizes
N <- 50000      # Total size of the covariate pool used for sampling
Ns <- 500       # Sample size for each external dataset
Nt <- 200       # Sample size for the target dataset
nt_all <- 1000  # Total sample size of the deployment dataset
nt <- 500       # Size of the calibration subset drawn from the deployment dataset

# Number of covariates
p <- 100

# Survival parameters
kappa_t <- 4    # Baseline hazard parameter for the target/deployment data
kappa_s <- 2    # Baseline hazard parameter for the external data

# Penalized Cox regression parameters
lambda <- "lambda.min"   # Alternative choice: "lambda.1se"
alpha <- 0.5               # Ridge penalty; alternative choices: 1 (lasso) or 0.5 (elastic net)

# Covariate shift parameters
ratio_scale_TeS <- 0.5   # Shift magnitude between deployment and external datasets
ratio_scale_TeTr <- 0.5  # Shift magnitude between deployment and target dataset

## ----generate, results = "hide"-----------------------------------------------
sim <- gen_surv_data(
  N = N,
  Ns = Ns,
  Nt = Nt,
  nt = nt_all,
  kappa_t = kappa_t,
  kappa_s = kappa_s,
  ratio_scale_TeS = ratio_scale_TeS,
  ratio_scale_TeTr = ratio_scale_TeTr,
  p = p,
  num_Xs = 3
)

## ----split--------------------------------------------------------------------
sample_calibrate <- sample(1:nt_all, nt)

data_calibrate <- sim$data_test[sample_calibrate, ]
data_evaldata  <- sim$data_test[-sample_calibrate, ]

# Target and external datasets
data_target <- sim$data_train
data_s <- sim$data_s

# Sizes of the four pieces of data
c(target = nrow(data_target),
  calibration = nrow(data_calibrate),
  evaluation = nrow(data_evaldata),
  external = sum(vapply(data_s, nrow, integer(1))))

## ----inputs-------------------------------------------------------------------
# Combine all external datasets
X_s_all <- do.call(rbind, lapply(1:3, function(i) {
  as.matrix(data_s[[i]][, 1:p])
}))

time_s_all <- do.call(c, lapply(1:3, function(i) data_s[[i]]$time))
status_s_all <- do.call(c, lapply(1:3, function(i) data_s[[i]]$status))
y_s_all <- Surv(time = time_s_all, event = status_s_all)

# Target dataset
X_target <- as.matrix(data_target[, 1:p])
y_target <- Surv(data_target$time, data_target$status)

# Deployment subsets
X_calibrate <- as.matrix(data_calibrate[, 1:p])
X_evaldata  <- as.matrix(data_evaldata[, 1:p])

# Predictor names
covariates <- colnames(data_target)[1:p]

## ----screen, warning = FALSE--------------------------------------------------
X_s_list <- lapply(data_s, function(d) as.matrix(d[, 1:p]))
y_s_list <- lapply(data_s, function(d) Surv(d$time, d$status))

screen <- detect_transferable(
  X_t = X_target, y_t = y_target,
  X_s_list = X_s_list, y_s_list = y_s_list,
  rule = "median", min_prop_positive = 0,
  alpha = alpha, V = 5, seed = 123
)
summary(screen)

## ----weights------------------------------------------------------------------
# External-to-deployment weights
w_source <- do.call(c, lapply(1:3, function(i) {
  calculate_weights(
    X_denominator = as.matrix(data_s[[i]][, 1:p]),
    X_numerator = X_calibrate,
    cv_lambda = FALSE,
    ncore = 1
  )$w
}))

# Target-to-deployment weights
weight_target <- calculate_weights(
  X_denominator = X_target,
  X_numerator = X_calibrate,
  cv_lambda = FALSE,
  ncore = 1
)
summary(weight_target)
w_target <- weight_target$w

## ----nontransfer-unw, eval=FALSE----------------------------------------------
# nontransfer_unw <- fit_coxrtl(
#   X_t = X_target, y_t = y_target,
#   transfer = FALSE, alpha = alpha, lambda = lambda,
#   type.measure = "deviance", seed = 123,
#   bootstrap = TRUE, nboot = 100
# )

## ----nontransfer-w, eval=FALSE------------------------------------------------
# nontransfer_w <- fit_coxrtl(
#   X_t = X_target, y_t = y_target,
#   w_t = w_target,
#   transfer = FALSE, alpha = alpha, lambda = lambda,
#   type.measure = "deviance", seed = 123,
#   bootstrap = TRUE, nboot = 100
# )

## ----transfer-unw, eval=FALSE-------------------------------------------------
# transfer_unw <- fit_coxrtl(
#   X_t = X_target, X_s = X_s_all, y_t = y_target, y_s = y_s_all,
#   transfer = TRUE, alpha1 = alpha, alpha2 = alpha, lambda = lambda,
#   type.measure = "deviance", seed = 123,
#   bootstrap = TRUE, nboot = 100
# )

## ----transfer-w---------------------------------------------------------------
transfer_w <- fit_coxrtl(
  X_t = X_target, X_s = X_s_all, y_t = y_target, y_s = y_s_all,
  w_t = w_target, w_s = w_source,
  transfer = TRUE, alpha1 = alpha, alpha2 = alpha, lambda = lambda,
  type.measure = "deviance", seed = 123,
  bootstrap = TRUE, nboot = 100
)
summary(transfer_w, head = 10)         # top coefficients with interval bounds
head(as.data.frame(transfer_w$ci))     # se / lower / upper for the

## ----predict------------------------------------------------------------------
preds <- predict_coxrtl(X = X_evaldata, transfer_w)$pred
risks <- predict_coxrtl(X = X_evaldata, transfer_w)$risk

## ----cindex-------------------------------------------------------------------
cindex <- survcomp::concordance.index(
    x = preds,
    surv.time = data_evaldata$time,
    surv.event = data_evaldata$status,
    method = "noether"
  )$c.index
cindex

## ----auc-readme---------------------------------------------------------------
eval_times <- as.numeric(quantile(data_evaldata$time[data_evaldata$status == 1], c(0.25, 0.5, 0.75)))
auc <- timeROC(
    T = data_evaldata$time,
    delta = data_evaldata$status,
    marker = preds,
    cause = 1,
    weighting = "marginal",
    times = eval_times
  )
auc

## ----test---------------------------------------------------------------------
preds_train <- predict_coxrtl(X = X_target, transfer_w)$pred

# Evaluation time points
times <- c(0.77, 1.47, 2.20, 2.93)

# Calculate Brier scores
briers <- predErr <- survAUC::predErr(
    Surv.rsp = Surv(
      data_target$time,
      data_target$status
    ),
    Surv.rsp.new = Surv(
      data_evaldata$time,
      data_evaldata$status
    ),
    lp = as.matrix(preds_train),
    lpnew = as.matrix(preds),
    times = times,
    type = "brier",
    int.type = "weighted"
  )
predErr

