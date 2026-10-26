#' CoxRTL: Deployment-Oriented Recalibrated Transfer Learning for Survival Prediction
#'
#' Implements deployment-oriented recalibrated transfer learning for survival
#' prediction.
#'
#' @keywords internal
"_PACKAGE"

# Imported wholesale so the exported helper functions behave exactly as the
# original top-level simulation scripts (which loaded these packages via
# library()).
#' @import MASS
#' @import survival
#' @import glmnet
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom stats coef nlm optim quantile rbinom runif rweibull sd
NULL

# `pos` and `rep_k` are the iterator variables of foreach() %dopar% loops in
# get_theta() and fit_coxrtl() (bootstrap replicates) respectively.
utils::globalVariables(c("pos", "rep_k"))
