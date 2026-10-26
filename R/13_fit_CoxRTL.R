# fit_coxrtl Function Documentation -------------------------------------------

#' Fit One Cox-Based Survival Model with Optional Confidence Intervals
#'
#' Fits penalized Cox models compared in the package workflow, 
#' optionally with coefficient-level confidence intervals. The 
#' models are selected by the logical argument \code{transfer}:
#' \itemize{
#'   \item \code{transfer = FALSE}: penalized Cox model trained on the target
#'   data only, i.e. without transfer learning
#'   (\code{\link{no_transfer_learning}}). Passing \code{w_t} reweights it
#'   toward the deployment population;
#'   \item \code{transfer = TRUE} (default): two-step transfer learning that
#'   borrows the external data (\code{\link{transfer_learning}}). Without
#'   \code{w_t} and \code{w_s} it uses unit weights and is not calibrated for
#'   covariate shift; supplying them reweights the target and external data
#'   toward the deployment population.
#' }
#'
#' Each call fits exactly \emph{one} model (argument \code{transfer}) and
#' returns its point estimate \code{beta}; \code{bootstrap = TRUE} additionally
#' returns nonparametric bootstrap confidence intervals in \code{ci} (with
#' \code{nboot} resamples), otherwise \code{ci} is \code{NULL}.
#'
#' @param X_t A numeric matrix of predictors for the target dataset.
#' @param X_s A numeric matrix of predictors for the \emph{pooled} external
#'   dataset(s). Pool the external datasets with \code{rbind()} before calling
#'   this function if several external sources are used. Required when
#'   \code{transfer = TRUE}; with \code{transfer = FALSE} it may be left
#'   \code{NULL}, since the target-only model does not use it.
#' @param y_t A survival response object for the target dataset, typically
#'   created using \code{\link[survival]{Surv}()}.
#' @param y_s A survival response object for the pooled external data, typically
#'   created using \code{\link[survival]{Surv}()}. Required when
#'   \code{transfer = TRUE}; otherwise it may be left \code{NULL} together with
#'   \code{X_s}.
#' @param w_t Optional numeric vector of density-ratio weights for the target
#'   dataset (see \code{\link{calculate_weights}}), reweighting it toward the
#'   deployment population. It is consumed by both models: with
#'   \code{transfer = FALSE} it weights the target-only fit, and with
#'   \code{transfer = TRUE} it weights both the pooled first step and the
#'   target-specific debiasing step. Leaving it \code{NULL} counts as unit
#'   weights.
#' @param w_s Optional numeric vector of density-ratio weights for the pooled
#'   external data, used when \code{transfer = TRUE} alone (whose first step
#'   pools the external data with the target data). Leaving it \code{NULL}
#'   counts as unit weights. It is ignored when \code{transfer = FALSE}, with a
#'   warning, because that model does not use the external data at all.
#' @param alpha Numeric. Elastic net mixing parameter shared by the requested
#'   model. \code{alpha = 0} is ridge, \code{alpha = 1} is lasso and
#'   \code{alpha = 0.5} is the elastic net. Required when
#'   \code{transfer = FALSE}; when \code{transfer = TRUE} it is only a fallback
#'   for \code{alpha1}/\code{alpha2}.
#' @param alpha1 Numeric, optional. Mixing parameter for the first-step combined
#'   model in \code{\link{transfer_learning}}. Defaults to \code{alpha}.
#' @param alpha2 Numeric, optional. Mixing parameter for the second-step
#'   target-specific debiasing model. Defaults to \code{alpha}.
#' @param lambda Character string or numeric value giving the \code{s} used to
#'   extract coefficients. It selects the penalty level of the target-only fit
#'   (\code{transfer = FALSE}); when \code{transfer = TRUE} it is only a
#'   fallback for \code{lambda1}/\code{lambda2}. Common choices are
#'   \code{"lambda.min"} (default) and \code{"lambda.1se"}.
#' @param lambda1 Character string or numeric, optional. Lambda rule for the
#'   first-step transfer-learning model. Defaults to \code{lambda}.
#' @param lambda2 Character string or numeric, optional. Lambda rule for the
#'   second-step transfer-learning model. Defaults to \code{lambda}.
#' @param type.measure Character string specifying the cross-validation loss
#'   used in \code{\link[glmnet]{cv.glmnet}}. For Cox models, common choices are
#'   \code{"deviance"} (default) and \code{"C"}.
#' @param seed Integer or NULL. Random seed passed to the fitted model for
#'   reproducibility. Default is NULL. With \code{bootstrap = TRUE} it
#'   additionally seeds the bootstrap resampling, so a given \code{seed} fully
#'   reproduces the intervals.
#' @param ncore Integer. Number of cores. When \code{bootstrap = FALSE} it is
#'   forwarded to the underlying
#'   \code{\link{no_transfer_learning}}/\code{\link{transfer_learning}} call to
#'   parallelize its cross-validation. When \code{bootstrap = TRUE} it instead
#'   parallelizes \emph{across bootstrap replicates} (the internal refits then
#'   run sequentially); \code{<= 1} (default) runs sequentially.
#' @param nfolds Integer. Number of cross-validation folds used by the underlying
#'   \code{\link[glmnet]{cv.glmnet}} call (both steps when
#'   \code{transfer = TRUE}). Default is 5.
#' @param standardize Logical. Whether to standardize the predictor columns
#'   before fitting. Passed to the underlying \code{cv.glmnet}/\code{glmnet}
#'   calls and to the fixed-lambda bootstrap refits, so the resampled
#'   coefficients are on the same scale as the point fit. Default is
#'   \code{FALSE}.
#' @param cox.ties Character. Method for handling tied survival times in the
#'   Cox partial likelihood: \code{"breslow"} (default) or \code{"efron"}.
#' @param transfer Logical. Whether to fit the transfer-learning model
#'   (\code{TRUE}, default) or the target-only model without transfer learning
#'   (\code{FALSE}). The target-only model does not use \code{X_s}/\code{y_s},
#'   so with \code{transfer = FALSE} it can be fitted from the target data
#'   alone.
#' @param bootstrap Logical. Whether to compute nonparametric bootstrap
#'   confidence intervals for the fitted coefficients (requires \code{nboot}).
#'   Default is \code{FALSE}, returning point estimates only.
#' @param nboot Integer. Number of bootstrap resamples. Required when
#'   \code{bootstrap = TRUE}.
#' @param conf.level Numeric in (0, 1). Confidence level of the bootstrap
#'   intervals. Default is 0.95.
#'
#' @return An object of class \code{"fit_coxrtl"}: a list with two components
#' \item{\code{beta}}{A named numeric vector of length \code{ncol(X_t)} with the
#' estimated coefficients of the requested model, named by covariate. This is
#' the point estimate in both shapes.}
#' \item{\code{ci}}{When \code{bootstrap = TRUE}, a list of three full-length,
#' named numeric vectors \code{se}, \code{lower} and \code{upper}: the standard
#' deviation across the bootstrap resamples and the \code{conf.level} quantile
#' interval of every coefficient. When \code{bootstrap = FALSE} (default),
#' \code{ci} is \code{NULL}.}
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' # small simulated example (see the vignette for the full-scale setup)
#' set.seed(123)
#' sim <- gen_surv_data(
#'   N = 5000, Ns = 500, Nt = 100, nt = 500,
#'   kappa_t = 1, kappa_s = 1,
#'   ratio_scale_TeS = 1, ratio_scale_TeTr = 1,
#'   p = 20, num_Xs = 3
#' )
#'
#' # split the deployment data into calibration and evaluation subsets
#' sample_cal <- sample(1:500, 250)
#' data_target    <- sim$data_train
#' data_calibrate <- sim$data_test[sample_cal, ]
#' data_evaldata  <- sim$data_test[-sample_cal, ]
#' data_s         <- sim$data_s
#'
#' p <- 20
#' X_t <- as.matrix(data_target[, 1:p])
#' y_t <- Surv(data_target$time, data_target$status)
#'
#' X_s <- do.call(rbind, lapply(data_s, function(d) as.matrix(d[, 1:p])))
#' time_s <- do.call(c, lapply(data_s, function(d) d$time))
#' status_s <- do.call(c, lapply(data_s, function(d) d$status))
#' y_s <- Surv(time_s, status_s)
#'
#' # density-ratio weights toward the deployment calibration data
#' w_t <- calculate_weights(
#'   X_denominator = X_t,
#'   X_numerator = as.matrix(data_calibrate[, 1:p]),
#'   cv_lambda = FALSE
#' )$w
#' w_s <- unlist(lapply(data_s, function(d)
#'   calculate_weights(
#'     X_denominator = as.matrix(d[, 1:p]),
#'     X_numerator = as.matrix(data_calibrate[, 1:p]),
#'     cv_lambda = FALSE
#'   )$w))
#'
#' # transfer model without weights: no covariate-shift recalibration
#' beta_unweighted <- fit_coxrtl(
#'   X_t = X_t, X_s = X_s, y_t = y_t, y_s = y_s,
#'   transfer = TRUE, alpha = 0, lambda = "lambda.min", seed = 123
#' )
#'
#' # the same model reweighted toward the deployment population
#' beta_transfer <- fit_coxrtl(
#'   X_t = X_t, X_s = X_s, y_t = y_t, y_s = y_s,
#'   w_t = w_t, w_s = w_s,
#'   transfer = TRUE, alpha = 0, lambda = "lambda.min", seed = 123
#' )
#'
#' # the target-only model needs no external data at all, and can also be
#' # reweighted toward the deployment population
#' beta_nontransfer <- fit_coxrtl(
#'   X_t = X_t, y_t = y_t,
#'   w_t = w_t,
#'   transfer = FALSE, alpha = 0, lambda = "lambda.min", seed = 123
#' )
#'
#' # reweighted again, now with bootstrap confidence intervals
#' fit_ci <- fit_coxrtl(
#'   X_t = X_t, X_s = X_s, y_t = y_t, y_s = y_s,
#'   w_t = w_t, w_s = w_s,
#'   transfer = TRUE, alpha = 1, lambda = "lambda.min", seed = 123,
#'   bootstrap = TRUE, nboot = 200
#' )
#' fit_ci$ci$lower       # lower bounds of the bootstrap intervals
#'}
#' @export
# --------------------------------------------------------------------------------

fit_coxrtl <- function(
    X_t,
    X_s = NULL,
    y_t,
    y_s = NULL,
    w_t = NULL,
    w_s = NULL,
    alpha = NULL,
    alpha1 = NULL,
    alpha2 = NULL,
    lambda = "lambda.min",
    lambda1 = NULL,
    lambda2 = NULL,
    type.measure = "deviance",
    seed = NULL,
    ncore = 1,
    nfolds = 5L,
    standardize = FALSE,
    cox.ties = "breslow",
    transfer = TRUE,
    bootstrap = FALSE,
    nboot = NULL,
    conf.level = 0.95
){

  if (!is.logical(transfer) || length(transfer) != 1L || is.na(transfer)) {
    stop("'transfer' must be a single logical (TRUE to fit the ",
         "transfer-learning model, FALSE to fit the target-only model).")
  }

  # Canonical name of the requested model: it drives the branching below and is
  # the label reported by print()/summary().
  model <- if (transfer) "transfer" else "nontransfer"

  if (!is.logical(bootstrap) || length(bootstrap) != 1L || is.na(bootstrap)) {
    stop("'bootstrap' must be a single logical (TRUE to compute bootstrap ",
         "confidence intervals).")
  }
  if (bootstrap) {
    if (!is.numeric(conf.level) || length(conf.level) != 1L ||
        is.na(conf.level) || conf.level <= 0 || conf.level >= 1) {
      stop("'conf.level' must be a single number in (0, 1).")
    }
    if (is.null(nboot) || !is.numeric(nboot) || length(nboot) != 1L ||
        is.na(nboot) || nboot < 1 || nboot != round(nboot)) {
      stop("'nboot' (number of bootstrap resamples) must be a single ",
           "integer >= 1 when bootstrap = TRUE.")
    }
    nboot <- as.integer(nboot)
  }

  # Resolve the transfer-learning tuning parameters from the shared values
  if (is.null(alpha1)) alpha1 <- alpha
  if (is.null(alpha2)) alpha2 <- alpha
  if (is.null(lambda1)) lambda1 <- lambda
  if (is.null(lambda2)) lambda2 <- lambda

  # Only the transfer model pools the external data, so 'X_s'/'y_s' are
  # optional overall but mandatory for it.
  if (model == "transfer" && (is.null(X_s) || is.null(y_s))) {
    stop("'X_s' and 'y_s' must both be supplied to fit the 'transfer' model. ",
         "Use model = \"nontransfer\" to fit on the target data alone.")
  }

  if (model == "nontransfer" && is.null(alpha)) {
    stop("'alpha' must be provided to fit the '", model, "' model.")
  }
  if (model == "transfer" && (is.null(alpha1) || is.null(alpha2))) {
    stop("'alpha' (or both 'alpha1' and 'alpha2') must be provided ",
         "for the transfer-learning model '", model, "'.")
  }

  # The two models are reweighted by different weights: the target-only model
  # consumes w_t alone, the transfer model consumes w_t and w_s. A w_s handed
  # to 'nontransfer' is unused, so it is reported (once) rather than silently
  # making the fit look reweighted when it is not.
  weighted <- if (model == "nontransfer") {
    !is.null(w_t)
  } else {
    !is.null(w_t) || !is.null(w_s)
  }

  if (model == "nontransfer" && !is.null(w_s)) {
    warning("'w_s' is ignored by the 'nontransfer' model, which fits the ",
            "target data only; use model = \"transfer\" to reweight the ",
            "external data.")
  }

  if (model == "transfer" && !weighted) {
    message("w_t and w_s not supplied: the transfer model is fitted with unit ",
            "weights and no covariate-shift reweighting.")
  }

  # With bootstrap = TRUE the requested cores parallelize the bootstrap
  # replicates, so the single full-data fit below runs sequentially.
  fit_ncore <- if (bootstrap) 1 else ncore

  # ------------------------------------------------------------
  # Step 1. Full-data fit of the requested single model (kept as a complete
  #         object so the selected lambda values are available for the CI step)
  # ------------------------------------------------------------

  fit <- switch(model,
    "nontransfer" = no_transfer_learning(
      X = X_t, y = y_t, alpha = alpha, lambda = lambda, w = w_t,
      seed = seed, type.measure = type.measure, ncore = fit_ncore,
      nfolds = nfolds, standardize = standardize, cox.ties = cox.ties),
    "transfer" = transfer_learning(
      X_t = X_t, X_s = X_s, y_t = y_t, y_s = y_s,
      w_t = w_t, w_s = w_s,
      alpha1 = alpha1, alpha2 = alpha2,
      lambda1 = lambda1, lambda2 = lambda2,
      seed = seed, type.measure = type.measure, ncore = fit_ncore,
      nfolds = nfolds, standardize = standardize, cox.ties = cox.ties)
  )

  beta <- fit$beta

  # ------------------------------------------------------------
  # Build the metadata reported by print()/summary() (model name, data sizes,
  # tuning, and the internal fitter's summary). No raw data is stored -- only
  # the fitted model's own output. It is attached below as the "fit_info"
  # attribute of the returned list(beta, ci) object.
  # ------------------------------------------------------------
  # 'X_s'/'y_s' are NULL for the target-only model, where nrow() would return
  # NULL rather than 0 and y_s[, 2] would error, so both are guarded here.
  fit_info <- list(
    model        = model,
    p            = ncol(X_t),
    n_target     = nrow(X_t),
    n_external   = if (is.null(X_s)) 0L else nrow(X_s),
    n_events_t   = sum(y_t[, 2]),
    n_events_s   = if (is.null(y_s)) 0L else sum(y_s[, 2]),
    weighted     = weighted,
    alpha        = alpha,
    alpha1       = alpha1,
    alpha2       = alpha2,
    type.measure = type.measure
  )
  # Store the complete summary the internal fitter returns besides beta, so
  # every quantity reported by no_transfer_learning()/transfer_learning() stays
  # available in fit_info: for the nontransfer model the resolved lambda, nzero
  # and deviance ratio; for the transfer model additionally the first-step
  # combined coefficients (coef_com), the second-step correction (delta) and the
  # per-step lambda / nzero / deviance-ratio statistics.
  fit_info$stats <- fit[setdiff(names(fit), "beta")]

  # The returned object is always a list(beta, ci) of class "fit_coxrtl": with
  # bootstrap = FALSE the intervals are simply NULL, so every fit exposes the
  # point estimate as $beta (and $ci = NULL). The metadata built above
  # (fit_info, bootstrap, conf.level, nboot) is carried as attributes of the
  # list; $beta itself is a plain named numeric vector.
  if (!bootstrap) {
    res <- list(beta = beta, ci = NULL)
    attr(res, "fit_info")  <- fit_info
    attr(res, "bootstrap") <- FALSE
    class(res) <- "fit_coxrtl"
    return(res)
  }

  # ------------------------------------------------------------
  # Step 2. Bootstrap confidence intervals for this single model
  # ------------------------------------------------------------

  # Reproducible replicate-specific seeds; drawn sequentially in the master so
  # that the result does not depend on the parallel scheduling of replicates.
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample(1:1e9, nboot, replace = TRUE)

  lam_num <- if (model == "nontransfer") {
    fit$lambda
  } else {
    c(fit$lambda1, fit$lambda2)
  }

  spec <- list(
    X_t     = X_t, y_t = y_t,
    X_s     = X_s, y_s = y_s,
    w_t     = w_t,
    w_s     = w_s,
    n_t     = nrow(X_t), n_s = if (is.null(X_s)) 0L else nrow(X_s),
    model   = model,
    alpha   = alpha, alpha1 = alpha1, alpha2 = alpha2,
    lam     = lam_num,
    standardize = standardize, cox.ties = cox.ties
  )

  backend <- start_parallel_backend(ncore)
  if (backend$owned) on.exit(stop_parallel_backend(backend), add = TRUE)
  `%par%` <- if (backend$active) foreach::`%dopar%` else foreach::`%do%`

  reps <- foreach(rep_k = 1:nboot, .export = "bootstrap_replicate") %par% {
    bootstrap_replicate(rep_k, seeds, spec)
  }

  n_failed <- sum(vapply(reps, is.null, logical(1)))
  if (n_failed > 0) {
    warning(sprintf("%d of %d bootstrap replicate(s) failed and were dropped.",
                    n_failed, nboot))
  }

  mats <- reps[!vapply(reps, is.null, logical(1))]
  if (length(mats) == 0) {
    stop("All bootstrap replicates failed for model '", model, "'.")
  }
  Bm <- do.call(cbind, mats)

  a <- 1 - conf.level
  ci <- list(
    se    = apply(Bm, 1, sd, na.rm = TRUE),
    lower = apply(Bm, 1, quantile, probs = a / 2, na.rm = TRUE),
    upper = apply(Bm, 1, quantile, probs = 1 - a / 2, na.rm = TRUE)
  )
  nm <- names(beta)
  ci <- lapply(ci, function(z) { names(z) <- nm; z })

  res <- list(beta = beta, ci = ci)
  attr(res, "fit_info")   <- fit_info
  attr(res, "bootstrap")  <- TRUE
  attr(res, "conf.level") <- conf.level
  attr(res, "nboot")      <- nboot
  class(res) <- "fit_coxrtl"
  res
}

# --------------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------------

#' Refit the Target-Only Model at a Fixed Lambda (Internal)
#'
#' Refits the unpenalized-path penalized Cox model of the \code{nontransfer}
#' model fitted by \code{\link{no_transfer_learning}} but \emph{without}
#' cross-validation, at a single fixed numeric \code{lambda}. Used inside each
#' bootstrap replicate of \code{\link{fit_coxrtl}}.
#'
#' @param X,y,alpha,w As in \code{\link{no_transfer_learning}} (a predictor
#'   matrix, a \code{Surv} response, the mixing parameter and optional weights).
#' @param lam A single numeric lambda value at which to refit.
#' @param standardize,cox.ties As in \code{\link{no_transfer_learning}}. Kept
#'   identical to the full-data fit so the refit is on the same scale.
#'
#' @return A named numeric vector of length \code{ncol(X)}.
#'
#' @keywords internal
# --------------------------------------------------------------------------------

refit_baseline_fixed <- function(X, y, alpha, lam, w = NULL,
                                 standardize = FALSE, cox.ties = "breslow") {
  fit <- if (is.null(w)) {
    glmnet(X, y, family = "cox", standardize = standardize, alpha = alpha,
           cox.ties = cox.ties, lambda = lam)
  } else {
    glmnet(X, y, family = "cox", standardize = standardize, alpha = alpha,
           cox.ties = cox.ties, lambda = lam, weights = w)
  }
  stats::coef(fit)[, 1]
}

#' Refit the Two-Step Transfer-Learning Model at Fixed Lambdas (Internal)
#'
#' Refits the two-step transfer-learning estimator of
#' \code{\link{transfer_learning}} \emph{without} cross-validation, at single
#' fixed numeric \code{lambda1}/\code{lambda2}. Used inside each bootstrap
#' replicate of \code{\link{fit_coxrtl}}.
#'
#' @param X_t,X_s,y_t,y_s,alpha1,alpha2,w_t,w_s As in
#'   \code{\link{transfer_learning}}.
#' @param lam1,lam2 Single numeric lambda values for the two steps.
#' @param standardize,cox.ties As in \code{\link{transfer_learning}}. Kept
#'   identical to the full-data fit so the refit is on the same scale.
#'
#' @return A named numeric vector of length \code{ncol(X_t)}.
#'
#' @keywords internal
# --------------------------------------------------------------------------------

refit_transfer_fixed <- function(X_t, y_t, X_s, y_s,
                                 alpha1, alpha2, lam1, lam2,
                                 w_t = NULL, w_s = NULL,
                                 standardize = FALSE, cox.ties = "breslow") {
  if (is.null(w_t)) w_t <- rep(1, nrow(X_t))
  if (is.null(w_s)) w_s <- rep(1, nrow(X_s))

  # Step 1: weighted penalized Cox on the pooled target + external data
  pen_com <- glmnet(
    x = rbind(X_t, X_s), y = c(y_t, y_s),
    standardize = standardize, alpha = alpha1, family = "cox",
    cox.ties = cox.ties, weights = c(w_t, w_s), lambda = lam1
  )
  coef_com <- stats::coef(pen_com)[, 1]

  # Step 2: target-specific offset-based debiasing
  offset <- as.matrix(X_t %*% coef_com)
  pen_d <- glmnet(
    x = X_t, y = y_t,
    standardize = standardize, alpha = alpha2, family = "cox",
    cox.ties = cox.ties, offset = offset, weights = w_t, lambda = lam2
  )
  delta <- stats::coef(pen_d)[, 1]

  coef_com + delta
}

#' One Bootstrap Replicate of a \code{fit_coxrtl} Model Fit (Internal)
#'
#' Resamples the target and external subjects with replacement and refits the
#' single model stored in \code{spec} at the full-data lambda values. Used
#' inside the bootstrap loop of \code{\link{fit_coxrtl}}; one replicate is
#' processed by one parallel worker when \code{ncore > 1}.
#'
#' @param b Integer index of the replicate (1-based).
#' @param seeds Integer vector of per-replicate seeds drawn in the master.
#' @param spec A list fully describing one replicate: matrices/response/weights
#'   of the target and external data, the requested model and its fixed numeric
#'   lambdas (see the body of \code{\link{fit_coxrtl}}).
#'
#' @return A named numeric vector of coefficients, or \code{NULL} if the
#'   replicate fails.
#'
#' @keywords internal
# --------------------------------------------------------------------------------

bootstrap_replicate <- function(b, seeds, spec) {
  tryCatch({
    set.seed(seeds[b])
    it <- sample.int(spec$n_t, spec$n_t, replace = TRUE)
    X_tb <- spec$X_t[it, , drop = FALSE]
    y_tb <- spec$y_t[it]
    w_tb <- if (is.null(spec$w_t)) NULL else spec$w_t[it]

    if (spec$n_s > 0L) {
      is_ <- sample.int(spec$n_s, spec$n_s, replace = TRUE)
      X_sb <- spec$X_s[is_, , drop = FALSE]
      y_sb <- spec$y_s[is_]
      w_sb <- if (is.null(spec$w_s)) NULL else spec$w_s[is_]
    } else {
      X_sb <- spec$X_s
      y_sb <- spec$y_s
      w_sb <- spec$w_s
    }

    m   <- spec$model
    lam <- spec$lam
    if (m == "nontransfer") {
      # w_tb is already NULL when the target data was not reweighted, and
      # refit_baseline_fixed() turns NULL into unit weights, so the resampled
      # fit mirrors the full-data fit in both weighting states.
      refit_baseline_fixed(X_tb, y_tb, spec$alpha, lam, w_tb,
                           standardize = spec$standardize,
                           cox.ties = spec$cox.ties)
    } else {
      # w_tb / w_sb are already NULL on whichever side was not reweighted, and
      # refit_transfer_fixed() turns NULL into unit weights, so the resampled
      # fit mirrors the full-data fit in every weighting combination.
      refit_transfer_fixed(X_tb, y_tb, X_sb, y_sb,
                           spec$alpha1, spec$alpha2, lam[1], lam[2],
                           w_tb, w_sb,
                           standardize = spec$standardize,
                           cox.ties = spec$cox.ties)
    }
  }, error = function(e) NULL)
}
