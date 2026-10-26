# start_parallel_backend Function Documentation ---------------------------------

#' Start a Local Parallel Backend for a Single Call
#'
#' Registers a local \code{doParallel} backend with \code{ncore} workers when
#' \code{ncore > 1} and no multi-worker backend is already registered, so that
#' the parallel calls made by \code{get_theta()}, \code{no_transfer_learning()}
#' and \code{transfer_learning()} stay self-contained within one function run.
#'
#' @param ncore Integer. Requested number of parallel workers. Values
#'   \code{<= 1} disable parallelism.
#'
#' @return A list describing the backend state:
#' \item{\code{active}}{Logical; whether the caller should run in parallel
#' (\code{ncore > 1}).}
#' \item{\code{owned}}{Logical; whether this call created the backend and must
#' shut it down again.}
#' \item{\code{cl}}{The created cluster object, present only when
#' \code{owned = TRUE}.}
#'
#' @details
#' When a multi-worker backend is already registered (for example by the user
#' via \code{doParallel::registerDoParallel()}), it is reused untouched and
#' \code{owned = FALSE}. Otherwise a cluster of \code{ncore} PSOCK workers is
#' created and registered for this call; \code{\link{stop_parallel_backend}}
#' must then be called (typically in \code{on.exit()}) to shut it down before
#' the function returns.
#'
#' @keywords internal
# --------------------------------------------------------------------------------

start_parallel_backend <- function(ncore) {
  use_parallel <- is.numeric(ncore) && length(ncore) == 1L &&
    !is.na(ncore) && ncore > 1
  if (!use_parallel) {
    return(list(active = FALSE, owned = FALSE))
  }
  if (foreach::getDoParWorkers() > 1) {
    # A parallel backend is already registered; reuse it and leave it alone.
    return(list(active = TRUE, owned = FALSE))
  }
  cl <- parallel::makeCluster(ncore)
  doParallel::registerDoParallel(cl)
  list(active = TRUE, owned = TRUE, cl = cl)
}

# stop_parallel_backend Function Documentation ---------------------------------

#' Stop a Local Parallel Backend Started by start_parallel_backend()
#'
#' Stops the cluster and resets \pkg{foreach} to the sequential backend when the
#' backend was created by \code{\link{start_parallel_backend}} (i.e.
#' \code{owned = TRUE}). For reused (non-owned) backends this is a no-op, leaving
#' the user's backend untouched.
#'
#' @param backend A list returned by \code{\link{start_parallel_backend}}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @details
#' Safe to call from \code{on.exit()}. Stopping is wrapped in \code{try()} so
#' that cleanup continues even if the cluster has already become unusable.
#'
#' @keywords internal
# --------------------------------------------------------------------------------

stop_parallel_backend <- function(backend) {
  if (isTRUE(backend$owned)) {
    try(parallel::stopCluster(backend$cl), silent = TRUE)
    foreach::registerDoSEQ()
  }
  invisible(NULL)
}
