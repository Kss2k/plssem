parseParallelArgs <- function(parallel, ncores) {
  if (is.null(parallel))
    parallel <- "no"

  parallel <- match.arg(
    tolower(parallel), c("no", "multicore", "multisession", "snow")
  )

  if (parallel == "snow")
    parallel <- "multisession"

  if (is.null(ncores) || !length(ncores) || !is.finite(ncores[[1L]]))
    ncores <- 1L

  ncores <- as.integer(ncores[[1L]])
  if (parallel == "no") workers <- 1L
  else                  workers <- ncores

  if (workers > 1L && parallel == "multicore" && .Platform$OS.type == "windows") {
    pls_msg_warn(paste0(
      "The `multicore` option is not supported on Windows.\n",
      "Falling back to `multisession`."
    ))

    parallel <- "multisession"
  }

  list(
    parallel = parallel,
    ncores   = ncores,
    workers  = workers
  )
}


plapply <- function(X,
                    FUN,
                    parallel = "no",
                    ncores   = 1L,
                    verbose  = FALSE,
                    iseed    = NULL,
                    label    = "Task",
                    packages = "plssem") {

  args <- parseParallelArgs(parallel, ncores)
  n    <- length(X)

  if (args$workers <= 1L) {
    if (!verbose) return(lapply(X, FUN))

    pb <- utils::txtProgressBar(
      min     = 0,
      max     = n,
      initial = 0,
      style   = 3,
      file    = stderr()
    )

    on.exit(close(pb), add = TRUE)

    return(lapply(seq_len(n), function(i) {
      tryCatch(
        utils::setTxtProgressBar(pb, i),
        error = \(e) pls_msg_warn(
          "Unable to update progress bar!", 
          "Message:", conditionMessage(e)
        )
      )

      FUN(X[[i]])
    }))
  }

  oldPlan <- future::plan()
  on.exit(future::plan(oldPlan), add = TRUE)

  if (args$parallel == "multicore") {
    future::plan(future::multicore, workers = args$workers)
  } else {
    future::plan(future::multisession, workers = args$workers)
  }

  # `future.seed = FALSE` warns whenever a worker touches the RNG, which all of
  # our tasks do; fall back to `TRUE` (proper, but non-reproducible, streams)
  future.seed <- if (is.null(iseed)) TRUE else iseed

  if (!verbose) {
    return(future.apply::future_lapply(
      X               = X,
      FUN             = FUN,
      future.seed     = future.seed,
      future.packages = packages
    ))
  }

  oldHandlers <- progressr::handlers()
  on.exit(progressr::handlers(oldHandlers), add = TRUE)
  progressr::handlers(progressr::handler_txtprogressbar(
    file  = stderr(),
    style = 3L
  ))

  progressr::with_progress({
    p <- progressr::progressor(along = seq_len(n))

    future.apply::future_lapply(
      X   = seq_len(n),
      FUN = function(i) {
        p(sprintf("%s %d/%d", label, i, n))
        FUN(X[[i]])
      },
      future.seed     = future.seed,
      future.packages = packages
    )
  })
}
