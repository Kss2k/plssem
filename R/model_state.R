initModelStatus <- function(tolerance, max.iter.0_5) {
  list(
    convergence    = FALSE,
    iterations     = 0L,
    iterations.0_5 = 0L,
    tolerance      = tolerance,
    max.iter.0_5   = max.iter.0_5,
    is.admissible  = TRUE,
    mcpls.update.args = NULL
  )
}


resetModelStatusLowerOrder <- function(model, hard.reset = FALSE) {
  model@status$convergence <- FALSE
  model@status$iterations.0_5 <- 0L

  if (hard.reset) {
    model@status$iterations <- 0L
    model@params$values.old <- NULL
  }

  model
}


updateModelInfo <- function(model, ...) {
  updates <- list(...)
  if (!length(updates))
    return(model)

  info <- model@info
  for (name in names(updates)) {
    info[[name]] <- updates[[name]]
  }

  info$estimator <- getEstimatorFromInfo(info)
  model@info <- info
  model
}


initModelParams <- function(model) {
  parnames <- getParamVecNames(model)
  k <- length(parnames)

  model@params <- list(
    names      = parnames,
    values     = rep(NA_real_, k),
    values.old = NULL,
    se         = rep(NA_real_, k),
    vcov       = NULL
  )

  model
}


initModelMcArgs <- function(min.iter,
                            max.iter,
                            mc.reps,
                            tol,
                            fixed.seed,
                            polyak.juditsky,
                            fn.args) {
  list(
    min.iter        = min.iter,
    max.iter        = max.iter,
    mc.reps         = mc.reps,
    tol             = tol,
    fixed.seed      = fixed.seed,
    polyak.juditsky = polyak.juditsky,
    fn.args         = fn.args,
    rng.seed        = NULL,
    p.start         = NULL
  )
}


initModelBootInfo <- function(bootstrap,
                              ncpus,
                              parallel,
                              R,
                              iseed,
                              optimize,
                              mc.boot.control) {
  list(
    bootstrap       = bootstrap,
    ncpus           = ncpus,
    parallel        = parallel,
    R               = R,
    iseed           = iseed,
    optimize        = optimize,
    mc.boot.control = mc.boot.control
  )
}


initModelInfo <- function(baseInfo,
                          parsed,
                          n,
                          ordered,
                          consistent,
                          verbose,
                          standardize,
                          reliabilities,
                          is.lower.order,
                          mc.args,
                          boot) {
  stopifnot(is.list(baseInfo), is.list(parsed))

  info <- baseInfo

  inds.x <- info$inds.x
  inds.y <- info$inds.y

  ordered.x <- intersect(inds.x, ordered)
  ordered.y <- intersect(inds.y, ordered)

  info$lme4.syntax    <- parsed$lme4.syntax
  info$is.mlm         <- parsed$is.mlm
  info$is.mcpls       <- parsed$is.mcpls
  info$is.probit      <- parsed$is.probit
  info$cluster        <- parsed$cluster
  info$consistent     <- consistent
  info$ordered        <- ordered
  info$ordered.x      <- ordered.x
  info$ordered.y      <- ordered.y
  info$intTermElems   <- parsed$intTermElems
  info$intTermNames   <- parsed$intTermNames
  info$is.nlin        <- parsed$is.nlin
  info$rng.seed       <- floor(stats::runif(1L, min = 0, max = 9999999))
  info$n              <- n
  info$estimator      <- getEstimatorFromInfo(info)
  info$verbose        <- verbose
  info$standardized   <- standardize
  info$reliabilities  <- reliabilities
  info$is.high.ord    <- FALSE
  info$is.lower.order <- isTRUE(is.lower.order)
  info$mc.args        <- mc.args
  info$boot           <- boot

  info
}
