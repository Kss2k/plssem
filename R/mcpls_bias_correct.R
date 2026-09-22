# One nested Newton step for the finite-sample bias of the MC-PLS estimator.
#
# MC-PLS solves f(X) = f(g(theta.hat)), i.e. theta.hat = b^-1(theta.hat*) with
# b = f o g the binding function. Both the curvature of b^-1 and the sampling
# variance of the auxiliary estimates theta.hat* are non-negligible in finite
# samples, so
#
#   E[b^-1(theta.hat*)] != b^-1(E[theta.hat*]),
#
# and the estimator carries a bias of order 1/n even though it is consistent.
# The bias is a smooth function of theta, and we can simulate from theta, so
# we estimate it by simulation and take one Newton step over it:
#
#   bias.hat = mean_b(theta.hat_b) - theta.hat,  theta.hat_b = estimator(g(theta.hat))
#   theta.bc = theta.hat - bias.hat = 2 * theta.hat - mean_b(theta.hat_b)
#
# This is the same Monte-Carlo consistency argument the estimator already
# makes about the PLS estimator, applied once more to MC-PLS itself. The inner
# loop mirrors `mcplsLoglik()`, except that each simulated data set is passed
# through the *whole* estimator (`estimatePLS()`) rather than through the
# auxiliary stage alone (`estimatePLS_Inner()`).
mcplsBiasCorrect <- function(object,
                             B           = 50L,
                             seed        = NULL,
                             reuse.start = FALSE,
                             clamp       = TRUE,
                             verbose     = interactive()) {
  combined <- combinedModel(object)

  pls_stopif(!is_mcpls(combined),
    "`mcpls_bias_correct()` is only implemented for MC-PLS models!"
  )

  B <- suppressWarnings(as.integer(B))
  pls_stopif(length(B) != 1L || is.na(B) || B < 2L,
    "`B` must be a single integer greater than 1."
  )

  pls_warnif(!isAdmissible(object),
    "The fitted model is inadmissible!",
    "The bias is estimated at an implausible parameter vector."
  )

  pls_warnif(B < 25L,
    sprintf("`B = %d` is small.", B),
    "The Monte-Carlo error of the correction falls as 1/sqrt(B);",
    "inspect the returned `se` before trusting the corrected estimates."
  )

  if (!is.null(seed)) {
    if (exists(".Random.seed")) {
      .Random.seed.orig <- .Random.seed
      on.exit({.Random.seed <<- .Random.seed.orig}, add = TRUE)
    }

    set.seed(seed)
  }

  # theta.hat. Only the parameters the root finder actually solves for are
  # corrected: residual variances and thresholds are functions of these, and
  # stepping them separately would break the constraints that tie them
  # together. They are recomputed by the estimator instead.
  par0  <- getFreeParamsTable(combined)
  free0 <- par0[par0$is.free, , drop = FALSE]
  nm0   <- getParNamesFromParTable(free0)
  theta <- stats::setNames(free0$est, nm = nm0)

  # g(): the estimator's own simulation scheme, at the observed sample size,
  # cutting at the estimated thresholds.
  parTable <- getParTableEstimates(
    combined, rm.tmp.ov = FALSE, clean.tmp.ind = FALSE
  )

  data <- modelData(object)
  vars <- colnames(data)
  n    <- NROW(data)

  if (isMLM(object)) {
    clusterSizes <- as.numeric(table(attr(data, "cluster")))
    clusterName  <- colnames(attr(data, "cluster"))

  } else {
    clusterSizes <- NULL
    clusterName  <- NULL

  }

  estimator <- combined@info$path.estimator
  mc.rescov <- object@info$mc.args$rescov

  use.full.rescov <- switch(mc.rescov,
    full    = TRUE,
    reduced = FALSE,
    auto    = estimator == "gls",
    pls_msg_stop("Unrecognized value for `mc.rescov` argument:", mc.rescov)
  )

  is.hi.ord <- isTRUE(combined@info$is.high.ord)

  # Each replicate is specified from scratch, exactly as a fresh `pls()` call
  # would, so that no state from the original fit leaks into it. Bootstrapping
  # and printing are switched off; everything else is inherited.
  spec.args <- object@info$spec.args

  pls_stopif(is.null(spec.args),
    "The fitted model does not carry the arguments needed to re-specify it.",
    "It was probably fitted with an older version of the package; refit it",
    "with `pls()` before calling `mcpls_bias_correct()`."
  )

  spec.args$verbose   <- FALSE
  spec.args$bootstrap <- FALSE

  # By default each replicate starts where a fresh `pls()` call would, at its
  # own auxiliary estimates, rather than at theta.hat. Starting at theta.hat
  # is faster, but it shrinks the estimated bias towards zero.
  p.start <- if (reuse.start) object@info$mc.args$p.start else NULL

  .f <- function(i, pb = NULL) {

    if (!is.null(pb)) {
      tryCatch(
        utils::setTxtProgressBar(pb, i),
        error = \(e) pls_msg_warn(
          "Unable to update progress bar!\nMessage:", conditionMessage(e)
        )
      )
    }

    sim <- simulateDataParTable(
      parTable     = parTable,
      N            = n,
      check.hi.ord = is.hi.ord,
      clusterSizes = clusterSizes,
      clusterName  = clusterName,
      full         = use.full.rescov,
      cut          = TRUE # cut from estimated thresholds
    )

    pls_stopif(!sim$is.admissible,
      "The estimated parameters do not imply an admissible data generating process!"
    )

    data.b <- sim$ov[vars]

    if (!is.null(sim$cluster))
      data.b <- cbind(data.b, sim$cluster)

    utils::capture.output(type = "message", { # capture real time output
      model.b <- do.call(specifyModel, c(list(data = data.b), spec.args))

      model.b <- suppressWarnings(estimatePLS(
        model   = model.b,
        # args passed onto mcpls
        verbose = FALSE,
        p.start = p.start
      ))
    })

    pls_stopif(!isAdmissible(model.b),
      "The replicate produced an inadmissible solution!"
    )

    par.b <- getFreeParamsTable(combinedModel(model.b))
    par.b[par.b$is.free, "est"]
  }

  X <- matrix(NA_real_, nrow = B, ncol = NROW(free0))
  colnames(X) <- nm0

  if (verbose) {
    pb <- utils::txtProgressBar(
      min     = 0,
      max     = B,
      initial = 0,
      style   = 3,
      file    = stderr()
    )

    on.exit(close(pb), add = TRUE)

  } else {
    pb <- NULL

  }

  for (i in seq_len(B)) {

    X[i,] <- tryCatch(.f(i = i, pb = pb), error = function(e) {
      pls_msg_warn(
        sprintf("Replicate %d failed! Message:", i), conditionMessage(e)
      )
      NA_real_
    })

  }

  complete   <- stats::complete.cases(X)
  X.complete <- X[complete, , drop = FALSE]
  ncomplete  <- sum(complete)

  pls_stopif(ncomplete < 2L,
    "Unable to estimate the bias of the MC-PLS estimator.",
    sprintf("Only %d usable replicate(s) out of %d remain.", ncomplete, B),
    "Increase `B` or inspect the warnings from the failed replicates."
  )

  pls_warnif(ncomplete < B,
    sprintf("%d of %d replicates were discarded.", B - ncomplete, B),
    "The remaining replicates are a selected sample, so the estimated bias",
    "is itself biased towards the admissible region."
  )

  # One Newton step over the simulated bias.
  bias      <- colMeans(X.complete) - theta
  se        <- apply(X.complete, MARGIN = 2L, FUN = stats::sd) / sqrt(ncomplete)
  corrected <- theta - bias

  if (clamp) {
    lower <- getMcLowerBounds(par0)
    upper <- getMcUpperBounds(par0)
    out   <- corrected < lower | corrected > upper

    pls_warnif(any(out),
      "Bias-corrected estimates outside the bounds of the parameter space",
      "were clamped:", paste(nm0[out], collapse = ", ")
    )

    corrected <- pmin(pmax(corrected, lower), upper)
  }

  list(
    est        = plssemVector(theta),
    bias       = plssemVector(bias),
    se         = plssemVector(se),
    est.bc     = plssemVector(corrected),
    replicates = X.complete,
    B          = ncomplete
  )
}
