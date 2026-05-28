bootstrap <- function(model,
                      zero.tol = 1e-10,
                      verbose  = model@info$verbose,
                      parallel = model@info$boot$parallel,
                      ncores   = model@info$boot$ncores,
                      R        = model@info$boot$R,
                      iseed    = model@info$boot$iseed) {

  if (is.null(parallel)) parallel <- "no"
  parallel <- match.arg(parallel, c("no", "multicore", "multisession", "snow"))
  if (parallel == "snow") parallel <- "multisession"

  data        <- model@data
  cluster     <- model@info$cluster
  is.probit   <- model@info$is.probit
  is.mcpls    <- model@info$is.mcpls
  ordered     <- model@info$ordered
  results     <- vector("list", R)
  mc.delta.se <- model@info$mc.args$delta.se

  # Check misspecified user arguments
  pls_warnif(
    parallel != "no" && ncores <= 1L,
    "The `boot.ncores` argument has to be larger than 1 for parallel\n",
    "bootstrapping to be enabled! Try `boot.ncores = 2`"
  )

  pls_warnif(
    parallel == "no" && ncores > 1L,
    'The `boot.parallel` must be "multicore" or "multisession" for parallel\n',
    'bootstrapping to be enabled! Try `boot.parallel="multisession"`'
  )

  baseModel <- model

  if (mc.delta.se && is.mcpls) {
    combinedModel(baseModel) <- NULL
    isMCPLS(baseModel, recursive = TRUE) <- FALSE
  }

  boot.optimize <- isTRUE(model@info$boot$optimize)
  mc.boot.control <- prepMCBootControl(
    boot.control  = model@info$boot$mc.boot.control,
    boot.optimize = boot.optimize,
    model         = model
  )

  combinedCoefs <- combinedModel(model)@params$values
  # parTemplate might be stale with mc.delta=TRUE...
  
  parTemplate <- stats::setNames(
    rep(NA_real_, length(combinedCoefs)),
    names(combinedCoefs)
  )

  formatBootPars <- function(par) {
    out <- parTemplate

    want <- names(combinedCoefs)
    have <- names(par)
    common <- intersect(have, want)

    if (length(common))
      out[common] <- par[common]

    out
  }

  P_START <- model@info$mc.args$p.start

  .bootf <- function(i) {
    tryCatch({
      sampleData <- resample(data, cluster = cluster)
      sampleS    <- getCorrMat(sampleData, ordered = ordered, probit = is.probit)
      model.b    <- baseModel

      model.b@data       <- sampleData
      model.b@matrices$S <- sampleS

      mc.args             <- model.b@info$mc.args
      boot.fixed.seed     <- mc.boot.control$fixed.seed
      boot.polyak         <- mc.boot.control$polyak.juditsky
      boot.pj.extrapolate <- mc.boot.control$pj.extrapolate
      boot.tol            <- mc.boot.control$tol
      boot.reps           <- mc.boot.control$mc.reps
      boot.p.start        <- if (mc.boot.control$reuse.p.start) P_START else NULL
      boot.verbose        <- mc.boot.control$verbose
      low.tol.penalty     <- mc.boot.control$low.tol.penalty

      utils::capture.output(type = "message", { # capture real time output

        model.b <- suppressWarnings(estimatePLS(
          model = model.b,
          # args passed onto mcpls
          fixed.seed      = boot.fixed.seed,
          verbose         = boot.verbose,
          polyak.juditsky = boot.polyak,
          pj.extrapolate  = boot.pj.extrapolate,
          tol             = boot.tol,
          mc.reps         = boot.reps,
          p.start         = boot.p.start
        ))

      })

      par <- formatBootPars(
        modelParams(combinedModel(model.b))$values
      )

      if (low.tol.penalty > 0 && is.mcpls && !mc.delta.se) {
        k <- length(par)
        par <- par + stats::rnorm(k, mean = 0, sd = low.tol.penalty)
      }

      attr(par, "id") <- i
      par

    }, error = \(e) {
      pls_msg_warn(
        paste0("Bootstrap replicate ", i, " failed: ", conditionMessage(e))
      )

      par <- parTemplate
      attr(par, "id") <- i
      par
    })
  }

  # iseed:
  # this is adapted from the lavaan package
  # - iseed is used for both serial and parallel
  # - if iseed is not set, iseed is generated + .Random.seed created/updated
  #     -> tmp.seed <- NA
  # - if iseed is set: don't touch .Random.seed (if it exists)
  #     -> tmp.seed <- .Random.seed (if it exists)
  #     -> tmp.seed <- NULL (if it does not exist)

  if (is.null(iseed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      stats::runif(1)
    }

    # identical(temp.seed, NA): Will not change .Random.seed in GlobalEnv
    temp.seed <- NA
    iseed <- stats::runif(1, 0, 999999999)

  } else {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      temp.seed <-
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    } else {
      # is.null(temp.seed): Will remove .Random.seed in GlobalEnv
      #                     if serial.
      #                     If parallel, .Random.seed will not be touched.
      temp.seed <- NULL
    }
  }

  if (verbose) pls_msg_note("Bootstrapping...")

  workers <- if (parallel == "no") 1L else ncores
  if (workers <= 1L) {
    set.seed(iseed)

    if (verbose) {
      pb <- utils::txtProgressBar(
        min     = 0,
        max     = R,
        initial = 0,
        style   = 3,
        file    = stderr()
      )

      on.exit(close(pb), add = TRUE)

      results <- lapply(seq_len(R), function(i) {
        tryCatch(
          utils::setTxtProgressBar(pb, i),
          error = \(e) pls_msg_warn(
            paste0("Unable to update progress bar!\nMessage: ", conditionMessage(e))
          )
        )

        .bootf(i)
      })

    } else {
      results <- lapply(seq_len(R), .bootf)

    }

  } else {
    oldPlan <- future::plan()
    on.exit(future::plan(oldPlan), add = TRUE)

    if (parallel == "multicore" && .Platform$OS.type == "windows") {
      pls_msg_warn(paste0(
        "The `boot.parallel = 'multicore'` option is not supported on Windows.\n",
        "Falling back to `boot.parallel = 'multisession'`."
      ))
      parallel <- "multisession"
    }

    if (parallel == "multicore") {
      future::plan(future::multicore, workers = workers)
    } else {
      future::plan(future::multisession, workers = workers)
    }

    if (verbose) {
      results <- progressr::with_progress({
        oldHandlers <- progressr::handlers()
        on.exit(progressr::handlers(oldHandlers), add = TRUE)
        progressr::handlers(progressr::handler_txtprogressbar(
          file = stderr(),
          style = 3L
        ))
        p <- progressr::progressor(along = seq_len(R))
        future.apply::future_lapply(
          X = seq_len(R),
          FUN = function(i) {
            p(sprintf("Bootstrap %d/%d", i, R))
            .bootf(i)
          },
          future.seed = iseed,
          future.packages = "plssem"
        )
      })

    } else {
      results <- future.apply::future_lapply(
        X = seq_len(R),
        FUN = .bootf,
        future.seed = iseed,
        future.packages = "plssem"
      )
    }
  }

  ids <- vapply(results, FUN.VALUE = integer(1L), FUN = \(x) attr(x, "id"))

  resultsMat <- do.call(rbind, results)
  rownames(resultsMat) <- ids
  colnames(resultsMat) <- names(combinedModel(model)@params$values)

  vcov <- stats::cov(resultsMat, use = "complete.obs")

  if (mc.delta.se) {
    combined <- combinedModel(model)
    params <- modelParams(combined)

    if (!is.null(params$Jacobian)) {
      Jacobian <- params$Jacobian
      pars <- intersect(colnames(Jacobian), colnames(vcov))

      vcov.sub <- vcov[pars, pars, drop = FALSE]
      J <- Jacobian[pars, pars, drop = FALSE]

      tryCatch({
        # Try to invert J
        J.inv <- tryCatch(
          solve(J),
          error = function(e) {
            pls_msg_warn("Jacobian is not positive definite!")
            MASS::ginv(J)
          }
        )

        vcov.mc <- J.inv %*% vcov.sub %*% t(J.inv)
        
        vcov[] <- 0
        vcov[pars, pars] <- vcov.mc[pars, pars]

      }, error = function(e) {
        pls_msg_warn(
          "Calculation of delta standard errors failed!",
          "Consider trying `mc.delta.se=FALSE` instead.",
          "Message:", conditionMessage(e)
        )

        vcov[] <<- NA_real_
      })
    }
  }

  se <- sqrt(diag(vcov))
  se[se <= zero.tol] <- NA_real_

  list(se = se, boot = resultsMat, vcov = vcov)
}


resample <- function(X, n.out = NROW(X), cluster = NULL, replace = TRUE) {
  if (is.null(cluster)) {
    idx <- sample(NROW(X), size = n.out, replace = replace)
    return(X[idx, , drop = FALSE])
  }

  pls_stopif(length(cluster) > 1L, "bootstrapping of multiple cluster variables is not implemented (yet)!")

  cluster.vals <- attr(X, "cluster")[, cluster, drop = TRUE]
  pls_stopif(NROW(cluster.vals) != NROW(X), "Cluster must be of same length as data!")

  clusters <- unique(cluster.vals)
  G <- length(clusters)

  clusters.sample <- sample(clusters, size = G, replace = replace)

  cluster.list <- lapply(clusters.sample, FUN = \(ci) X[cluster.vals==ci, , drop=FALSE])
  # create new (fresh) cluster indices
  # if we sample the same cluster twice, we want the model to treat them
  # as different clusters
  indices.list <- lapply(
    X = seq_along(cluster.list),
    FUN = \(idx) matrix(idx, nrow = NROW(cluster.list[[idx]]), ncol = 1L,
                        dimnames = list(NULL, cluster))
  )

  Y <- do.call(rbind, cluster.list)
  attr(Y, "cluster") <- do.call(rbind, indices.list)

  Y
}


CONFIG_BOOT_CONTROL <- list(
  min.iter        = 10L,
  max.iter        = 250L,
  reuse.p.start   = TRUE,
  fixed.seed      = TRUE,
  polyak.juditsky = TRUE,
  pj.extrapolate  = TRUE,
  tol.mult        = 2,
  mc.reps.mult    = 0.5,
  verbose         = FALSE,
  low.tol.penalty = 0
)


prepMCBootControl <- function(boot.control, boot.optimize, model) {
  mc.args <- model@info$mc.args

  if (!boot.optimize) {
    # Overwrite fields by arguments in mc.args
    # I.e., we do not try to optimize the algorithm
    # differently from the main procedure

    boot.control$min.iter        <- mc.args$min.iter
    boot.control$max.iter        <- mc.args$max.iter
    boot.control$polyak.juditsky <- mc.args$polyak.juditsky
    boot.control$pj.extrapolate  <- mc.args$pj.extrapolate
    boot.control$reuse.p.start   <- FALSE
    boot.control$fixed.seed      <- mc.args$fixed.seed
    boot.control$tol             <- mc.args$tol
    boot.control$mc.reps         <- mc.args$mc.reps

  }

  if (is.null(boot.control))
    boot.control <- list()

  if (is.null(boot.control$min.iter))
    boot.control$min.iter <- CONFIG_BOOT_CONTROL$min.iter

  if (is.null(boot.control$max.iter))
    boot.control$max.iter <- CONFIG_BOOT_CONTROL$max.iter

  if (is.null(boot.control$polyak.juditsky))
    boot.control$polyak.juditsky <- CONFIG_BOOT_CONTROL$polyak.juditsky

  if (is.null(boot.control$pj.extrapolate))
    boot.control$pj.extrapolate <- CONFIG_BOOT_CONTROL$pj.extrapolate

  if (is.null(boot.control$reuse.p.start))
    boot.control$reuse.p.start <- CONFIG_BOOT_CONTROL$reuse.p.start

  if (is.null(boot.control$fixed.seed))
    boot.control$fixed.seed <- CONFIG_BOOT_CONTROL$fixed.seed

  if (is.null(boot.control$tol))
    boot.control$tol <- CONFIG_BOOT_CONTROL$tol.mult * mc.args$tol

  if (is.null(boot.control$mc.reps))
    boot.control$mc.reps <- floor(CONFIG_BOOT_CONTROL$mc.reps.mult * mc.args$mc.reps)

  if (is.null(boot.control$verbose))
    boot.control$verbose <- CONFIG_BOOT_CONTROL$verbose

  if (is.null(boot.control$low.tol.penalty)) {
    # Lowering the tolerance while using a warm start lowers the variance
    # We adress this by adding a noise penalty to the parameters.
    # If you're reading this I hope you like magic numbers...
    boot.control$low.tol.penalty <- max(0)#, pi * (boot.control$tol - mc.args$tol))
  }

  boot.control
}
