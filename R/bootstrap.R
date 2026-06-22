bootstrap <- function(model,
                      zero.tol = 1e-10,
                      verbose  = model@info$verbose,
                      parallel = model@info$boot$parallel,
                      ncores   = model@info$boot$ncores,
                      R        = model@info$boot$R,
                      iseed    = model@info$boot$iseed,
                      drop.inadmissible = isTRUE(model@info$boot$drop.inadmissible)) {

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

  # The base fit may itself be inadmissible. Reset the admissibility status on
  # the template such that the bootstraps don't inherit the inadmissibility flag
  isAdmissible(baseModel, recursive = TRUE) <- TRUE

  boot.optimize <- isTRUE(model@info$boot$optimize)
  mc.boot.control <- prepMCBootControl(
    boot.control  = model@info$boot$mc.boot.control,
    boot.optimize = boot.optimize,
    model         = model
  )

  combinedCoefs <- combinedModel(model)@params$values

  probsTemplate <- model@thresholdStruct@proportions
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

      model.b@thresholdStruct <- updateThresholds(updateProportions(
        thr  = baseModel@thresholdStruct,
        data = sampleData
      ))

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

      par <- formatBootPars(modelParams(combinedModel(model.b))$values)
      probs <- model.b@thresholdStruct@proportions

      if (low.tol.penalty > 0 && is.mcpls && !mc.delta.se) {
        k <- length(par)
        par <- par + stats::rnorm(k, mean = 0, sd = low.tol.penalty)
      }

      inadmissible <- !isAdmissible(model.b)

      if (inadmissible && drop.inadmissible) {
        # drop the replicate by treating it like a failed fit; the all-NA row
        # is excluded from the (co-)variances via `use = "complete.obs"` below
        out <- c(parTemplate, probsTemplate + NA_real_)
      } else {
        out <- c(par, probs)
      }

      attr(out, "id") <- i
      attr(out, "inadmissible") <- inadmissible
      out

    }, error = \(e) {
      pls_msg_warn(
        paste0("Bootstrap replicate ", i, " failed: ", conditionMessage(e))
      )

      out <- c(parTemplate, probsTemplate + NA_real_)
      attr(out, "id") <- i
      attr(out, "inadmissible") <- NA # errored before admissibility is known
      out
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

  inadmissible <- vapply(
    results, FUN.VALUE = logical(1L),
    FUN = \(x) isTRUE(attr(x, "inadmissible"))
  )
  n.inadmissible <- sum(inadmissible)

  if (n.inadmissible) {
    verb <- if (drop.inadmissible) "Dropped" else "Kept"
    pls_msg_warn(sprintf(
      "%s %d (out of %d) bootstrap replicate(s) with inadmissible solutions.",
      verb, n.inadmissible, R
    ))
  }

  pls_warnif(drop.inadmissible && (R - n.inadmissible) < 2L,
    "Fewer than 2 admissible bootstrap replicates remain after dropping",
    "inadmissible solutions; standard errors cannot be estimated reliably."
  )

  resultsMat <- do.call(rbind, results)
  rownames(resultsMat) <- ids
  colnames(resultsMat) <- c(names(combinedModel(model)@params$values), names(probsTemplate))

  vcov.joint <- stats::cov(resultsMat, use = "complete.obs")
  par.names <- names(combinedModel(model)@params$values)
  vcov <- vcov.joint[par.names, par.names, drop = FALSE]

  if (mc.delta.se) {
    combined <- combinedModel(model)
    params <- modelParams(combined)

    if (!is.null(params$Jacobian0)) {
      Jacobian0 <- params$Jacobian0
      pars.free <- intersect(colnames(Jacobian0), colnames(vcov.joint))

      J0 <- Jacobian0[pars.free, pars.free, drop = FALSE]

      tryCatch({
        J0.inv <- invertMcJacobian(J0)

        if (!is.null(params$Jacobian1)) {
          Jacobian1 <- params$Jacobian1

          pars.all <- intersect(rownames(Jacobian1), colnames(vcov))
          missing0 <- setdiff(colnames(vcov), rownames(Jacobian1))
          missing1 <- setdiff(pars.free, colnames(Jacobian1))

          pls_warnif(length(missing0),
            "Missing entries in Jacobian for some parameters:",
            paste0(missing0, collapse = ", ")
          )

          pls_stopif(length(missing1), # this really shouldn't happen!
            "Missing entries in Jacobian for some free parameters:",
            paste0(missing1, collapse = ", ")
          )

          J1 <- Jacobian1[pars.all, pars.free, drop = FALSE]
          D.par <- J1 %*% J0.inv
          prob.names <- intersect(names(probsTemplate), colnames(vcov.joint))

          if (length(prob.names)) {
            # Implicit delta method:
            # dp = J0^-1 da - J0^-1 Jp dc
            # dy = J1 dp + Gp dc
            Jp <- params$JacobianProbs0[pars.free, prob.names, drop = FALSE]
            Gp <- params$JacobianProbs1[pars.all, prob.names, drop = FALSE]
            D.probs <- Gp - D.par %*% Jp
            D <- cbind(D.par, D.probs)
            vcov.sub <- vcov.joint[
              c(pars.free, prob.names), c(pars.free, prob.names), drop = FALSE
            ]

          } else {
            D <- D.par
            vcov.sub <- vcov.joint[pars.free, pars.free, drop = FALSE]
          }

          vcov.mc.full <- D %*% vcov.sub %*% t(D)

          vcov[] <- 0
          vcov[pars.all, pars.all] <- vcov.mc.full[pars.all, pars.all]

        } else {
          # Just use standard errors for free parameters
          vcov.sub <- vcov.joint[pars.free, pars.free, drop = FALSE]
          vcov.mc.free <- J0.inv %*% vcov.sub %*% t(J0.inv)
          vcov[] <- 0
          vcov[pars.free, pars.free] <- vcov.mc.free[pars.free, pars.free]

        }

      }, error = function(e) {
        pls_msg_warn(
          "Calculation of delta standard errors failed!",
          "Consider trying `mc.delta.se=FALSE` instead.",
          "Message:", conditionMessage(e)
        )

        vcov[] <<- NA_real_
      })

      # The (co-)variances structure of the interaction terms, is often
      # not a function of the free model parameters (though some are).
      # E.g., the variance of a quadratic terms should analytically always be 2,
      # and X~X:Z should always be zero. The delta-method should therefore not
      # add any information here. However, it might seem like it, due to sampling error.
      # Here we just remove these from the vcov
      pars   <- rownames(vcov)
      is.cov <- grepl("~~", pars)
      is.int <- grepl(":", pars)

      split  <- stringr::str_split_fixed(
        pars[is.cov & is.int], pattern = "~~", n = 2
      )

      if (NROW(split)) {
        split.sub <- split[
          !grepl("~", split[,1L]) & !grepl("~", split[,2L]) & # remove random effect variances
          !is.na(split[,1L])      & !is.na(split[,2L]), , drop = FALSE
        ]

        if (NROW(split.sub)) {
          rm <- paste0(split.sub[,1L], "~~", split.sub[,2L])
          vcov[rm,] <- 0
          vcov[,rm] <- 0
        }
      }

    }
  }

  se <- sqrt(diag(vcov))
  se[se <= zero.tol] <- NA_real_

  list(se = se, boot = resultsMat[, par.names, drop = FALSE], vcov = vcov)
}


invertMcJacobian <- function(J, rcond.tol = 1e-10) {
  condition <- tryCatch(rcond(J), error = function(e) 0)

  if (!is.finite(condition) || condition < rcond.tol) {
    pls_msg_warn(
      "MC-PLS root-equation Jacobian is ill-conditioned; ",
      "using a generalized inverse for delta-method standard errors."
    )
    return(MASS::ginv(J))
  }

  tryCatch(
    solve(J),
    error = function(e) {
      pls_msg_warn(
        "Unable to invert the MC-PLS root-equation Jacobian; ",
        "using a generalized inverse for delta-method standard errors."
      )
      MASS::ginv(J)
    }
  )
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
