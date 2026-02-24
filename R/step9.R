estimatePLS_Step9 <- function(model) {
  # Handle ordered thresholds
  ordered    <- model$info$ordered
  ordered.x  <- model$info$ordered.x
  ordered.y  <- model$info$ordered.y
  is.probit  <- model$info$is.probit
  is.cexp    <- model$info$is.cexp
  is.nlin    <- model$info$is.nlin
  mc.reps    <- model$info$mc.reps
  etas       <- model$info$etas
  indsLvs    <- model$info$indsLvs
  is.ordered <- length(ordered) > 0

  model$status$finished <- TRUE
  model$status$iterations.0_9 <- model$status$iterations.0_9 + 1L

  if (model$status$iterations.0_9 >= model$status$max.iter.0_9) {
    warning("Maximum number of 0 through 9 iterations reached!")
    model$status$convergence <- FALSE
    model$status$finished    <- TRUE
    return(model)
  }

  if (is.ordered && is.probit && !is.nlin) {
    for (ord in ordered) {
      tau <- getThresholdsFromQuantiles(X = model$data, variable = ord)
      model$params$values <- c(model$params$values, tau)
    }

    model$params$se <- rep(NA_real_, length(model$params$values))

  } else if (length(ordered) && is.cexp) {
    params.old <- model$params$values.old
    params.new <- model$params$values
    X          <- as.data.frame(model$data)

    if (!is.null(params.old)) {
      paths.new <- params.new[!grepl("\\||=~|~~|~1", names(params.new))]
      paths.old <- params.old[!grepl("\\||=~|~~|~1", names(params.old))]

      eps <- mean(abs(paths.new - paths.old))

      converged <- eps <= model$status$tolerance
      failed    <- is.na(eps)

      if (failed) {
        warning("NAs in estimated model parameters!")
        model$status$convergence <- FALSE
        model$status$finished    <- TRUE

        return(model)

      } else if (converged) {
        model$status$convergence <- TRUE
        model$status$finished    <- TRUE

      } else {
        model$status$convergence <- FALSE
        model$status$finished    <- FALSE
      }

    } else {
      model$status$finished <- FALSE
    }

    parTable <- getParTableEstimates(model = model, rm.tmp = FALSE)
    parTable <- addColonPI_ParTable(parTable, model = model)
    parTable <- parTable[!(parTable$op == "=~" & grepl(":", parTable$lhs)), , drop = FALSE]

    sim.ov <- simulateDataParTable(
      parTable = parTable,
      N        = mc.reps,
      seed     = model$info$rng.seed
    )$ov

    # In theory we only need to update the thresholds for indicators of
    # endogenous variables. For now we just update everything, as it's
    # a cleaner way of getting the threshold parameters (though a bit slower)
    # It's also cleaner when bootstrapping the model
    for (ord.x in ordered.x) {
      rescaled <- rescaleOrderedVariableAnalytic(
        name = ord.x, data = X
      )

      model$params$values <- c(model$params$values, rescaled$thresholds)
      X[[ord.x]] <- rescaled$values
    }

    if (USE_NON_LINEAR_PROBIT_CORR_MAT) {
      for (eta in etas) {
        inds.eta <- indsLvs[[eta]]
        ord.eta  <- intersect(ordered, inds.eta)

        if (!length(ord.eta))
          next

        sim.cont <- sim.ov[, inds.eta, drop = FALSE]
        sim.ord  <- sim.ov[, inds.eta, drop = FALSE]
        colnames(sim.cont) <- paste0(".as_continous__", inds.eta)

        for (ord.y in ord.eta) {
          rescaled <- rescaleOrderedVariableMonteCarlo(
            name = ord.y, data = X, sim.ov = sim.ov
          )

          model$params$values <- c(model$params$values, rescaled$thresholds)
          X[[ord.y]] <- rescaled$values
          sim.ord[, ord.y] <- rescaled$sim.y.rescaled
        }

        sim.eta <- cbind(sim.cont, sim.ord)
        S.eta <- stats::cor(sim.eta)

        model$matrices$probit2cont[[eta]] <- S.eta
      }

    } else {
      for (ord.y in ordered.y) {
        rescaled <- rescaleOrderedVariableMonteCarlo(
          name = ord.y, data = X, sim.ov = sim.ov
        )
        model$params$values <- c(model$params$values, rescaled$thresholds)
        X[[ord.y]] <- rescaled$values
      }
    }

    model$params$se <- rep(NA_real_, length(model$params$values))

    model$params$values.old <- model$params$values
    model$data <- as.matrix(X)

    if (!is.probit) {
      # if we have a probit model, we will just get
      # the same polychoric matrix as before
      # So there is no point in updating the model...
      model$matrices$S <- getCorrMat(X, ordered = ordered, probit = is.probit)
    }
  }

  model
}
