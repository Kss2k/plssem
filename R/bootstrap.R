bootstrap <- function(model,
                      zero.tol = 1e-10,
                      verbose  = model@info$verbose,
                      parallel = model@info$boot$parallel,
                      ncpus    = model@info$boot$ncpus,
                      R        = model@info$boot$R,
                      iseed    = model@info$boot$iseed) {

  # The logic here follows the lavaan package
  is.mc <- is.snow <- FALSE
  if (is.null(parallel)) parallel <- "no"
  parallel <- match.arg(parallel, c("no", "multicore", "snow"))

  data      <- model@data
  cluster   <- model@info$cluster
  is.probit <- model@info$is.probit
  is.mcpls  <- model@info$is.mcpls
  ordered   <- model@info$ordered
  results   <- vector("list", R)
  verbose   <- verbose && !is.snow

  model.base <- model

  boot.optimize <- isTRUE(model@info$boot$optimize)
  mc.boot.control <- prepMCBootControl(
    boot.control  = model@info$boot$mc.boot.control,
    boot.optimize = boot.optimize,
    model         = model
  )

  progress <- function(i) {
    available <- getOption("width")
    available <- if (is.null(available)) 80L else available
    len.out   <- max(40L, min(available - 20L, 60L))

    finished <- floor((i / R) * len.out)
    left     <- len.out - finished

    printf(
      paste0("\rBootstrap [%i] |", strrep("=", finished), strrep(" ", left), "|"), R
    )
  }

  na.par <- stats::setNames(rep(NA_real_, length(model@params$values)),
                            names(model@params$values))

  P_START <- model@info$mc.args$p.start

  .bootf <- function(i) {
    if (verbose) progress(i)

    tryCatch({
      sampleData <- resample(data, cluster = cluster)
      sampleS    <- getCorrMat(sampleData, ordered = ordered, probit = is.probit)
      model.b    <- model.base

      model.b@data       <- sampleData
      model.b@matrices$S <- sampleS

      mc.args         <- model.b@info$mc.args
      boot.fixed.seed <- mc.boot.control$fixed.seed
      boot.polyak     <- mc.boot.control$polyak.juditsky
      boot.tol        <- mc.boot.control$tol
      boot.reps       <- mc.boot.control$mc.reps
      boot.p.start    <- if (mc.boot.control$reuse.p.start) P_START else NULL
      boot.verbose    <- mc.boot.control$verbose
      low.tol.penalty <- mc.boot.control$low.tol.penalty

      utils::capture.output(type = "message", { # capture real time output

        model.b <- suppressWarnings(estimatePLS(
          model = model.b,
          # args passed onto mcpls
          fixed.seed      = boot.fixed.seed,
          verbose         = boot.verbose,
          polyak.juditsky = boot.polyak,
          tol             = boot.tol,
          mc.reps         = boot.reps,
          p.start         = boot.p.start
        ))

      })

      par <- combinedModel(model.b)@params$values

      if (low.tol.penalty > 0 && is.mcpls) {
        k <- length(par)
        par <- par + stats::rnorm(k, mean = 0, sd = low.tol.penalty)
      }

      attr(par, "id") <- i
      par

    }, error = \(e) {
      warning("Bootstrap replicate ", i, " failed: ", conditionMessage(e))

      par <- na.par
      attr(par, "id") <- i
      par
    })
  }


  if (parallel != "no" && ncpus > 1L) {
    switch(
      parallel,
      multicore = { is.mc   <- .Platform$OS.type != "windows" },
      snow      = { is.snow <- TRUE                           },
                  { ncpus   <- 1L                             }
    )

    loadNamespace("parallel") # before recording seed!
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

  if (!(ncpus > 1L && (is.mc || is.snow))) { # Only for serial
    set.seed(iseed)
  }

  # this is adapted from the boot function in package boot
  if (ncpus > 1L && (is.mc || is.snow)) {

    if (is.mc) {
      RNGkind("L'Ecuyer-CMRG") # to allow for reproducible results
      set.seed(iseed)
      results <- parallel::mclapply(seq_len(R), .bootf, mc.cores = ncpus)

    } else if (is.snow) {
      cl <- tryCatch(
        parallel::makePSOCKcluster(rep("localhost", ncpus)),
        error = function(e) {
          warning2(
            "Failed to start PSOCK cluster; falling back to serial bootstrap.\n",
            "Message: ", conditionMessage(e)
          )
          NULL
        }
      )

      if (is.null(cl)) {
        results <- lapply(seq_len(R), .bootf)
      } else {
        on.exit(parallel::stopCluster(cl), add = TRUE)

      # Ensure workers have the package namespace loaded so internal helper
      # functions referenced by `.bootf` resolve correctly.
      # When running via devtools::test() the package is loaded from source
      # but not installed — workers must load from the same source tree.
      pkg_path <- if (requireNamespace("pkgload", quietly = TRUE)) {
        tryCatch(pkgload::package_file(), error = function(e) NULL)
      } else {
        NULL
      }
      if (!is.null(pkg_path)) {
        parallel::clusterCall(cl, function(p) pkgload::load_all(p, quiet = TRUE), pkg_path)
      } else {
        parallel::clusterEvalQ(cl, loadNamespace("plssem"))
      }

      # No need for
      # if(RNGkind()[1L] == "L'Ecuyer-CMRG")
      # clusterSetRNGStream() always calls `RNGkind("L'Ecuyer-CMRG")`
      parallel::clusterSetRNGStream(cl, iseed = iseed)
        results <- parallel::parLapply(cl, seq_len(R), .bootf)

      }

    }

  } else {
    results <- lapply(seq_len(R), .bootf)

  }

  if (verbose) {
    progress(R)
    cat("\n")
  }

  ids <- vapply(results, FUN.VALUE = integer(1L), FUN = \(x) attr(x, "id"))

  resultsMat <- do.call(rbind, results)
  rownames(resultsMat) <- ids
  colnames(resultsMat) <- names(combinedModel(model)@params$values)

  vcov <- stats::cov(resultsMat, use = "complete.obs")
  se <- sqrt(diag(vcov))
  se[se <= zero.tol] <- NA_real_

  list(se = se, boot = resultsMat, vcov = vcov)
}


resample <- function(X, n.out = NROW(X), cluster = NULL, replace = TRUE) {
  if (is.null(cluster)) {
    idx <- sample(NROW(X), size = n.out, replace = replace)
    return(X[idx, , drop = FALSE])
  }

  stopif(length(cluster) > 1L, "bootstrapping of multiple cluster variables is not implemented (yet)!")

  cluster.vals <- attr(X, "cluster")[, cluster, drop = TRUE]
  stopif(NROW(cluster.vals) != NROW(X), "Cluster must be of same length as data!")

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
