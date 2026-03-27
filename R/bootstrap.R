bootstrap <- function(model,
                      zero.tol = 1e-10,
                      verbose  = model$info$verbose,
                      parallel = model$info$boot$parallel,
                      ncpus    = model$info$boot$ncpus,
                      R        = model$info$boot$R,
                      iseed    = model$info$boot$iseed) {

  # The logic here follows the lavaan package
  is.mc <- is.snow <- FALSE
  if (is.null(parallel)) parallel <- "no"
  parallel <- match.arg(parallel, c("no", "multicore", "snow"))

  data      <- model$data
  cluster   <- model$info$cluster
  is.probit <- model$info$is.probit
  ordered   <- model$info$ordered
  results   <- vector("list", R)
  verbose   <- verbose && !is.snow

  progress <- function(i) {
    available <- getOption("width")
    available <- if (is.null(available)) 80L else available
    len.out   <- max(40L, min(available - 20L, 60L))

    finished <- floor((i / R) * len.out)
    left     <- len.out - finished

    printf(
      paste0("\rBoostrap [%i] |", strrep("=", finished), strrep(" ", left), "|"), R
    )
  }

  .bootf <- function(i) {
    if (verbose) progress(i)

    sampleData       <- resample(data, cluster = cluster)
    model$matrices$S <- getCorrMat(sampleData, ordered = ordered, probit = is.probit)
    model$data       <- sampleData

    utils::capture.output(type = "message", { # capture real time output
      model <- suppressWarnings(estimatePLS(
        model = model,
        # args passed onto mcpls
        fixed.seed = TRUE,
        verbose    = FALSE
      ))
    })

    par <- model$params$values
    attr(par, "id") <- i

    par 
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
  # this is adapted from the Lavaan package
  # - iseed is used for both serial and parallel
  # - if iseed is not set, iseed is generated + .Random.seed created/updated
  #     -> tmp.seed <- NA
  # - if iseed is set: don't touch .Random.seed (if it exists)
  #     -> tmp.seed <- .Random.seed (if it exists)
  #     -> tmp.seed <- NULL (if it does not exist)

  if (is.null(iseed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      runif(1)
    }

    # identical(temp.seed, NA): Will not change .Random.seed in GlobalEnv
    temp.seed <- NA
    iseed <- runif(1, 0, 999999999)

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

      cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Ensure workers have the package namespace loaded so internal helper
      # functions referenced by `.bootf` resolve correctly.
      parallel::clusterEvalQ(cl, loadNamespace("plssem"))

      # No need for
      # if(RNGkind()[1L] == "L'Ecuyer-CMRG")
      # clusterSetRNGStream() always calls `RNGkind("L'Ecuyer-CMRG")`
      parallel::clusterSetRNGStream(cl, iseed = iseed)
      results <- parallel::parLapply(cl, seq_len(R), .bootf)

    }

  } else {
    results <- lapply(seq_len(R), .bootf)

  }

  if (verbose) {
    progress(R)
    cat("\n")
  }

  resultsMat <- do.call(rbind, results) 
  names(resultsMat) <- names(model$params$values)

  se <- apply(resultsMat, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)
  se[se <= zero.tol] <- NA_real_

  list(se = se, boot = resultsMat)
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

