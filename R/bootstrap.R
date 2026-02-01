bootstrap <- function(model, R = 50L, zero.tol = 1e-10) {
  data      <- model$data
  cluster   <- model$info$cluster
  is.probit <- model$info$is.probit
  ordered   <- model$info$ordered
  results   <- vector("list", R)

  cli::cli_progress_bar(
    name  = sprintf("bootstrap[%i]", R),
    type  = "iterator",
    total = R
  )

  for (i in seq_len(R)) {
    cli::cli_progress_update()

    sampleData       <- resample(data, cluster = cluster)
    model$matrices$S <- getCorrMat(sampleData, ordered = ordered, probit = is.probit)
    model$data       <- sampleData

    stats::capture.output(type = "message", { # capture real time output
      model <- suppressWarnings(estimatePLS(model))
    })

    results[[i]] <- model$params$values
  }

  cli::cli_progress_done()

  suppressWarnings({ # TODO: Fix mismatching thresholds in bootstrapping
    resultsMat <- do.call(rbind, results) 
  })
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

