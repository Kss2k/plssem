bootstrap <- function(model, R = 50L, zero.tol = 1e-10) {
  data    <- model$data 
  cluster <- model$info$cluster
  results <- vector("list", R)

  cli::cli_progress_bar(
    name  = "bootstrap",
    type  = "iterator",
    total = R
  )

  for (i in seq_len(R)) {
    cli::cli_progress_update()

    sampleData       <- resample(data, cluster = cluster)
    model$matrices$S <- stats::cov(as.data.frame(sampleData), use = "pairwise.complete.obs")
    model$data       <- sampleData

    model <- estimatePLS(model)
    results[[i]] <- model$params$values
  }

  cli::cli_progress_done()

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
  stopif(NROW(cluster.vals) != NROW(X), "Cluster must be of same lenght as data!")

  clusters <- unique(cluster.vals)
  G <- length(clusters)

  clusters.sample <- sample(clusters, size = G, replace = replace)

  .boot <- function(Z)
    do.call(rbind, lapply(clusters.sample, FUN = \(ci) Z[cluster.vals==ci, , drop=FALSE]))

  Y <- .boot(X)
  attr(Y, "cluster") <- .boot(attr(X, "cluster"))

  Y
}

