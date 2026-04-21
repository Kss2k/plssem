kNN_ImputeMissing <- function(data, k = 5, ordered = NULL) {
  X <- as.matrix(data)
  complete <- complete.cases(X)

  if (all(complete))
    return(as.data.frame(X))

  Y <- X[complete,, drop=FALSE]
  X.std <- apply(X, MARGIN = 2, FUN = standardizeAtomic)
  Y.std <- X.std[complete,, drop=FALSE]

  Z <- X[!complete,, drop=FALSE]
  Z.std <- X.std[!complete,, drop=FALSE]

  Zpattern <- !is.na(Z)
  patterns <- unique(Zpattern)

  if (length(ordered)) {
    ncat <- apply(X[,ordered,drop=FALSE], MARGIN = 2L, FUN = \(x) length(uniqueComplete(x)))
    nominal <- ordered[ncat == 2]
    ordered <- setdiff(ordered, nominal)
  } else {
    nominal <- NULL
  }

  if (length(ordered) && k %% 2 == 0) {
    warning2("k should be odd when imputing ordered variables!\n",
             "Using k=", k + 1, " instead of k=", k, "!")
    k <- k + 1
  }

  stopif(NROW(Y) < k, "kNN imputation expects that there exists at least k\n",
         "(k=", k, ") complete observations, ", NROW(Y), " were found!")


  cols <- colnames(X)

  for (i in seq_len(NROW(patterns))) {
    p.i    <- patterns[i, ]

    if (!any(p.i)) {
      warning2("Unable to impute values for rows with all missing values.")
      next
    }

    rows.z <- apply(Zpattern, MARGIN = 1L, FUN = \(row) all(row == p.i))

    Y.std.i <- Y.std[, p.i, drop = FALSE]
    Y.i     <- Y[, p.i, drop = FALSE]

    Z.std.i <- Z.std[rows.z, p.i, drop = FALSE]
    Z.i     <- Z[rows.z, p.i, drop = FALSE]

    neighbours <- FNN::knnx.index(data = Y.std.i, query = Z.std.i, k = k)

    for (variable in cols[!p.i]) {
      if      (variable %in% nominal) .f <- calcMode 
      else if (variable %in% ordered) .f <- median
      else                            .f <- mean

      y <- Y[,variable]
      z <- vapply(
        X = seq_len(NROW(Z.i)),
        FUN.VALUE = numeric(1L),
        FUN = \(j) .f(y[neighbours[j, ]])
      )

      Z[rows.z, variable] <- z
    }
  }

  X[!complete, ] <- Z
 
  as.data.frame(X)
}


calcMode <- function(x) {
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}
