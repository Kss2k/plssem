plsCor <- function(data, ordered = NULL) {
  
  vars <- colnames(data)
  ord  <- vars %in% ordered
  k    <- length(vars)

  S <- matrix(
    NaN, nrow = k, ncol = k,
    dimnames = list(vars, vars)
  )
  diag(S) <- 1

  for (i in seq_len(k)) {
    x <- data[,i, drop = TRUE]

    if (ord[[i]]) y <- as.integer(as.ordered(x))
    else          y <- as.numeric(y)

    data[,i] <- y
  }

  for (i in seq_len(k)) {
    ord.i <- ord[[i]]
    x <- data[,i, drop = TRUE]

    for (j in seq_len(i-1)) {
      ord.j <- ord[[j]]
      y <- data[,j, drop = TRUE]

      if (!ord.i && !ord.j) {
        rho <- cor(x, y)
      } else if (ord.i && ord.j) {
        rho <- tryCatch(plsPolychor(x, y), error  = \(e) NA)
      } else {
        pls_msg_stop("Not implemented yet!")
      }

      S[i, j] <- S[j, i] <- rho
    }
  }

  S
}
