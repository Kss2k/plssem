plsCor <- function(data, ordered = NULL) {
  if (!is.data.frame(data))
    data <- data.frame(data)

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
    else          y <- as.numeric(x)

    data[,i] <- y
  }

  trycor <- function(.f, type) {
    function(...) {
      tryCatch(.f(...), error = function(e) {
        pls_msg_warn(
          "Estimation of", type, "correlation failed!",
          "Message:", conditionMessage(e)
        )
        NaN
      })
    }
  }

  hasNA <- any(!stats::complete.cases(data))

  pls_warnif(hasNA,
    "Found missing values in data!",
    "Using pairwise complete observations..."
  )

  polychor   <- trycor(plsPolychor, "polychoric")
  pearson    <- trycor(stats::cor, "Pearson")
  polyserial <- trycor(plsPolyserial, "polyserial")

  for (i in seq_len(k)) {
    ord.i <- ord[[i]]
    x <- data[,i, drop = TRUE]

    for (j in seq_len(i-1)) {
      ord.j <- ord[[j]]
      y <- data[,j, drop = TRUE]

      if (hasNA) {
        .x <- x
        ok <- !is.na(x) & !is.na(y)
        x  <- x[ok]
        y  <- y[ok]
      }

      if      (!ord.i && !ord.j) rho <- pearson(x, y)
      else if ( ord.i &&  ord.j) rho <- polychor(x, y)
      else if ( ord.i && !ord.j) rho <- polyserial(y, x) 
      else if (!ord.i &&  ord.j) rho <- polyserial(x, y) 
      else pls_msg_stop( # Should be impossible
        "Unexpected condition in `plsCor()`, this is likely a bug!"
      )

      if (hasNA) {
        x <- .x # restore in case we removed observations
      }

      S[i, j] <- S[j, i] <- rho
    }
  }

  S
}
