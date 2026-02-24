simulateDataParTable <- function(parTable, N = 1e5, seed = NULL) {
  if (!is.null(seed) && exists(".Random.seed")) .Random.seed.orig <- .Random.seed
  else                                          .Random.seed.orig <- NULL

  on.exit({
    if (!is.null(.Random.seed.orig)) .Random.seed <<- .Random.seed.orig
  })

  if (!is.null(seed))
    set.seed(seed)

  xis     <- getXis(parTable)
  etas    <- getSortedEtas(parTable)
  lvs     <- getLVs(parTable)
  indsLVs <- getIndsLVs(parTable, lVs = lvs)
  ovs     <- getOVs(parTable)

  Psi.x <- diag(1, length(xis))
  dimnames(Psi.x) <- list(xis, xis)

  for (i in seq_along(xis)) {
    xi.i <- xis[[i]]

    for (j in seq_len(i - 1)) {
      xi.j <- xis[[j]]

      p.ij <- parTable[
        (parTable$rhs == xi.i & parTable$op == "~~" & parTable$lhs == xi.j) |
        (parTable$rhs == xi.j & parTable$op == "~~" & parTable$lhs == xi.i),
        "est"
      ]

      Psi.x[i, j] <- Psi.x[j, i] <- p.ij
    }
  }

  Xi <- mvtnorm::rmvnorm(n = N, mean = rep(0, length(xis)), sigma = Psi.x)
  Xi <- as.data.frame(scale(Xi))
  colnames(Xi) <- xis

  undefIntTerms <- getIntTerms(parTable)
  elemsIntTerms <- stringr::str_split(undefIntTerms, pattern = ":")
  names(elemsIntTerms) <- undefIntTerms

  for (eta in etas) {

    for (intTerm in undefIntTerms) {
      elems <- elemsIntTerms[[intTerm]]

      if (all(elems %in% colnames(Xi))) {
        Xi[[intTerm]] <- multiplyIndicatorsCpp(Xi[elems])
        Xi[[intTerm]] <- Xi[[intTerm]] - mean(Xi[[intTerm]])

        undefIntTerms <- setdiff(undefIntTerms, intTerm)
      }
    }

    predRows <- parTable[parTable$lhs == eta & parTable$op == "~", , drop = FALSE]
    
    vals <- numeric(N)

    for (i in seq_len(NROW(predRows))) {
      row  <- predRows[i, ]
      beta <- row$est
      pred <- row$rhs

      vals <- vals + beta * Xi[[pred]]
    }

    resvar <- max(1 - var(vals), 0)
    vals <- vals + rnorm(N, mean = 0, sd = sqrt(resvar))

    Xi[[eta]] <- vals
  }

  Inds <- list()

  for (lv in lvs) {
    for (ind in indsLVs[[lv]]) {
      lambda <- parTable[parTable$lhs == lv &
                         parTable$op == "=~" &
                         parTable$rhs == ind, "est"]
      epsilon <- max(1 - lambda^2, 0)

      vals <- lambda * Xi[[lv]] + rnorm(N, mean = 0, sd = sqrt(epsilon))
      vals <- (vals - mean(vals)) / sd(vals)


      Inds[[ind]] <- vals
    }
  }

  Inds <- as.data.frame(Inds)
  All  <- cbind(Xi, Inds)
  Lv   <- All[lvs]
  Ov   <- All[ovs]

  list(
    all = All,
    ov  = Ov,
    lv  = Lv
  )
}
