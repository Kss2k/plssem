#' @export
pls_predict <- function(object,
                        approach = c("earliest", "direct"),
                        newdata = NULL,
                        std.ord.exp = FALSE,
                        benchmark = "R2", # pearson and tetrachoric correlations
                        ...) {
  approach <- match.arg(tolower(approach), c("earliest", "direct"))
  benchmark <- match.arg(tolower(benchmark), c("r2"))

  W <- object$fit$fitWeights
  L <- object$fit$fitLambda

  ordered <- object$info$ordered
  outerX <- getOuterDataMatrices(object, std.ord.exp = std.ord.exp)
  X.cont <- outerX$X.cont
  X.ord  <- outerX$X.ord

  Y <- X.cont %*% W

  if (approach == "earliest") {
    parTable <- getParTableEstimates(object, rm.tmp = FALSE)
    xis      <- getXis(parTable)
    etas     <- getSortedEtas(parTable)

    undefIntTerms <- getIntTerms(parTable)
    elemsIntTerms <- stringr::str_split(undefIntTerms, pattern = ":")
    names(elemsIntTerms) <- undefIntTerms

    Y[,c(etas, undefIntTerms)] <- NA_real_

    for (eta in etas) {

      for (intTerm in undefIntTerms) {
        elems <- elemsIntTerms[[intTerm]]

        if (all(elems %in% colnames(Y))) {
          vals <- multiplyIndicatorsCpp(Y[,elems])
          Y[,intTerm] <- vals - mean(vals)

          undefIntTerms <- setdiff(undefIntTerms, intTerm)
        }
      }

      cond <- parTable$lhs == eta & parTable$op == "~"
      predRows <- parTable[cond, , drop = FALSE]

      vals <- numeric(NROW(Y))

      for (i in seq_len(NROW(predRows))) {
        row  <- predRows[i, ]
        beta <- row$est
        pred <- row$rhs

        vals <- vals + beta * Y[,pred]
      }

      Y[,eta] <- vals
    }
  }

  X.cont.pred <- Y %*% t(L)

  if (length(ordered)) {
    Tau <- outerX$Tau
    X.ord.pred  <- X.cont.pred

    for (ord in ordered) {
      breaks <- c(-Inf, sort(Tau[[ord]]), Inf)
      X.ord.pred[,ord] <- cut(X.ord.pred[,ord], breaks = breaks, labels = FALSE)
    }

  } else {
    X.ord.pred <- NULL

  }

  if (benchmark == "r2") {
    
    getr2 <- function(variable) {
      x <- X.cont.pred[,variable]

      if (variable %in% ordered) {
        y <- X.ord[,variable]
        r <- tetracor(x = x, y = y)
      } else {
        y <- X.cont[,variable]
        r <- cor(x = x, y = y)
      }

      r^2
    }

    benchmarked <- vapply(
      X = colnames(X.cont.pred),
      FUN.VALUE = numeric(1L),
      FUN = getr2
    )

  }

  out <- list(
    Y           = plssemMatrix(Y),
    X.cont      = plssemMatrix(X.cont),
    X.cont.pred = plssemMatrix(X.cont.pred),
    X.ord       = plssemMatrix(X.ord),      # Might be NULL
    X.ord.pred  = plssemMatrix(X.ord.pred), # Might be NULL
    benchmark   = plssemVector(benchmarked)
  )

  class(out) <- "PlsSemPredict"
  out
}


print.PlsSemPredict <- function(x, ...) {
  # Summarize prediction object
}

assignScoresOrdinalNormal <- function(x, std.ord.exp = FALSE, probs = NULL) {
  x.i <- reindex(x)

  # observed category proportions
  freq <- as.numeric(table(x))
  K    <- length(freq)
  n    <- sum(freq)

  stopif(K < 2, "Need at least 2 ordered categories for variable!")

  # cumulative probs for interior thresholds (K-1 of them)
  if (is.null(probs)) probs <- cumsum(freq)[-K] / n  # drop last; sums to 1
  C <- sort(probs[probs < 1])

  # thresholds on standard-normal scale
  tauf <- stats::qnorm(C)
  tau <- c(-Inf, tauf, Inf) # length K-1

  # truncated-normal means for each category interval (a_i, b_i]
  # mu_i = (phi(a_i) - phi(b_i)) / (Phi(b_i) - Phi(a_i))
  a <- tau[1:K]
  b <- tau[2:(K+1)]
  Phi_a <- stats::pnorm(a)
  Phi_b <- stats::pnorm(b)
  phi_a <- stats::dnorm(a)
  phi_b <- stats::dnorm(b)
  denom <- pmax(Phi_b - Phi_a, eps)
  mu    <- (phi_a - phi_b) / denom

  # map each observed category to its conditional mean
  x.out <- rep(NA_real_, length(x.i))
  for (i in seq_len(K)) x.out[x.i == i] <- mu[i]

  # optional standardization of the mapped scores (sample standardization)
  if (std.ord.exp)
    x.out <- standardizeAtomic(x.out)

  # labels for interior thresholds
  labels.t   <- paste0(name, "|t", seq_len(K - 1))
  thresholds <- stats::setNames(tauf, labels.t)

  attr(x.out, "thresholds") <- thresholds
  x.out
}


assignScoresOrdinalMonteCarlo <- function(x, y, y.i, std.ord.exp = FALSE,
                                          probs = NULL) {
  # x = observed categories, y = continous monte-carlo sampled values,
  # y.i = simulated categories

  x.i  <- reindex(x) # re-index
  y.i  <- reindex(y.i) # re-index
  freq <- as.numeric(table(y.i))
  K    <- length(freq)
  n    <- sum(freq)

  stopif(K < 2, "Need at least 2 ordered categories for variable!")

  # map each observed category to its conditional mean
  x.out <- rep(NA_real_, length(x.i))
  for (i in seq_len(K)) {
    x.out[x.i == i] <- mean(y[y.i == i])
  }

  # optional standardization of the mapped scores (sample standardization)
  if (std.ord.exp)
    x.out <- standardizeAtomic(x.out)

  # cumulative probs for interior thresholds (K-1 of them)
  if (is.null(probs)) probs <- cumsum(freq)[-K] / n  # drop last; sums to 1
  C <- sort(probs[probs < 1])

  # thresholds on standard-normal scale
  tauf <- collapse::fquantile(x, probs = C)
  tau <- c(-Inf, tauf, Inf) # length K-1

  # labels for interior thresholds
  labels.t   <- paste0("y|t", seq_len(K - 1))
  thresholds <- stats::setNames(tauf, labels.t)

  attr(x.out, "thresholds") <- thresholds

  x.out
}


getOuterDataMatrices <- function(model, newdata = NULL, std.ord.exp = FALSE) {
  ordered  <- model$info$ordered
  is.mcpls <- model$info$is.mcpls
  olddata  <- model$data
  ordered  <- intersect(colnames(olddata), ordered)
 
  if (!is.null(newdata)) {
    missing <- setdiff(colnames(olddata), colnames(newdata))
    stopif(length(missing), "Missing variables in `newdata`!\n",
           "Missing: ", paste0(missing, collapse = ", "))

    newdata <- as.matrix(as.data.frame(newdata)[colnames(olddata)])

  } else {
    newdata <- olddata

  }

  if (!length(ordered))
    return(list(X.cont = newdata, X.ord = NULL)) # no need to assign scores

  newdata.cont <- newdata
  newdata.ord  <- apply(newdata, MARGIN = 2, reindex)

  PROBS <- getPROBS(data = olddata, ordered = ordered)
  Tau   <- stats::setNames(vector("list", length(ordered)), nm = ordered)

  if (is.mcpls) {
    sim.ov.cont <- model$matrices$sim.ov.cont
    sim.ov.ord  <- model$matrices$sim.ov.ord

    for (ord in ordered) {
      y <- assignScoresOrdinalMonteCarlo(
        x = newdata.cont[,ord], y = sim.ov.cont[[ord]],
        y.i = sim.ov.ord[[ord]], std.ord.exp = std.ord.exp,
        probs = PROBS[[ord]]
      )

      newdata.cont[,ord] <- y
      Tau[[ord]] <- attr(y, "thresholds")
    }

  } else {

    for (ord in ordered) {
      y <- assignScoresOrdinalNormal(
        x = newdata.cont[,ord], std.ord.exp = std.ord.exp, probs = PROBS[[ord]]
      )

      newdata.cont[,ord] <- y
      Tau[[ord]] <- attr(y, "thresholds")
    }

  }

  list(
    X.cont = newdata.cont,
    X.ord  = newdata.ord,
    Tau    = Tau
  )
}


pls_construct_scores <- function(object, ...) {
  predicted <- pls_predict(object, ...)
  predicted$Y
}


reindex <- function(x) as.integer(as.ordered(x))
