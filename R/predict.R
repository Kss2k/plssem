#' Predict from a fitted PLS-SEM model
#'
#' @param object A fitted \code{plssem} model.
#' @param approach Prediction approach. If \code{approach = "earliest"} (default), then
#'   only indicators of exogenous benchmark.vars are used for prediction. If
#'   \code{approach = "direct"}, then all indicators are used.
#' @param benchmark.vars What predictions should be benchmarked? If \code{benchmark.vars = "endog"}
#'   (default), preidction benchmarks are applied to indicators of endogenous benchmark.vars.
#'   If \code{benchmark.vars = "exog"}, preidction benchmarks are applied to indicators
#'   of exogenous benchmark.vars. If \code{benchmark.vars = "all"}, preidction benchmarks are applied to
#'   all of the indicators in the model.
#' @param newdata Optional new data matrix/data frame.
#' @param std.ord.exp Logical; standardize ordinal expectation scores.
#' @param benchmark Benchmark type(s). Either length 1 (recycled) or one entry
#'   per indicator (optionally named). Supported: \code{"r2"}, \code{"rmse"},
#'   \code{"mae"}, \code{"q2_predict"}, \code{"acc"}, \code{"ord_mae"}.
#' @param ... Additional arguments passed to internal helpers.
#' @return A \code{PlsSemPredict} object with matrices and benchmark results.
#'
#' @export
pls_predict <- function(object,
                        approach = c("earliest", "direct"),
                        newdata = NULL,
                        std.ord.exp = FALSE,
                        benchmark = "R2",
                        benchmark.vars = c("endog", "exog", "all"),
                        ...) {
  # TODO:
  #  1. Allow the user to pass only indicators of exogenous variables, if
  #     approach='earliest'.
  #  2. Allow the user to pass ordinal variables with a subset of categories
  #  3. Use standardization parameters from the training data, don't just
  #     standardize the test data directly.

  approach <- match.arg(tolower(approach), c("earliest", "direct"))
  benchmark.vars <- match.arg(tolower(benchmark.vars), c("endog", "exog", "all"))

  W <- object$fit$fitWeights
  L <- object$fit$fitLambda

  ordered <- object$info$ordered
  outerX <- getOuterDataMatrices(object, newdata = newdata,
                                 std.ord.exp = std.ord.exp)
  X.cont <- outerX$X.cont
  X.ord  <- outerX$X.ord
  ordered <- intersect(ordered, colnames(X.cont))

  Y <- X.cont %*% W

  if (approach == "earliest") {
    parTable <- getParTableEstimates(object, rm.tmp = FALSE)
    xis      <- getXis(parTable)
    etas     <- getSortedEtas(parTable)

    undefIntTerms <- getIntTerms(parTable)
    elemsIntTerms <- stringr::str_split(undefIntTerms, pattern = ":")
    names(elemsIntTerms) <- undefIntTerms

    Y.sub <- as.data.frame(Y)[xis]

    for (eta in etas) {

      for (intTerm in undefIntTerms) {
        elems <- elemsIntTerms[[intTerm]]

        if (all(elems %in% colnames(Y.sub))) {
          vals <- multiplyIndicatorsCpp(Y.sub[elems])
          Y.sub[[intTerm]] <- vals - mean(vals)

          undefIntTerms <- setdiff(undefIntTerms, intTerm)
        }
      }

      cond <- parTable$lhs == eta & parTable$op == "~"
      predRows <- parTable[cond, , drop = FALSE]

      vals <- numeric(NROW(Y.sub))

      for (i in seq_len(NROW(predRows))) {
        row  <- predRows[i, ]
        beta <- row$est
        pred <- row$rhs

        vals <- vals + beta * Y.sub[[pred]]
      }

      Y.sub[[eta]] <- vals
    }

    Y <- as.matrix(Y.sub)
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

  Y           <- plssemMatrix(Y, is.public = TRUE)
  X.cont      <- plssemMatrix(X.cont, is.public = TRUE)
  X.cont.pred <- plssemMatrix(X.cont.pred, is.public = TRUE)
  X.ord       <- plssemMatrix(X.ord, is.public = TRUE)
  X.ord.pred  <- plssemMatrix(X.ord.pred, is.public = TRUE)
  ordered     <- unique(c(ordered, stringr::str_remove_all(ordered, TEMP_OV_PREFIX)))
  all.vars    <- colnames(X.cont.pred)
  benchmarked <- NULL

  inds.x <- stringr::str_remove_all(object$info$inds.x, TEMP_OV_PREFIX)
  inds.y <- stringr::str_remove_all(object$info$inds.y, TEMP_OV_PREFIX)

  pred.vars <- switch(benchmark.vars,
    all = all.vars,
    exog = intersect(all.vars, inds.x),
    endog = intersect(all.vars, inds.y),
    stop2("Unrecognize value for benchmark.vars argument: ", benchmark.vars)
  )

  if (!is.null(benchmark)) {
    allowedBenchmarks <- c("r2", "rmse", "mae", "acc", "ord_mae", "q2_predict")
    benchmark <- .normalizeBenchmarkType(benchmark)

    if (length(benchmark) == 1L) {
      benchmarkByVar <- stats::setNames(rep(benchmark, length(pred.vars)), pred.vars)

    } else if (is.null(names(benchmark)) || all(names(benchmark) == "")) {
      stopif(length(benchmark) != length(pred.vars),
             "`benchmark` must have length 1 or one entry per variable.")
      benchmarkByVar <- stats::setNames(benchmark, pred.vars)

    } else {
      missing <- setdiff(pred.vars, names(benchmark))
      extra <- setdiff(names(benchmark), pred.vars)

      stopif(length(missing),
             "Missing benchmark type(s) for variable(s): ",
             paste0(missing, collapse = ", "))
      warnif(length(extra),
             "Ignoring benchmark type(s) for unknown variable(s): ",
             paste0(extra, collapse = ", "))

      benchmarkByVar <- benchmark[pred.vars]
      benchmarkByVar <- stats::setNames(as.character(benchmarkByVar), pred.vars)
    }

    invalid <- setdiff(unique(benchmarkByVar), allowedBenchmarks)
    stopif(length(invalid),
           "Invalid `benchmark` type(s): ", paste0(invalid, collapse = ", "),
           "\nAllowed: ", paste0(allowedBenchmarks, collapse = ", "))

    trainOuterX <- getOuterDataMatrices(object, newdata = object$data,
                                        std.ord.exp = std.ord.exp)
    trainMean <- colMeans(
      plssemMatrix(trainOuterX$X.cont, is.public = TRUE), na.rm = TRUE
    )

    values <- mapply(
      FUN = .benchmark,
      variable = pred.vars,
      benchmark = unname(benchmarkByVar),
      X.cont = list(X.cont),
      X.cont.pred = list(X.cont.pred),
      X.ord = list(X.ord),
      X.ord.pred = list(X.ord.pred),
      ordered = list(ordered),
      trainMean = list(trainMean),
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )

    benchmarked <- data.frame(
      variable = pred.vars,
      metric = unname(benchmarkByVar),
      value = as.numeric(values),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    Y           = Y,
    X.cont      = X.cont,
    X.cont.pred = X.cont.pred,
    X.ord       = X.ord,      # Might be NULL
    X.ord.pred  = X.ord.pred, # Might be NULL
    benchmark   = benchmarked
  )

  class(out) <- "PlsSemPredict"
  out
}


#' Print a \code{PlsSemPredict} object
#'
#' @param x A \code{PlsSemPredict} object.
#' @param ... Additional arguments for compatibility with the generic.
#' @return The input object, invisibly.
#'
#' @export
print.PlsSemPredict <- function(x, ...) {
  fields <- names(x)
  hasOrd <- !is.null(x$X.ord) || !is.null(x$X.ord.pred)

  printf("PlsSemPredict object\n")
  printf("Available fields: %s\n\n", paste0("$", fields, collapse = ", "))

  printHead <- function(name, mat) {
    printf("%s [%i x %i] (head)\n", name, nrow(mat), ncol(mat))
    print(plssemMatrix(utils::head(mat)))
    cat("\n")
  }

  if (!is.null(x$Y))           printHead("Y", x$Y)
  if (!is.null(x$X.cont))      printHead("X.cont", x$X.cont)
  if (!is.null(x$X.cont.pred)) printHead("X.cont.pred", x$X.cont.pred)

  if (hasOrd) {
    if (!is.null(x$X.ord))      printHead("X.ord", x$X.ord)
    if (!is.null(x$X.ord.pred)) printHead("X.ord.pred", x$X.ord.pred)
  }

  if (!is.null(x$benchmark)) {
    if (is.data.frame(x$benchmark) &&
        all(c("variable", "metric", "value") %in% names(x$benchmark))) {
      bm <- x$benchmark
      metrics <- unique(bm$metric)

      printf("Benchmark summary\n")
      for (m in metrics) {
        vals <- bm$value[bm$metric == m]
        vals <- vals[is.finite(vals)]

        if (!length(vals)) {
          printf("%s: <empty>\n", m)
        } else {
          stats <- c(mean = mean(vals), median = stats::median(vals),
                     min = min(vals), max = max(vals))
          fmt <- formatNumeric(stats, digits = 3L)
          printf("%s: n=%i, mean=%s, median=%s, min=%s, max=%s\n",
                 m, length(vals), fmt["mean"], fmt["median"],
                 fmt["min"], fmt["max"])
        }

        bm.m <- bm[bm$metric == m, c("variable", "value"), drop = FALSE]
        bm.m$value <- formatNumeric(bm.m$value, digits = 3L)
        rownames(bm.m) <- NULL
        print(bm.m, row.names = FALSE)
      }

      cat("\n")

    } else {
      bench <- as.numeric(x$benchmark)
      bench <- bench[is.finite(bench)]

      printf("Benchmark")
      if (!length(bench)) {
        printf(": <empty>\n\n")
      } else {
        stats <- c(mean = mean(bench), median = stats::median(bench),
                   min = min(bench), max = max(bench))
        fmt <- formatNumeric(stats, digits = 3L)
        printf(": n=%i, mean=%s, median=%s, min=%s, max=%s\n\n",
               length(bench), fmt["mean"], fmt["median"],
               fmt["min"], fmt["max"])
      }
    }
  }

  invisible(x)
}


.benchmark <- function(variable, benchmark,
                       X.cont, X.cont.pred,
                       X.ord = NULL, X.ord.pred = NULL,
                       ordered = character(0), trainMean = NULL) {
  type <- benchmark
  xObs <- X.cont[,variable]
  xPred <- X.cont.pred[,variable]

  switch(
    type,
    r2 = .bm_r2(
      variable = variable,
      xObs     = xObs,
      xPred    = xPred,
      ordered  = ordered,
      yOrd     = if (is.null(X.ord)) NULL else X.ord[,variable]
    ),
    rmse = .bm_rmse(xObs = xObs, xPred = xPred),
    mae  = .bm_mae(xObs = xObs, xPred = xPred),
    q2_predict = .bm_q2_predict(
      variable = variable,
      xObs = xObs,
      xPred = xPred,
      trainMean = trainMean
    ),
    acc = {
      .bm_check_ord_only(type = type, variable = variable, ordered = ordered)
      .bm_check_ord_mats(type = type, X.ord = X.ord, X.ord.pred = X.ord.pred)
      .bm_acc(yObs = X.ord[,variable], yPred = X.ord.pred[,variable])
    },
    ord_mae = {
      .bm_check_ord_only(type = type, variable = variable, ordered = ordered)
      .bm_check_ord_mats(type = type, X.ord = X.ord, X.ord.pred = X.ord.pred)
      .bm_ord_mae(yObs = X.ord[,variable], yPred = X.ord.pred[,variable])
    },
    stop("Unhandled benchmark type: ", type, call. = FALSE)
  )
}

assignScoresOrdinalNormal <- function(x, std.ord.exp = FALSE, probs = NULL, eps = 1e-5) {
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
  labels.t   <- paste0("y|t", seq_len(K - 1))
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
  tauf <- collapse::fquantile(y, probs = C)
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
    nm.o <- colnames(olddata)
    nm.n <- colnames(newdata)

    is.tmp <- grepl(TEMP_OV_PREFIX, nm.o)

    if (any(is.tmp)) {
      tmpReplacements <- stats::setNames(
        nm.o[is.tmp], stringr::str_remove_all(nm.o[is.tmp], TEMP_OV_PREFIX)
      )

      keys <- intersect(nm.n, names(tmpReplacements))
      nm.n <- stats::setNames(nm.n, nm = nm.n)
      nm.n[keys] <- tmpReplacements[keys]

      colnames(newdata) <- nm.n
    }

    missing <- setdiff(colnames(olddata), colnames(newdata))
    stopif(length(missing), "Missing variables in `newdata`!\n",
           "Missing: ", paste0(missing, collapse = ", "))

    newdata.df <- as.data.frame(newdata)[colnames(olddata)]
    is.ord <- vapply(newdata.df, FUN.VALUE = logical(1L), FUN = is.ordered)
    newdata.df[is.ord] <- lapply(newdata.df[is.ord], reindex)

    if (model$info$standardized)
      newdata <- Rfast::standardise(as.matrix(newdata.df))
    else
      newdata <- as.matrix(newdata.df)

  } else {
    newdata <- olddata

  }

  if (!length(ordered))
    return(list(X.cont = newdata, X.ord = NULL)) # no need to assign scores

  # Check that we have the same number of categories in newdata and olddata
  for (ord in ordered) {
    ncatOld <- length(uniqueComplete(olddata[,ord]))
    ncatNew <- length(uniqueComplete(newdata[,ord]))
    
    stopif(ncatNew < ncatOld,
           "There are fewer categories for ", ord, " in the test data (", 
           ncatNew, "),\nthan in the training data (", ncatOld, ")!")

    stopif(ncatNew > ncatOld,
           "There are more categories for ", ord, " in the test data (", 
           ncatNew, "),\nthan in the training data (", ncatOld, ")!")
  }

  newdata.cont <- newdata
  newdata.ord  <- apply(newdata, MARGIN = 2, reindex)
  dimnames(newdata.ord) <- dimnames(newdata)

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


#' Construct latent variable scores
#'
#' Convenience wrapper around [pls_predict()] returning only the predicted
#' latent scores matrix.
#'
#' @param object A fitted \code{plssem} model.
#' @param ... Passed to [pls_predict()].
#' @return A \code{PlsSemMatrix} of predicted latent scores.
#'
#' @export
pls_construct_scores <- function(object, ...) {
  predicted <- pls_predict(object, ...)
  predicted$Y
}


.normalizeBenchmarkType <- function(x) {
  x <- tolower(x)
  x[x %in% c("ordmae", "ord-mae", "ordinalmae", "ordinal_mae")] <- "ord_mae"
  x[x %in% c("q2", "q2predict", "q2-predict", "q2_predict")] <- "q2_predict"
  x
}


.bm_check_ord_only <- function(type, variable, ordered) {
  stopif(!(variable %in% ordered),
         "Benchmark `", type, "` is only available for ordered variables: ",
         variable)
}


.bm_check_ord_mats <- function(type, X.ord, X.ord.pred) {
  stopif(is.null(X.ord) || is.null(X.ord.pred),
         "Ordinal matrices are NULL; cannot compute `", type, "`.")
}


.bm_r2 <- function(variable, xObs, xPred, ordered, yOrd = NULL) {
  if (variable %in% ordered) {
    stopif(is.null(yOrd), "Missing observed ordinal values for: ", variable)
    r <- tryCatchNA(tetracor(x = xPred, y = yOrd))
  } else {
    r <- tryCatchNA(stats::cor(x = xPred, y = xObs, use = "complete.obs"))
  }
  r^2
}


.bm_rmse <- function(xObs, xPred) {
  sqrt(mean((xPred - xObs)^2, na.rm = TRUE))
}


.bm_mae <- function(xObs, xPred) {
  mean(abs(xPred - xObs), na.rm = TRUE)
}


.bm_q2_predict <- function(variable, xObs, xPred, trainMean) {
  stopif(is.null(trainMean), "Missing training mean; cannot compute Q2_predict.")
  sse <- sum((xPred - xObs)^2, na.rm = TRUE)
  sst <- sum((xObs - trainMean[variable])^2, na.rm = TRUE)
  if (!is.finite(sst) || sst <= 0) return(NA_real_)
  1 - sse / sst
}


.bm_acc <- function(yObs, yPred) {
  mean(yPred == yObs, na.rm = TRUE)
}


.bm_ord_mae <- function(yObs, yPred) {
  mean(abs(yPred - yObs), na.rm = TRUE)
}


reindex <- function(x) as.integer(as.ordered(x))
