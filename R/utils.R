getPLS_Data <- function(data,
                        indicators,
                        cluster = NULL,
                        consistent = TRUE,
                        standardize = TRUE,
                        ordered = NULL,
                        is.probit = NULL,
                        is.cexp = NULL) {
  vars <- c(indicators, cluster)
  missing <- !vars %in% colnames(data)

  if (any(missing))
    stop("Missing variables: ", paste0(vars[missing], collapse = ", "))

  data <- as.data.frame(data)[vars]

  missingCases <- !stats::complete.cases(data)
  if (any(missingCases)) {
    warning("Removing missing data list wise for factor scores.\n",
            "Removing missing data pair wise in covariance matrix.\n",
            "TODO: Add multiple imputation")
  }
 
  if (standardize) {
    data <- standardizeDataFrame(
      data    = data,
      cluster = cluster
    )
  }

  S <- getCorrMat(data[indicators], probit = is.probit, ordered = ordered)
  X <- as.matrix(data[indicators])

  if (is.cexp) for (ord in intersect(indicators, ordered)) {
    # Should serve as a better starting point than not correcting at all...
    X[,ord] <- rescaleOrderedVariableAnalytic(ord, data = data)$values
  }

  if (!is.null(cluster)) {
    if (!is.character(cluster))
      stop("`cluster` must be a character string, if lme4.syntax is provided!")

    attr(X, "cluster") <- data[, cluster, drop = FALSE]
  }

  list(X = X, S = S)
}


getNonZeroElems <- function(x) {
  as.vector(x[!is.na(x) & x != 0])
}


standardizeDataFrame <- function(data, cluster = NULL, subset = colnames(data)) {
  subset <- setdiff(subset, cluster)

  Z <- lapply(data[subset], FUN = standardizeAtomic)
  data[subset] <- Z

  data
}


standardizeAtomic <- function(x) {
  (x - mean(x)) / stats::sd(x)
}


getPathCoefs <- function(y, X, C) {
  # y: dependent factor 
  # X: independent factors
  # C: correlation matrix
  solve(C[X, X]) %*% C[X, y]
}



weightsProdInds <- function(wx, wy) {
  combos <- as.data.frame(expand.grid(wx, wy))
  colnames(combos) <- c("wx", "wy")
  w <- apply(combos, MARGIN = 1, FUN = function(row) 
        row[[1]] * row[[2]])

  if (!is.null(names(wx)) && !is.null(names(wy))) {
    comboNames <- as.data.frame(expand.grid(names(wx), names(wy)))
    colnames(comboNames) <- c("wx", "wy")
    names(w) <- apply(comboNames, MARGIN = 1, FUN = function(row) 
                      paste0(row[[1]], ":", row[[2]]))
  }

  w
}


diagPartitioned <- function(X, Y) {
  out <- rbind(cbind(X, matrix(0, nrow = nrow(X), ncol = ncol(Y))),
               cbind(matrix(0, nrow = nrow(Y), ncol = ncol(X)), Y))
  colnames(out) <- c(colnames(X), colnames(Y))
  rownames(out) <- c(rownames(X), rownames(Y)) 
  out
}


diag2 <- function(X) {
  Y <- diag(diag(X))
  dimnames(Y) <- dimnames(X)
  Y
}


printf <- function(...) {
  cat(sprintf(...))
  utils::flush.console()
}


cov2cor <- function(vcov) {
  if (is.null(vcov))
    return(NULL)

  sd <- sqrt(abs(diag(vcov))) # use `abs()`, in case some variances are negative

  if (length(sd) == 1L) D <- matrix(1 / sd, nrow = 1L, ncol = 1L)
  else                  D <- diag(1 / sd)

  structure(D %*% vcov %*% D, dimnames = dimnames(vcov))
}


stopif <- function(cond, ...) {
  if (isTRUE(cond)) stop(..., call. = FALSE)
}


warnif <- function(cond, ...) {
  if (isTRUE(cond)) warning(..., call. = FALSE)
}


getCorrMat <- function(data, probit = FALSE, ordered = NULL) {
  if (probit) getPolyCorr(data, ordered = ordered)
  else        getPearsonCorr(data)
}


getPearsonCorr <- function(data) {
  stats::cor(as.data.frame(data), use = "pairwise.complete.obs")
}


getPolyCorr <- function(data, ordered = NULL) {
  data <- as.data.frame(data)
  lavaan::lavCor(data, ordered = ordered)
}


getThresholdsFromQuantiles <- function(X, variable) {
  x   <- X[, variable]
  pct <- table(x) / length(x)
  tau <- stats::qnorm(cumsum(pct)[-length(pct)])
  lab <- paste0(variable, "|t", seq_along(tau))

  stats::setNames(tau, nm = lab)
}


getCorrMatsProbit2cont <- function(data, ordered, lvs, selectLambda) {
  probit2cont <- list()
  inds        <- rownames(selectLambda)

  if (length(ordered)) for (lv in lvs) {
    inds.lv <- inds[selectLambda[,lv]]
    ord.lv  <- intersect(inds.lv, ordered)

    if (!length(ord.lv)) # just skip
      next

    X.cont <- data[, inds.lv, drop = FALSE]
    X.ord  <- data[, inds.lv, drop = FALSE]

    colnames(X.cont) <- paste0(".as_continous__", inds.lv)
    X <- cbind(X.cont, X.ord)

    suppressWarnings({
      # Suppress this:
      #> Warning Message:
      #> lavaan->lav_samplestats_step2():  
      #> correlation between variables x1 and .as_continous__x1 
      #> is (nearly) 1.0 
      S.lv <- getCorrMat(X, probit = TRUE, ordered = ord.lv)
    })

    probit2cont[[lv]] <- S.lv
  }

  probit2cont 
}


formatNumeric <- function(x, digits = 3, scientific = FALSE,
                          justify = "right", width = NULL) {
  digits_fmt <- if (is.finite(digits)) max(0L, as.integer(digits)) else 3L
  digits_fmt_fmt <- max(1L, digits_fmt)
  if (is.numeric(x)) {
    x_round <- round(x, digits_fmt)
    format(x_round, nsmall = digits_fmt, digits = digits_fmt_fmt,
           trim = FALSE, justify = justify, scientific = scientific,
           width = width)
  } else {
    format(x, trim = FALSE, justify = justify, scientific = scientific,
           width = width)
  }
}


getOrderedResidualCorrection <- function(lvs, indsLvs, ordered, X) {
  res.ord  <- NULL
  res.cont <- NULL

  for (lv in lvs) {
    inds <- indsLvs[[lv]]
    inds.ov <- intersect(inds, ordered)

    syntax <- paste0(lv, "=~", paste0(inds.ov, collapse="+"))
    fit.c <- lavaan::cfa(syntax, X)
    fit.o <- lavaan::cfa(syntax, X, ordered = inds.ov)
  
    res.ord.lv  <- 1 - lavaan::lavInspect(fit.o, what = "r2")
    res.cont.lv <- 1 - lavaan::lavInspect(fit.c, what = "r2")

    new <- setdiff(names(res.ord.lv), names(res.ord))
    res.ord  <- c(res.ord, res.ord.lv[new])
    res.cont <- c(res.cont, res.cont.lv[new])
  }

  correction <- res.ord / res.cont
  correction[is.na(correction)] <- 1L

  list(
    res.ord = res.ord,
    res.cont = res.cont,
    correction = correction
  )
}


getIntTerms <- function(parTable) {
  unique(parTable[grepl(":", parTable$rhs), "rhs"])
}
