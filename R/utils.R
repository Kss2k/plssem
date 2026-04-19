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
    # "TODO: Add multiple imputation"
    message("Removing missing data using list wise deletion...")
    data <- data[!missingCases, , drop = FALSE]
  }
 
  if (standardize) {
    data <- standardizeDataFrame(
      data    = data,
      cluster = cluster
    )
  }

  S <- getCorrMat(data[indicators], probit = is.probit, ordered = ordered)
  X <- as.matrix(data[indicators])

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


warning2 <- function(...) {
  warning(..., call. = FALSE)
}


stop2 <- function(...) {
  stop(..., call. = FALSE)
}


stopif <- function(cond, ...) {
  if (isTRUE(cond)) stop2(...)
}


warnif <- function(cond, ...) {
  if (isTRUE(cond)) warning2(...)
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


tetracor <- function(x, y) {
  # x is continous, y is ordinal
  X <- data.frame(x, y)
  lavaan::lavCor(X, ordered = "y")[1, 2]
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


getIntTerms <- function(parTable) {
  unique(parTable[grepl(":", parTable$rhs), "rhs"])
}


quickdf <- function(l) {
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  l
}


tryCatchNA <- function(expr) {
  tryCatch(expr, error = \(e) NA_real_)  
}


getReflectiveLVs <- function(parTable) {
  unique(parTable[parTable$op == "=~", "lhs"])
}


getFormativeLVs <- function(parTable) {
  unique(parTable[parTable$op == "<~", "lhs"])
}


getLVs <- function(parTable) {
  # reflective <- getReflectiveLVs(parTable)
  # formative  <- getFormativeLVs(parTable)
  # Don't get reflective and formative constructs seperately,
  # as we want to keep the sorting in the partable
  unique(parTable[parTable$op %in% c("<~", "=~"), "lhs"])
}


getIndicators <- function(parTable, observed = TRUE, op = c("=~", "<~")) {
  indicators <- unique(parTable[!grepl(":", parTable$rhs) &
                                parTable$op %in% op, "rhs"])

  if (observed) indicators <- indicators[!indicators %in% getLVs(parTable)]
  indicators
}


getReflectiveIndicators <- function(..., op = "=~") {
  getIndicators(..., op = op)
}


getOVs <- function(parTable) {
  lVs    <- getLVs(parTable)
  select <- parTable$op %in% c("=~", "~", "~~", "<~")
  vars   <- unique(c(parTable$lhs[select], parTable$rhs[select]))

  vars[!vars %in% lVs & !grepl(":", vars)]
}


getStructVars <- function(parTable) {
  struct <- parTable[parTable$op == "~", , drop = FALSE]
  unique(c(struct$lhs, struct$rhs))
}


getStructOVs <- function(parTable) {
  intersect(getStructVars(parTable), getOVs(parTable))
}


getEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  lVs <- getLVs(parTable)

  cond.lhs <- parTable$op == "~"
  cond.rhs <- parTable$op %in% c("=~", "<~") & parTable$rhs %in% lVs

  if (isLV) cond.lhs <- cond.lhs & parTable$lhs %in% lVs

  etas.lhs <- parTable[cond.lhs, "lhs"]
  etas.rhs <- parTable[cond.rhs, "rhs"]

  etas <- unique(c(etas.rhs, etas.lhs))
  stopif(checkAny && !length(etas), "No etas found")

  etas
}


getSortedEtas <- function(parTable, isLV = FALSE, checkAny = TRUE) {
  unsortedEtas <- getEtas(parTable, isLV = isLV, checkAny = checkAny)

  cond1 <- parTable$op == "~"
  cond2 <- parTable$op %in% c("=~", "<~") & parTable$rhs %in% unsortedEtas

  structExprs <- parTable[cond1, , drop = FALSE]
  measrExprs  <- parTable[cond2, , drop = FALSE]

  if (NROW(measrExprs)) {
    measr2struct <- measrExprs
    measr2struct$lhs <- measrExprs$rhs
    measr2struct$op  <- "~"
    measr2struct$rhs <- measrExprs$lhs

    structExprs <- rbind(structExprs, measr2struct)
  }

  sortedEtas  <- character(0L)

  while (length(sortedEtas) < length(unsortedEtas) && nrow(structExprs) > 0) {
    stopif(all(unique(structExprs$lhs) %in% structExprs$rhs), "Model is non-recursive")

    for (i in seq_len(nrow(structExprs))) {
      if ((eta <- structExprs[i, "lhs"]) %in% structExprs$rhs) next

      sortedEtas  <- c(eta, sortedEtas)
      structExprs <- structExprs[!grepl(eta, structExprs$lhs), ]
      break
    }
  }

  if (!all(sortedEtas %in% unsortedEtas) &&
      length(sortedEtas) != length(unsortedEtas)) {
      warning("unable to sort etas")
      return(unsortedEtas)
  }

  sortedEtas
}


getXis <- function(parTable, etas = NULL, isLV = TRUE, checkAny = TRUE) {
  if (is.null(etas)) etas <- getEtas(parTable, isLV = isLV)
  # add all LVs which are not etas
  xis <- unique(parTable[parTable$op %in% c("=~", "<~") & !parTable$lhs %in% etas, "lhs"])

  if (!isLV) { # add any other variabels found in structural expressions
    xis <- unique(c(xis, parTable[parTable$op == "~" &
                                  !parTable$rhs %in% etas, "rhs"]))
  }

  xis <- xis[!grepl(":", xis)] # remove interaction terms

  stopif(checkAny && !length(xis), "No xis found")
  xis
}


getIndsLVs <- function(parTable, lVs, isOV = FALSE, ovs = NULL) {
  if (!length(lVs)) return(NULL)

  measr <- parTable[parTable$op %in% c("=~", "<~") & parTable$lhs %in% lVs, ]
  stopif(!NROW(measr), "No measurement expressions found, for", lVs)

  if (isOV) .f <- \(lV) measr[measr$lhs == lV & measr$rhs %in% ovs, "rhs"]
  else      .f <- \(lV) measr[measr$lhs == lV, "rhs"]

  lapplyNamed(lVs, FUN = .f, names = lVs)
}


lapplyNamed <- function(X, FUN, ..., names = X) {
  structure(lapply(X, FUN, ...), names = names)
}


tr <- function(X) {
  sum(diag(X))
}


uniqueComplete <- function(x) {
  unique(x[!is.na(x)])
}
