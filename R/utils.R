getNonZeroElems <- function(x) {
  as.vector(x[!is.na(x) & x != 0])
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
  cond1 <- grepl(":", parTable$rhs)
  cond2 <- !grepl("\\(|\\)\\|", parTable$rhs)
  unique(parTable[cond1 & cond2, "rhs"])
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


getCovOnlyVars <- function(parTable) {
  covRows <- parTable[parTable$op == "~~", , drop = FALSE]
  if (!NROW(covRows)) return(NULL)

  covVars <- unique(c(covRows$lhs, covRows$rhs))
  inds    <- getIndicators(parTable)

  setdiff(covVars, inds)
}


getStructOVs <- function(parTable) {
  struct <- getStructVars(parTable)
  cov    <- getCovOnlyVars(parTable)

  intersect(union(struct, cov), getOVs(parTable))
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


getHigherOrderLVs <- function(parTable) {
  lVs                  <- getLVs(parTable)
  isHigherOrder        <- logical(length(lVs))
  names(isHigherOrder) <- lVs

  for (lV in lVs) {
    inds <- parTable[parTable$lhs == lV & parTable$op %in% c("<~", "=~"), "rhs"] |>
      stringr::str_split(pattern = ":") |> unlist()

    if (any(inds %in% lVs)) isHigherOrder[[lV]] <- TRUE
  }

  if (!any(isHigherOrder)) NULL else lVs[isHigherOrder]
}


getParTableFromParNames <- function(parnames) {
  OP <- "~~|=~|~1|~|<~"
  op <- stringr::str_extract(parnames, pattern = OP)
  lr <- stringr::str_split_fixed(parnames, pattern = OP, n = 2)

  lhs <- lr[, 1]
  rhs <- lr[, 2]
  op[is.na(op)] <- "~"

  list(lhs = lhs, op = op, rhs = rhs)
}


namedListUnion <- function(x, y) {
  if      (is.null(y)) return(x)
  else if (is.null(x)) return(y)

  new <- setdiff(names(y), names(x))
  x[new] <- y[new]
  x
}
