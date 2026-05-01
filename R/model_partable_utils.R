getIntTerms <- function(parTable) {
  cond1 <- grepl(":", parTable$rhs)
  cond2 <- !grepl("\\(|\\)\\|", parTable$rhs)
  unique(parTable[cond1 & cond2, "rhs"])
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
  indicators <- unique(parTable[parTable$op %in% op, "rhs"])

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
      structExprs <- structExprs[structExprs$lhs != eta, , drop = FALSE]
      break
    }
  }

  if (!all(unsortedEtas %in% sortedEtas) ||
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


getRandomEffectLabels <- function(parTable) {
  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  rlhs <- lhs[grepl("~", lhs) & op == "~~"]
  rrhs <- rhs[grepl("~", rhs) & op == "~~"]

  union(rlhs, rrhs)
}


getIndsLVs <- function(parTable, lVs, isOV = FALSE, ovs = NULL) {
  if (!length(lVs)) return(NULL)

  measr <- parTable[parTable$op %in% c("=~", "<~") & parTable$lhs %in% lVs, ]
  stopif(!NROW(measr), "No measurement expressions found, for", lVs)

  if (isOV) {
    if (is.null(ovs)) ovs <- getOVs(parTable)
    measr <- measr[measr$rhs %in% ovs, , drop = FALSE]
  }

  split <- split(measr$rhs, measr$lhs)
  split[vapply(split, FUN.VALUE = integer(1L), FUN = length) > 0L]
}


getHigherOrderLVs <- function(parTable) {
  lvs <- getLVs(parTable)
  cond <- parTable$op %in% c("=~", "<~") & parTable$rhs %in% lvs
  unique(parTable[cond, "lhs"])
}


getParTableFromParNames <- function(parnames) {
  split <- splitParameterNames(parnames)
  data.frame(lhs = split$lhs, op = split$op, rhs = split$rhs, stringsAsFactors = FALSE)
}
