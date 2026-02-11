TEMP_OV_PREFIX <- ".TEMP_OV__"


parseModelArguments <- function(syntax,
                                data,
                                pi.match = NULL,
                                pi.match.recycle = NULL,
                                ordered = NULL,
                                probit = NULL) {
  stopif(length(syntax) > 1L || !is.character(syntax),
         "`syntax` must be a string of length 1!")

  parTable <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  data     <- as.data.frame(data)

  intTermNames <- unique(parTable[grepl(":", parTable$rhs), "rhs"])
  intTermElems <- stringr::str_split(intTermNames, pattern = ":")
  names(intTermElems)  <- intTermNames
  is.nlin <- length(intTermElems) > 0L

  # Int Terms
  # Check for observed (structural) variables
  structovs <- getStructOVs(parTable)
  ovs       <- getOVs(parTable)

  # Remove any (x + z + ... + y | cluster1 + cluster2 + ... + cluster3) expressions
  structovs <- structovs[!grepl("\\(|\\)", structovs)]
  ovs       <- ovs[!grepl("\\(|\\)", ovs)]

  vars       <- intersect(ovs, colnames(data))
  is.ordered <- vapply(data[vars], FUN.VALUE = logical(1L), FUN = is.ordered)
  ordered    <- intersect(union(ordered, vars[is.ordered]), vars)

  for (ord in ordered)
    data[[ord]] <- as.integer(as.ordered(data[[ord]]))

  missing <- setdiff(ovs, colnames(data))
  stopif(length(missing), "Missing observed variables in data:\n  ",
         paste(missing, collapse = ", "))

  for (ov in structovs) {
    tmp.ov <- paste0(TEMP_OV_PREFIX, ov)

    if (length(ordered))
      ordered[ordered==ov] <- tmp.ov

    data[[tmp.ov]] <- data[[ov]]
    parTable <- rbind(
      parTable,
      data.frame(lhs = ov, op = "=~", rhs = tmp.ov, mod = "1")
    )
  }

  syntax <- parTableToSyntax(parTable)

  isMultilevel <- grepl("\\(.*\\|.*\\)", parTable$rhs) & parTable$op == "~"
  if (any(isMultilevel)) {
    multilevelEtas <- unique(parTable[isMultilevel, "lhs"])
    lme4.syntax <- character(0L)

    for (eta in multilevelEtas) {
      rhs <- parTable[parTable$op == "~" & parTable$lhs == eta, "rhs"]

      lme4.syntax <- c(
        lme4.syntax,
        sprintf("%s~%s", eta, paste0(rhs, collapse = "+"))
      )
    }

    cluster <- getClusterFromMultilevelStrings(parTable[isMultilevel, "rhs"])
    parTable.pls <- parTable[!isMultilevel, , drop = FALSE]

  } else {
    parTable.pls <- parTable
    cluster <- NULL
    lme4.syntax <- NULL
  }

  if (is.null(probit))
    probit <- length(ordered) > 0L

  is.probit <- probit && !is.nlin
  is.cexp   <- probit && is.nlin

  list(
    syntax       = syntax,
    data         = data,
    parTable.pls = parTable.pls,
    parTable.all = parTable,
    cluster      = cluster,
    lme4.syntax  = lme4.syntax,
    intTermElems = intTermElems,
    intTermNames = intTermNames,
    is.nlin      = is.nlin,
    ordered      = ordered,
    is.probit    = is.probit,
    is.cexp      = is.cexp
  )
}


getClusterFromMultilevelStrings <- function(strings) {
  split <- stringr::str_remove_all(strings, pattern = "\\(|\\)") |>
    stringr::str_split_fixed(pattern = "\\|", n = 2L)

  unique(unlist(stringr::str_split(split[,2L], pattern = "\\+")))
}
