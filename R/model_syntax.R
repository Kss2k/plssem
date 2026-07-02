plsParseModelSyntax <- function(syntax) {
  parTable <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  parTable$start <- parseStartModifiers(parTable$mod)
  parTable
}


parseModelArguments <- function(parTable,
                                data,
                                ordered = NULL,
                                probit = NULL,
                                mcpls = FALSE,
                                mc.fast.lmer = NULL,
                                consistent = TRUE,
                                is.lower.order = FALSE,
                                strict = TRUE) {
  # make sure we're working with a data.frame
  data <- asDataFrame(data)

  # Check for interation terms
  checkLhsIntTerms(parTable) # make sure they are all independent

  intTermNames <- getIntTerms(parTable)
  intTermElems <- stats::setNames(
    stringr::str_split(intTermNames, pattern = ":"), nm = intTermNames
  )

  is.nlin <- length(intTermElems) > 0L

  # Check for MIMIC blocks
  # If we have constructs which have a mix of "<~" and "=~" we redefine
  # lv=~ov -> ov~lv. We add a suffix such that we can clean it up later
  mode.c <- intersect(
    getFormativeLVs(parTable),
    getReflectiveLVs(parTable)
  )

  for (lv in mode.c) {
    # what needs to be redefined?
    idxs <- which(parTable$lhs == lv & parTable$op == "=~")

    for (idx in idxs) {
      ov <- parTable[idx, "rhs"]
      tmp <- paste0(ov, TEMP_MIMIC_SUFFIX)

      if (ov %in% colnames(data))
        data[[tmp]] <- data[[ov]]

      if (ov %in% ordered)
        ordered <- c(ordered, tmp)

      # redefine measurement model
      parTable[idx, "lhs"] <- tmp
      parTable[idx, "op"]  <- "~"
      parTable[idx, "rhs"] <- lv
    }
  }

  # Check for dupliacted indicators
  isind <- parTable$op %in% MOPS
  isstr <- parTable$op == "~"

  indicators <- parTable[isind, "rhs"]

  # A variable is structural if it is connected by a regression path (`~`) or
  # explicitly declared as a path-less structural node (e.g. `x ~ 1`). A `~~`
  # no longer promotes a variable to a structural variable; it only frees a
  # (co)variance among existing nodes. This means covariances may be specified
  # within the measurement model without renaming the indicators involved.
  declvars   <- getDeclaredStructVars(parTable)
  structvars <- unique(c(parTable[isstr, "rhs"], declvars))

  # A duplicated indicator can occur under two circumstances:
  #   1. It's an indicator which is part of two constructs
  #   2. It's an indicator which is used both as a structural variable
  #      and as an indicator. This is quite likely in higher order models.
  #      In these cases we rename the indicators, to avoid
  #      name clashes.
  dupMsr <- unique(indicators[duplicated(indicators)])     # Case 1
  dupStr <- unique(indicators[indicators %in% structvars]) # Case 2
  dupAll <- unique(c(dupMsr, dupStr))

  for (ind in dupAll) {

    cond <- isind & parTable$rhs == ind
    k    <- sum(cond)

    newnames <- paste0(
      parTable[cond, "rhs"], TEMP_IND_SUFFIX, seq_len(k)
    )
    parTable[cond, "rhs"] <- newnames

    if (!is.null(ordered) && ind %in% ordered)
      ordered <- union(ordered, newnames)

    for (newname in newnames)
      data[[newname]] <- data[[ind]]

  }

  # Check for observed (structural) variables
  structovs <- getStructOVs(parTable)
  ovs       <- getOVs(parTable)

  # Remove any (x + z + ... + y | cluster1 + cluster2 + ... + cluster3) expressions
  structovs <- structovs[!grepl("\\(|\\)", structovs)]
  ovs       <- ovs[!grepl("\\(|\\)", ovs)]

  vars       <- intersect(ovs, colnames(data))
  data       <- checkAndFixDTypesPLS_Data(data, check = vars)
  is.ordered <- vapply(data[vars], FUN.VALUE = logical(1L), FUN = is.ordered)
  ordered    <- intersect(union(ordered, vars[is.ordered]), vars)

  for (ord in ordered)
    data[[ord]] <- as.integer(as.ordered(data[[ord]]))

  missing <- setdiff(ovs, colnames(data))
  pls_stopif(length(missing), paste0("Missing observed variables in data:\n  ",
             paste(missing, collapse = ", ")))

  for (ov in structovs) {
    tmp.ov <- paste0(TEMP_OV_PREFIX, ov)

    if (length(ordered))
      ordered[ordered==ov] <- tmp.ov

    data[[tmp.ov]] <- data[[ov]]

    # Replace (potentially) existing measurement expressions
    # parTable[parTable$op %in% MOPS &
    #          parTable$rhs == ov, "rhs"] <- tmp.ov

    # Add measurment equation
    parTable <- rbind(
      parTable,
      data.frame(lhs = ov, op = "<~", rhs = tmp.ov, mod = "", start = NA)
    )
  }

  # `ADDITIONAL_STRUCT_VAR_OP` (`~1`) rows only *declare* structural nodes; they
  # are not estimable equations. The structural side (GLS) reads them straight
  # off the raw parTable, and observed declarations have, by now, been turned
  # into measurement equations above, so drop them before the PLS measurement
  # model.
  parTable <- parTable[parTable$op != ADDITIONAL_STRUCT_VAR_OP, , drop = FALSE]

  # Recompile syntax
  syntax <- parTableToSyntax(parTable)

  # Check for multilevel/mixed effects
  isMultilevel <- grepl("\\(.*\\|.*\\)", parTable$rhs) & parTable$op %in% c("~", "~~")
  if (any(isMultilevel)) {
    multilevelEtas <- unique(parTable[isMultilevel & parTable$op == "~", "lhs"])
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

  COALLESCE <- \(x, y) if (is.null(x)) isTRUE(y) else isTRUE(x)

  has.ord            <- length(ordered) > 0L
  is.mlm             <- length(lme4.syntax) > 0L
  use.mcpls.default  <- (has.ord && (is.nlin || is.lower.order)) || is.mlm
  is.mcpls           <- COALLESCE(mcpls, use.mcpls.default)
  is.probit          <- COALLESCE(probit, has.ord && !is.mcpls)
  mc.fast.lmer       <- COALLESCE(mc.fast.lmer, is.mcpls)

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
    is.mcpls     = is.mcpls,
    is.mlm       = is.mlm,
    mc.fast.lmer = mc.fast.lmer,
    consistent   = consistent && (!is.mcpls || is.mlm) # Don't use consistency correction
                                                       # for single-level MC-PLS models
  )
}


getClusterFromMultilevelStrings <- function(strings) {
  split <- stringr::str_remove_all(strings, pattern = "\\(|\\)") |>
    stringr::str_split_fixed(pattern = "\\|", n = 2L)

  unique(unlist(stringr::str_split(split[,2L], pattern = "\\+")))
}


parseStartModifiers <- function(mod) {
  pattern <- "^start\\((.*)\\)$"

  idx   <- which(grepl(pattern, mod))
  start <- rep(NA_real_, length(mod))

  vals <- stringr::str_extract(mod[idx], pattern = pattern, group = 1L)
  vals <- suppressWarnings(as.numeric(vals))

  if (anyNA(vals)) {
    bad <- mod[idx][is.na(vals)]

    pls_msg_warn(
      "Could not parse start values for some `start()` modifiers!",
      "Failed to parse:", paste0(bad, collapse = ", ")
    )
  }

  start[idx] <- vals
  start
}
