TEMP_OV_PREFIX <- ".TEMP_OV__"
TEMP_IND_SUFFIX <- "__TEMP_IND_"

TEMP_OV_PREFIX_PATTERN <- paste0("^", TEMP_OV_PREFIX)
TEMP_IND_SUFFIX_PATTERN <- paste0(TEMP_IND_SUFFIX, "([0-9]+)$")

parseModelArguments <- function(parTable,
                                data,
                                ordered = NULL,
                                probit = NULL,
                                mcpls = FALSE,
                                mc.fast.lmer = NULL,
                                consistent = TRUE,
                                is.lower.order = FALSE) {
  data <- as.data.frame(data)

  # Check for interation terms
  intTermNames <- getIntTerms(parTable)
  intTermElems <- stringr::str_split(intTermNames, pattern = ":")
  names(intTermElems)  <- intTermNames
  is.nlin <- length(intTermElems) > 0L

  # Check for dupliacted indicators
  isind <- parTable$op %in% MOPS
  isstr <- parTable$op == "~"
  iscov <- parTable$op == "~~"

  indicators <- parTable[isind, "rhs"]
  covvars    <- unique(c(parTable[iscov, "lhs"], parTable[iscov, "rhs"]))
  structvars <- unique(c(parTable[isstr, "rhs"], covvars))

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
  stopif(length(missing), "Missing observed variables in data:\n  ",
         paste(missing, collapse = ", "))

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
      data.frame(lhs = ov, op = "<~", rhs = tmp.ov, mod = "")
    )
  }

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

  has.ord <- length(ordered) > 0L
  is.mcpls      <- if (is.null(mcpls)) has.ord && (is.nlin || is.lower.order) else mcpls
  is.probit     <- if (is.null(probit)) has.ord && !is.mcpls else probit
  is.mc.fast.lmer <- if (is.null(mc.fast.lmer)) is.mcpls else isTRUE(mc.fast.lmer)
  is.mlm    <- length(lme4.syntax) > 0L

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
    is.mcpls        = is.mcpls,
    is.mc.fast.lmer = is.mc.fast.lmer,
    is.mlm          = is.mlm,
    consistent   = consistent && (!is.mcpls || is.mlm) # Don't use consistency correction
                                                       # for single-level MC-PLS models
  )
}


getClusterFromMultilevelStrings <- function(strings) {
  split <- stringr::str_remove_all(strings, pattern = "\\(|\\)") |>
    stringr::str_split_fixed(pattern = "\\|", n = 2L)

  unique(unlist(stringr::str_split(split[,2L], pattern = "\\+")))
}
