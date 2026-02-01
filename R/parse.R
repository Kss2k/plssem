TEMP_OV_PREFIX <- ".TEMP_OV__"


parseModelArguments <- function(syntax, data) {
  parTable <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  data     <- as.data.frame(data)

  # Check for observed (structural) variables
  structovs <- getStructOVs(parTable)
  ovs       <- getOVs(parTable)

  structovs <- structovs[!grepl("\\(|\\)", structovs)]
  ovs       <- ovs[!grepl("\\(|\\)", ovs)]

  missing <- setdiff(ovs, colnames(data))
  stopif(length(missing), "Missing observed variables in data:\n  ",
         paste(missing, collapse = ", "))

  for (ov in structovs) {
    tmp.ov <- paste0(TEMP_OV_PREFIX, ov)

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

  list(
    syntax       = syntax,
    data         = data,
    parTable.pls = parTable.pls,
    parTable.all = parTable,
    cluster      = cluster,
    lme4.syntax  = lme4.syntax
  )
}


getClusterFromMultilevelStrings <- function(strings) {
  split <- stringr::str_remove_all(strings, pattern = "\\(|\\)") |>
    stringr::str_split_fixed(pattern = "\\|", n = 2L)

  unique(unlist(stringr::str_split(split[,2L], pattern = "\\+")))
}
