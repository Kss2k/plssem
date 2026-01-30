TEMP_OV_PREFIX <- ".TEMP_OV__"


parseModelArguments <- function(syntax, data) {
  parTable <- modsem::modsemify(syntax)
  data     <- as.data.frame(data)

  # Check for observed (structural) variables
  structovs <- getStructOVs(parTable)
  ovs       <- getOVs(parTable)

  missing <- setdiff(ovs, colnames(data))
  stopif(length(missing), "Missing observed variables in data:\n  ",
         paste(missing, collapse = ", "))

  if (!length(structovs))
    return(list(syntax = syntax, data = data, parTable = parTable))

  for (ov in structovs) {
    tmp.ov <- paste0(TEMP_OV_PREFIX, ov)

    data[[tmp.ov]] <- data[[ov]]
    parTable <- rbind(
      parTable,
      data.frame(lhs = ov, op = "=~", rhs = tmp.ov, mod = "1")
    )
  }

  syntax <- parTableToSyntax(parTable)

  list(syntax = syntax, data = data, parTable = parTable)
}
