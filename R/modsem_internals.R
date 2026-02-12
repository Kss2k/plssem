.internalModsemAPI <- function(nm) {
  # Function for calling internal modsem objects
  # using ::: directly leads to an R CMD check note

  eval(parse(text=paste0("modsem:::", nm)))
}


parTableToSyntax      <- .internalModsemAPI("parTableToSyntax")
getOVs                <- .internalModsemAPI("getOVs")
getLVs                <- .internalModsemAPI("getLVs")
getEtas               <- .internalModsemAPI("getEtas")
getStructOVs          <- .internalModsemAPI("getStructOVs")
getParTableLabels     <- .internalModsemAPI("getParTableLabels")
multiplyIndicatorsCpp <- .internalModsemAPI("multiplyIndicatorsCpp")
simulateDataParTable  <- .internalModsemAPI("simulateDataParTable")
createProdInds        <- .internalModsemAPI("createProdInds")
allignLhsRhs          <- .internalModsemAPI("allignLhsRhs")
getIndicators         <- .internalModsemAPI("getIndicators")
