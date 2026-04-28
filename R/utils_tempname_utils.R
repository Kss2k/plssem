removeTempOvPrefix <- function(x) {
  stringr::str_remove(x, pattern = TEMP_OV_PREFIX_PATTERN)
}


removeTempIndSuffix <- function(x) {
  stringr::str_remove(x, pattern = TEMP_IND_SUFFIX_PATTERN)
}


hasTempOvPrefix <- function(x) {
  startsWith(x, prefix = TEMP_OV_PREFIX)
}


hasTempIndSuffix <- function(x) {
  grepl(TEMP_IND_SUFFIX_PATTERN, x)
}


removeTempAffixes <- function(x) {
  x |> removeTempOvPrefix() |> removeTempIndSuffix()
}

