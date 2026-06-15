TEMP_OV_PREFIX <- ".TEMP_OV__"
TEMP_IND_SUFFIX <- "__TEMP_IND_"
TEMP_MIMIC_SUFFIX <- "__TEMP_MIMIC"

TEMP_OV_PREFIX_PATTERN    <- TEMP_OV_PREFIX # plain for now
TEMP_IND_SUFFIX_PATTERN   <- paste0(TEMP_IND_SUFFIX, "([0-9]+)")
TEMP_MIMIC_SUFFIX_PATTERN <- TEMP_MIMIC_SUFFIX # plain for now


removeTempOvPrefix <- function(x) {
  stringr::str_remove_all(x, pattern = TEMP_OV_PREFIX_PATTERN)
}


removeTempIndSuffix <- function(x) {
  stringr::str_remove_all(x, pattern = TEMP_IND_SUFFIX_PATTERN)
}


removeTempMimicSuffix <- function(x) {
  stringr::str_remove_all(x, pattern = TEMP_MIMIC_SUFFIX)
}


hasTempOvPrefix <- function(x) {
  grepl(TEMP_OV_PREFIX_PATTERN, x)
}


hasTempIndSuffix <- function(x) {
  grepl(TEMP_IND_SUFFIX_PATTERN, x)
}


hasTempMimicSuffix <- function(x) {
  grepl(TEMP_MIMIC_SUFFIX_PATTERN, x)
}


removeTempAffixes <- function(x) {
  x |>
    removeTempOvPrefix() |>
    removeTempIndSuffix() |>
    removeTempMimicSuffix()
}


hasTempAffixes <- function(x) {
  hasTempOvPrefix(x) |
  hasTempIndSuffix(x) |
  hasTempMimicSuffix(x)
}
