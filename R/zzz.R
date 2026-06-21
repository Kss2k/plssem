PKG_INFO <- rlang::env(
  version = NULL,
  debug   = FALSE
)


UNICODE_MSG_STRINGS <- list(
  # SimDesign.RobbinsMonro()
  SimDesign.RobbinsMonro0 = "\rItertion: %i; Max change in E(\u03b8) = %.3f",
  SimDesign.RobbinsMonro1 = "\rIteration: %i; Max change in \u03b8 = %.3f",
  # limitQ()
  limitQ0 = "Reliability for %s is larger than 1! Q\u00B2(%s) = %.2f",
  limitQ1 = "Reliability for %s is close to or smaller than zero! Q\u00B2(%s) = %.2g"
)


ASCII_MSG_STRINGS <- list(
  # SimDesign.RobbinsMonro()
  SimDesign.RobbinsMonro0 = "\rItertion: %i; Max change in E(p) = %.3f",
  SimDesign.RobbinsMonro1 = "\rIteration: %i; Max change in p = %.3f",
  # limitQ()
  limitQ0 = "Reliability for %s is larger than 1! Q^2(%s) = %.2f",
  limitQ1 = "Reliability for %s is close to or smaller than zero! Q^2(%s) = %.2g"
)


MSG_STRINGS <- rlang::env(
  unicode = FALSE,
  strings = ASCII_MSG_STRINGS
)


getPackageVersion <- function(pkgname) {
  tryCatch({
    read.dcf(file = system.file("DESCRIPTION", package = pkgname),
             fields = "Version")
  }, error = function(e) {
    pls_msg_warn("Failed to get package version")
    "??" # replace this with a hard-coded value?
  })
}


.onLoad <- function(libname, pkgname) {
  PKG_INFO$version <- getPackageVersion(pkgname)

  utf8 <- tryCatch(isTRUE(l10n_info()[["UTF-8"]]), error = \(e) FALSE)
  MSG_STRINGS$unicode <- utf8

  if (utf8) MSG_STRINGS$strings <- UNICODE_MSG_STRINGS
  else      MSG_STRINGS$strings <- ASCII_MSG_STRINGS
}


.onAttach <- function(libname, pkgname) {
  version <- getPackageVersion(pkgname)
  message <- sprintf("This is %s (%s). Please report any bugs!", pkgname, version)

  packageStartupMessage(message)
}
