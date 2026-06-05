printf <- function(...) {
  cat(sprintf(...))
  utils::flush.console()
}


messagef <- function(...) {
  message(sprintf(...), appendLF = FALSE)
  utils::flush.console()
}


quickdf <- function(l) {
  class(l) <- "data.frame"
  attr(l, "row.names") <- .set_row_names(length(l[[1]]))
  l
}


tryCatchNA <- function(expr) {
  tryCatch(expr, error = \(e) NA_real_)
}


lapplyNamed <- function(...) {
  lapply(...) |> stats::setNames(names(...))
}


uniqueComplete <- function(x) {
  unique(x[stats::complete.cases(x)])
}


namedListUnion <- function(x, y) {
  out <- x
  out[names(y)] <- y
  out
}


emptyNamedList <- function(nm) {
  stats::setNames(vector("list", length(nm)), nm = nm)
}


expandVcov <- function(vcov, labels) {
  if (!NROW(vcov) || !NCOL(vcov)) {
    return(matrix(
      0, nrow = length(labels), ncol = length(labels),
      dimnames = list(labels, labels)
    ))
  }

  labels.vcov <- colnames(vcov)
  labels.vv <- intersect(labels, labels.vcov)
  labels.zz <- setdiff(labels, labels.vcov)

  m <- length(labels.vv)
  k <- length(labels.zz)

  Vvv <- vcov[labels.vv, labels.vv]
  Vzz <- matrix(0, nrow = k, ncol = k, dimnames = list(labels.zz, labels.zz))
  Vvz <- matrix(0, nrow = m, ncol = k, dimnames = list(labels.vv, labels.zz))

  V <- rbind(cbind(Vvv, Vvz), cbind(t(Vvz), Vzz))
  V[labels, labels] # sort
}


rename <- function(.X, ...) {
  newNames <- list(...)
  oldNames <- names(newNames)

  for (old in oldNames) {
    names(.X)[names(.X) == old] <- newNames[[old]]
  }

  .X
}
