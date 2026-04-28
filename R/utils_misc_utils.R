printf <- function(...) {
  cat(sprintf(...))
  utils::flush.console()
}


warning2 <- function(...) {
  warning(..., call. = FALSE)
}


stop2 <- function(...) {
  stop(..., call. = FALSE)
}


stopif <- function(cond, ...) {
  if (isTRUE(cond)) stop2(...)
}


warnif <- function(cond, ...) {
  if (isTRUE(cond)) warning2(...)
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
  unique(x[complete.cases(x)])
}


namedListUnion <- function(x, y) {
  out <- x
  out[names(y)] <- y
  out
}

