printf <- function(...) {
  cat(sprintf(...))
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

