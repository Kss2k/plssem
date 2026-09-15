getNonZeroElems <- function(x) {
  as.vector(x[!is.na(x) & x != 0])
}


getOlsPathCoefs <- function(y, X, C) {
  # y: dependent factor
  # X: independent factors
  # C: correlation matrix
  # `solve(A, b)` avoids forming the explicit inverse.
  solve(C[X, X, drop = FALSE], C[X, y])
}


weightsProdInds <- function(wx, wy) {
  combos <- as.data.frame(expand.grid(wx, wy))
  colnames(combos) <- c("wx", "wy")
  w <- apply(combos, MARGIN = 1, FUN = function(row)
        row[[1]] * row[[2]])

  if (!is.null(names(wx)) && !is.null(names(wy))) {
    comboNames <- as.data.frame(expand.grid(names(wx), names(wy)))
    colnames(comboNames) <- c("wx", "wy")
    names(w) <- apply(comboNames, MARGIN = 1, FUN = function(row)
                      paste0(row[[1]], ":", row[[2]]))
  }

  w
}


diagPartitioned <- function(X, Y) {
  out <- rbind(cbind(X, matrix(0, nrow = nrow(X), ncol = ncol(Y))),
               cbind(matrix(0, nrow = nrow(Y), ncol = ncol(X)), Y))
  colnames(out) <- c(colnames(X), colnames(Y))
  rownames(out) <- c(rownames(X), rownames(Y))
  out
}


diag2 <- function(X) {
  if (NROW(X) <= 1L)
    return(X[1, 1, drop=FALSE])

  n <- NROW(X)
  Y <- matrix(0, n, n, dimnames = dimnames(X))
  Y[1L + 0L:(n - 1L) * (n + 1L)] <- X[1L + 0L:(n - 1L) * (n + 1L)]
  Y
}


tr <- function(x) {
  sum(diag(x))
}


isPositiveDefinite <- function(X, tol = 1e-8) {
  # A Cholesky factorisation succeeds exactly for positive-definite matrices and
  # is far cheaper than a full eigendecomposition. `pivot = TRUE` reports the
  # rank instead of erroring, and the pivoted diagonal gives the same tolerance
  # test as the smallest eigenvalue.
  if (anyNA(X)) return(FALSE)
  # `chol(pivot = TRUE)` warns rather than errors on a non-PD matrix, so no
  # condition handler is needed -- and handlers are expensive relative to a
  # factorisation this small.
  R <- suppressWarnings(chol.default(X, pivot = TRUE))
  attr(R, "rank") == NCOL(X) && min(diag(R)^2) > tol
}


gradSymMat <- function(gX) {
  # gX is the gradient of a matrix X, assuming the matrix X is non-symmetric
  # here we return the gradient of X assuming X is symmetric
  2 * gX - diag(diag(gX))
}
