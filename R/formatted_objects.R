plssemMatrix <- function(mat, symmetric = isSymmetric(mat), is.public = FALSE) {
  if (is.null(mat)) return(mat)

  if (is.public) {
    rn <- rownames(mat)
    cn <- colnames(mat)

    if (!is.null(rn)) {
      isTempRn <- startsWith(rn, prefix = TEMP_OV_PREFIX)
      rnClean <- stringr::str_remove_all(rn, pattern = TEMP_OV_PREFIX)
      isDupTempRn <- rnClean[isTempRn] %in% rnClean[!isTempRn]

      keepRn <- rep(TRUE, length(rn))
      keepRn[isTempRn][isDupTempRn] <- FALSE
      mat <- mat[keepRn, , drop = FALSE]
    }

    if (!is.null(cn)) {
      isTempCn <- startsWith(cn, prefix = TEMP_OV_PREFIX)
      cnClean <- stringr::str_remove_all(cn, pattern = TEMP_OV_PREFIX)
      isDupTempCn <- cnClean[isTempCn] %in% cnClean[!isTempCn]

      keepCn <- rep(TRUE, length(cn))
      keepCn[isTempCn][isDupTempCn] <- FALSE
      mat <- mat[ , keepCn, drop = FALSE]
  }

  class(mat) <- unique(c("PlsSemMatrix", class(mat)))

  if (symmetric)
    class(mat) <- unique(c("PlsSemSymmetricMatrix", class(mat)))

  mat
}


#' @export
print.PlsSemSymmetricMatrix <- function(x, digits = 3, sep = " ", ...) {
  if (nrow(x) == ncol(x)) x[upper.tri(x)] <- NA

  print.PlsSemMatrix(x, digits = digits, sep = sep, ...)
}


#' @export
print.PlsSemMatrix <- function(x, digits = 3, shift = 0L, ...) {
  y <- matrix(formatNumeric(x, digits = digits), nrow = nrow(x),
              ncol = ncol(x), dimnames = dimnames(x))

  y <- unclass(y)
  # Remove NAs
  y[is.na(x)] <- ""

  if (!is.null(colnames(x))) {
    colnames(y) <- abbreviate(colnames(x), minlength = digits + 3L) }
  if (shift > 0L) {
    empty.string <- rep(strrep(x = " ", times = shift), times = nrow(x))

    if (!is.null(rownames(x))) rownames(y) <- paste0(empty.string, rownames(x))
    else                       rownames(y) <- empty.string
  }

  print(y, ..., quote = FALSE, right = TRUE)

  invisible(x)
}


plssemParTable <- function(parTable, is.public = FALSE) {
  if (is.null(parTable)) return(parTable)

  rownames(parTable) <- NULL
  class(parTable) <- unique(c("PlsSemParTable", class(parTable)))
  parTable
}


plssemVector <- function(vec, is.public = FALSE) {
  if (is.null(vec)) return(vec)

  if (is.public) {
    nm <- names(vec)
    isTemp <- startsWith(nm, prefix = TEMP_OV_PREFIX)

    if (any(isTemp)) {
      clean <- stringr::str_remove_all(nm, pattern = TEMP_OV_PREFIX)

      isDupTemp <- clean[isTemp] %in% clean[!isTemp]

      keep <- rep(TRUE, length(nm))
      keep[isTemp][isDupTemp] <- FALSE

      vec <- stats::setNames(vec[keep], nm = clean[keep])
    }
  }

  class(vec) <- unique(c("PlsSemVector", class(vec)))

  vec
}


#' @export
print.PlsSemVector <- function(x, digits = 3L, ...) {
  y <- formatNumeric(x, digits = digits)
  print(y, quote = FALSE)
}


#' @export
print.PlsSemParTable <- function(x, nd = 3L, ...) {
  row.names <- rownames(x)
  y <- lapply(x, \(x) if (is.numeric(x)) round(x, nd) else x) |>
    as.data.frame()

  rownames(y) <- row.names

  print(y, ...)

  invisible(x)
}


formatNumeric <- function(x, digits = 3, scientific = FALSE,
                          justify = "right", width = NULL) {
  digits_fmt <- if (is.finite(digits)) max(0L, as.integer(digits)) else 3L
  digits_fmt_fmt <- max(1L, digits_fmt)
  if (is.numeric(x)) {
    x_round <- round(x, digits_fmt)
    format(x_round, nsmall = digits_fmt, digits = digits_fmt_fmt,
           trim = FALSE, justify = justify, scientific = scientific,
           width = width)
  } else {
    format(x, trim = FALSE, justify = justify, scientific = scientific,
           width = width)
  }
}

