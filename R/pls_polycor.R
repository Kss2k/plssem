rawcor <- function(x, y) {
  if (!is.numeric(x)) x <- as.numeric(x)
  if (!is.numeric(y)) y <- as.numeric(y)
  cor(x, y)
}


plsPolychor <- function(x, y,
                        control = list(),
                        maxrho =.999,
                        start = rawcor(x, y)) {
  tab <- table(x, y)
  zerorows <- apply(tab, 1, \(x) all(x == 0))
  zerocols <- apply(tab, 2, \(x) all(x == 0))

  zr <- sum(zerorows)
  zc <- sum(zerocols)

  pls_warnif(zr > 0, paste0(
    zr, " row", suffix <- if(zr == 1) "" else "s",
    " with zero marginal", suffix," removed"
  ))

  pls_warnif(zc > 0, paste0(
    zc, " column", suffix <- if(zc == 1) "" else "s",
    " with zero marginal", suffix, " removed"
  ))

  tab <- tab[!zerorows, ,drop=FALSE]  
  tab <- tab[, !zerocols, drop=FALSE] 

  r <- nrow(tab)
  c <- ncol(tab)

  if (r < 2) {
    pls_msg_warn("the table has fewer than 2 rows")
    return(NA)
  }

  if (c < 2) {
    pls_msg_warn("the table has fewer than 2 columns")
    return(NA)
  }

  n <- sum(tab)
  rc <- qnorm(cumsum(rowSums(tab))/n)[-r]
  cc <- qnorm(cumsum(colSums(tab))/n)[-c]

  row.cuts <- c(-Inf, rc, Inf)
  col.cuts <- c(-Inf, cc, Inf)

  ridx <- matrix(seq_len(r), nrow = r, ncol = c, byrow = FALSE)
  cidx <- matrix(seq_len(c), nrow = r, ncol = c, byrow = TRUE)
  keep <- which(tab > 0)

  t <- tab[keep]
  i <- ridx[keep]
  j <- cidx[keep]

  cache_rho <- NA_real_
  cache_P <- NULL
  cache_G <- NULL

  compute <- function(rho) {
    if (!is.na(cache_rho) && identical(rho, cache_rho)) {
      return(list(P = cache_P, G = cache_G))
    }

    lower_x <- row.cuts[i]
    lower_y <- col.cuts[j]
    upper_x <- row.cuts[i + 1L]
    upper_y <- col.cuts[j + 1L]

    P <- pbinorm(
      lower_x = lower_x, lower_y = lower_y,
      upper_x = upper_x, upper_y = upper_y,
      rho = rho
    )

    G <- (
      dbinorm(upper_x, upper_y, rho = rho, force_zero = TRUE) -
      dbinorm(lower_x, upper_y, rho = rho, force_zero = TRUE) -
      dbinorm(upper_x, lower_y, rho = rho, force_zero = TRUE) +
      dbinorm(lower_x, lower_y, rho = rho, force_zero = TRUE)
    )

    cache_rho <<- rho
    cache_P <<- P
    cache_G <<- G

    list(P = P, G = G)
  }

  objective <- function(rho) {
    vals <- compute(rho)
    -sum(t * log(vals$P))
  }

  gradient <- function(rho) {
    vals <- compute(rho)
    -sum(t * vals$G / vals$P)
  }

  opt <- suppressWarnings(nlminb(
    objective = objective,
    gradient = gradient,
    start = start,
    lower = -abs(maxrho), upper = abs(maxrho)
  ))

  opt$par
}


binBvn <- function(rho, rc, cc) {  
  row.cuts <- c(-Inf, rc, Inf)
  col.cuts <- c(-Inf, cc, Inf)

  r <- length(row.cuts) - 1
  c <- length(col.cuts) - 1

  idx <- expand.grid(seq_len(r), seq_len(c))
  i <- idx[[1L]]
  j <- idx[[2L]]

  p <- pbinorm(
    lower_x = row.cuts[i], lower_y = col.cuts[j],
    upper_x = row.cuts[i+1], upper_y = col.cuts[j+1],
    rho = rho
  )

  matrix(p, nrow = r, ncol = c)
}


gradBinvBvn <- function(rho, rc, cc) {
  # gradient with respect to rho
  row.cuts <- c(-Inf, rc, Inf)
  col.cuts <- c(-Inf, cc, Inf)

  r <- length(row.cuts) - 1
  c <- length(col.cuts) - 1

  idx <- expand.grid(seq_len(r), seq_len(c))
  i <- idx[[1L]]
  j <- idx[[2L]]

  lower_x <- row.cuts[i]
  lower_y <- col.cuts[j]
  upper_x <- row.cuts[i+1]
  upper_y <- col.cuts[j+1]

  g <- (
    dbinorm(upper_x, upper_y, rho = rho, force_zero = TRUE) -
    dbinorm(lower_x, upper_y, rho = rho, force_zero = TRUE) -
    dbinorm(upper_x, lower_y, rho = rho, force_zero = TRUE) +
    dbinorm(lower_x, lower_y, rho = rho, force_zero = TRUE)
  )

  matrix(g, nrow = r, ncol = c)
}


pbinorm <- function(upper_x = NULL, upper_y = NULL, rho = 0.0,
                    lower_x = -Inf, lower_y = -Inf, check = FALSE) {

  n <- length(upper_x)
  stopifnot(length(upper_y) == n)
  if (n > 1L) {
    if (length(rho) == 1L) {
      rho <- rep(rho, n)
    }
    if (length(lower_x) == 1L) {
      lower_x <- rep(lower_x, n)
    }
    if (length(lower_y) == 1L) {
      lower_y <- rep(lower_y, n)
    }
  }

  upper_only <- all(lower_x == -Inf & lower_y == -Inf)
  if (upper_only) {
    # pbivnorm::pbivnorm does not handle +/-Inf well...
    upper_x[upper_x == +Inf] <- exp(10)
    upper_y[upper_y == +Inf] <- exp(10)
    upper_x[upper_x == -Inf] <- -exp(10)
    upper_y[upper_y == -Inf] <- -exp(10)
    res <- pbivnorm::pbivnorm(upper_x, upper_y, rho = rho)
  } else {
    # pbivnorm::pbivnorm does not handle +/-Inf well...
    upper_x[upper_x == +Inf] <- exp(10) # better pnorm?
    upper_y[upper_y == +Inf] <- exp(10)
    lower_x[lower_x == -Inf] <- -exp(10)
    lower_y[lower_y == -Inf] <- -exp(10)
    res <- pbivnorm::pbivnorm(upper_x, upper_y, rho = rho) -
      pbivnorm::pbivnorm(lower_x, upper_y, rho = rho) -
      pbivnorm::pbivnorm(upper_x, lower_y, rho = rho) +
      pbivnorm::pbivnorm(lower_x, lower_y, rho = rho)
  }

  res
}


dbinorm <- function(u, v, rho, force_zero = FALSE) {
  # dirty hack to handle extreme large values for rho
  # note that u, v, and rho are vectorized!
  rho_limit <- 0.9999
  abs_rho <- abs(rho)
  idx <- which(abs_rho > rho_limit)
  if (length(idx) > 0L) {
    rho[idx] <- sign(rho[idx]) * rho_limit
  }

  r <- 1 - rho * rho
  out <- 1 / (2 * pi * sqrt(r)) *
       exp(-0.5 * (u * u - 2 * rho * u * v + v * v) / r)

  # if abs(u) or abs(v) are very large (say, >10), set result equal
  # to exactly zero
  idx <- which(abs(u) > 10 | abs(v) > 10)
  if (length(idx) > 0L && force_zero) {
    out[idx] <- 0
  }

  out
}
