plsPolychor <- function(x, y,
                        control = list(),
                        maxrho =.999,
                        start = rawcor(x, y)) {
  freq <- fastIntTab(x, y)
  zerorows <- rowSums(freq) == 0
  zerocols <- colSums(freq) == 0

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

  freq <- freq[!zerorows, ,drop=FALSE]  
  freq <- freq[, !zerocols, drop=FALSE] 

  r <- nrow(freq)
  c <- ncol(freq)

  if (r < 2) {
    pls_msg_warn("the cross table has fewer than 2 rows")
    return(NA)
  }

  if (c < 2) {
    pls_msg_warn("the cross table has fewer than 2 columns")
    return(NA)
  }

  n <- sum(freq)
  rc <- qnorm(cumsum(rowSums(freq))/n)[-r]
  cc <- qnorm(cumsum(colSums(freq))/n)[-c]
  kx <- length(rc)
  ky <- length(cc)

  # We can ignore computing a corner probability if none of the four
  # adjacent cells is nonzero.
  nzero <- freq > 0 & !is.na(freq)
  nz00  <- rbind(cbind(nzero, FALSE), FALSE)   # top-left corner of each nonzero cell
  nz10  <- rbind(cbind(FALSE, nzero), FALSE)   # top-right corner of each nonzero cell
  nz01  <- rbind(FALSE, cbind(nzero, FALSE))   # bottom-left corner of each nonzero cell
  nz11  <- rbind(FALSE, cbind(FALSE, nzero))   # bottom-right corner of each nonzero cell
  keep  <- nz11 | nz10 | nz01 | nz00

  # Keep only interior corners between finite thresholds.
  innerCorners <- keep[seq_len(kx) + 1, seq_len(ky) + 1, drop = FALSE]
  outerCorners <- unname(rbind(
    FALSE, cbind(FALSE, innerCorners, FALSE), FALSE
  ))

  keep.inner.idx <- which(innerCorners)
  keep.outer.idx <- which(outerCorners)
  keep.freq.idx  <- which(nzero)

  pcorners <- unname(rbind(
    0,
    cbind(0, matrix(NA, nrow = kx, ncol = ky), pnorm(rc)),
    c(0, pnorm(cc), 1)
  ))

  gcorners <- unname(rbind(
    0, cbind(0, matrix(NA, nrow = kx, ncol = ky), 0), 0)
  )
    
  cache.rho  <- NA_real_ # for now
  P          <- NULL
  G          <- NULL
  t          <- freq[keep.freq.idx]
  upper.x    <- rep(rc, times = length(cc))[keep.inner.idx]
  upper.y    <- rep(cc, each = length(rc))[keep.inner.idx]
  nr         <- nrow(pcorners)

  freq.row <- ((keep.freq.idx - 1L) %% r) + 1L
  freq.col <- ((keep.freq.idx - 1L) %/% r) + 1L

  # Precompute corners for each nonzero frequency cell.
  idx11 <- freq.row + 1L + freq.col * nr
  idx10 <- freq.row + 1L + (freq.col - 1L) * nr
  idx01 <- freq.row + freq.col * nr
  idx00 <- freq.row + (freq.col - 1L) * nr
  
  updateCache <- function(rho) {
    if (!is.na(cache.rho) && identical(rho, cache.rho))
      return(list(P = P, G = G))
  
    cache.rho <<- rho
  
    # Get probabilities for corners
    pcorners[keep.outer.idx] <- pbivnorm::pbivnorm(
      x   = upper.x,
      y   = upper.y,
      rho = rho
    )

    P <<- pcorners[idx11] - pcorners[idx10] -
      pcorners[idx01] + pcorners[idx00]

    # Get densities for corners
    gcorners[keep.outer.idx] <- dbinorm(
      u   = upper.x,
      v   = upper.y,
      rho = rho,
      force.zero = TRUE # numerical stability
    )
    
    # Get (truncated) densites from corners (for gradient)
    G <<- gcorners[idx11] - gcorners[idx10] -
      gcorners[idx01] + gcorners[idx00]

    list(G = G, P = P)
  }

  plsPolycorObjective <- function(rho) {
    cache <- updateCache(rho = rho)
    -sum(t * log(cache$P), na.rm = TRUE)
  }

  plsPolycorGradient <- function(rho) {
    cache <- updateCache(rho = rho)
    -sum(t * cache$G / cache$P, na.rm = TRUE)
  }
  
  # try 1
  optim <- suppressWarnings(nlminb(
    objective = plsPolycorObjective,
    gradient = plsPolycorGradient,
    start = start, control = control,
    lower = -abs(maxrho), upper = abs(maxrho)
  ))

  # try 2
  if (optim$convergence != 0L) {
    # try again, with different starting value
    optim <- suppressWarnings(nlminb(
      objective = plsPolycorObjective,
      gradient = plsPolycorGradient,
      start = 0.0, control = control,
      lower = -abs(maxrho), upper = abs(maxrho)
    ))
  }

  # check convergence
  pls_warnif(optim$convergence != 0L,
    "estimation of polychoric correlation did not converge!"
  )

  optim$par
}


dbinorm <- function(u, v, rho, force.zero = FALSE, rho.lim = 0.9999) {
  # dirty hack to handle extreme large values for rho
  # note that u, v, and rho are vectorized!
  abs.rho <- abs(rho)
  idx <- which(abs.rho > rho.lim)
  if (length(idx) > 0L)
    rho[idx] <- sign(rho[idx]) * rho.lim

  r <- 1 - rho * rho
  out <- 1 / (2 * pi * sqrt(r)) *
       exp(-0.5 * (u * u - 2 * rho * u * v + v * v) / r)

  # if abs(u) or abs(v) are very large (say, >10), set result equal
  # to exactly zero
  idx <- which(abs(u) > 10 | abs(v) > 10)
  if (length(idx) > 0L && force.zero)
    out[idx] <- 0

  out
}


rawcor <- function(x, y) {
  if (!is.numeric(x)) x <- as.numeric(x)
  if (!is.numeric(y)) y <- as.numeric(y)
  cor(x, y)
}


fastIntTab <- function(x, y = NULL) {
  if (is.null(y)) {
    ok <- !is.na(x)
    x <- as.integer(x[ok])

    nr <- max(x)
    return(tabulate(x, nbins = nr))
  }

  ok <- !is.na(x) & !is.na(y)
  x <- as.integer(x[ok])
  y <- as.integer(y[ok])

  nr <- max(x)
  nc <- max(y)

  matrix(
    tabulate(x + (y - 1L) * nr, nbins = nr * nc),
    nrow = nr,
    ncol = nc
  )
}
