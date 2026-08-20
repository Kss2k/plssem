plsMapOrderedToExpectations <- function(z) {
  x <- as.integer(as.ordered(z))

  freq <- table(x)
  pct  <- cumsum(freq) / sum(freq)
  tau  <- c(-Inf, stats::qnorm(pct[-length(pct)]), Inf)
  k    <- length(freq)

  y <- numeric(length(z))

  for (i in seq_len(k)) {
    a <- tau[i]
    b <- tau[i+1]
    
    pdf.a <- stats::dnorm(a)
    pdf.b <- stats::dnorm(b)
    cdf.a <- stats::pnorm(a)
    cdf.b <- stats::pnorm(b)
    
    y[x == i] <- (pdf.a - pdf.b) / (cdf.b - cdf.a)
  }

  y
}
