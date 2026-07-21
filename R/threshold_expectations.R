plsMapOrderedToExpectations <- function(z) {
  x <- as.integer(as.ordered(z))

  freq <- table(x)
  pct  <- cumsum(freq) / sum(freq)
  tau  <- c(-Inf, qnorm(pct[-length(pct)]), Inf)
  k    <- length(freq)

  y <- numeric(length(z))

  for (i in seq_len(k)) {
    a <- tau[i]
    b <- tau[i+1]
    
    pdf.a <- dnorm(a)
    pdf.b <- dnorm(b)
    cdf.a <- pnorm(a)
    cdf.b <- pnorm(b)
    
    y[x == i] <- (pdf.a - pdf.b) / (cdf.b - cdf.a)
  }

  y
}
