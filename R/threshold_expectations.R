plsMapOrderedToExpectations <- function(z) {
  x <- as.integer(as.ordered(z))

  freq <- table(x)
  pct  <- freq[-length(freq)] / sum(freq)
  k    <- length(freq)
  tau  <- c(-Inf, qnorm(pct), Inf)

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
