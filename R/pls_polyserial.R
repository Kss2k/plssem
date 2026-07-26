plsPolyserial <- function(x, y,
                          control = list(),
                          maxrho =.999,
                          start = rawcor(x, y)) {
  if (!is.integer(y))
    y <- as.integer(y)
  if (!is.numeric(x))
    x <- as.numeric(x)

  freq <- fastIntTab(y) # y is ordinal

  pls_stopif(length(freq) <= 1,
    "Ordinal variable must have at least two (observed) categories!"
  )
  stopifnot(length(x) == length(y))
  
  n <- length(y)
  thr <- c(-Inf, qnorm(cumsum(freq)/n)[-length(freq)], Inf)
  tau0 <- thr[y]
  tau1 <- thr[y+1]
  finite.lower <- is.finite(tau0)
  finite.upper <- is.finite(tau1)

  cache.rho  <- NA_real_ # for now
  z          <- (x - mean(x)) / sd(x)
  logy       <- NULL
  lower      <- NULL
  upper      <- NULL
  plower     <- NULL
  pupper     <- NULL
  prob       <- NULL
  zrlower    <- NULL
  zrupper    <- NULL
  sigma      <- NULL

  updateCache <- function(rho) {
    if (!is.na(cache.rho) && identical(rho, cache.rho))
      return(invisible(NULL))
  
    cache.rho <<- rho
    ey    <- z * rho
    var   <- 1 - rho * rho
    sigma <<- sqrt(var)
    lower <<- tau0 - ey
    upper <<- tau1 - ey

    plower <<- pnorm(lower, sd = sigma)
    pupper <<- pnorm(upper, sd = sigma)
    prob   <<- pupper - plower

    logy <<- log(prob)

    rho.var <- rho / var
    zrlower <<- z - rho.var * lower
    zrupper <<- z - rho.var * upper

    invisible(NULL)
  }

  plsPolyserialObjective <- function(rho) {
    updateCache(rho)
    -sum(logy)
  }

  plsPolyserialGradient <- function(rho) {
    updateCache(rho)

    lowerTerm <- dnorm(lower, 0, sigma) * zrlower
    upperTerm <- dnorm(upper, 0, sigma) * zrupper

    lowerTerm[!finite.lower] <- 0
    upperTerm[!finite.upper] <- 0

    -sum((lowerTerm - upperTerm) / prob)
  }
  
  # try 1
  optim <- suppressWarnings(nlminb(
    objective = plsPolyserialObjective,
    gradient = plsPolyserialGradient,
    start = start, control = control,
    lower = -abs(maxrho), upper = abs(maxrho)
  ))

  # try 2
  if (optim$convergence != 0L) {
    # try again, with different starting value
    optim <- suppressWarnings(nlminb(
      objective = plsPolyserialObjective,
      gradient = plsPolyserialGradient,
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
