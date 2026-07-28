plsPolyserial <- function(x, y,
                          control = list(),
                          maxrho =.999,
                          start = NULL) {
  if (!is.integer(y))
    y <- as.integer(as.ordered(y))
  if (!is.numeric(x))
    x <- as.numeric(x)

  freq <- fastIntTab(y) # y is ordinal

  pls_stopif(length(freq) <= 1,
    "Ordinal variable must have at least two (observed) categories!"
  )
  stopifnot(length(x) == length(y))

  n <- length(y)
  thr.inner <- stats::qnorm(cumsum(freq)/n)[-length(freq)]
  thr <- c(-Inf, thr.inner, Inf)
  tau0 <- thr[y]
  tau1 <- thr[y+1]
  finite.lower <- is.finite(tau0)
  finite.upper <- is.finite(tau1)

  cache.rho  <- NA_real_ # for now
  z          <- (x - mean(x)) / stats::sd(x)
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

    plower <<- stats::pnorm(lower, sd = sigma)
    pupper <<- stats::pnorm(upper, sd = sigma)
    prob   <<- pupper - plower

    logy <<- log(prob)

    rho.var <- rho / var
    zrlower <<- z - rho.var * lower
    zrupper <<- z - rho.var * upper

    invisible(NULL)
  }

  plsPolyserialObjective <- function(rho) {
    if (!is.finite(rho))
      return(NaN)

    updateCache(rho)
    -sum(logy)
  }

  plsPolyserialGradient <- function(rho) {
    if (!is.finite(rho))
      return(NaN)

    updateCache(rho)

    lowerTerm <- stats::dnorm(lower, 0, sigma) * zrlower
    upperTerm <- stats::dnorm(upper, 0, sigma) * zrupper

    lowerTerm[!finite.lower] <- 0
    upperTerm[!finite.upper] <- 0

    -sum((lowerTerm - upperTerm) / prob)
  }

  if (is.null(start)) {
    # Starting values from Olsson 1982 eq 38
    cor.xy <- stats::cor(x, y)
    sd.y   <- stats::sd(y) * sqrt((n - 1) / n)
    start  <- cor.xy * sd.y / sum(stats::dnorm(thr.inner))

    if (!is.finite(start) || abs(start) > maxrho)
      start <- cor.xy
  }

  if (!is.finite(start) || abs(start) > maxrho)
    start <- 0.0

  # try 1
  optim <- .nlminb(
    objective = plsPolyserialObjective,
    gradient = plsPolyserialGradient,
    start = start, control = control,
    lower = -abs(maxrho), upper = abs(maxrho)
  )

  # try 2
  if (optim$convergence != 0L) {
    # try again, with different starting value
    retry <- .nlminb(
      objective = plsPolyserialObjective,
      gradient = plsPolyserialGradient,
      start = 0.0, control = control,
      lower = -abs(maxrho), upper = abs(maxrho)
    )
    if (!is.na(retry$par))
      optim <- retry
  }

  # check convergence
  pls_warnif(optim$convergence != 0L,
    "estimation of polyserial correlation did not converge!",
    "Message:", optim$message
  )

  if (is.na(optim$par)) start else optim$par
}
