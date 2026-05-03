robbinsMonro1951 <- function(p, f, tol, min.iter, max.iter, verbose,
                             polyak.juditsky = FALSE, pj.extrapolate = TRUE, fn.args, ...) {
  # Wrapper for SimDesign::RobbinsMonro/SimDesign.RobbinsMonro

  args.required <- list(
      p               = p,
      f               = f,
      tol             = tol,
      miniter         = min.iter,
      maxiter         = max.iter,
      verbose         = verbose,
      Polyak_Juditsky = polyak.juditsky
  )

  args <- c(args.required, fn.args, list(...))

  # mcfit <- do.call(SimDesign::RobbinsMonro, args)
  mcfit <- do.call(SimDesign.RobbinsMonro, args)
  if (verbose) cat("\n")

  if (polyak.juditsky && pj.extrapolate)
    mcfit$root <- getConvergencePoints(mcfit$history)

  mcfit
}


# The code below is adapted from the `SimDesign` package
# https://github.com/philchalmers/SimDesign/blob/main/R/RobbinsMonro.R
# The SimDesign package depends on the `qs2` package, which currently
# 04.25.2026 has some installation issues on r-release-macos-x86_64
# We only use the RobbinsMonro function, which does not depend on the
# `qs2` package. In the future we might revert back to `SimDesign::RobbinsMonro`.
SimDesign.RobbinsMonro <- function(f, p, ...,
                                   Polyak_Juditsky = FALSE,
                                   maxiter = 500L, miniter = 100L, k = 3L,
                                   tol = .00001, verbose = interactive(),
                                   fn.a = function(iter, a = 1, b = 1/2, c = 0, ...)
                                     a / (iter + c)^b) {
  if(maxiter < miniter) maxiter <- miniter
  history <- rbind(p, matrix(NA, nrow=maxiter, ncol=length(p)))
  k.succ <- 0
  pbar_last <- pbar <- p

  for(i in 1L:maxiter){
    a <- fn.a(iter=i, ...)
    fp <- f(p)
    p <- p - a * fp
    history[i + 1L, ] <- p
    change <- max(abs(history[i,]-p))
    if(Polyak_Juditsky){
      pbar_last <- pbar
      pbar <- PK_average(history)
      change <- max(abs(pbar_last - pbar))
    }
    if(verbose){
      if(Polyak_Juditsky)
        cat(sprintf("\rIter: %i; Max change in E(p) = %.3f",
                    i, change))
      else
        cat(sprintf("\rIter: %i; Max change in p = %.3f",
                    i, change))
      utils::flush.console()
    }
    if(i > miniter && all(change < tol)){
      k.succ <- k.succ + 1L
      if(k.succ == k) break
    } else k.succ <- 0L
  }

  converged <- i < maxiter
  history <- history[0L:i + 1L, , drop=FALSE]
  ret <- list(iter=i, root=if(Polyak_Juditsky) pbar else p,
              terminated_early=converged,
              history=history, Polyak_Juditsky=Polyak_Juditsky)
  class(ret) <- 'RM'
  ret
}


PK_average <- function(history) {
  ret <- colMeans(history, na.rm=TRUE)
  matrix(ret, ncol=ncol(history))
}


getConvergencePoint <- function(y, t = seq_along(y)) {
  c.aitken <- aitkenAccelerate(y)

  tryCatch({
    fit <- stats::nls(
      y ~ c + a * exp(-k * t),
      start = list(
        c = mean(utils::tail(y, 3)),
        a = y[1] - mean(utils::tail(y, 3)),
        k = 0.1
      ),
      algorithm = "port",
      lower = c(c = -Inf, a = -Inf, k = 0)
    )

    c.nls <- coef(fit)[["c"]]
    span  <- diff(range(y))
    k.fit <- coef(fit)[["k"]]

    # Fall back to Aitken if nls converged to an implausible solution:
    # wild extrapolation beyond the observed range, near-zero decay rate
    # (c poorly identified), or large disagreement with Aitken
    bad.range <- c.nls < min(y) - span || c.nls > max(y) + span
    bad.k     <- k.fit < sqrt(.Machine$double.eps)
    bad.agree <- abs(c.nls - c.aitken) > 0.5 * span

    if (bad.range || bad.k || bad.agree) c.aitken else c.nls
  }, error = \(e) c.aitken)
}


# Aitken's delta^2 sequence acceleration applied to all triplets in the history,
# returning the median over valid estimates. More robust than a single nls()
# call when trajectories are noisy: no tuning needed, no optimisation.
# The median naturally discards wild estimates from transient-phase triplets,
# so no tail-window selection is needed.
aitkenAccelerate <- function(y) {
  n <- length(y)
  if (n < 3L) return(mean(y, na.rm = TRUE))

  ests <- vapply(X = seq_len(n - 2L), FUN = function(i) {

    p0 <- y[i]
    p1 <- y[i + 1L]
    p2 <- y[i + 2L]

    denom <- p2 - 2 * p1 + p0
    if (abs(denom) < .Machine$double.eps^0.5)
      return(NA_real_)

    p0 - (p1 - p0)^2 / denom

  }, FUN.VALUE = numeric(1L))

  valid <- ests[is.finite(ests)]
  if (length(valid) == 0L)
    return(mean(y, na.rm = TRUE))

  stats::median(valid)
}


getConvergencePoints <- function(history) {
  history <- history[
    stats::complete.cases(history), , drop = FALSE
  ]

  apply(X = history, MARGIN = 2, FUN = getConvergencePoint)
}
