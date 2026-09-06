# The code below was originally adapted from the `SimDesign` package
# https://github.com/philchalmers/SimDesign/blob/main/R/RobbinsMonro.R
# The SimDesign package depends on the `qs2` package, which currently
# 04.25.2026 has some installation issues on r-release-macos-x86_64
# We only use the RobbinsMonro function, which does not depend on the
# `qs2` package. In the future we might revert back to `SimDesign::RobbinsMonro`.
robbinsMonro1951 <- function(p,
                             f,
                             ...,
                             k = 3,
                             tol = 0.00001,
                             min.iter = 100L,
                             max.iter = 500L,
                             lower = NULL,
                             upper = NULL,
                             verbose = interactive(),
                             polyak.juditsky = FALSE,
                             pj.extrapolate = TRUE,
                             fn.a = \(iter, a = 1, b = 1/2, c = 0, ...) a / (iter + c)^b,
                             fn.args = list(),
                             k.dyn.bound = 5,
                             diverge.factor = 10,
                             diverge.min.iter = 10L,
                             diverge.ema.decay = 0.9,
                             diag.secant = FALSE,
                             ds.eps = 1e-8,
                             ds.clip = 10) {

  if (max.iter < min.iter)
    max.iter <- min.iter

  history.p <- history.f <- rbind(p, matrix(NA, nrow=max.iter, ncol=length(p)))
  k.succ <- 0
  pbar.last <- pbar <- p
  diverged <- FALSE

  # Track the best (smallest-residual) point seen so far so a run that starts
  # converging normally but later diverges (e.g. a non-monotone `f()`) can be caught
  best.p           <- p
  best.resid.norm  <- Inf

  # Divergence is judged on an EMA-smoothed residual norm.
  resid.ema      <- NULL
  best.resid.ema <- Inf

  # Diagonal secant method, which can flip a coordinate's step
  # direction automatically when its local relationship is non-monotone.
  p.prev  <- NULL
  fp.prev <- NULL

  if (is.null(lower)) lower <- rep(-Inf, length(p))
  if (is.null(upper)) upper <- rep( Inf, length(p))

  for (i in seq_len(max.iter)) {
    a  <- do.call(fn.a, c(list(iter = i), fn.args))
    fp <- f(p, ...)

    resid.norm <- sqrt(sum(fp^2))
    if (resid.norm < best.resid.norm) {
      best.resid.norm <- resid.norm
      best.p          <- p
    }

    if (is.null(resid.ema)) {
      resid.ema <- resid.norm
    } else {
      resid.ema <- diverge.ema.decay * resid.ema + 
        (1 - diverge.ema.decay) * resid.norm
    }

    if (resid.ema < best.resid.ema)
      best.resid.ema <- resid.ema

    if (diag.secant && !is.null(p.prev)) {
      delta.p  <- p - p.prev
      delta.fp <- fp - fp.prev
      a.ds     <- delta.p / delta.fp

      # fall back to the standard step if the secant is degenerate
      bad <- !is.finite(a.ds) | abs(delta.fp) < ds.eps
      a.ds[bad] <- a

      # keep the secant step from taking a wild jump
      a.max <- ds.clip * a
      a.ds  <- pmax(pmin(a.ds, a.max), -a.max)

      step <- a.ds * fp

    } else {
      step <- a * fp # standard step, or no secant yet on the first iteration
    }

    if (diag.secant) {
      p.prev  <- p
      fp.prev <- fp
    }

    p <- p - step

    # Does f() return (dynamic) bounduaries?
    lower.fp <- attr(fp, "lower")
    upper.fp <- attr(fp, "upper")

    # Combine bounds
    lower.i <- if (!is.null(lower.fp)) pmax(lower, lower.fp) else lower
    upper.i <- if (!is.null(upper.fp)) pmin(upper, upper.fp) else upper

    # clamp
    p[p < lower.i] <- lower.i[p < lower.i]
    p[p > upper.i] <- upper.i[p > upper.i]

    history.p[i + 1L, ] <- p
    history.f[i + 1L, ] <- fp
    change <- max(abs(history.p[i,]-p))

    if (polyak.juditsky) {
      pbar.last <- pbar
      pbar <- pjRunningAverage(history.p)
      change <- max(abs(pbar.last - pbar))
    }

    if (verbose) {
      if (polyak.juditsky) {
        msg <- MSG_STRINGS$strings$SimDesign.RobbinsMonro0
        messagef(msg, i, change)
      } else {
        msg <- MSG_STRINGS$strings$SimDesign.RobbinsMonro1
        messagef(msg, i, change)
      }
    }

    if (i > min.iter && all(change < tol)) {
      k.succ <- k.succ + 1L
      if (k.succ == k) break
    } else k.succ <- 0L

    if (i > diverge.min.iter && resid.ema > diverge.factor * best.resid.ema) {
      diverged <- TRUE
      break
    }
  }

  converged <- i < max.iter && !diverged
  history.p <- history.p[0L:i + 1L, , drop=FALSE]
  history.f <- history.f[0L:i + 1L, , drop=FALSE]

  if (verbose) messagef("\n")

  if (polyak.juditsky) {

    # Get approximation, irrespective of the dynamic bounds
    if (pj.extrapolate) {
      pbar <- pbar0 <- getConvergencePoints(
        history = history.p, lower = lower, upper = upper
      )
    } else {
      pbar0 <- pbar
    }

    # Check dynamic bounds
    has.dyn.bound <- !is.null(lower.fp) || !is.null(upper.fp)

    if (has.dyn.bound) for (.i in seq_len(k.dyn.bound)) {
      fp <- f(pbar, ...)

      # Get dynamic bounds
      lower.fp <- attr(fp, "lower")
      upper.fp <- attr(fp, "upper")

      # Combine bounds
      lower.i <- if (!is.null(lower.fp)) pmax(lower, lower.fp) else lower
      upper.i <- if (!is.null(upper.fp)) pmin(upper, upper.fp) else upper

      if (all(pbar >= lower.i & pbar <= upper.i))
        break # nothing needs to be done

      # see if any previously discarded values in pbar0 are
      # within the current bounds
      inside0 <- pbar0 >= lower.i & pbar0 <= upper.i
      pbar[inside0] <- pbar0[inside0]

      # clamp and try again
      pbar[pbar < lower.i] <- lower.i[pbar < lower.i]
      pbar[pbar > upper.i] <- upper.i[pbar > upper.i]
    }
  }

  # if we have a divergent solution, we should set p to best.p
  # only relevant if polyak.juditsky=FALSE
  if (diverged)
    p <- best.p

  ret <- list(
    iter            = i,
    root            = if (polyak.juditsky) pbar else p,
    history.p       = history.p,
    history.f       = history.f,
    lower           = lower.i,
    upper           = upper.i,
    converged       = converged,
    polyak.juditsky = polyak.juditsky,
    diverged        = diverged,
    resid.norm      = best.resid.norm,
    best.p          = best.p
  )

  ret
}


pjRunningAverage <- function(history) {
  colMeans(history, na.rm=TRUE)
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


getConvergencePoints <- function(history, lower = NULL, upper = NULL) {
  history <- history[
    stats::complete.cases(history), , drop = FALSE
  ]

  pbar <- apply(X = history, MARGIN = 2, FUN = getConvergencePoint)

  if (!is.null(lower))
    pbar[pbar < lower] <- lower[pbar < lower]

  if (!is.null(upper))
    pbar[pbar > upper] <- upper[pbar > upper]

  pbar
}
