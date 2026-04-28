robbinsMonro1951 <- function(p, f, tol, min.iter, max.iter, verbose,
                             polyak.juditsky, fn.args, ...) {
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
  t <- sum(rowSums(!is.na(history)) > 0L)
  ret <- colSums(history[1:(t - 1L), , drop=FALSE], na.rm=TRUE) / t
  matrix(ret, ncol=ncol(history))
}
