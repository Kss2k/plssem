rescaleOrderedVariableMonteCarlo <- function(name,
                                             data,
                                             sim.ov,
                                             smooth_eps = 0,
                                             eps_expand = 1e-12,
                                             standardize = TRUE) {
  x   <- as.ordered(data[[name]])
  x.i <- as.integer(x)
  y   <- sim.ov[, name]

  # standardize for scale stability (no distributional assumption)
  y_mean <- mean(y, na.rm = TRUE)
  y_sd   <- stats::sd(y,  na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
  y <- (y - y_mean) / y_sd

  # observed category proportions (optional Laplace smoothing)
  tab <- as.numeric(table(x))
  K   <- length(tab)
  n   <- sum(tab)
  if (smooth_eps > 0) {
    p_obs <- (tab + smooth_eps) / (n + K * smooth_eps)
  } else {
    p_obs <- tab / n
  }

  # empirical CDF from data; quantile cutpoints on simulated y
  cdf_obs <- cumsum(p_obs[-length(p_obs)]) # drop last, we know that it sums to 1
  q <- stats::quantile(y, probs = cdf_obs, names = FALSE, type = 7)
  q <- c(-Inf, q, Inf) # [0, ..., last(p_obs)]

  # compute conditional means in each [q[i], q[i+1]) (last bin inclusive)
  mu <- numeric(K)
  for (i in seq_len(K)) {
    lo <- q[i]
    hi <- q[i + 1]

    # left-closed/right-open, last bin inclusive
    if (i < K) {
      idx <- (y >= lo & y < hi)
    } else {
      idx <- (y >= lo & y <= hi)
    }

    # guard for degenerate bins (ties / finite-N quantile artifacts)
    if (!any(idx)) {
      # expand a hair around the cutpoints to catch ties
      if (i < K) {
        idx <- (y >= (lo - eps_expand) & y < (hi + eps_expand))
      } else {
        idx <- (y >= (lo - eps_expand) & y <= (hi + eps_expand))
      }
      # still empty? fallback to midpoint
      if (!any(idx)) mu[i] <- 0.5 * (lo + hi)
    }

    if (any(idx)) mu[i] <- mean(y[idx])
  }

  # exact centering to remove residual drift (weighted by observed category probs)
  mu <- mu - sum(p_obs * mu, na.rm = TRUE)

  # map each observed category to its conditional mean
  x.out <- rep(NA_real_, length(x))
  for (i in seq_len(K))
    x.out[x.i == i] <- mu[i]

  if (standardize)
    x.out <- standardizeAtomic(x.out)

  labels.t <- paste0(name, "|t", seq_len(sum(is.finite(q))))
  thresholds <- stats::setNames(q[is.finite(q)], labels.t)

  list(values = x.out, thresholds = thresholds)
}


rescaleOrderedVariableAnalytic <- function(name,
                                           data,
                                           smooth_eps = 0,
                                           standardize = TRUE) {
  x   <- as.ordered(data[[name]])
  x.i <- as.integer(x)

  # observed category proportions (optional Laplace smoothing)
  tab <- as.numeric(table(x))
  K   <- length(tab)
  n   <- sum(tab)
  if (K < 2) stop("Need at least 2 ordered categories for '", name, "'.")

  if (smooth_eps > 0) {
    p_hat <- (tab + smooth_eps) / (n + K * smooth_eps)
  } else {
    p_hat <- tab / n
  }

  # cumulative probs for interior thresholds (K-1 of them)
  C <- cumsum(p_hat)[-K]  # drop last; sums to 1
  # numeric guards
  eps <- .Machine$double.eps^0.5
  C <- pmin(pmax(C, eps), 1 - eps)

  # thresholds on standard-normal scale
  t_int <- stats::qnorm(C)                 # length K-1
  q     <- c(-Inf, t_int, Inf)             # length K+1 for interval bounds

  # truncated-normal means for each category interval (a_i, b_i]
  # mu_i = (phi(a_i) - phi(b_i)) / (Phi(b_i) - Phi(a_i))
  a <- q[1:K]
  b <- q[2:(K+1)]
  Phi_a <- stats::pnorm(a)
  Phi_b <- stats::pnorm(b)
  phi_a <- stats::dnorm(a)
  phi_b <- stats::dnorm(b)
  denom <- pmax(Phi_b - Phi_a, eps)
  mu    <- (phi_a - phi_b) / denom

  # exact centering to remove small drift from smoothing/finite-n
  mu <- mu - sum(p_hat * mu)

  # map each observed category to its conditional mean
  x.out <- rep(NA_real_, length(x.i))
  for (i in seq_len(K)) x.out[x.i == i] <- mu[i]

  # optional standardization of the mapped scores (sample standardization)
  if (standardize)
    x.out <- standardizeAtomic(x.out)

  # labels for interior thresholds
  labels.t   <- paste0(name, "|t", seq_len(K - 1))
  thresholds <- stats::setNames(t_int, labels.t)

  list(values = x.out, thresholds = thresholds)
}
