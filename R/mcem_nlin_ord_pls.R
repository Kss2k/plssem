#' @export
mcem_nlin_ord_pls <- function(syntax, data, ordered = NULL, max.iter.mc = 100, mc.reps = 1e4, rng.seed = 2983472,
                              tol = 1e-4, ...) {
  data <- as.data.frame(data)
  orderedData <- colnames(data)[sapply(data, FUN = is.ordered)]
  data[orderedData] <- lapply(data[orderedData], FUN = as.integer)

  fit0 <- pls(syntax=syntax, data=data, ...)

  info <- fit0$info
  data <- as.data.frame(fit0$data)
  vars <- colnames(data)
  ordered <- intersect(vars, union(orderedData, ordered))

  PROBS <- stats::setNames(vector("list", length(ordered)), nm = ordered)

  for (ord in ordered) {
    freq <- table(data[[ord]])
    pct  <- cumsum(freq) / sum(freq)
    PROBS[[ord]] <- pct[-length(pct)]
  }

  par0 <- getFreeParamsTable(parameter_estimates(fit0))
  par1 <- par0[c("lhs", "op", "rhs", "est")]

  for (iter in seq_len(max.iter.mc)) {
    temperature <- 1 - (iter - 1) / max.iter.mc # this could be improved

    sim <- simulateDataParTable(par1, N = mc.reps, seed = rng.seed)
    sim.ov <- sim$ov

    for (ord in ordered)
      sim.ov[,ord] <- ordinalize(sim.ov[,ord], probs = PROBS[[ord]])

    fit2 <- pls(syntax=syntax, data=sim.ov, ...)
    par2 <- getFreeParamsTable(parameter_estimates(fit2))


    diff <- par2$est - par0$est
    objective <- mean(abs(diff))
    par1$est <- par1$est - temperature * diff

    printf("\rIter: %i, diff: %.2g", iter, objective)

    if (objective < tol) {
      printf("\nSolution converged!\n")
      break
    }
  }

  if (iter == max.iter.mc)
    printf("\nSolution dit not converge!")

  par1
}


ordinalize <- function(x, probs) {
  probs <- sort(probs[probs < 1]) # skip 100%, we analytically know it to be Inf
  breaks <- c(-Inf, quantile(x, probs = probs), Inf)

  cut(x, breaks = breaks, labels = FALSE)
}


getFreeParamsTable <- function(parTable) {
  lhs <- parTable$lhs
  op  <- parTable$op
  rhs <- parTable$rhs

  cond1 <- !(lhs == rhs & op == "~~")
  cond2 <- !((grepl(":", lhs) | grepl(":", rhs)) & op == "~~")

  out <- parTable[cond1 & cond2, , drop = FALSE]
  attr(out, "cond") <- cond1 & cond2

  out
}
