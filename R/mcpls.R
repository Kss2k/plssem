#' @export
mcpls <- function(
  syntax,
  data,
  ordered = NULL,
  consistent = FALSE, # we get consistent estimates as a natural by-product
  max.iter.mc = 100,
  mc.reps = 20000,
  rng.seed = NULL,
  tol = 1e-3,
  miniter = 25,
  maxiter = 250,
  ...
) {
  data <- as.data.frame(data)
  orderedData <- colnames(data)[sapply(data, FUN = is.ordered)]
  data[orderedData] <- lapply(data[orderedData], FUN = as.integer)

  fit0 <- pls(syntax = syntax, data = data, consistent = consistent, ...)
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

  .f <- function(p) {
    par1$est <- p
    
    sim <- simulateDataParTable(par1, N = mc.reps, seed = rng.seed)
    sim.ov <- sim$ov
    for (ord in ordered)
      sim.ov[, ord] <- ordinalize(sim.ov[, ord], probs = PROBS[[ord]])

    fit0$data <- scale(sim.ov[vars])
    fit0$matrices$S <- cov(fit0$data)

    fit2 <- estimatePLS(fit0)
    par2 <- getFreeParamsTable(getParTableEstimates(fit2))
 
    par2$est - par0$est
  }


  par1$est <- SimDesign::RobbinsMonro(p = par1$est,
                                      f = .f, tol = tol,
                                      miniter = miniter,
                                      maxiter = maxiter)$root
  cat("\n")

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
