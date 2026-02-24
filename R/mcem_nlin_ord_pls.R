mcem_nlin_ord_pls <- function(syntax, data, ordered = NULL, max.iter.mc = 100, mc.reps = 1e4, rng.seed = 2983472, ...) {
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

  par0 <- parameter_estimates(fit0)
  par1 <- par0[c("lhs", "op", "rhs", "est")]

  par1$est.last  <- par1$est - 0.05
  par1$est.new   <- par1$est + 0.05
  par1$delta.raw <- 0.1

  for (iter in seq_len(max.iter.mc)) {
    temperature <- 1 - sqrt((iter - 1) / max.iter.mc) # this could be improved
    par1$delta.rel <- (par1$est.new - par1$est.last) / par1$delta.raw

    sim <- simulateDataParTable(par1, N = mc.reps, seed = rng.seed)
    sim.ov <- sim$OV[[1L]]

    for (ord in ordered)
      sim.ov[,ord] <- ordinalize(sim.ov[,ord], probs = PROBS[[ord]])

    fit2 <- pls(syntax=syntax, data=sim.ov, ...)
    par2 <- parameter_estimates(fit2)


    diff <- par2$est - par0$est
    par1$delta.raw <- diff / par1$delta.rel
    
    zero <- par1$delta.raw == 0 | par1$delta.rel == 0
    par1$delta.raw[zero] <- 0
    par1$delta.rel[zero] <- 1

    par1$est.last  <- par1$est.new
    par1$est.new   <- par2$est
    par1$est       <- par1$est - temperature * diff
    browser()
  }
}


ordinalize <- function(x, probs) {
  probs <- sort(probs[probs < 1]) # skip 100%, we analytically know it to be Inf
  breaks <- c(-Inf, quantile(x, probs = probs), Inf)

  cut(x, breaks = breaks, labels = FALSE)
}
