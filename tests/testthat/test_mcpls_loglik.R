devtools::load_all()

sim.m <- '
  X =~ 0.8 * x1 + 0.6 * x2 + 0.7 * x3
  Z =~ 0.8 * z1 + 0.6 * z2 + 0.7 * z3
  Y =~ 0.8 * y1 + 0.6 * y2 + 0.7 * y3
  Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
  X ~~ 0.2 * Z
'

m.start.low <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + start(0) * X:Z
'

m.start.high <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + start(0.5) * X:Z
'

m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z #+ X:Z
'

# Continous data
set.seed(79624)
parSim <- modsem::modsemify(sim.m)
parSim$est <- as.numeric(parSim$mod)

R <- 400
N <- c(150, 200, 300, 1000, 2000)
pars <- c("Y~X", "Y~Z", "Y~X:Z")
thr <- -1.04

results <- NULL
for (i in seq_len(R)) {
  printf("Iter: %4d/%4d", i, R)
  results.i <- NULL

  sim.cont <- simulateDataParTable(parSim, N = max(N))$ov
  sim.ord <- as.data.frame(lapply(sim.cont, \(x) as.integer(x > thr)))

  for (n in N) {
    printf(", n=%d", n)
    par <- coef(pls(m, data = sim.cont[1:n,], consistent = TRUE))
    results.ij <- data.frame(n = n, iter = i, par = names(par), est = unname(par))
    results.i <- rbind(results.i, results.ij)
  }

  printf(".\n")
  results <- rbind(results, results.i)
}

library(dplyr)
results |>
  filter(par=="Y~Z") |>
  group_by(n) |>
  summarize(mu = mean(est), med = median(est), se = sd(est), ci.lower = -1.96 * se + mu, ci.upper = 1.96 * se + mu)
