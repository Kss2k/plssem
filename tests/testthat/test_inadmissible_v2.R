devtools::load_all()

sim.m <- '
  X =~ 0.5 * x1 + 0.5 * x2 + 0.5 * x3
  Z =~ 0.5 * z1 + 0.5 * z2 + 0.5 * z3
  Y =~ 0.5 * y1 + 0.5 * y2 + 0.5 * y3
  Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
  X ~~ 0.2 * Z
'

parTable <- modsem::modsemify(sim.m)
parTableSim <- parTableClean <- parTable

parTableClean$mod <- ""
parTableSim$est <- as.numeric(parTableSim$mod)

mod.m <- modsem:::parTableToSyntax(parTableClean)

# Continous data
set.seed(2341039)
sim.cont <- simulateDataParTable(parTableSim, N = 200)$ov

# Ordinal data
sim.ord <- as.data.frame(lapply(sim.cont, \(x) as.integer(x > 0)))

fit <- pls(mod.m, data = sim.ord, ordered = colnames(sim.ord))
testthat::expect_true(is_admissible(fit))
