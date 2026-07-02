devtools::load_all()

sim.m <- '
  X =~ 0.5 * x1 + 0.5 * x2 + 0.5 * x3
  Z =~ 0.5 * z1 + 0.5 * z2 + 0.5 * z3
  Y =~ 0.5 * y1 + 0.5 * y2 + 0.5 * y3
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

# Continous data
set.seed(79624)
parSim <- modsem::modsemify(sim.m)
parSim$est <- as.numeric(parSim$mod)
sim.cont <- simulateDataParTable(parSim, N = 200)$ov

# Ordinal data
sim.ord <- as.data.frame(lapply(sim.cont, \(x) as.integer(x > 0)))

fit.start.low <- pls(m.start.low, data = sim.ord, ordered = colnames(sim.ord))
fit.start.high <- pls(m.start.high, data = sim.ord, ordered = colnames(sim.ord))

mcplsLoglik(fit.start.low)
mcplsLoglik(fit.start.high)
