devtools::load_all()

sim.m <- '
  X =~ 0.5 * x1 + 0.5 * x2 + 0.5 * x3
  Z =~ 0.5 * z1 + 0.5 * z2 + 0.5 * z3
  Y =~ 0.5 * y1 + 0.5 * y2 + 0.5 * y3
  Y  ~ 0.4 *  X + 0.5 *  Z + 0.3 * X:Z
  X ~~ 0.2 * Z
'

m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'

# Ordinal data
set.seed(79624)
parSim <- modsem::modsemify(sim.m)
parSim$est <- as.numeric(parSim$mod)
sim.cont <- simulateDataParTable(parSim, N = 200)$ov
sim.ord  <- as.data.frame(lapply(sim.cont, \(x) as.integer(x > 0)))

fit <- pls(m, data = sim.ord, ordered = colnames(sim.ord))
bc  <- suppressWarnings(
  mcpls_bias_correct(fit, B = 5, seed = 123, verbose = FALSE, clamp = FALSE))

testthat::expect_equal(names(bc$est), names(bc$est.bc))
testthat::expect_true(bc$B >= 2 && bc$B <= 5)
testthat::expect_equal(NROW(bc$replicates), bc$B)

# The Newton step itself
testthat::expect_equal(as.numeric(bc$est.bc), as.numeric(bc$est) - as.numeric(bc$bias))
testthat::expect_equal(as.numeric(bc$bias),
                       unname(colMeans(bc$replicates)) - as.numeric(bc$est))

# Only the free parameters are corrected; thresholds and residual variances
# are recomputed by the estimator instead.
testthat::expect_true(all(!grepl("\\|", names(bc$est))))

# The correction must not depend on state left behind by a previous call:
# the same seed reproduces the same replicates.
bc2 <- suppressWarnings(
  mcpls_bias_correct(fit, B = 5, seed = 123, verbose = FALSE, clamp = FALSE))
testthat::expect_equal(as.numeric(bc$bias), as.numeric(bc2$bias))

testthat::expect_error(mcpls_bias_correct(fit, B = 1), "`B` must be")
