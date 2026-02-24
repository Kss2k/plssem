devtools::load_all()

m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z + X:X
'

# summary(modsem::modsem(m, oneIntOrdered, "lms"), standardized = TRUE)
fit <- pls(m, modsem::oneInt, bootstrap = TRUE, sample = 50)
summary(fit)


fit <- pls(m, oneIntOrdered, bootstrap = FALSE)
summary(fit)


fit <- mcem_nlin_ord_pls(m, oneIntOrdered, mc.reps = 1e4)
