devtools::load_all()

m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z + X:X
'

# summary(modsem::modsem(m, oneIntOrdered, "lms"), standardized = TRUE)
fit <- pls(m, modsem::oneInt, bootstrap = TRUE, boot.R = 100,
           boot.parallel = "snow", boot.ncpus = 2)
summary(fit)



fit <- pls(m, oneIntOrdered, bootstrap = TRUE, boot.R = 500,
           boot.parallel = "multicore", boot.ncpus = 4,
           boot.optimize = TRUE)
summary(fit)
