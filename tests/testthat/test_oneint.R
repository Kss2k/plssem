devtools::load_all()
library(cSEM)

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

mcsem <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X.Z + X.X
'

fit.csem <- csem(.model = mcsem, .data = modsem::oneInt, .disattenuate = TRUE,
               .resample_method = "none")
summarize(fit.csem)
calculateR2(fit.csem)
