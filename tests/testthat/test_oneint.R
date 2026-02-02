devtools::load_all()
library(cSEM)

m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

fit <- pls(m, modsem::oneInt, bootstrap = TRUE)
summary(fit)


fit <- pls(m, oneIntOrdered, bootstrap = FALSE)
summary(fit)

mcsem <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X.Z
'

summarize(csem(.model = mcsem, .data = oneIntOrdered, .disattenuate = TRUE))
