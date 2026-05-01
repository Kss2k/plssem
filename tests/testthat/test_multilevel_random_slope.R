devtools::load_all()

syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  W =~ w1 + w2 + w3
  Y ~ X + Z + (1 + X + Z | cluster)
  W ~ X + Z + (1 + X + Z | cluster)
"

testthat::expect_no_error({
  fit.c <- pls(syntax, data = randomSlopes,
             consistent = TRUE, bootstrap = TRUE, mcpls = TRUE,
               mc.fixed.seed = TRUE)
  summary(fit.c)
})

library(stringr)
library(lme4)

lmer(y1 ~ x1 + z1 + (x1 + z1 | cluster), data = randomSlopes,
     control = lmerControl(calc.derivs = FALSE,
                           optCtrl = list(xtol_abs = 1e-2, ftol_abs = 1e-2,
                                          max_eval = 1, max_iter = 1, maxit = 1,
                                          maxfun = 0)))

syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  W =~ w1 + w2 + w3
  Y ~ X + Z + (1 + X + Z | cluster)
  W ~ X + Z + (1 + X + Z | cluster)
"

testthat::expect_no_error({
  fit.o <- pls(syntax, data = randomSlopesOrdered,
             consistent = TRUE, bootstrap = TRUE, mcpls = TRUE)
  summary(fit.o)
})

# With interaction terms

syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  W =~ w1 + w2 + w3
  Y ~ X + Z + X:Z + (1 + X + Z + X:Z | cluster)
  W ~ X + Z + X:Z + (1 + X + Z + X:Z | cluster)
"

testthat::expect_no_error({
  fit.c <- pls(syntax, data = randomSlopes,
               consistent = TRUE, bootstrap = FALSE)
  summary(fit.c)
})
