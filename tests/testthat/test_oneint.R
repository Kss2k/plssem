devtools::load_all()

m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z + X:X
'

testthat::expect_no_error({
  fit <- pls(m, modsem::oneInt, bootstrap = TRUE, boot.R = 100,
             boot.parallel = "snow", boot.ncpus = 2)
  summary(fit)
})


testthat::expect_no_error({
  fit <- pls(m, oneIntOrdered, bootstrap = TRUE, boot.R = 50,
             boot.parallel = "multicore", boot.ncpus = 4)
  summary(fit)
})

pls_predict(fit, approach = "earliest")
