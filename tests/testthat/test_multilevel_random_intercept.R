devtools::load_all()

syntax <- '
  f =~ y1 + y2 + y3
  f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
'

testthat::expect_no_error({
  fit <- pls(syntax, data = randomIntercepts, bootstrap = FALSE, boot.R = 50)
  summary(fit)
})


syntax <- '
  f =~ y1 + y2 + y3
  f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
'

testthat::expect_no_error({
  fit <- pls(syntax, data = randomInterceptsOrdered, bootstrap = TRUE, boot.R = 50)
  summary(fit)
})
