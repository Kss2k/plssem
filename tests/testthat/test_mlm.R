devtools::load_all()

# Standardized estimates of the loadings and regressions of `fit.pls`, and the
# corresponding (standardized) estimates of `fit.lav`, matched by parameter and level
compareStdSolutions <- function(fit.pls, fit.lav, levels = 1:2) {
  pt.pls <- parameter_estimates(fit.pls)
  pt.pls <- pt.pls[pt.pls$op %in% c("=~", "~") & pt.pls$level %in% levels, ]

  pt.lav <- lavaan::parameterEstimates(fit.lav, standardized = TRUE)

  key <- \(pt) paste(pt$lhs, pt$op, pt$rhs, pt$level)

  data.frame(
    par    = paste(pt.pls$lhs, pt.pls$op, pt.pls$rhs),
    level  = pt.pls$level,
    pls    = pt.pls$est,
    lavaan = pt.lav$std.all[match(key(pt.pls), key(pt.lav))]
  )
}

model <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'

set.seed(23124)
testthat::expect_no_error({
  fit.pls <- pls(model, data = randomSlopes, cluster = "cluster", bootstrap = TRUE)
  fit.lav <- lavaan::sem(model, data = randomSlopes, cluster = "cluster")
})

testthat::test_that("pls() and lavaan give similar standardized estimates (two-level model)", {
  est <- compareStdSolutions(fit.pls, fit.lav)
  testthat::expect_false(anyNA(est$lavaan))

  # level 1 has many more observations than level 2, so we use a stricter
  # tolerance at level 1 (the standard errors are approx. 0.01-0.03 at level 1,
  # and 0.04-0.09 at level 2)
  diff <- abs(est$pls - est$lavaan)
  testthat::expect_lt(max(diff[est$level == 1]), 0.02)
  testthat::expect_lt(max(diff[est$level == 2]), 0.05)
})

modelr <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ rv("s1")*x1 + rv("s2")*x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
        # the random slopes are latent variables at the between level
        s1 + s2 ~ w1 + w2
'

set.seed(23984)
testthat::expect_no_error({
  fit.pls <- pls(modelr, data = randomSlopes, cluster = "cluster", bootstrap = TRUE)
  fit.lav <- lavaan::sem(modelr, data = randomSlopes, cluster = "cluster")
}) 

testthat::test_that("pls() and lavaan give similar standardized estimates (random slopes)", {
  # Only level 2 is compared. In lavaan the fixed effects of the random slopes
  # are the intercepts of s1 and s2 (`fw ~ x1` and `fw ~ x2` are fixed to zero),
  # and the standardized estimates at level 1 are based on a different variance
  # of fw, such that the level 1 estimates aren't comparable
  est <- compareStdSolutions(fit.pls, fit.lav, levels = 2)
  testthat::expect_false(anyNA(est$lavaan))
  testthat::expect_true(all(c("s1 ~ w1", "s2 ~ w2") %in% est$par))

  testthat::expect_lt(max(abs(est$pls - est$lavaan)), 0.05)
})
