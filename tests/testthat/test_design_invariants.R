devtools::load_all()

testthat::test_that("fit C and matrices C stay in sync", {
  testthat::skip_if_not_installed("modsem")

  syntax <- "
    X =~ x1 + x2 + x3
    Z =~ z1 + z2 + z3
    Y =~ y1 + y2 + y3

    Y ~ X + Z
  "

  fit_pls <- pls(syntax, data = modsem::oneInt, consistent = FALSE)
  testthat::expect_equal(fit_pls@matrices$C, fit_pls@fit$fitC)

  fit_plsc <- pls(syntax, data = modsem::oneInt, consistent = TRUE)
  testthat::expect_equal(fit_plsc@matrices$C, fit_plsc@fit$fitC)
})


testthat::test_that("fit measures uses MCPLS flag via is_mcpls", {
  testthat::skip_if_not_installed("modsem")

  syntax <- "
    X =~ x1 + x2 + x3
    Z =~ z1 + z2 + z3
    Y =~ y1 + y2 + y3

    Y ~ X + Z + X:Z
  "

  fit <- pls(syntax, data = modsem::oneInt, mcpls = FALSE)
  testthat::expect_false(is_mcpls(fit))
  testthat::expect_type(fit_measures(fit)$chisq, "double")
})

