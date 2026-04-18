testthat::test_that("pls_predict works for continuous models", {
  set.seed(1)

  data <- randomIntercepts[1:300, , drop = FALSE]

  syntax <- "
    f =~ y1 + y2 + y3
    f ~ x1 + x2 + x3 + w1 + w2
  "

  fit <- pls(syntax, data = data, bootstrap = FALSE)

  for (approach in c("earliest", "direct")) {
    for (bm in c("r2", "rmse", "mae", "q2", "q2_predict")) {
      pred <- pls_predict(fit, approach = approach, benchmark = bm)

      testthat::expect_s3_class(pred, "PlsSemPredict")
      testthat::expect_true(is.null(pred$X.ord))
      testthat::expect_true(is.null(pred$X.ord.pred))

      testthat::expect_true(is.data.frame(pred$benchmark))
      testthat::expect_true(all(c("variable", "metric", "value") %in% names(pred$benchmark)))
      testthat::expect_equal(nrow(pred$benchmark), ncol(pred$X.cont.pred))
      testthat::expect_true(all(pred$benchmark$metric %in% c("r2", "rmse", "mae", "q2_predict")))
      testthat::expect_type(pred$benchmark$value, "double")
    }
  }

  pred <- pls_predict(fit, approach = "direct", benchmark = "rmse")
  testthat::expect_output(print(pred), "PlsSemPredict object")

  testthat::expect_error(
    pls_predict(fit, benchmark = "acc"),
    "only available for ordered variables"
  )
  testthat::expect_error(
    pls_predict(fit, benchmark = "ord_mae"),
    "only available for ordered variables"
  )
})


testthat::test_that("pls_predict works for ordinal (non-MCPLS) models", {
  set.seed(1)

  data <- randomInterceptsOrdered[1:300, , drop = FALSE]
  ordered <- colnames(data)

  syntax <- "
    f =~ y1 + y2 + y3
    f ~ x1 + x2 + x3 + w1 + w2
  "

  fit <- pls(syntax, data = data, ordered = ordered, bootstrap = FALSE, mcpls = FALSE)

  for (approach in c("earliest", "direct")) {
    for (bm in c("r2", "rmse", "mae", "q2_predict", "acc", "ord_mae")) {
      if (bm == "r2") {
        pred <- suppressWarnings(pls_predict(fit, approach = approach, benchmark = bm))
      } else {
        pred <- pls_predict(fit, approach = approach, benchmark = bm)
      }

      testthat::expect_s3_class(pred, "PlsSemPredict")
      testthat::expect_true(!is.null(pred$X.ord))
      testthat::expect_true(!is.null(pred$X.ord.pred))

      testthat::expect_true(is.data.frame(pred$benchmark))
      testthat::expect_equal(nrow(pred$benchmark), ncol(pred$X.cont.pred))
      testthat::expect_true(all(pred$benchmark$metric == bm))
    }
  }

  pred <- suppressWarnings(pls_predict(
    fit,
    approach = "earliest",
    benchmark = c(
      y1 = "acc",
      y2 = "ordmae",
      y3 = "r2",
      x1 = "rmse",
      x2 = "mae",
      x3 = "q2",
      w1 = "rmse",
      w2 = "rmse"
    )
  ))

  testthat::expect_setequal(
    unique(pred$benchmark$metric),
    c("acc", "ord_mae", "r2", "rmse", "mae", "q2_predict")
  )
})


testthat::test_that("pls_predict works for MC-PLS ordinal interaction models", {
  set.seed(1)

  data <- oneIntOrdered[1:300, , drop = FALSE]
  ordered <- colnames(data)

  syntax <- "
    X =~ x1 + x2 + x3
    Z =~ z1 + z2 + z3
    Y =~ y1 + y2 + y3
    Y ~ X + Z + X:Z
  "

  fit <- pls(
    syntax,
    data = data,
    ordered = ordered,
    bootstrap = FALSE,
    mc.min.iter = 2L,
    mc.max.iter = 10L,
    mc.reps = 200L,
    mc.tol = 1e-2,
    mc.fixed.seed = TRUE,
    verbose = FALSE
  )

  bench <- c(
    x1 = "r2",
    x2 = "rmse",
    x3 = "mae",
    z1 = "q2",
    z2 = "acc",
    z3 = "ord_mae",
    y1 = "rmse",
    y2 = "rmse",
    y3 = "rmse"
  )

  for (approach in c("earliest", "direct")) {
    pred <- suppressWarnings(pls_predict(fit, approach = approach, benchmark = bench))

    testthat::expect_s3_class(pred, "PlsSemPredict")
    testthat::expect_true(!is.null(pred$X.ord))
    testthat::expect_true(!is.null(pred$X.ord.pred))

    testthat::expect_true(is.data.frame(pred$benchmark))
    testthat::expect_equal(nrow(pred$benchmark), ncol(pred$X.cont.pred))
    testthat::expect_setequal(
      unique(pred$benchmark$metric),
      c("r2", "rmse", "mae", "q2_predict", "acc", "ord_mae")
    )
  }
})
