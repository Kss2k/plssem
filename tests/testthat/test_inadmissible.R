devtools::load_all()


tryCatch(
  setwd("test_data"),
  error=\(...) tryCatch(
    setwd("tests/testthat/test_data"),
    error=\(...) setwd(".") # last resort
  )
)


m1 <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
'


d1 <- readRDS("inadmissible_data_1.rds")
f1 <- pls(m1, data = d1, ordered = colnames(d1))
testthat::expect_true(isAdmissible(f1))


d2 <- readRDS("inadmissible_data_2.rds")
f2 <- pls(m1, data = d2, ordered = colnames(d2))
testthat::expect_true(isAdmissible(f2))


d3 <- readRDS("inadmissible_data_3.rds")
f3 <- pls(m1, data = d3, ordered = colnames(d3))
testthat::expect_true(isAdmissible(f3))
