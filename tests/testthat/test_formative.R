devtools::load_all()

testthat::expect_no_error({
  m <- "
    X <~ x1 + x2 + x3
    Z <~ z1 + z2 + z3
  "

  pt <- modsem::modsemify(m)
})
