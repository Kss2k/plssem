devtools::load_all()

m <- "
  X <~ x1 + x2 + x3
  Z <~ z1 + z2 + z3
"

pt <- modsem::modsemify(m)
