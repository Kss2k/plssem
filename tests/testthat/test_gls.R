devtools::load_all()
library(modsem)

m <- '
  int1 ~ att1 + sn1 + pbc1
  b1 ~ int1 + pbc1
'

mod <- createGlsModel(
  modsemify(m),
  data.cov = cov(TPB)
)

glsCalcSigmaHat(glsFillModel(mod, mod$info$start))
fn <- \(x) glsObjective(glsFillModel(mod, x))

nlminb(
  start = mod$info$start,
  objective = fn,
  control = list(iter.max = 4)
)
