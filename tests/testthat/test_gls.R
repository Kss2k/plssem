devtools::load_all()
library(modsem)

m <- '
  int1 ~ att1 + sn1 + pbc1
  b1 ~ int1 + pbc1
'

mod <- createGlsModel(
  modsemify(m),
  data.cov = NULL #cov(TPB)
)

glsModelCovMatrix(mod) <- cov(TPB)
fn <- \(x) glsObjective(glsFillModel(mod, x))
gr <- \(x) glsGradient(glsFillModel(mod, x))

glsEstimateModel(mod, cov(TPB))
opt <- nlminb(
  start = mod@info$start,
  objective = fn,
  gradient = gr
)

glsFillModel(mod, opt$par)
