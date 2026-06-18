devtools::load_all()
library(modsem)
library(lavaan)

m <- '
  int1 ~ att1 + sn1 + pbc1
  b1 ~ att1 + sn1 + pbc1 + pbc1:sn1
  int1 ~~ b1
'

mod <- GlsPathModel(
  modsemify(m),
  data.cov = NULL #cov(TPB)
)


# Compare with gls estimator in lavaan
addRevCov <- function(pt) {
  cov <- pt[pt$op == "~~" & pt$lhs != pt$rhs, ,drop=FALSE]
  rev <- cov
  rev$lhs <- cov$rhs
  rev$rhs <- cov$lhs
  rbind(pt, rev)
}

X <- TPB
X$`pbc1:sn1` <- X$pbc1 * X$sn1

ppls <- addRevCov(glsEstimateParameters(mod, cov(X))@parTable)
plav <- addRevCov(lavaan::parameterEstimates(lavaan::sem(m, TPB)))

testthat::expect_equal(nrow(ppls), nrow(plav))
pall <- merge(ppls, plav, by = c("lhs", "op", "rhs"))
testthat::expect_equal(pall$est.x, pall$est.y, tol = 6e-4)

tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ ATT + SN + PBC + PBC:SN

  BEH ~~ INT
'

testthat::expect_no_error({
  fit <- pls(tpb, TPB)
  summary(fit)

  fit <- pls(tpb, TPB, mcpls=TRUE, bootstrap=TRUE)
  summary(fit)
})
