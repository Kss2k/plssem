# `plssem` <img src="man/figures/plssem.png" alt="Logo" align = "right" height="139" class="logo">
<!-- badges: start -->
[![Tests](https://github.com/kss2k/plssem/actions/workflows/tests.yml/badge.svg)](https://github.com/kss2k/plssem/actions/workflows/tests.yml)
[![R-CMD-check](https://github.com/kss2k/plssem/actions/workflows/checks.yml/badge.svg)](https://github.com/kss2k/plssem/actions/workflows/checks.yml)
[![PKGDOWN-Build](https://github.com/kss2k/plssem/actions/workflows/pkgdown-build.yml/badge.svg)](https://github.com/kss2k/plssem/actions/workflows/pkgdown-build.yml)
[![CRAN](https://www.r-pkg.org/badges/version/plssem)](https://cran.r-project.org/package=plssem)
<!-- [![](https://cranlogs.r-pkg.org/badges/grand-total/plssem)](https://cran.r-project.org/package=plssem) -->
<!-- badges: end -->

The goal of the [`plssem`](https://kss2k.github.io/plssem/) package is to allow the
estimation of complex Structural Equation Models (SEMs) using the
PLS-SEM framework. This package expands the PLS-SEM (and PLSc-SEM) framework
to handle categorical data, non-linear models, and multilevel structures, using
[Monte-Carlo Consistent Partial Least Squares Structural Equation Modelling (MC-PLSc-SEM)](https://osf.io/preprints/psyarxiv/fwzj6_v1)

[`plssem`](https://kss2k.github.io/plssem/) is currently under development. The end goal is to allow the consistent estimation
of non-linear multilevel SEMs with ordinal and categorical data, using the MC-PLSc-SEM framework.

## Installation

The package can be downloaded from `CRAN`.

```r
install.packages("plssem")
```

The development version of the package can be installed from `GitHub`.

```r
devtools::install_github("kss2k/plssem")
```

## Examples
### Linear Model with Continuous Data
```r
library(plssem)
library(modsem)

tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
'

fit <- pls(tpb, TPB, bootstrap = TRUE)
summary(fit)
```

### Linear Model with Ordered Data
```r
tpb <- '
# Outer Model (Based on Hagger et al., 2007)
  ATT =~ att1 + att2 + att3 + att4 + att5
  SN =~ sn1 + sn2
  PBC =~ pbc1 + pbc2 + pbc3
  INT =~ int1 + int2 + int3
  BEH =~ b1 + b2

# Inner Model (Based on Steinmetz et al., 2011)
  INT ~ ATT + SN + PBC
  BEH ~ INT + PBC
'

fit <- pls(tpb, TPB_Ordered, bootstrap = TRUE)
summary(fit)
```

### Multilevel Random Slopes Model with Continuous Data
```r
syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  W =~ w1 + w2 + w3
  Y ~ X + Z + (1 + X + Z | cluster)
  W ~ X + Z + (1 + X + Z | cluster)
"

fit <- pls(syntax, data = randomSlopes, bootstrap = TRUE)
summary(fit)
```

### Multilevel Random Slopes Model with Ordered Data
```r
syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  W =~ w1 + w2 + w3
  Y ~ X + Z + (1 + X + Z | cluster)
  W ~ X + Z + (1 + X + Z | cluster)
"

fit <- pls(syntax, data = randomSlopesOrdered, bootstrap = TRUE)
summary(fit)
```

### Multilevel Random Intercepts Model with Continuous Data
```r
syntax <- '
  f =~ y1 + y2 + y3
  f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
'

fit <- pls(syntax, data = randomIntercepts, bootstrap = TRUE)
summary(fit)
```

### Multilevel Random Intercepts Model with Ordered Data
```r
syntax <- '
  f =~ y1 + y2 + y3
  f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
'

fit <- pls(syntax, data = randomInterceptsOrdered, bootstrap = TRUE)
summary(fit)
```

### Interaction Model with Continuous Data
```r
m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

fit <- pls(m, modsem::oneInt, bootstrap = TRUE)
summary(fit)
```

### Interaction Model with Ordered Data
```r
m <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3

  Y ~ X + Z + X:Z
'

fit <- pls(m, oneIntOrdered, bootstrap = TRUE)
summary(fit)
```

### Multilevel Random Slopes Interaction Model with Continuous Data
```r
syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  W =~ w1 + w2 + w3
  Y ~ X + Z + X:Z + (1 + X + Z + X:Z | cluster)
  W ~ X + Z + X:Z + (1 + X + Z + X:Z | cluster)
"

fit <- pls(m, randomSlopes, bootstrap = TRUE)
summary(fit)
```
