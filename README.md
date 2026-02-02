# `plssem`

The goal of the `plssem` package is to allow the 
estimation of complex Structural Equation Models (SEMs) using the 
PLS-SEM framework. This package expands the PLS-SEM (and PLSc-SEM) framework
to handle categorical data, non-linear models, and multilevel structures.

`plssem` is currently under development. The end goal is to allow the estimation
of non-linear multilevel PLS-SEM (and PLSc-SEM) models with ordinal/categorical data. 

## Installation
The package currently needs to be installed from `GitHub`. Currently, it
also depends on the latest development version of `modsem`.

```r
devtools::install_github("kss2k/modsem")
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

## TODO
1. Fix mismatching thresholds in bootstrapping (`R/bootstrap.R`, line 30)
