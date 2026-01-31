# `plssem`

The goal of the `plssem` package is to allow the 
estimation of complex Structural Equation Models (SEMs) using the 
PLS-SEM framework. This package expands the PLS-SEM (and PLSc-SEM) framework
to handle categorical data, non-linear models, and multilevel structures.

`plssem` is currently under development. The end goal is to allow the estimation
of non-linear multilevel PLS-SEM (and PLSc-SEM) models with ordinal/categorical data. 

## Installation
The package currently needs to be installed from `GitHub`.

```r
devtools::install_github("kss2k/plssem")
```

## Example
### Linear Model
```r
library(plssem)

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

fit <- pls(tpb, modsem::TPB, consistent = TRUE, 
           bootstrap = TRUE, sample = 500)
summary(fit)
```

### Multilevel Model
```r
pls.syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z
"

lme4.syntax <- "Y ~ X + Z + (1 + X + Z | cluster)"

fit <- pls(pls.syntax, data = randomSlopes, lme4.syntax = lme4.syntax,
           cluster = "cluster", consistent = TRUE)
```

## TODO

1. Add handling of ordinal data
2. Add handling of multilevel models
3. Add handling of missing data (multiple imputation)
4. Add handling of ordinal data in multilevel models
5. Add handling of interaction effects
6. Add handling of ordinal data with interaction effects
7. Add `parameter_estimates()` function.
8. Add `coef()` function.
9. Add `vcov()` function.
10. Finish `summary()` function.
    - Add fit measures
    - Add `R^2` and `R^2-adj`
