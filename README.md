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
