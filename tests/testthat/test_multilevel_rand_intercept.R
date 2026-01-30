devtools::load_all()
library(lavaan)
library(lme4)

lmer.syntax <- "f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)"

model.cfa <- 'f =~ y1 + y2 + y3'

model <- '
    level: 1
        f ~ x1 + x2 + x3
    level: 2
        f ~ w1 + w2
'

fit.cfa <- cfa(model.cfa, Demo.twolevel)
Z <- as.data.frame(round(lavPredict(fit.cfa, append = TRUE, method = "bartlett"), 5))
X <- round(as.matrix(Demo.twolevel), 5)
Y <- as.data.frame(apply(merge(X, Z), MARGIN = 2L, FUN = standardizeAtomic))

fit.lmer <- lmer(lmer.syntax, data = Y)
#> Fixed Effects:
#> (Intercept)           x1           x2           x3           w1           w2  
#>   -0.008186     0.329428     0.267710     0.137230     0.095767     0.082591  
fit.sem <- sem(model, Y, cluster = "cluster")

summary(fit.lmer)
summary(fit.sem, standardized = TRUE)


model.pls <- '
f =~ y1 + y2 + y3
f ~ x1 + x2 + x3 + w1 + w2
'

fit.pls <- pls(model.pls, data = Demo.twolevel, consistent = TRUE, lme4.syntax = lmer.syntax,
               cluster = "cluster")
fit.pls
#> $fixef$f
#>          f~1         f~x1         f~x2         f~x3         f~w1         f~w2 
#> -0.009932809  0.382335745  0.307879834  0.156605698  0.106674833  0.100416354 
