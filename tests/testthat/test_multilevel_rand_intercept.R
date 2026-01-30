devtools::load_all()

lmer.syntax <- "f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)"

model.pls <- '
f =~ y1 + y2 + y3
f ~ x1 + x2 + x3 + w1 + w2
'

fit.pls <- pls(model.pls, data = randomIntercepts, consistent = TRUE, lme4.syntax = lmer.syntax,
               cluster = "cluster")
fit.pls
#> $fixef$f
#>          f~1         f~x1         f~x2         f~x3         f~w1         f~w2 
#> -0.009932809  0.382335745  0.307879834  0.156605698  0.106674833  0.100416354 
