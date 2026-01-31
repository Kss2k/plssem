devtools::load_all()

lmer.syntax <- "f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)"

model.pls <- '
f =~ y1 + y2 + y3
f ~ x1 + x2 + x3 + w1 + w2
'

fit.pls <- pls(model.pls, data = randomIntercepts, consistent = TRUE, lme4.syntax = lmer.syntax,
               cluster = "cluster", bootstrap = TRUE)
summary(fit.pls)
