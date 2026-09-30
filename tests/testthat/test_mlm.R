devtools::load_all()

model <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'

fit.pls <- mpls(model, data = randomSlopes, cluster = "cluster")
fit.pls

fit.lav <- lavaan::sem(model, data = randomSlopes, cluster = "cluster")
lavaan::summary(fit.lav, standardized = TRUE)
lavaan::lavInspect(fit.lav, "icc")

modelr <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ rv("s1")*x1 + rv("s2")*x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
        # the random slopes are latent variables at the between level
        s1 + s2 ~ w1 + w2
'

set.seed(23984)
fit.pls <- mpls(modelr, data = randomSlopes, cluster = "cluster")
fit.pls

fit.lav <- lavaan::sem(model = modelr, data = randomSlopes, cluster = "cluster")
lavaan::summary(fit.lav, standardized = TRUE)
lavaan::lavInspect(fit.lav, "icc")
