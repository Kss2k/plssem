devtools::load_all()

model <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'

fit.pls <- pls(model, data = randomSlopes, cluster = "cluster", bootstrap = TRUE, approach.weights = "pca",
               boot.R = 500, mc.delta.jacobian.k = 5, boot.ncores = 4, level2.cov = "muml", small.sample = FALSE, mc.reps = 20000)
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
fit.pls <- pls(modelr, data = randomSlopesOrdered,
               ordered = colnames(randomSlopesOrdered),
               cluster = "cluster", boot.R = 5000, boot.ncores = 5,
               mc.delta.jacobian.k = 5, boot.parallel = "multisession", level2.cov = "muml",
               bootstrap = TRUE, approach.weights = "pca")
fit.pls
  

fit.lav <- lavaan::sem(model = modelr, data = randomSlopes, cluster = "cluster")
lavaan::summary(fit.lav, standardized = TRUE)
lavaan::lavInspect(fit.lav, "icc")
