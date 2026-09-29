devtools::load_all()

model <- '
    level: 1
        X1 <~ x1
        X2 <~ x2
        X3 <~ x3
        fw =~ y1 + y2 + y3
        fw ~ X1 + X2 + X3
    level: 2
        W1 <~ w1
        W2 <~ w2
        fb =~ y1 + y2 + y3
        fb ~ W1 + W2
'

fit <- mpls(model, data = lavaan::Demo.twolevel, cluster = "cluster")

fit.lav <- lavaan::sem(model, data = lavaan::Demo.twolevel, cluster = "cluster")
lavaan::summary(fit.lav, standardized = TRUE)
lavaan::lavInspect(fit.lav, "icc")
