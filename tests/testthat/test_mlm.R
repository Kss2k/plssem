devtools::load_all()

model <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ X1 + X2 + X3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ W1 + W2
'

fit <- mpls(model, data = lavaan::Demo.twolevel, cluster = "cluster")
