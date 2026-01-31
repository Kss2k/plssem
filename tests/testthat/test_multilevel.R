devtools::load_all()

pls.syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  W =~ w1 + w2 + w3
  Y ~ X + Z
  W ~ X + Z
"

lme4.syntax <- "
Y ~ X + Z + (1 + X + Z | cluster)
W ~ X + Z + (1 + X + Z | cluster)
"

fit <- pls(pls.syntax, data = randomSlopes, lme4.syntax = lme4.syntax,
           cluster = "cluster", consistent = TRUE, bootstrap = TRUE)
fit
