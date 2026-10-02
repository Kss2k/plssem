library(plssem)
library(modsem)

cleanTmpCols <- function(X) {
  colnames(X) <- stringr::str_remove_all(
    colnames(X), pattern = "\\.tmp$"
  )
  X
}

set.seed(2308257)
level.1 <- '
    x1 <~ 1 * x1.tmp
    x2 <~ 1 * x2.tmp
    x3 <~ 1 * x3.tmp

    fw =~ 0.7 * y1 + 0.6 * y2 + 0.8 * y3
    # fw ~ rv("s1")*x1 + rv("s2")*x2 + x3

    fw ~ 0.5 * x1 + 0.4 * x2 + 0.2 * x3 # fixed effects
    fw ~ 0.2 * s1:x1 + 0.1 * s2:x2      # random effects
'

level.2 <- '
    w1 <~ 1 * w1.tmp
    w2 <~ 1 * w2.tmp
    s1 <~ 1 * s1.tmp
    s2 <~ 1 * s2.tmp

    fb =~ 0.9 * y1 + 0.8 * y2 + 0.9 * y3
    fb ~ 0.2 * w1 + 0.15 * w2

    # the random slopes are latent variables at the between level
    s1 ~ 0.2 * w1 + 0.4 * w2
    s2 ~ 0.1 * w1 + 0.2 * w2
'

par.l1 <- modsemify(level.1)
par.l1$est <- as.numeric(par.l1$mod)

par.l2 <- modsemify(level.2)
par.l2$est <- as.numeric(par.l2$mod)

clusters <- 200
clusterSize <- 15
cluster.idx <- rep(seq_len(clusters), each = clusterSize)

sim.l2 <- plssem:::simulateDataParTable(
  parTable = par.l2,
  N = clusters
)

ov.l2 <- cleanTmpCols(sim.l2$ov)
ov.l2.full <- ov.l2[cluster.idx,,drop=FALSE]

sim.l1 <- plssem:::simulateDataParTable(
  parTable  = par.l1,
  N         = clusters * clusterSize,
  exogenous = data.frame(
    s1 = ov.l2.full[,"s1",drop = TRUE],
    s2 = ov.l2.full[,"s2",drop = TRUE]
  )
)

ov.l1 <- cleanTmpCols(sim.l1$ov)

icc <- list(y1 = 0.333, y2 = 0.200, y3 = 0.500)

randomSlopes <- data.frame(
  y1 = sqrt(1 - icc$y1) * ov.l1$y1 + sqrt(icc$y1) * ov.l2.full$y1,
  y2 = sqrt(1 - icc$y2) * ov.l1$y2 + sqrt(icc$y2) * ov.l2.full$y2,
  y3 = sqrt(1 - icc$y3) * ov.l1$y3 + sqrt(icc$y3) * ov.l2.full$y3,
  x1 = ov.l1$x1,
  x2 = ov.l1$x2,
  x3 = ov.l1$x3,
  w1 = ov.l2.full$w1,
  w2 = ov.l2.full$w2,
  cluster = cluster.idx
)

save(randomSlopes, file = "data/randomSlopes.rda")
