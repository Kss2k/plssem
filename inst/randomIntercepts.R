library(mvtnorm)
library(lme4)

std1 <- \(x) (x - mean(x)) / sd(x)

model <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'

var_f      <- 1
var_ov     <- 1

beta_f      <- 1.4
zeta_fw     <- 0.4
zeta_fb     <- 0.2
gamma_f_x1  <- 0.3
gamma_f_x2  <- 0.2
gamma_f_x3  <- 0.1
gamma_f_w1  <- 0.15
gamma_f_w2  <- 0.1

cov_x1_x2 <- 0.1
cov_x1_x3 <- 0.0
cov_x2_x3 <- 0.1

cov_w1_w2 <- 0.4
var_beta_f <- 0.7

lambda_1   <- 1
lambda_2   <- .7
lambda_3   <- .8

epsilon    <- 0.2
beta_1     <- 1.2
beta_2     <- 0.8
beta_3     <- 1.5
n          <- 10000
k          <- 400 # number of categories in cluster

sim_data <- function(N = n, K = k) {
  residual <- function(epsilon) rnorm(N, sd = sqrt(epsilon))


  create_ind <- function(lv, beta, lambda, epsilon) {
    beta + lambda * lv + residual(epsilon)
  }

  SXI <- diag(rep(var_ov, 3))
  SXI[1, 2] <- SXI[2, 1] <- cov_x1_x2
  SXI[1, 3] <- SXI[3, 1] <- cov_x1_x3
  SXI[2, 3] <- SXI[3, 2] <- cov_x2_x3

  XI <- rmvnorm(N, sigma = SXI)
  
  x1 <- XI[, 1]
  x2 <- XI[, 2]
  x3 <- XI[, 3]
  
  cluster <- sample(K, N, replace = TRUE)

  x1 <- x1 - ave(x1, cluster)
  x2 <- x2 - ave(x2, cluster)
  x3 <- x3 - ave(x3, cluster)
  w1 <- rnorm(K, mean = 0, sd = sqrt(var_ov))[cluster]
  w2 <- rnorm(K, mean = 0, sd = sqrt(var_ov))[cluster]
 
  fw <- beta_f + gamma_f_x1 * x1 + gamma_f_x2 * x2 + gamma_f_x3 * x3 + residual(zeta_fw)

  fb <- rep(0, N)
  for (i in unique(cluster)) {
    cond <- cluster == i

     
    dbeta <- rnorm(1L, mean = 0, sd = sqrt(var_beta_f))
    fb[cond] <- fb[cond] + dbeta + gamma_f_w1 * w1[cond] + gamma_f_w2 * w2[cond] + rnorm(1L, mean = 0, sqrt(zeta_fb))
  }

  f <- fw + fb

  y1 <- lambda_1 * f + residual(epsilon)
  y2 <- lambda_2 * f + residual(epsilon)
  y3 <- lambda_3 * f + residual(epsilon)
  
  data <- data.frame(
    y1, y2, y3, x1, x2, x3, w1, w2,
    cluster
  )

  data.lv <- data.frame(
    x1, x2, x3, w1, w2, f,
    cluster
  )

  data.lv.std <- as.data.frame(lapply(data.lv, FUN = standardizeAtomic))

  cat("Unstandardized results:\n")
  print(
    lmer('f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)', data = data.lv)
  )
  #> Fixed Effects:
  #> (Intercept)           x1           x2  
  #>     1.47350      0.29505      0.20201  
  #>          x3           w1           w2  
  #>     0.09237      0.15636      0.11234  
  
  cat("Standardized results:\n")
  print(
    lmer('f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)', data = data.lv.std)
  )
  #> Fixed Effects:
  #> (Intercept)           x1           x2  
  #>    -0.01216      0.23596      0.16155  
  #>          x3           w1           w2  
  #>     0.07467      0.12640      0.09235  
  data
}


set.seed(2308257)
randomIntercepts <- sim_data()


save(randomIntercepts, file = "data/randomIntercepts.Rda")
