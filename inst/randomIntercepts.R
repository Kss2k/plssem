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

beta_f      <- 0
gamma_f_x1  <- 0.3
gamma_f_x2  <- 0.2
gamma_f_x3  <- 0.1
gamma_f_w1  <- 0.15
gamma_f_w2  <- 0.1

cov_x1_x2 <- 0.2
cov_x1_x3 <- 0.3
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
n          <- 2500
k          <- 5 # number of categories in cluster

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

  w1 <- rnorm(K, mean = 0, sd = sqrt(var_ov))[cluster]
  w2 <- rnorm(K, mean = 0, sd = sqrt(var_ov))[cluster]
 
  browser()
  f <- gamma_f_x1 * x1 + gamma_f_x2 * x2 + gamma_f_x3 * x3 +
    gamma_f_w1 * w1 + gamma_f_w2 * w2

  for (i in unique(cluster)) {
    cond <- cluster == i
    dbeta <- beta_f + rnorm(1L, mean = 0, sd = sqrt(var_beta_f))
    
    f[cond] <- f[cond] + dbeta
  }

  zeta_f <- 1 - var(f) # empirical solution
  f <- f + residual(zeta_f)

  
  data <- data.frame(
    x1, x2, x3, w1, w2,
    cluster
  )

  data.lv <- data.frame(
    x1, x2, x3, w1, w2, f,
    cluster
  )

  cat("Unstandardized results:\n")
  print(
    lmer('f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)', data = data.lv)
  )
  # 
  # cat("Standardized results:\n")
  # print(
  #   lmer('Y ~ X + Z + (1 + X + Z | cluster)',
  #        data = data.frame(X=std1(X), Z=std1(Z), Y=std1(Y), cluster))
  # )
  data
}


set.seed(2308257)
randomIntercepts <- sim_data()

save(randomSlopes, file = "data/randomSlopes.Rda")
