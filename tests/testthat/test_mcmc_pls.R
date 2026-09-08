devtools::load_all()

m <- '
  X =~ l * x1 + l * x2 + l * x3
  Z =~ l * z1 + l * z2 + l * z3
  Y =~ l * y1 + l * y2 + l * y3

  Y ~ b1 * X + b2 * Z + b3 * X:Z + b4 * X:X
  # l  :~ dnorm(mean = 0.8, sd = 0.1)
  # b1 :~ dnorm(mean = 0.4, sd = 0.05)
  # b2 :~ dnorm(mean = 0.4, sd = 0.05)
  # b3 :~ dnorm(mean = 0.45, sd = 0.05)
  # b4 :~ dnorm(mean = 0, sd = 0.005)
'

if (FALSE) { # don't run on GitHub
  set.seed(23942)
  fit <- mcmc_pls(m, modsem::oneInt, boot.R = 500, warmup = 2000, iter = 4000, sampler = "Metropolis-Hastings",
   parallel = "multisession", chains = 2)
  round(apply(fit, MARGIN = 2, FUN = mean), 3)
  round(apply(fit, MARGIN = 2, FUN = sd), 3)


  fit <- mcmc_pls(m, oneIntOrdered, ordered = colnames(oneIntOrdered),
                  boot.R = 500, warmup = 2000, iter = 4000, sampler = "Metropolis-Hastings")
  round(apply(fit[[1]], MARGIN = 2, FUN = mean), 3)
  round(apply(fit[[1]], MARGIN = 2, FUN = sd), 3)
}
