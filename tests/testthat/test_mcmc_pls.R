devtools::load_all()

m <- '
  X =~ l * x1 + l * x2 + l * x3
  Z =~ l * z1 + l * z2 + l * z3
  Y =~ l * y1 + l * y2 + l * y3

  Y ~ "dnorm(.4, .1)" * X + "dnorm(.35, .1)" * Z + "dnorm(.45, .1)" * X:Z + "dnorm(0, .005)" * X:X
  l :~ dnorm(mean = 0.8, sd = 0.3)
'

if (FALSE) { # don't run on GitHub
  set.seed(23942)
  fit <- mcmc_pls(
    m, modsem::oneInt, boot.R = 500, warmup = 2000, iter = 4000, sampler = "Metropolis-Hastings",
    parallel = "multisession", chains = 2
  )

  round(apply(fit$samples, MARGIN = 2, FUN = mean), 3)
  round(apply(fit$samples, MARGIN = 2, FUN = sd), 3)


  fit <- mcmc_pls(m, oneIntOrdered, ordered = colnames(oneIntOrdered),
                  boot.R = 500, warmup = 2000, iter = 4000, sampler = "Metropolis-Hastings")
  round(apply(fit$samples, MARGIN = 2, FUN = mean), 3)
  round(apply(fit$samples, MARGIN = 2, FUN = sd), 3)
}
