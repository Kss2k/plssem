devtools::load_all()

m <- '
  X =~ load * x1 + load * x2 + load * x3
  Z =~ load * z1 + load * z2 + load * z3
  Y =~ load * y1 + load * y2 + load * y3

  Y ~ "dnorm(.4, .1)" * X +
     "dnorm(.35, .1)" * Z +
     "dnorm(.45, .1)" * X:Z +
     "dnorm(0, .005)" * X:X

  load :~ dnorm(.8, .5)
'

if (FALSE) { # don't run on GitHub
  set.seed(23942)
  fit <- bpls(
    m, modsem::oneInt, boot.R = 500, warmup = 200, iter = 400, sampler = "Metropolis-Hastings",
    parallel = "multisession", chains = 2
  )

  round(apply(fit$samples, MARGIN = 2, FUN = mean), 3)
  round(apply(fit$samples, MARGIN = 2, FUN = sd), 3)


  fit <- bpls(m, oneIntOrdered, ordered = colnames(oneIntOrdered),
                  boot.R = 5000, warmup = 5000, iter = 10000, sampler = "Metropolis-Hastings",
                  parallel = "multicore", chains = 3)
  round(apply(fit$samples, MARGIN = 2, FUN = mean), 3)
  round(apply(fit$samples, MARGIN = 2, FUN = sd), 3)

  
}
