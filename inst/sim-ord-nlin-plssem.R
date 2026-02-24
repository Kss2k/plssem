set.seed(2308257)

library(mvtnorm)
library(modsem)
library(tidyr)
library(ggplot2)
library(plssem)
library(cSEM)
library(dplyr)

# ──────────────────────────────────────────────────────────────────────────────
# Setup Simulation
# ──────────────────────────────────────────────────────────────────────────────


model <- '
X =~ x1 + x2 + x3
Z =~ z1 + z2 + z3
Y =~ y1 + y2 + y3

Y ~ X + Z + X:Z #+ X:X + Z:Z
'

corr_X_Z  <- 0.2 


gamma_Y_X  <- 0.4
gamma_Y_Z  <- 0.5
gamma_Y_XZ <- 0.3
gamma_Y_ZZ <- 0 # exclude for now
gamma_Y_XX <- 0 # exclude for now

# cov(xy,uv) = E(x)E(u)cov(y,v) + 
#              E(x)E(v)cov(y,u) +
#              E(y)E(u)cov(x,v) +
#              E(y)E(v)cov(x,u) + 
#              cov(x,u)cov(y,v) +
#              cov(x,v)cov(y,u)
#
# E(x) = E(y) = E(u) = E(v) = 0
# 
# cov(xy,uv) = cov(x,u)cov(y,v) +
#              cov(x,v)cov(y,u)

cov_XZ_XX <- 2*corr_X_Z
var_XX <- 2
var_ZZ <- 2
cov_XZ_ZZ <- 2*corr_X_Z
cov_XX_ZZ <- 2*corr_X_Z^2
var_XZ    <- 1 + corr_X_Z^2

mean_XZ <- corr_X_Z
mean_XX <- 1
mean_ZZ <- 1

lambda_1   <- 0.90
lambda_2   <- 0.80
lambda_3   <- 0.85

n          <- 1000

alpha_Y <- - (gamma_Y_XZ * mean_XZ + gamma_Y_XX * mean_XX + gamma_Y_ZZ * mean_ZZ) # we want E[Y] = 0
zeta_Y  <- 1 - (
  gamma_Y_X^2 + gamma_Y_Z^2 + gamma_Y_XZ^2*var_XZ + gamma_Y_XX^2*var_XX + gamma_Y_ZZ^2*var_ZZ +
  2 * gamma_Y_X*gamma_Y_Z*corr_X_Z +
  2 * gamma_Y_XX*gamma_Y_XZ*cov_XZ_XX +
  2 * gamma_Y_ZZ*gamma_Y_XZ*cov_XZ_ZZ +
  2 * gamma_Y_XX*gamma_Y_ZZ*cov_XX_ZZ
)

lhs <- "Y"
rhs <- c("X", "Z", "X:Z", "X:X", "Z:Z")
parTable.true <- data.frame(
  lhs = lhs,
  op = "~",
  rhs = rhs,
  par = sprintf("%s~%s", lhs, rhs),
  est.true = c(gamma_Y_X, gamma_Y_Z, gamma_Y_XZ, gamma_Y_XX, gamma_Y_ZZ)
)


sim_cont_data <- function(N = n[1]) {
  residual <- function(epsilon) rnorm(N, sd = sqrt(epsilon))

  create_ind <- function(lv, lambda) {
    epsilon <- 1 - lambda^2
    lambda * lv + residual(epsilon)
  }

  SXI <- matrix(c(1, corr_X_Z,
                  corr_X_Z, 1), nrow = 2)
  XI <- rmvnorm(N, sigma = SXI)

  X <- XI[, 1]
  Z <- XI[, 2]

  Y <- 
    gamma_Y_X * X + 
    gamma_Y_Z * Z + 
    gamma_Y_XZ * X * Z +
    gamma_Y_XX * X * X +
    gamma_Y_ZZ * Z * Z +
    residual(zeta_Y)
 
  x1 <- create_ind(X, lambda_1)
  x2 <- create_ind(X, lambda_2)
  x3 <- create_ind(X, lambda_2)
  
  z1 <- create_ind(Z, lambda_1)
  z2 <- create_ind(Z, lambda_2)
  z3 <- create_ind(Z, lambda_2)
  
  y1 <- create_ind(Y, lambda_1)
  y2 <- create_ind(Y, lambda_2)
  y3 <- create_ind(Y, lambda_2)
  
  data <- data.frame(
    x1, x2, x3,
    z1, z2, z3,
    y1, y2, y3
  )

  # return(var(Y))
  data
}


rthreshold <- function(k, offset = runif(1, min = -0.7, max = 0.7), sigma = 0.4) {
  t <- seq_len(k) - mean(seq_len(k)) + offset
  t + runif(k, min = -sigma, max = sigma)
}


# Based on Rhemtulla et al., 2012 and Schubert et al., 2018
list_thresholds <- list(
  Uneven = list(
    `2` = \() rthreshold(1),
    `3` = \() rthreshold(2),
    `4` = \() rthreshold(3),
    `5` = \() rthreshold(4),
    `6` = \() rthreshold(5),
    `7` = \() rthreshold(6)
  ),
  Symmetric = list(
  `2` = c( 0.00),
  `3` = c(-0.83,  0.83),
  `4` = c(-1.25,  0.00,  1.25),
  `5` = c(-1.50, -0.50,  0.50, 1.50),
  `6` = c(-1.60, -0.83,  0.00, 0.83, 1.60),
  `7` = c(-1.79, -1.07, -0.36, 0.36, 1.07, 1.79)
  ),
  Moderate = list(
    `2` = c( 0.36),
    `3` = c(-0.50,  0.76),
    `4` = c(-0.31,  0.79,  1.66),
    `5` = c(-0.70,  0.39,  1.16,  2.05),
    `6` = c(-1.05,  0.08,  0.81,  1.44,  2.33),
    `7` = c(-1.43, -0.43,  0.38,  0.94,  1.44,  2.54)
  ),
  Extreme = list(
    `2` = c( 1.04),
    `3` = c( 0.58,  1.13),
    `4` = c( 0.28,  0.71,  1.23),
    `5` = c( 0.05,  0.44,  0.84,  1.34),
    `6` = c(-0.13,  0.25,  0.61,  0.99,  1.48),
    `7` = c(-0.25,  0.13,  0.47,  0.81,  1.18,  1.64)
  ),
  Alt.Mod = list(
    `2` = c(-0.36),
    `3` = c(-0.76,  0.50),
    `4` = c(-1.66, -0.79,  0.31),
    `5` = c(-2.05, -1.16, -0.39,  0.70),
    `6` = c(-2.33, -1.44, -0.81, -0.08,  1.05),
    `7` = c(-2.54, -1.44, -0.94, -0.38,  0.43,  1.43)
  ),
  Ald.Ext = list(
    `2` = c(-1.04),
    `3` = c(-1.13, -0.58),
    `4` = c(-1.23, -0.71, -0.28),
    `5` = c(-1.34, -0.84, -0.44, -0.05),
    `6` = c(-1.48, -0.99, -0.61, -0.25,  0.13),
    `7` = c(-1.64, -1.18, -0.81, -0.47, -0.13,  0.25)
  )
)


cut_data <- function(data, thr, choose = NULL) {
  standardize <- \(x) (x - mean(x)) / sd(x)

  if (is.null(choose))
    choose <- colnames(data)

  means <- data

  if (is.function(thr))
    thr <- thr()

  for (i in seq_along(choose)) {
    var <- choose[[i]]
    x <- data[[var]]
    breaks <- c(-Inf, thr, Inf)
    y <- cut(standardize(x), breaks = breaks)
    y <- as.integer(as.ordered(as.integer(y)))

    data[[var]] <- y
    z <- rep(NA_real_, length(y))

    for (v in unique(y))
        z[y == v] <- mean(x[y == v])

    means[[var]] <- z
  }

  attr(data, "means") <- as.data.frame(means)
  data
}


R <- 200L


get_output <- function(expr,
                       parfun = parameter_estimates,
                       method = NA,
                       id = NA,
                       cond = "",
                       ncat = "",
                       ...) {
  start <- Sys.time()
  fit <- tryCatch(eval(expr), 
                  error = \(e) {
                    warning(sprintf("%s (%d) failed!, message:\n %s",
                                    method, id, e))
                    NULL
                  })
  end <- Sys.time()
  elapsed <- end - start
 
  out <- list(
    fit = fit,
    elapsed = elapsed,
    method = method,
    id = id,
    cond = cond,
    ncat = as.integer(ncat),
    pars = if (!is.null(fit)) parfun(fit, ...) else NULL
  )
  
  class(out) <- "simoutput"
  out
}


print.simoutput <- function(x, ...) {
  cat(sprintf("ID: %i, Method: %s, Elapsed: %s Cond: %s, NCAT: %d\n", 
              x$id, x$method, capture.output(x$elapsed), x$cond, x$ncat))
  if (!is.null(x$pars)) print(x$pars[x$pars$op == "~", ])
  else cat("<FAILED>: NULL")
}


print_sep <-  \() cat(strrep("─", options("width")[[1]]), "\n")


parameter_estimates.cSEMResults <- function(object, ...) {
  paths <- cSEM::summarize(object)$Estimates$Path_estimates
  par <- stringr::str_remove_all(paths$Name, " ")
  par <- stringr::str_replace_all(par, "\\.", ":")
  lhs <- stringr::str_split_i(par, "~", i = 1L)
  rhs <- stringr::str_split_i(par, "~", i = 2L)
  op  <- "~"
  est <- paths$Estimate
  se  <- paths$Std_err
  pvalue <- paths$p_value
  z      <- est / se
  ci.lower <- est - 1.96 * se
  ci.upper <- est + 1.96 * se

  plssem:::plssemParTable(data.frame(
    lhs = lhs, op = op, rhs = rhs, est = est,
    se = se, z = z, pvalue = pvalue,
    ci.lower = ci.lower, ci.upper = ci.upper, label = ""
  ))
}

# ──────────────────────────────────────────────────────────────────────────────
# Run Simulation
# ──────────────────────────────────────────────────────────────────────────────
id <- 0
total <- R * sum(sapply(list_thresholds, length))
results <- NULL

for (cond in names(list_thresholds)) {
  categories <- names(list_thresholds[[cond]])

  for (ncat in categories) {
    thresholds <- list_thresholds[[cond]][[ncat]]

    results_sub <- vector("list", length=R)
    for (i in seq_len(R)) {
      id <- id + 1

      if (!is.null(results_sub[[i]])) {
        message(sprintf("Skipping iteration %i, as it has already been run...", i))
        next
      }

      data_cont_i <- sim_cont_data()
      data_cat_i  <- cut_data(data = data_cont_i, thr = thresholds)
      ordered <- colnames(data_cat_i)

      print_sep()
      cat(sprintf("Iteration %d/%d:\n", id, total)) 
      print_sep()
      
      fitted_i <- list(
        plsc.ord = get_output(
          expr = suppressMessages(pls(model, data = data_cat_i, ordered = ordered,
                                      consistent.probit = TRUE)),
          method = "OrdPLSc", id = id, cond = cond, ncat = ncat
        ),

        plsc = get_output(
          expr = suppressMessages(pls(model, data = data_cat_i)),
          method = "PLSc", id = id, cond = cond, ncat = ncat
        ),
        
        pls.ord.c = get_output(
          expr = suppressMessages(pls(model, data = data_cat_i, ordered = ordered, consistent = FALSE)),
          method = "OrdPLS", id = id, cond = cond, ncat = ncat
        ),
        
        pls = get_output(
          expr = suppressMessages(pls(model, data = data_cat_i, consistent = FALSE)),
          method = "PLS", id = id, cond = cond, ncat = ncat
        ),

        pls.ord.csem = get_output(
          expr = csem(.model = stringr::str_replace_all(model, ":", "."),
                      .data = as.data.frame(lapply(data_cat_i, as.ordered))),
          method = "cSEM (OrdPLS+2SMM)", id = id, cond = cond, ncat = ncat
        ),

        pls.ord.mcem= get_output(
          expr = mcem_nlin_ord_pls(model, data_cat_i, ordered = ordered, mc.reps = 1e4),
          method = "OrdPLSc (MCEM)", id = i, cond = cond, ncat = ncat,
          parfun = \(x) x
        )

        #,

        # lms.cont = get_output( # to have an unbiased reference point
        #   expr = modsem(model, data = data_cont_i, method = "lms"),
        #   method = "LMS-cont", id = i, cond = cond, ncat = ncat,
        #   parfun = standardized_estimates
        # )
      ) 
      
      print(fitted_i)
      
      results_sub[[i]] <- fitted_i
    }

    results <- c(results, results_sub)
  }
}


# ──────────────────────────────────────────────────────────────────────────────
# Parameter Estimates
# ──────────────────────────────────────────────────────────────────────────────
resd <- NULL
cols <- c("id", "method", "lhs", "op", "rhs", 
          "par", "est", "se", "pvalue", "elapsed",
          "ncat", "cond")


failed.pars <- parTable.true
failed.pars$est <- NA

failed.pars <- data.frame(
  lhs = lhs,
  op = "~",
  rhs = rhs,
  par = sprintf("%s~%s", lhs, rhs),
  est = NA,
  se = NA,
  pvalue = NA,
  method = NA,
  id = NA,
  elapsed = NA,
  ncat = NA,
  cond = NA
)


results_dfs <- vector("list", length(results))
for (id in seq_along(results)) {
  results_id <- results[[id]]
  cat(sprintf("\r%d/%d", id, length(results)))

  resid <- NULL
  for (fitted in results_id) {
    if (is.null(fitted$pars)) resi <- failed.pars  
    else                      resi <- fitted$pars
    resi$id <- fitted$id

    colnames(resi) <- stringr::str_replace_all(
      string = colnames(resi),
      pattern = c(std.error = "se", p.value = "pvalue")
    )
    
    resi$method <- fitted$method
    resi$ncat   <- fitted$ncat
    resi$cond   <- fitted$cond
    resi$par <- paste0(resi$lhs,  resi$op, resi$rhs)
    resi$elapsed <- fitted$elapsed
  
    resid <- rbind(resid, resi[cols])
  }

  results_dfs[[id]] <- resid
}
    
resd <- do.call("rbind", results_dfs)
cat("\n")

resd$par <- paste0(resd$lhs,  resd$op, resd$rhs)
resd <- merge(resd, parTable.true, all.x = TRUE)
resd$bias <- resd$est - resd$est.true



# Mehtods to look at
methods <- unique(resd$method)
plot_results <- function(compare = methods, param = "Y~X:Z") {
  resd |>
    filter(op == "~" & method %in% compare & par == param) |>
    group_by(par, method, cond, ncat) |>
    dplyr::summarize(bias = mean(bias, na.rm = TRUE)) |>
    ggplot(aes(x = method, y = bias, colour = method, fill = method)) +
    geom_col(alpha = 0.2) +
    facet_grid(rows = vars(cond), cols = vars(ncat), scales = "fixed") +
    ggtitle(sprintf("Bias for %s by method", param)) +
    ylab("Bias") +
    xlab("Method") + 
    theme_bw()
}


table_results <- function(compare = methods, alpha = 0.05) {
  resd |>
    filter(op == "~" & method %in% compare) |> 
    group_by(par, method) |>
    summarize(est.mean = mean(est, na.rm = TRUE),
              se.obs = sd(est, na.rm = TRUE),
              se.exp = mean(se, na.rm = TRUE),
              percent.sig = sum(pvalue < alpha, na.rm = TRUE) / 
                            sum(!is.na(pvalue))) |>
    print(n=Inf)
}


plot_results(param = "Y~X:Z")
plot_results(param = "Y~Z")
plot_results(param = "Y~X")
table_results()

saveRDS(resd, sprintf("results-sim-ord-nlin-csem-%s.rds", Sys.time()))
