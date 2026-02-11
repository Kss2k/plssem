library(mvtnorm)
library(lme4)
library(plssem)
library(dplyr)
library(ggplot2)
library(plssem)
library(dplyr)

std1 <- \(x) (x - mean(x)) / sd(x)

model <- '
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  
  Y ~ X + Z + (1 + X + Z | cluster)
'

var_x     <- 1
var_z     <- 1

beta_y      <- 1.4
gamma_y_x  <- 0.3
gamma_y_z  <- 0.2

cov_x_z <- 0.2

var_beta_y <- 0.2
var_gamma_y_x <- 0.15
var_gamma_y_z <- 0.05
cov_gamma_x_z <- 0.02

lambda_1   <- 1
lambda_2   <- .7
lambda_3   <- .8

epsilon    <- 0.2
beta_1     <- 1.2
beta_2     <- 0.8
beta_3     <- 1.5
n          <- 10000
k          <- 400 # number of categories in cluster


var_fixed <-
  gamma_y_x^2 * var_x +
  gamma_y_z^2 * var_z +
  2 * gamma_y_x * gamma_y_z * cov_x_z


var_random <-
  var_beta_y +
  var_gamma_y_x * var_x +
  var_gamma_y_z * var_z +
  2 * cov_gamma_x_z * cov_x_z

var_y_proj <- var_fixed + var_random
zeta_y <- 1 - var_y_proj


lhs <- c(rep("Y", 2), "Y~1", "Y~X", "Y~Z", "Y~X")
rhs <- c("X", "Z", "Y~1", "Y~X", "Y~Z", "Y~Z")
op  <- c(rep("~", 2), rep("~~", 4))
parTable.true <- data.frame(
  lhs = lhs,
  op  = op,
  rhs = rhs,
  par = sprintf("%s%s%s", lhs, op, rhs),
  est.true = c(gamma_y_x, gamma_y_z, var_beta_y, var_gamma_y_x, var_gamma_y_z, cov_gamma_x_z)
)


sim_cont_data <- function(N = n, K = k) {
  residual <- function(epsilon) rnorm(N, sd = sqrt(epsilon))


  create_ind <- function(lv, beta, lambda, epsilon) {
    beta + lambda * lv + residual(epsilon)
  }

  SXI <- diag(rep(var_ov, 2))
  SXI[1, 2] <- SXI[2, 1] <- cov_x_z
  XI <- rmvnorm(N, sigma = SXI)
  
  x <- XI[, 1]
  z <- XI[, 2]
  
  cluster <- sample(K, N, replace = TRUE)
  Kc <- length(unique(cluster))

  # Y = XB * Zu + e
  # Y = dependent variable
  # B = fixed effects
  # X = design matrix of B
  # Z = design matrix of u
  # u = random effects
  # e = residual

  # X = [1, x1, x2, x3, w1, w2]
  # Z denotes within and between level categories
  # 1...k columns define the random intercepts additional columns contain
  # values of indpendent variables with random slopes

  X <- matrix(c(rep(1, N), x, z), nrow = N, ncol = 3, byrow = FALSE)
  Z <- matrix(0L, nrow = N, ncol = Kc*3)

  
  for (ki in seq_len(Kc)) {
    Z[cluster==ki, ki] <- 1L
    Z[cluster==ki, Kc + ki] <- x[cluster==ki]
    Z[cluster==ki, 2 * Kc + ki] <- z[cluster==ki]
  }

  Beta <- c(beta_y, gamma_y_x, gamma_y_z)
  G <- diag(c(var_beta_y, var_gamma_y_x, var_gamma_y_z))
  G[2, 3] <- G[3, 2] <- cov_gamma_x_z

  U <- rmvnorm(Kc, mean = c(0, 0, 0), sigma = G)
  u <- as.vector(U)
  e <- rnorm(N, mean = 0, sd = sqrt(zeta_y))

  y <- X %*% Beta + Z %*% u + e

  y1 <- lambda_1 * y + residual(epsilon)
  y2 <- lambda_2 * y + residual(epsilon)
  y3 <- lambda_3 * y + residual(epsilon)

  x1 <- lambda_1 * x + residual(epsilon)
  x2 <- lambda_2 * x + residual(epsilon)
  x3 <- lambda_3 * x + residual(epsilon)

  z1 <- lambda_1 * z + residual(epsilon)
  z2 <- lambda_2 * z + residual(epsilon)
  z3 <- lambda_3 * z + residual(epsilon)
  
  data.ov <- data.frame(
    x1, x2, x3,
    z1, z2, z3,
    y1, y2, y3,
    cluster
  )

  data.ov
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
                       save.fit = FALSE, # save memory...
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
    fit = if (save.fit) fit else NULL,
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
  if (!is.null(x$pars)) print(x$pars[x$pars$op == "~" | grepl("~", x$pars$lhs), ])
  else cat("<FAILED>: NULL")
}


print_sep <-  \() cat(strrep("─", options("width")[[1]]), "\n")


parameter_estimates.cSEMResults <- function(object, ...) {
  paths <- summarize(object)$Estimates$Path_estimates
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
      ordered <- setdiff(colnames(data_cont_i), "cluster")
      data_cat_i  <- cut_data(data = data_cont_i, thr = thresholds,
                              choose = ordered)

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
        
        pls.ord = get_output(
          expr = suppressMessages(pls(model, data = data_cat_i, ordered = ordered, consistent = FALSE)),
          method = "OrdPLS", id = id, cond = cond, ncat = ncat
        ),
        
        pls = get_output(
          expr = suppressMessages(pls(model, data = data_cat_i, consistent = FALSE)),
          method = "PLS", id = id, cond = cond, ncat = ncat
        ),

        plsc.cont = get_output(
          expr = suppressMessages(pls(model, data = data_cont_i)),
          method = "PLSc (cont)", id = id, cond = cond, ncat = ncat
        ),

        pls.cont = get_output(
          expr = suppressMessages(pls(model, data = data_cont_i, consistent = FALSE)),
          method = "PLS (cont)", id = id, cond = cond, ncat = ncat,
        )
      ) 
      
      print(fitted_i)
      
      results_sub[[i]] <- fitted_i
    }

    results <- c(results, results_sub)
  }
}

results <- c(results, results_sub)


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
    filter(method %in% compare & par == param) |>
    group_by(par, method, cond, ncat) |>
    summarize(bias = mean(bias)) |>
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


plot_results(param = "f~x1")
plot_results(param = "f~x2")
plot_results(param = "f~x3")
plot_results(param = "f~w1")
plot_results(param = "f~w2")
plot_results(param = "f~1~~f~1")
table_results()

saveRDS(resd, sprintf("results-random-intercepts-%s.rds", Sys.time()))
