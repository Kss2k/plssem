set.seed(2308257)

library(mvtnorm)
library(modsem)
library(tidyr)
library(ggplot2)
library(plssem)
library(cSEM)
library(dplyr)

setwd("C:/Users/slupp/Desktop/CODE/R/repos/plssem")

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




get_output <- function(expr,
                       parfun = parameter_estimates,
                       method = NA,
                       id = NA,
                       cond = "",
                       ncat = "",
                       seed = NULL,
                       write = TRUE,
                       ...) {
  start <- Sys.time()

  if (!is.null(seed))
    set.seed(seed)

  fit <- tryCatch(eval(expr),
                  error = \(e) {
                    warning(sprintf("%s (%d) failed!, message:\n %s",
                                    method, id, e))
                    NULL
                  })
  end <- Sys.time()
  elapsed <- end - start

  out <- list(
    fit = NULL, #fit,
    elapsed = elapsed,
    method = method,
    id = id,
    cond = cond,
    ncat = as.integer(ncat),
    pars = if (!is.null(fit)) parfun(fit, ...) else NULL,
    seed = seed
  )

  if (write)
    writeOutputToFile(out)

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
  estimates <- cSEM::summarize(object)$Estimates

  vals <- rbind(
    estimates$Loading_estimates,
    estimates$Path_estimates
  )

  par <- stringr::str_remove_all(vals$Name, " ")
  par <- stringr::str_replace_all(par, "\\.", ":")
  lhs <- stringr::str_split_i(par, "=~|~", i = 1L)
  rhs <- stringr::str_split_i(par, "=~|~", i = 2L)
  op  <- stringr::str_extract(par, "=~|~")
  est <- vals$Estimate
  se  <- vals$Std_err
  pvalue <- vals$p_value
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
files <- dir("inst/subset/")
results <- do.call(c,
  lapply(paste0("inst/subset/",  files), readRDS)
)

drop <- vapply(results, FUN.VALUE = logical(1L), FUN = is.null)
results <- results[!drop]

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

    resi[setdiff(cols, colnames(resi))] <- NA
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
    dplyr::summarize(
      se = sd(bias),
      bias = mean(bias, na.rm = TRUE),
      ci.lower = bias - se * qnorm(0.975),
      ci.upper = bias + se * qnorm(0.975)
    ) |>
    ggplot(aes(x = method, y = bias, ymin = ci.lower,
               ymax = ci.upper, colour = method, fill = method)) +
    geom_col(alpha = 0.2) +
    geom_errorbar(width = 0.7) +
    facet_grid(rows = vars(cond), cols = vars(ncat), scales = "fixed") +
    ggtitle(sprintf("Bias for %s by method with 95%s CI", param, "%")) +
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

