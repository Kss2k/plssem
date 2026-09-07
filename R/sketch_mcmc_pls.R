PRIOR_OP <- ":~"
PRIOR_OP_ALIAS <- "==" # override this operator, as it's not used for anything
                       # this can be fixed properly in modsem later


mcmc_pls <- function(syntax,
                     data,
                     ...,
                     priors = list(),
                     Q = list(
                       r = \(x, s) as.vector(mvtnorm::rmvnorm(n = 1, mean = x, sigma = s)),
                       d = \(x, y, s, log = TRUE) mvtnorm::dmvnorm(matrix(y, nrow = 1), mean = x, sigma = s, log = log)
                     ),
                     chains = 1L,
                     iter = 2000,
                     warmup = floor(iter / 2),
                     N = 20000,
                     acceptance.rate = \(d) 0.234 + 0.21 / d, # acceptance ratio by the number of dimensions in a block
                     sampler = c("Metropolis-Hastings", "Gibbs"),
                     sample.thresholds = FALSE, # currently too difficult to sample...
                     # capture
                     bootstrap = NULL,
                     consistent = FALSE,
                     mcpls = FALSE,
                     probit = FALSE) {

  sampler <- match.arg(
    tolower(sampler), c("metropolis-hastings", "gibbs")
  )

  # have priors been specified in the syntax?
  syntax <- stringr::str_replace_all(syntax, stringr::fixed(PRIOR_OP), PRIOR_OP_ALIAS)
  input <- modsem::modsemify(syntax, parentheses.as.string = TRUE)
  input[input$op == PRIOR_OP_ALIAS, "op"] <- PRIOR_OP

  inputModel  <- input[input$op != PRIOR_OP, , drop = FALSE]
  inputPriors <- input[input$op == PRIOR_OP, , drop = FALSE]

  if (NROW(inputPriors)) {
    priorsSyntax <- lapply(
      stats::setNames(
        paste0("(function(x) x |> ", inputPriors$rhs, ")"),
        nm = inputPriors$lhs
      ),
      FUN = \(expr) eval(parse(text=expr))
    )

  } else {
    priorsSyntax <- NULL

  }

  fit0 <- pls(
    syntax = parTableToSyntax(input),
    data   = data,
    bootstrap = TRUE, # neccessary
    consistent = consistent,
    mcpls = mcpls,
    probit = probit,
    ...
  )

  ordered <- combinedModel(fit0)@info$ordered
  thresholdStruct0 <- combinedModel(fit0)@thresholdStruct
  data <- modelData(fit0)
  vars <- colnames(data)

  if (sample.thresholds) exclude.pars <- c("~1", ":=")
  else                   exclude.pars <- c("~1", ":=", "|")

  parTable <- getFreeParamsTable(fit0, exclude = exclude.pars)
  parTable$par <- paste0(parTable$lhs, parTable$op, parTable$rhs)
  pars <- parTable[parTable$is.free, "par"]

  vcov <- vcov(fit0, use.labels = FALSE)[pars, pars, drop = FALSE]
  coef <- coef(fit0, use.labels = FALSE)[pars]

  labs <- stats::setNames(
    names(coef(fit0, use.labels = TRUE)),
    names(coef(fit0, use.labels = FALSE))
  )

  missingPriors <- setdiff(pars, names(priors))
  for (missing in missingPriors) {
    lab <- tryCatch(labs[[missing]], error = \(.e) NULL)

    if (!is.null(lab) && lab %in% names(priorsSyntax))
      priors[[missing]] <- priorsSyntax[[lab]]
    else
      priors[[missing]] <- \(...) 1 # flat prior
  }

  lower <- NULL
  upper <- NULL

  L <- function(x) {
    parTablex <- parTable[c("lhs", "op", "rhs", "est", "is.free")]
    parTablex[parTablex$is.free, "est"] <- x

    sim <- simulateDataParTable(
      parTable = parTablex,
      N = N,
      cut = sample.thresholds
    )
   
    sim.ov <- sim$ov

    if (!sample.thresholds && length(ordered)) {
      sim.ov  <- ordinalizeDataFrame(
        df = sim.ov, thresholdStruct = thresholdStruct0
      )
    }

    lower <<- sim$lower
    upper <<- sim$upper

    fit.sim <- fit0
    Y <- Rfast::standardise(as.matrix(sim.ov[vars]))
    S <- Rfast::cova(Y)

    if (!is.null(sim$cluster))
      attr(Y, "cluster") <- sim$cluster

    # Update observed-data (lowest-order) model input
    modelData(fit.sim)  <- Y
    indCorrMatrix(fit.sim) <- S

    if (sample.thresholds && any(parTablex$op == "|")) {
      fit.sim@thresholdStruct <- ThresholdStruct(
        data = Y,
        ordered = ordered
      )
    }

    fit.y <- estimatePLS_Inner(fit.sim)

    y <- coef(fit.y, use.labels=FALSE)[pars]
    mvtnorm::dmvnorm(matrix(y, nrow = 1), mean = coef, sigma = vcov, log = TRUE)
  }

  P <- function(x) {
    p <- numeric(length(x))
    for (i in seq_along(x)) {
      par <- pars[[i]]
      p[[i]] <- priors[[par]](x[[i]])
    }
    log(p)
  }

  # for (chain in chains) {
  pls_stopif(warmup >= iter, "warmup must be less than iter!")
  x  <- coef
  S0 <- vcov
  S  <- S0

  if (sampler == "metropolis-hastings") {
    blocks <- list(seq_along(x))

  } else {
    split <- as.data.frame(splitParameterNames(pars))
    op <- split$op

    oblocks <- lapply(ordered, function(ord) {
      which(split$lhs == ord & op == "|")
    })

    mblock <- list(which(split$op == "=~"))
    pblock <- list(which(split$op == "~"))
    cblock <- list(setdiff( # covariances
      seq_along(x),
      c(unlist(oblocks), unlist(mblock), unlist(pblock))
    ))

    blocks <- c(oblocks, mblock, pblock, cblock)
  }


  # One log proposal-SD multiplier per block
  # If S is already a reasonable estimate of the target covariance,
  # 2.38 / sqrt(d) is a good starting point.
  log.scale <- vapply(
    X = blocks,
    FUN.VALUE = numeric(1),
    FUN = function(idx) {
      d <- length(idx)
      if (d > 0) log(2.38 / sqrt(d)) else NA_real_
    }
  )

  # Number of adaptation steps for each block
  adapt.n <- integer(length(blocks))

  rejections <- numeric(length(blocks))
  acceptances <- numeric(length(blocks))
 
  SAMPLES <- matrix(NA, nrow = iter, ncol = length(x), dimnames = list(NULL, pars))
  for (i in seq_len(iter)) {
    mode <- if (i > warmup) "sampling" else "warmup"
    cat(sprintf(
      "Iter %d/%d, mode: %s, acceptances: %d, rejections: %d, acceptance rate: %.3f, step = %.2g...\n",
      i, iter, mode, sum(acceptances), sum(rejections),
      sum(acceptances)/(sum(acceptances)+sum(rejections)),
      suppressWarnings(mean(scale))
    ))

    if (i %% 10 == 0) {
      # Update S
      wpct <- warmup / iter
      n    <- iter - warmup
      n1   <- floor(i * wpct)
      n0   <- max(0, n - n1)
      sub  <- tail(SAMPLES[seq_len(i), , drop = FALSE], n = n1)

      if (NROW(sub) > 10) {
        S1 <- stats::cov(sub, use = "complete.obs")
        S <- ((n0 - 1) * S0 + (n1 * 1) * S1) / (n0 + n1 - 2)
      }
    }

    x.star <- x

    for (block in seq_along(blocks)) {
      idx <- blocks[[block]]
      d <- length(idx) # dimension

      if (d <= 0)
        next # nothing to do

      # sample proposal
      scale <- exp(log.scale[[block]])
      S.star.b <- scale^2 * S[idx, idx, drop = FALSE]
      x.star[idx] <- Q$r(x = x[idx], s = S.star.b)

      # Constrain by the last iterations lower and upper bounds
      if (!is.null(upper)) x.star[idx] <- pmin(x.star[idx], upper[idx])
      if (!is.null(lower)) x.star[idx] <- pmax(x.star[idx], lower[idx])

      q.star <- Q$d(x = x[idx], y = x.star[idx], s = S.star.b, log = TRUE)
      q.x    <- Q$d(x = x.star[idx], y = x[idx], s = S.star.b, log = TRUE)

      Lx <- L(x)
      Px <- P(x)
      Lx.star <- L(x.star)
      Px.star <- P(x.star)

      a <- min(log(1), (q.x - q.star) + (Lx.star + Px.star) - (Lx + Px)) # q.x/q.star = 1 for symmetric distributions
      k <- log(runif(1, min = 0, max = 1))

      accept <- k<=a && !is.na(k<=a)
      rejections[[block]] <- rejections[[block]] + as.integer(!accept)
      acceptances[[block]] <- acceptances[[block]] + as.integer(accept)
    
      adapt.n[[block]] <- adapt.n[[block]] + 1L

      # Robbins-Monro learning rate
      eta <- 0.5 / (10 + adapt.n[[block]])^0.6

      target <- acceptance.rate(d)
      log.scale[[block]] <- log.scale[[block]] + eta * (as.numeric(accept) - target)

      if (accept)
        x[idx] <- x.star[idx]
    }

    SAMPLES[i,] <- x
  }

  # }
  SAMPLES[(warmup+1):iter,,drop=FALSE]
}
