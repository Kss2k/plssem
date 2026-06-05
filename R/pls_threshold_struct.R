ThresholdStruct <- function(data, ordered = NULL, digits = 5L) {
  if (!length(ordered))
    return(methods::new("ThresholdStruct"))

  indices     <- emptyNamedList(ordered)
  levels      <- emptyNamedList(ordered)
  thresholds  <- NULL
  proportions <- NULL

  for (ord in ordered) {
    # round to avoid bad floating point comparisons
    x <- round(data[, ord, drop = TRUE], digits)

    # get values
    levels.x      <- sort(unique(x))
    freq.x        <- table(factor(x, levels = levels.x))
    pct.x         <- cumsum(freq.x) / sum(freq.x)
    proportions.x <- unname(pct.x[-length(pct.x)])
    thresholds.x  <- stats::qnorm(proportions.x)   

    # get indices
    idx.x <- seq_along(thresholds.x) + length(thresholds)

    # set names
    names(proportions.x) <- paste0(ord, "|P", seq_along(proportions.x))
    names(thresholds.x)  <- paste0(ord, "|t", seq_along(thresholds.x))

    # set indices
    indices[[ord]] <- idx.x
    levels[[ord]]  <- levels.x

    # set values
    thresholds  <- c(thresholds, thresholds.x)
    proportions <- c(proportions, proportions.x)
  }

  methods::new("ThresholdStruct",
    ordered     = ordered,
    indices     = indices,
    thresholds  = thresholds,
    proportions = proportions,
    levels      = levels,
    digits      = digits
  )
}


updateProportions <- function(thr, data) {
  if (!length(thr@ordered))
    return(thr)

  for (ord in thr@ordered) {
    # round to avoid bad floating point comparisons
    x <- round(data[, ord, drop = TRUE], thr@digits)

    freq.x        <- table(factor(x, levels = thr@levels[[ord]]))
    pct.x         <- cumsum(freq.x) / sum(freq.x)
    proportions.x <- pct.x[-length(pct.x)]

    # TODO: handle missing levels
    thr@proportions[thr@indices[[ord]]] <- proportions.x
  }

  thr
}


updateThresholds <- function(thr, sim.cont = NULL, zero.tol = .Machine$double.eps^0.5) {
  if (!length(thr@ordered))
    return(thr)

  probs <- thr@proportions
  probs <- pmin(pmax(probs, zero.tol), 1 - zero.tol)

  for (ord in thr@ordered) {
    idx <- thr@indices[[ord]]
    probs.x <- probs[idx]

    if (!is.null(sim.cont)) {
      sim.x <- sim.cont[, ord]
      thr.x <- collapse::fquantile(sim.x, probs = probs.x)
    } else {
      thr.x <- stats::qnorm(probs.x)
    }

    thr@thresholds[idx] <- thr.x
  }

  thr
}
