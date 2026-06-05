#' Unstandardized Parameter Estimates
#'
#' Transform parameter estimates from a fitted PLS-SEM model to observed- and
#' latent-variable scales. Variables not selected through \code{unstandardized}
#' remain on their standardized scales.
#'
#' @param model A fitted \code{PlsModel} object.
#' @param unstandardized Character vector naming variables to unstandardize, or
#'   one of \code{"all"}, \code{"ov"}, or \code{"lv"}.
#' @param se Character string selecting delta-method standard errors
#'   (\code{"delta"}) or no standard errors (\code{"none"}).
#' @param scale.uncertainty Should scale uncertainty be included?
#'   defaults to \code{FALSE}.
#' @param eps Positive numeric finite-difference step used for the delta-method
#'   Jacobian.
#' @param zero.tol Non-negative numeric tolerance below which standard errors
#'   are returned as missing.
#' @param rm.tmp.ov Logical; whether rows involving temporary observed
#'   variables should be removed from the returned parameter table.
#' @param clean.tmp.ind Logical; whether rows involving temporary indicators
#'   should be cleaned from the returned parameter table.
#'
#' @return A \code{PlsSemParTable} containing transformed estimates and (when
#'   requested) delta-method standard errors. The transformed covariance matrix
#'   is stored in the \code{"vcov"} attribute.
#'
#' @examples
#' \dontrun{
#' tpb <- '
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT <~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH <~ b1 + b2
#'
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC + INT:PBC
#' '
#'
#' fit <- pls(tpb, modsem::TPB, bootstrap = TRUE, boot.R = 50)
#' unstandardized_estimates(fit)
#' }
#' @export
unstandardized_estimates <- function(model, unstandardized = "all",
                                     se = c("delta", "none"),
                                     scale.uncertainty = FALSE,
                                     eps = 1e-4, zero.tol = 1e-10,
                                     rm.tmp.ov = TRUE, clean.tmp.ind = TRUE) {
  se <- tolower(se)
  se <- match.arg(se)

  combined <- combinedModel(model)
  fit      <- modelFit(combined)
  info     <- modelInfo(combined)

  lvs      <- info$lvs.linear
  inds.a   <- info$inds.a
  inds.b   <- info$inds.b
  ovs      <- union(inds.a, inds.b)
  ordered  <- intersect(info$ordered, ovs)
  cont     <- setdiff(ovs, ordered)
  scale    <- info$scale
  resvar   <- diag(fit$fitTheta)
  n        <- NROW(modelData(combined))

  # do we have a multilevel/mixed-effects model?
  pls_stopif(isMLM(combined),
    "`unstandardized_estimates()` is not available for",
    "mixed-effects/multilevel models (yet)!"
  )

  if (length(unstandardized) == 1L) {

    uvars <- switch(
      tolower(unstandardized),
      all = c(lvs, ovs),
      ov  = ovs,
      lv  = lvs,
      unstandardized
    )

  } else {
    uvars <- unstandardized

  }

  # check for potential temp_ov/temp_ind
  # relevant for all but `unstandardized="all"`
  if (length(setdiff(ovs, uvars))) {
    tmp0 <- paste0(TEMP_OV_PREFIX, uvars)
    tmp1 <- paste0(uvars, TEMP_IND_SUFFIX)
    tmp2 <- paste0(tmp0, TEMP_IND_SUFFIX)
    tmp  <- c(tmp0, tmp1, tmp2)

    uvars <- union(uvars, intersect(tmp, ovs))
  }

  # validate
  missing <- setdiff(uvars, c(ovs, lvs))
  pls_stopif(length(missing),
    "Unknown variables in `unstandardized`! Variables:",
    paste0(missing, collapse = ", ")
  )

  parTable <- getParTableEstimates(
    combined, rm.tmp.ov = FALSE, clean.tmp.ind = FALSE
  )

  # get target sds for observed variables, and placeholders for lvs
  sds <- stats::setNames(
    rep(1, length(ovs) + length(lvs)), nm = c(ovs, lvs)
  )

  # for continuous (unstandardized) variables we can just use the observed sds
  for (x in intersect(cont, uvars)) {
    sds[[x]] <- tryCatch(scale[[x]], error = \(...) NA)
  }

  # for ordered variables we try to fix the residual variance to 1
  # (for indicators which have residual variances)
  for (x in intersect(ordered, uvars)) {
    if (x %in% inds.b) sds[[x]] <- 1 # best guess
    else sds[[x]] <- tryCatch(1 / sqrt(resvar[[x]]), error = \(...) NA)
  }

  # if we don't know, we assume they are standard normal
  sds[is.na(sds)] <- 1

  parTableUstd <- unstandardizedEstimatesInternal(
    parTable      = parTable,
    unstandardize = uvars,
    lvs           = lvs,
    inds.a        = inds.a,
    inds.b        = inds.b,
    sds           = sds
  )

  # full vcov
  vcov <- modelParams(combined)$vcov
  if (se == "delta" && NROW(vcov) && NCOL(vcov)) {
    # V[pars]
    pars <- paste0(parTable$lhs, parTable$op, parTable$rhs)
    Vpar <- expandVcov(vcov, labels = pars)

    # V[vars]
    # use approximated std.errors
    sdsu <- sds[uvars]
    sdss <- sds[setdiff(names(sds), uvars)]

    if (scale.uncertainty) {
      # add uncertainty estimates to the scale of the variables
      Vsds <- diag(sdsu^2 / (2 * (n - 1)), nrow = length(sdsu))
      dimnames(Vsds) <- list(names(sdsu), names(sdsu))

      # V[vars,pars]
      V <- diagPartitioned(Vpar, Vsds)

    } else {
      # else we just consider the uncertainty in the parameters
      V <- Vpar

    }

    # jacobian for delta method standard errors
    J <- matrix(
      0, nrow = length(pars), ncol = NCOL(V),
      dimnames = list(pars, colnames(V))
    )

    # if scale.uncertainty=FALSE i never reaches any of the values in
    # idx.var, and thus it stays constant.
    idx.par <- seq_along(pars)
    idx.var <- seq_along(sdsu) + length(pars)

    sds0  <- c(sdsu, sdss)
    par0  <- parTable
    x0    <- c(parTable$est, sdsu)
    y0    <- parTableUstd$est

    for (i in seq_len(NROW(V))) {
      # get templates
      pari   <- par0
      sdsi   <- sds0
      xi     <- x0
      xi[i]  <- x0[i] + eps

      # fill in templats
      pari$est                 <- xi[idx.par]
      sdsi[seq_along(idx.var)] <- xi[idx.var]

      # get estimated solution
      parTableUstd.i <- unstandardizedEstimatesInternal(
        parTable      = pari,
        unstandardize = uvars,
        lvs           = lvs,
        inds.a        = inds.a,
        inds.b        = inds.b,
        sds           = sdsi
      )

      yi <- parTableUstd.i$est

      J[,i] <- (yi - y0) / eps
    }

    Vstd <- J %*% V %*% t(J)

    se <- sqrt(diag(Vstd))
    se[se <= zero.tol] <- NA_real_

    parTableUstd$se <- se

  } else {
    parTableUstd$se <- NA_real_
    Vstd <- NULL
  }

  # add z-stats and format
  parTableOut <- addZStatsParTable(parTableUstd)

  if (rm.tmp.ov)
    parTableOut <- removeTempOvRowsParTable(parTableOut)

  if (clean.tmp.ind)
    parTableOut <- cleanTempIndRowsParTable(parTableOut)

  parTableOut <- plssemParTable(parTableOut)

  # add vcov?
  if (!is.null(Vstd))
    attr(parTableOut, "vcov") <- plssemMatrix(Vstd, is.public = TRUE)

  parTableOut
}


unstandardizedEstimatesInternal <- function(parTable,
                                            unstandardize,
                                            lvs, inds.a, inds.b,
                                            sds) {
  # keep original parTable for restoring the higher order measurement model
  uinds.a  <- intersect(inds.a, unstandardize)
  uinds.b  <- intersect(inds.b, unstandardize)
  ulvs     <- intersect(lvs,    unstandardize)

  # Indicators Mode A and Mode B
  for (ind in c(uinds.a, uinds.b)) {
    idxl  <- which(parTable$op == "=~" & parTable$rhs == ind)
    idxc  <- which(parTable$op == "<~" & parTable$rhs == ind)
    idxvl <- which(parTable$op == "~~" & parTable$lhs == ind)
    idxvr <- which(parTable$op == "~~" & parTable$rhs == ind)
    idxt  <- which(parTable$op == "|"  & parTable$lhs == ind)

    target <- sds[[ind]]

    parTable[idxl,  "est"] <- parTable[idxl,  "est"] * target # loadings
    parTable[idxc,  "est"] <- parTable[idxc,  "est"] / target # weights
    parTable[idxvl, "est"] <- parTable[idxvl, "est"] * target
    parTable[idxvr, "est"] <- parTable[idxvr, "est"] * target
    parTable[idxt,  "est"] <- parTable[idxt,  "est"] * target
  }

  # Latent Variables (Mode A) and Composites (Mode B)
  for (lv in ulvs) {

    # Measurement model
    idxma <- which(parTable$op == "=~" & parTable$lhs == lv)
    idxmb <- which(parTable$op == "<~" & parTable$lhs == lv)

    pls_stopif(length(idxma) && length(idxmb),
      "Did not expect latent variable to be both of mode A and B!",
      "Variable:", paste0("`", lv, "`")
    )

    if (length(idxma)) { # MODE A
      # fix first loading to 1
      loadings <- parTable[idxma, "est", drop = TRUE]
      parTable[idxma, "est"] <- loadings / loadings[[1L]]

      first.ind <- parTable[idxma[[1L]], "rhs", drop = TRUE]
      first.res <- parTable[
        parTable$lhs == first.ind &
        parTable$op  == "~~"      &
        parTable$rhs == first.ind, "est", drop = TRUE
      ][[1L]]
      first.var <- sds[[first.ind]]^2

      target <- sqrt(first.var - first.res)

    } else if (length(idxmb)) { # MODE B

      # fix first loading to 1
      weights <- parTable[idxmb, "est", drop = TRUE]

      uweights <- weights / weights[[1L]]
      names(uweights) <- parTable[idxmb, "rhs", drop = TRUE]

      parTable[idxmb, "est"] <- unname(uweights)

      k <- length(uweights)
      S <- matrix(
        0, nrow = k, ncol = k,
        dimnames = list(names(uweights), names(uweights))
      )

      for (i in seq_len(k)) {
        lhs <- names(uweights)[[i]]

        for (j in seq_len(i)) {
          rhs <- names(uweights)[[j]]

          p <- parTable[
            (parTable$lhs == lhs & parTable$op == "~~" & parTable$rhs == rhs) |
            (parTable$lhs == rhs & parTable$op == "~~" & parTable$rhs == lhs),
            "est", drop = TRUE
          ]

          S[i, j] <- S[j, i] <- p
        }
      }

      target <- c(sqrt(t(uweights) %*% S %*% uweights))
    }

    sds[[lv]] <- target

    # (co-)variances
    idxvl <- which(parTable$op == "~~" & parTable$lhs == lv)
    idxvr <- which(parTable$op == "~~" & parTable$rhs == lv)
    parTable[idxvl, "est"] <- parTable[idxvl, "est"] * target
    parTable[idxvr, "est"] <- parTable[idxvr, "est"] * target

    # outgoing paths
    idxp0 <- which(parTable$rhs == lv & parTable$op == "~")
    parTable[idxp0, "est"] <- parTable[idxp0, "est"] / target

    # incoming paths
    idxp1 <- which(parTable$lhs == lv & parTable$op == "~")
    parTable[idxp1, "est"] <- parTable[idxp1, "est"] * target
  }

  intTerms <- getIntTerms(parTable)

  for (intTerm in intTerms) {
    # first we rescale the interaction term
    # note that we've already rescaled the lhs side, we just need to
    # update the rhs side, we have
    # std(b3) =     b3  * sd(x) *  sd(z) / sd(y)
    #     b3  = std(b3) * sd(y) / (sd(x) * sd(z))
    # at this point we have b3 = std(b3) * sd(y), and thus divide
    # by sd(x) * sd(z).

    elems <- stringr::str_split_1(intTerm, pattern = ":")
    idxp  <- which(parTable$rhs == intTerm & parTable$op == "~")
    idxvl <- which(parTable$lhs == intTerm & parTable$op == "~~")
    idxvr <- which(parTable$rhs == intTerm & parTable$op == "~~")
    target <- prod(sds[elems])

    parTable[idxp,  "est"] <- parTable[idxp,  "est"] / target
    parTable[idxvl, "est"] <- parTable[idxvl, "est"] * target
    parTable[idxvr, "est"] <- parTable[idxvr, "est"] * target
  }

  parTable
}
