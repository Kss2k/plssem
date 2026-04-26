getConsistentCorrMat <- function(model, Q) {
  lvs          <- model@info$lvs.linear
  intTermElems <- model@info$intTermElems
  intTermNames <- model@info$intTermNames
  C            <- model@matrices$C
  selectFrom   <- model@matrices$select$nlinFrom

  for (i in seq_along(lvs)) {
    for (j in seq_len(i - 1)) {
      x <- lvs[[i]]
      y <- lvs[[j]]
      C[x, y] <- C[y, x] <- C[x, y] / sqrt(Q[[x]]^2 * Q[[y]]^2)
    }
  }

  if (!model@info$is.nlin)
    return(C)

  for (i in seq_along(intTermElems)) {
    xz.i <- intTermNames[[i]]

    for (j in seq_len(i)) {
      xz.j <- intTermNames[[j]]
      pij  <- f2(xz.i, xz.j, selectFrom, .Q = Q, .H = model@factorScores)
      C[xz.i, xz.j] <- C[xz.j, xz.i] <- pij
    }

    for (y in lvs) {
      pij <- f2(xz.i, y, selectFrom, .Q = Q, .H = model@factorScores)
      C[xz.i, y] <- C[y, xz.i] <- pij
    }
  }

  C
}


getConsistentLoadings <- function(model, Q) {
  lvs     <- model@info$lvs.linear
  indsLvs <- model@info$indsLvs
  lambda  <- model@matrices$lambda
  modes   <- model@info$modes

  for (lv in lvs) {
    wq   <- lambda[indsLvs[[lv]], lv]
    mode <- modes[[lv]]

    lq <- switch(mode,
      A = wq %*% (Q[[lv]] / t(wq) %*% wq),
      B = wq,
      NA_real_
    )

    lambda[indsLvs[[lv]], lv] <- lq
  }

  lambda
}


getConstructQualities <- function(model) {
  lambda  <- model@matrices$lambda
  lvs     <- model@info$lvs.linear
  modes   <- model@info$modes
  indsLvs <- model@info$indsLvs
  rel     <- model@info$reliabilities
  S       <- model@matrices$S
  isHigherOrderComposite <- model@info$isHigherOrderComposite

  Q        <- vector("numeric", length(lvs))
  names(Q) <- lvs

  for (lv in lvs) {
    inds.lv <- indsLvs[[lv]]

    if (lv %in% names(rel)) {
      Q[[lv]] <- sqrt(rel[[lv]])
      next

    } else if (isHigherOrderComposite[[lv]] && !is.null(rel)) {
      w.i   <- lambda[inds.lv, lv, drop = FALSE]
      S.ii  <- S[inds.lv, inds.lv, drop = FALSE]
      w.i.t <- t(w.i)

      idx         <- intersect(inds.lv, names(rel))
      diag(S.ii)[idx] <- rel[idx]
      Q[[lv]] <- sqrt(c(w.i.t %*% S.ii %*% w.i))
      next

    } else if (length(inds.lv) <= 1L || modes[[lv]] != "A") {
      Q[[lv]] <- 1
      next
    }

    w.i   <- lambda[inds.lv, lv, drop = FALSE]
    S.ii  <- S[inds.lv, inds.lv, drop = FALSE]
    w.i.t <- t(w.i)

    c.i <- sqrt(
      (w.i.t %*% (S.ii - diag2(S.ii)) %*% w.i) /
      (w.i.t %*% (w.i %*% w.i.t - diag2(w.i %*% w.i.t)) %*% w.i)
    )

    Q[[lv]] <- sqrt((w.i.t %*% w.i)^2 * c.i^2)
  }

  Q
}


needsAdditionalDistributionalAssumptions <- function(terms) {
  elems   <- stringr::str_split(terms, pattern = ":")
  k       <- vapply(elems, FUN.VALUE = integer(1L), FUN = length)
  quadcube <- vapply(elems, FUN.VALUE = integer(1L), FUN = \(x) any(duplicated(x)))
  any(k | quadcube)
}


# getConsistentCorrelationMatrix — work in progress, currently disabled.
# See nlin_djikstra_correction.R for the partial implementation.
