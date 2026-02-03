
getConsistenCorrMat <- function(model, Q) {
  lVs          <- model$info$lVs.linear
  intTermElems <- model$info$intTermElems
  intTermNames <- model$info$intTermNames
  C            <- model$matrices$C
 
  for (i in seq_along(lVs)) {
    for (j in seq_len(i - 1)) {
      x <- lVs[[i]]
      y <- lVs[[j]]
      C[x, y] <- C[y, x] <- C[x, y] / sqrt(Q[[x]]^2 * Q[[y]]^2)
    }
  }

  selectFrom <- model$matrices$nlinSelectFrom
  for (i in seq_along(intTermElems)) {
    xz.i <- intTermNames[[i]]

    for (j in seq_len(i)) {
      xz.j <- intTermNames[[j]]
      C[xz.i, xz.j] <- f2(xz.i, xz.j, selectFrom, .Q = Q, .H = model$factorScores)
    }
    

    for (y in lVs) {
      xz.j <- intTermNames[[j]]
      C[xz.i, y] <- f2(xz.i, y, selectFrom, .Q = Q, .H = model$factorScores)
    }

  }

  C
}


getConsistentLoadings <- function(model, Q) {
  lVs     <- model$info$lVs.linear
  indsLvs <- model$info$indsLvs
  lambda  <- model$matrices$lambda

  for (lV in lVs) {
    wq <- lambda[indsLvs[[lV]], lV]
    lambda[indsLvs[[lV]], lV] <- wq %*% (Q[[lV]] / t(wq) %*% wq)
  }

  lambda
}


getConstructQualities <- function(model) {
  gamma   <- model$matrices$gamma
  lVs     <- model$info$lVs.linear
  indsLvs <- model$info$indsLvs
  lambda  <- model$matrices$lambda
  # S = sample covariance matrix 
  S <- model$matrices$S
  Q <- vector("numeric", length(lVs))
  names(Q) <- lVs

  for (lv in lVs) {
    inds.lv <- indsLvs[[lv]]
    
    if (length(inds.lv) <= 1L) {
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

    Q.i2 <- (w.i.t %*% w.i)^2 * c.i^2
    Q.i  <- sqrt(Q.i2)

    Q[[lv]] <- Q.i
  }

  Q
}


needsAdditionalDistributionalAssumptions <- function(terms) {
  elems <- stringr::str_split(terms, pattern = ":")
 
  k <- vapply(elems, FUN.VALUE = integer(1L), FUN = length)
  quadcube <- vapply(elems, FUN.VALUE = integer(1L), FUN = \(x) any(duplicated(x)))

  any(k | quadcube)
}


# This didn't quite work...
# Maybe I'll pick it up later...
# getConsistentCorrelationMatrix <- function(Q, C, X) {
#   CQ <- C
#   CQ[TRUE] <- NA
# 
#   E <- stats::setNames(rep(NA, NCOL(CQ)), nm = colnames(CQ))
# 
#   vars     <- colnames(C)
#   intterms <- vars[grepl(":", vars)]
#   linvars  <- setdiff(vars, intterms)
# 
#   # assumeNormal <- needsAdditionalDistributionalAssumptions(intterms)
# 
#   getConsistentCorrelation <- function(i, j) {
#     # i = name (with class) of ith variable
#     # j = name (with class) of jth variable
#     # Q = consistensies
# 
#     if (!is.na(CQ[i, j]))
#       return(CQ[i, j])
# 
#     elems.i <- stringr::str_split_1(i, pattern = ":")
#     elems.j <- stringr::str_split_1(j, pattern = ":")
#     elems   <- c(elems.i, elems.j)
# 
#     k   <- length(elems)
#     k.i <- length(elems.i)
#     k.j <- length(elems.j)
# 
#     if (k <= 2L && i == j) {
#       return(CQ[i, j] <<- CQ[j, i] <<- 1)
# 
#     } else if (k <= 2L) { # no interaction terms
#       cov.i.j <- C[i, j] / (Q[i] * Q[j])
#       return((CQ[i, j] <<- CQ[j, i] <<- cov.i.j))
#     }
# 
#     # From here on all equations are based on
#     # Djikstra & Schermelleh-Engel (2014)
#     # https://link.springer.com/article/10.1007/s11336-013-9370-0
# 
#     stopif(k > 6L || length(unique(elems)) > 4L,
#            "4-way interactions / quartic terms and beyond are not allowed!")
# 
#     freq <- sort(table(elems), decreasing = TRUE) # The frequency of each variable decides the correction
#     vars <- names(freq)
# 
#     w1 <- which.max(freq)
#     w2 <- which.max(freq[-w1])
#     w3 <- which.max(freq[-c(w1, w2)])
#     w4 <- which.max(freq[-c(w1, w2, w3)])
# 
#     m1 <- freq[w1]
#     m2 <- freq[-w1][w2]
#     m3 <- freq[-c(w1, w2)][w3]
#     m4 <- freq[-c(w1, w2, w3)][w4]
#    
#     if (m1 == 1L) {
#       # in the case of no duplicates the correction is quite simple
#       # E.prod <- matrixStats::rowProds(X[, elems, drop = FALSE])
#       # if (assumeNormal)
#       #   return((CQ[i, j] <<- CQ[j, i] <<- 0))
# 
#       E.prod <- C[i, j]
#       Q.prod <- prod(Q[elems])
#       cov.i.j <- E.prod / Q.prod
#       return((CQ[i, j] <<- CQ[j, i] <<- cov.i.j))
#     }
# 
#     # We need to identify the expected means of products before centering.
#     # This is so we can map from E[x*z] to Cov(x, z). Since we're working with
#     # centered variables, it is actually just the consistent covariance.
#     # While it's not efficient, we just recursively call this function
#     # Without correcting we end up returning E[i*j] = E[i]*E[j] + Cov(i, j)
#     # Thus we have to subtract the offset E[i]*E[j], recursively getting the
#     # expectations of the centered variables we get E[i] = cov(i.i, i.j) and so on
#     getExpectedFromElems <- function(l, elems.l) {
#       k.l <- length(elems.l)
# 
#       if (k.l <= 1L)
#         return(0)
#       
#       E.l <- E[[l]]
#       if (!is.na(E.l))
#         return(E.l)
# 
#       li <- elems.l[1L]
#       lj <- paste0(elems.l[-1L], collapse = ":")
#      
#       # E[l] = cov(l.i, l.j)
#       E.l <- getConsistentCorrelation(li, lj)
#       (E[l] <<- E.l)
#     }
#    
#     E.i <- getExpectedFromElems(i, elems.i)
#     E.j <- getExpectedFromElems(j, elems.j)
#     offset <- E.i * E.j
# 
#     #==========================================================================#
#     if (k == 3L)  {
#     #==========================================================================#
# 
#       # E[x^3]
#       if (m1 == 3L) {
#         # Here we have to make some normality assumptions
#         # leading us to the conclustion that E[x^3] = 0
#         return(CQ[i, j] <<- CQ[j, i] <<- 0)
# 
#       # E[x^2 * z]
#       } else if (m1 == 2L && m2 == 1L) {
# 
#         # if (assumeNormal)
#         #   return((CQ[i, j] <<- CQ[j, i] <<- 0))
# 
#         # E[x.c^2 * z.c] = Q.x^2 * Q.z * E[x^2 * z]
#         # E[x^2 * z] = E[x.c^2 * z.c] / (Q.x^2 * Q.z)
#         x <- vars[1L]
#         z <- vars[2L]
# 
#         Qx <- Q[x]
#         Qz <- Q[z]
# 
#         C.xc2.zc <- C[i, j]
#         # E.xc2.zc <- C.xc2.zc + offset
#         E.xc2.zc <- mean(X[,x]^2*X[,z])
# 
#         # We often have some offset caused by the model centering the
#         # interaction terms. Which is not accounted for in E[x * z * ... * y]
#         # C[x, z, ..., y] = E[x * z * ... * y] - offset
# 
#         browser()
#         cov.x2.z <- E.xc2.zc / (Qx^2 * Qz) - offset
#         return((CQ[i, j] <<- CQ[j, i] <<- cov.x2.z))
# 
#       # E[x*z*y]
#       } else {
#         # Should already handled
#         stop("Unexpected combo (m1=1 & m2=2)?")
# 
#       }
#         
#     #==========================================================================#
#     } else if (k == 4L) {
#     #==========================================================================#
# 
#       # E[x^4]
#       if (m1 == 4L) {
#         # Here we just have to assume normality
#         return((CQ[i, j] <<- CQ[j, i] <<- 2))
# 
#       # E[x^3 * z] = E[x^2 * xz]
#       } else if (m1 == 3L && m2 == 1L) {
#         x <- vars[1L]
#         z <- vars[2L]
# 
#         # if (assumeNormal) {
#         #   p.ij <- getConsistentCorrelation(x, z)
#         #   return((CQ[i, j] <<- CQ[j, i] <<- 3 * p.ij))
#         # }
# 
#         # E[x.c^3 * z.c] = Q.x^3 * Q.z * E[x^3 * z] + 3 * E[x.c * z.c] * (1 - Q.x^2)
#         # E[x.c^3 * z.c] - 3 * E[x.c * z.c] * (1 - Q.x^2) = Q.x^3 * Q.z * E[x^3 * z] 
#         # (E[x.c^3 * z.c] - 3 * E[x.c * z.c] * (1 - Q.x^2)) / (Q.x^2 * Q.z) = E[x^3 * z] 
#         # E[x^3 * z] = (E[x.c^3 * z.c] - 3 * E[x.c * z.c] * (1 - Q.x^2)) / (Q.x^2 * Q.z)
# 
#         E.xc.zc <- C[x, z] # mean(X[,x], X[,z])
#         # E.xc3.zc <- C[i, j] + offset
#         E.xc3.zc <- mean(X[,x]^3 * X[,z])
# 
#         Qx <- Q[x]
#         Qz <- Q[z]
# 
#         E.x3.z <- (E.xc3.zc - 3 * E.xc.zc * (1 - Qx^2)) / (Qx^2 * Qz)
# 
#         cov.x3.z <- E.x3.z - offset
#         return((CQ[i, j] <<- CQ[j, i] <<- cov.x3.z))
# 
#       # E[x^2 * z^2] = E[xx * zz] = E[xz * xz]
#       } else if (m1 == 2L && m2 == 2L) {
#         # UNDER NORMALITY THIS IS DIFFERENT
#           # Under normality we would get the same as above
#           # But we don't have to assume normality here.
#           # We could include a flag here. If we have a quadratic term
#           # it would be consistent to assume normality and just return 0
#         # E[x.c^2 * z.c^2] = Q.x^2 * Q.z^2 * (E[x^2 * z^2] - 1) + 1
#         # E[x.c^2 * z.c^2] - 1 = Q.x^2 * Q.z^2 * (E[x^2 * z^2] - 1)
#         # (E[x.c^2 * z.c^2] - 1) / (Q.x^2 * Q.z^2) = E[x^2 * z^2] - 1
#         # (E[x.c^2 * z.c^2] - 1) / (Q.x^2 * Q.z^2) + 1 = E[x^2 * z^2]
#         # E[x^2 * z^2] = (E[x.c^2 * z.c^2] - 1) / (Q.x^2 * Q.z^2) + 1
# 
#         x <- vars[1L]
#         z <- vars[2L]
# 
#         # E.xc2.zc2 <- C[i, j] + offset
#         E.xc2.zc2 <- mean(X[,x]^2 * X[,z]^2)
# 
#         Qx <- Q[x]
#         Qz <- Q[z]
# 
#         E.x2.z2 <- (E.xc2.zc2 - 1) / (Qx * Qz) + 1
#         cov.x2.z2 <- E.x2.z2 - offset
#         return((CQ[i, j] <<- CQ[j, i] <<- cov.x2.z2))
# 
#       # E[x^2 * z * y] = E[xx * z * y] = E[xz * x * y]
#       } else if (m1 == 2L && m2 == 1L) {
#         # E[x.c^2 * z.c * y.c] = Q.x^2 * Q.z * Q.y * E[x^2 * z * y] + p(x,z) * Q.z * Q.y * (1 - Q.x^2)
#         # E[x.c^2 * z.c * y.c] - p(x,z) * Q.z * Q.y * (1 - Q.x^2) = Q.x^2 * Q.z * Q.y * E[x^2 * z * y] 
#         # (E[x.c^2 * z.c * y.c] - p(x,z) * Q.z * Q.y * (1 - Q.x^2)) / (Q.x^2 * Q.z * Q.y) = E[x^2 * z * y] 
#         # E[x^2 * z * y] = (E[x.c^2 * z.c * y.c] - p(z,y) * Q.z * Q.y * (1 - Q.x^2)) / (Q.x^2 * Q.z * Q.y)
#         
#         x <- vars[1L]
#         z <- vars[2L]
#         y <- vars[3L]
# 
#         # E.xc2.zc.yc <- C[i, j] + offset
#         E.xc2.zc.yc <- mean(X[,x]^2 * X[,z] * X[,y])
# 
#         Qx <- Q[x]
#         Qz <- Q[z]
#         Qy <- Q[y]
# 
#         cov.z.y <- getConsistentCorrelation(z, y)
#         E.x2.z.y <- (E.xc2.zc.yc - cov.z.y * Qz * Qy * (1 - Qx^2)) / (Qx^2 * Qz * Qy)
# 
#         cov.x2.z.y <- E.x2.z.y - offset
#         return((CQ[i, j] <<- CQ[j, i] <<- cov.x2.z.y))
# 
#       } else {
#         browser()
#         stop("Unexpected term!")
#       }
# 
#     #==========================================================================#
#     } else if (k == 5L) {
#     #==========================================================================#
#         warning("DEBUG: k=5 not implemented yet!")
# 
#     #==========================================================================#
#     } else if (k == 6L) {
#     #==========================================================================#
#         warning("DEBUG: k=6 not implemented yet!")
#     }
# 
#   }
# 
#   diag(CQ[linvars, linvars]) <- 1
# 
#   for (i in seq_along(linvars)) {
#     x <- vars[[i]]
# 
#     for (j in seq_len(i-1L)) {
#       y <- vars[[j]]
# 
#       CQ[x, y] <- CQ[y, x] <- C[x, y] / sqrt(Q[[x]]^2 * Q[[y]]^2)
#     }
#   }
# 
#   for (xz in intterms) for (y in linvars)
#     getConsistentCorrelation(xz, y)
#   
#   for (i in seq_along(intterms)) {
#     xz.i <- intterms[[i]]
# 
#     for (j in seq_len(i)) {
#       xz.j <- intterms[[j]]
# 
#       getConsistentCorrelation(xz.i, xz.j)
#     }
#   }
# 
#   CQ
# }

