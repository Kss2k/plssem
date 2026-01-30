#' randomSlopes
#'
#' @name randomSlopes
#' @docType data
#' @description A simulated dataset.
#' @examples
#'
#' pls.syntax <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   Y ~ X + Z + X:Z
#' "
#' 
#' lme4.syntax <- "Y ~ X + Z + (1 + X + Z | cluster)"
#'
#' fit <- pls(pls.syntax, data = randomSlopes, lme4.syntax = lme4.syntax,
#'            cluster = "cluster")
#' summary(fit)
NULL


#' randomIntercepts
#'
#' @name randomIntercepts
#' @docType data
#' @description A simulated dataset.
#' @examples
#'
#' pls.syntax <- "
#'   f =~ y1 + y2 + y3
#'   f ~ x1 + x2 + x3 + w1 + w2
#' "
#' 
#' lme4.syntax <- "f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)"
#' 
#' 
#' fit <- pls(model.pls, data = randomIntercepts,
#'            lme4.syntax = lmer.syntax, cluster = "cluster")
#' summary(fit)
NULL
