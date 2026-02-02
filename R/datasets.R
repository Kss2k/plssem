#' randomSlopes
#'
#' @name randomSlopes
#' @docType data
#' @description A simulated dataset.
#' syntax <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   W =~ w1 + w2 + w3
#'   Y ~ X + Z + (1 + X + Z | cluster)
#'   W ~ X + Z + (1 + X + Z | cluster)
#' "
#' 
#' fit <- pls(syntax, data = randomSlopes)
#' fit
NULL


#' randomIntercepts
#'
#' @name randomIntercepts
#' @docType data
#' @description A simulated dataset.
#' @examples
#'
#' syntax <- '
#'   f =~ y1 + y2 + y3
#'   f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
#' '
#' 
#' fit <- pls(syntax, data = randomIntercepts)
#' summary(fit)
NULL


#' randomInterceptsOrdered
#'
#' @name randomInterceptsOrdered
#' @docType data
#' @description A simulated dataset.
#' @examples
#'
#' syntax <- '
#'   f =~ y1 + y2 + y3
#'   f ~ x1 + x2 + x3 + w1 + w2 + (1 | cluster)
#' '
#' 
#' fit <- pls(syntax, data = randomInterceptsOrdered)
#' summary(fit)
NULL


#' randomSlopesOrdered
#'
#' @name randomSlopesOrdered
#' @docType data
#' @description A simulated dataset.
#' @examples
#'
#' syntax <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   W =~ w1 + w2 + w3
#'   Y ~ X + Z + (1 + X + Z | cluster)
#'   W ~ X + Z + (1 + X + Z | cluster)
#' "
#' 
#' fit <- pls(syntax, data = randomSlopesOrdered)
#' fit
#' summary(fit)
NULL


#' TPB_Ordered 
#'
#' @name TPB_Ordered
#' @docType data
#' @description A simulated dataset.
#' @examples
#' 
#' tpb <- ' 
#' # Outer Model (Based on Hagger et al., 2007)
#'   ATT =~ att1 + att2 + att3 + att4 + att5
#'   SN =~ sn1 + sn2
#'   PBC =~ pbc1 + pbc2 + pbc3
#'   INT =~ int1 + int2 + int3
#'   BEH =~ b1 + b2
#' 
#' # Inner Model (Based on Steinmetz et al., 2011)
#'   INT ~ ATT + SN + PBC
#'   BEH ~ INT + PBC 
#' '
#' 
#' fit <- pls(tpb, TPB_Ordered)
#' summary(fit)
NULL


#' oneIntOrdered 
#'
#' @name oneIntOrdered 
#' @docType data
#' @description A simulated dataset.
#' @examples
#' 
#' m <- '
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#' 
#'   Y ~ X + Z + X:Z
#' '
#' 
#' fit <- pls(m, oneIntOrdered)
#' summary(fit)
NULL
