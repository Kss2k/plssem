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


#' Titanic Passenger Survival Data Set.
#'
#' @name titanic
#' @docType data
#'
#' @description
#' This dataset has been re-packaged for convenience from
#' https://github.com/paulhendricks/titanic
#'
#' \describe{
#' \item{PassengerId}{Passenger ID}
#' \item{Survived}{Passenger Survival Indicator}
#' \item{Pclass}{Passenger Class}
#' \item{Name}{Name}
#' \item{Sex}{Sex}
#' \item{Age}{Age}
#' \item{SibSp}{Number of Siblings/Spouses Aboard}
#' \item{Parch}{Number of Parents/Children Aboard}
#' \item{Ticket}{Ticket Number}
#' \item{Fare}{Passenger Fare}
#' \item{Cabin}{Cabin}
#' \item{Embarked}{Port of Embarkation}
#' \item{Female}{Dummy variable for \code{Sex="female"}}
#' }
#'
#' @format A data frame with 1309 rows and 12 variables:
#' @source https://www.kaggle.com/c/titanic/data
#'
#' @examples
#'
#' fit <- pls("Survived ~ Age + Female + Age:Female",
#'            data = titanic, ordered = "Survived")
#' pls_predict(fit, benchmark = "acc")
NULL
