calcImpliedConstructCorrMat <- function(model, saturated = FALSE) {
  fit <- model$fit
  M   <- model$matrices

  if (saturated)
    return(M$C)

  Gamma  <- t(fit$fitStructural)
  Lambda <- fit$fitLambda
  I      <- diag(NCOL(Gamma))
  Phi    <- fit$fitCov

  B <- I - Gamma
  B.inv <- solve(B)

  B.inv %*% Phi %*% t(B.inv)
}


calcImpliedIndicatorCorrMat <- function(model, saturated = FALSE) {
  fit    <- model$fit
  M      <- model$matrices
  mode.b <- model$info$mode.b

  S      <- M$S
  Phi    <- calcImpliedConstructCorrMat(model, saturated = saturated)
  Lambda <- fit$fitLambda
  Theta  <- diag2(S) - diag2(Lambda %*% t(Lambda))

  Implied <- Lambda %*% Phi %*% t(Lambda) + Theta

  # Correct blocks for formative constructs
  if (length(mode.b)) {
    MeasurementB <- M$select$lambda[ ,mode.b, drop=FALSE]
    BlocksB <- (MeasurementB %*% t(MeasurementB)) == 1
    Implied[BlocksB] <- S[BlocksB]
  }

  Implied
}


calcSRMR <- function(model, saturated = FALSE, diagonal = TRUE) {

  Expected <- calcImpliedIndicatorCorrMat(model, saturated = saturated)
  Observed <- model$matrices$S

  Diff <- cov2cor(Expected) - cov2cor(Observed)
  sqrt(mean(Diff[lower.tri(Diff, diag = diagonal)]^2))
}


calcChisq <- function(model, saturated = FALSE) {
  O <- model$matrices$S
  E <- calcImpliedIndicatorCorrMat(model, saturated = saturated)
  N <- NROW(model$data)
  p <- NCOL(O)
  Einv <- solve(E)

  as.vector(
    (N - 1) * (tr(O %*% Einv) - log(det(O %*% Einv)) - p)
  )
}
