#' Compute the residuals
#'
#' @param vY vector, y the outcome
#' @param mX matrix, X the predictors
#' @param vBeta vector, the parameters beta
#' @return vResiduals, vector of residuals
#' @examples
#' residuals(vY, mX, vBeta)
residuals <- function(vY, mX, vBeta) {
  vResiduals <- vY - mX %*% vBeta
  return(vResiduals)
}
#' Loss function
#'
#' @param vResiduals vector, the residuals
#' @param vBeta vector, the parameters beta
#' @param dLambda double, the lambda parameter
#' @param dAlpha double, the alpha parameter
#' @return dLoss, double the value of the loss function
#' @examples
#' loss(vResiduals, vBeta, dLambda, dAlpha)
loss <- function(vResiduals, vBeta, dLambda, dAlpha) {
  
  vBeta <- as.matrix(vBeta)
  
  # init
  iN <- nrow(vResiduals)
  dConstr <- dLambda*((1 - dAlpha)/(2*t(vBeta) %*% vBeta) + dAlpha * norm(vBeta, type = "1"))
  
  # compute the loss function
  dLoss <- (2*iN)^(-1) * (t(vResiduals) %*% vResiduals) + dConstr
  
  return(dLoss)
}

#' Change of loss function
#'
#' @param vBeta0 vector, last iterations beta
#' @param vBeta1 vector, this iterations beta
#' @param vResiduals0 vector, residual vector
#' @param vResiduals1, vector, residual vector
#' @param dLambda double, the lambda parameter
#' @param dAlpha double, the alpha parameter
#' @return dLoss, double the value of the loss function
#' @examples
#' loss_change(vBeta0, vBeta1, vResiduals, dLambda, dAlpha)
loss_change <- function(vBeta_0, vBeta_1, vResiduals0, vResiduals1, dLambda, dAlpha) {
  
  # compute the loss
  dLoss0 <- loss(vResiduals0, vBeta_0, dLambda, dAlpha)
  dLoss1 <- loss(vResiduals1, vBeta_1, dLambda, dAlpha)
  
  # compute the change
  dChange <- (dLoss0 - dLoss1)/dLoss0
  
  return(dChange)
}
#' Get the D_jj element, vectorized
#'
#' @param vBeta vector, vector of betas
#' @param dEps double, precision variable
#' @return vD, vector, diagonal elements of D
#' @examples
#' computeDjj(vBeta, dEps)
computeDjj <- Vectorize(function(vBeta, dEps) {
  dD <- 1 / (max(c(abs(vBeta), dEps)))
  return(dD)
}, vectorize.args = "vBeta")
#' Create the D matrix
#'
#' @param vBeta vector, vector of betas
#' @param dEps double, precision variables
#' @return mD, matrix, matrix D with d_jj on diagonal
#' @examples
#' getD(vBeta, dEps)
getD <- function(vBeta, dEps) {
  vDiag <- computeDjj(vBeta, dEps)
  mD <- diag(vDiag)
  return(mD)
}
#' Create the A matrix
#'
#' @param mA1 matrix, pulled out value of A matrix
#' @param mD matrix, the matrix D
#' @param dAlpha double, alpha parameter
#' @param dLambda double, lambda parameter
#' @return mD, matrix, matrix D with d_jj on diagonal
#' @examples
#' getA(mX, mD, dAlpha, dLambda)
getA <- function(mA1, mD, dAlpha, dLambda) {
  mA <- mA1 + dLambda * dAlpha * mD
  return(mA)
  
}
#' Calculate the update for Beta
#'
#' @param mX matrix, the predictor's data
#' @param vY vector, the vector of the outcome variable
#' @param mA matrix, A matrix
#' @return vBetaUpdate, vector, the updated betas
#' @examples
#' updateBeta(mX, vY, mA)
updateBeta <- function(mX, vY, mA) {
  
  # init
  iN <- nrow(mX)
  
  # calculate
  vBetaUpdate <- solve(mA, (1/iN) * (t(mX) %*% vY))
  
  return(vBetaUpdate)
  
}
#' Perform one iteration of the MM
#'
#' @param mX matrix, the predictor's data
#' @param vY vector, the vector of the outcome variable
#' @param vBeta vector, the betas
#' @param dEps double, the precision epsilon
#' @param dAlpha double, the alpha parameter
#' @param dLambda double, the lambda parameter
#' @return vBetaUpdate, vector, the updated betas
#' @examples
#' itElasticNetMM(mX, vY, vBeta, dEps, dAlpha, dLambda)
itElasticNetMM <- function(mX, vY, vBeta, dEps, dAlpha, dLambda) {
  
  # perform one iteration
  mD <- getD(vBeta, dEps)
  mA <- getA(mX, mD, dAlpha, dLambda)
  vBetaUpdate <- updateBeta(mX, vY, mA)
  
  return(vBetaUpdate)
}
#' Run the full MM
#' 
#' @param mX matrix, the predictor's data
#' @param vY vector, the vector of the outcome variable
#' @param vBeta vector, the betas
#' @param dEps double, the precision epsilon
#' @param dAlpha double, the alpha parameter
#' @param dLambda double, the lambda parameter
ElasticNetMM <- function(mX, vY, vBeta, dEps, dAlpha, dLambda) {
  # loop objects
  iP <- ncol(mX)
  iN <- nrow(mX)
  ik <- 1
  dLossChange <- 0
  vBeta0 <- runif(ncol(mX))
  dEps <- 10e-6
  mXtX <- t(mX) %*% mX
  
  # Pull out part one of A
  dA1 <- (1/iN) * mXtX + (dLambda * (1 - dAlpha)) * diag(iP)
  
  while ((ik == 1) | (dLossChange > dEps)) {
    
    # update counter
    ik <- ik + 1
    
    # perform one update
    mD <- getD(vBeta0, dEps)
    mA <- getA(dA1, mD, dAlpha, dLambda)
    vBeta1 <- updateBeta(mX, vY, mA)
    
    # get residuals
    vResiduals0 <- residuals(vY, mX, vBeta0)
    vResiduals1 <- residuals(vY, mX, vBeta1)
    
    # compute loss change
    dLossChange <- loss_change(vBeta0, vBeta1, vResiduals0, vResiduals1, dLambda, dAlpha)
    
    # output
    print(paste0("Iteration: ", ik, " \n"))
    print(paste0("Loss Change: ", dLossChange, " \n"))
    
    # update the beta
    vBeta0 <- vBeta1
  }
  
  print("Beta Estimate")
  print(vBeta0)
  return(vBeta0)
}

