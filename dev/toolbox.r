' Compute the residuals
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
ElasticNetMM <- function(mX, vY, dEps, dAlpha, dLambda) {
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
    mD <- getD(vBeta0, 10e-8)
    mA <- getA(dA1, mD, dAlpha, dLambda)
    vBeta1 <- updateBeta(mX, vY, mA)
    
    # get residuals
    vResiduals0 <- residuals(vY, mX, vBeta0)
    vResiduals1 <- residuals(vY, mX, vBeta1)
    
    # compute loss change
    dLossChange <- loss_change(vBeta0, vBeta1, vResiduals0, vResiduals1, dLambda, dAlpha)
    
    # output
    # print(paste0("Iteration: ", ik, " \n"))
    # print(paste0("Loss Change: ", dLossChange, " \n"))
    # 
    # update the beta
    vBeta0 <- vBeta1
  }
  
  # print("Beta Estimate")
  # print(vBeta0)
  return(vBeta0)
}
#' MAPE
#' 
#' @param vBeta_A vector, estimates of model A
#' @param vBeta_B vector, estimates of model B
#' @return MAPE, vector - main absolute percentage error
MAPE <- function(vBeta_A, vBeta_B) {
  return(abs((vBeta_A - vBeta_B)/vBeta_A) / 100)
}
#' MAE
#' 
#' @param vBeta_A vector, estimates of model A
#' @param vBeta_B vector, estimates of model B
#' @return MAE, vector - main absolute percentage error
MAE <- function(vBeta_A, vBeta_B) {
  return(abs(vBeta_A - vBeta_B))
}
#' Create a dataframe matching up estimates
#' 
#' @param vBeta_A vector, estimates of model A
#' @param vBeta_B vector, estimates of model B
#' @param vColNames, vector of length 3, col names to use
#' @return dfCompareBetaTable, df
CompareEstimates <- function(vBeta_A, vBeta_B, vColNames = c("GLMNET", "MM", "Predictor")) {
  # compare estimates with each other
  dfCompareBeta <- cbind(vBeta_A, vBeta_B)
  vNameCoef <- rownames(dfCompareBeta)
  
  dfCompareBetaTable <- dfCompareBeta %>% as_tibble()
  dfCompareBetaTable$Predictor <- vNameCoef
  colnames(dfCompareBetaTable) <- vColNames
  
  return(dfCompareBetaTable)
}


#' Function to CV k_fold to look for optimal lambda
#' 
#' @param mX matrix, the predictor's data
#' @param vY vector, the vector of the outcome variable
#' @param nfolds integer, number of folds
#' @param vBeta vector, the betas
#' @param dEps double, the precision epsilon
#' @param dAlpha double, the alpha parameter
#' @param lLambda list, list of lambda parameter
#' @return mRSS, matrix of Residual sum of squares, each column is for folds, each row is for each lambda value
k_fold_lambda <- function(mX, vy, nfolds, vBeta, dEps, dAlpha,lLambda){
  #Shuffle the data
  set.seed(321)
  data <-cbind(vy,mX)
  data<-data[sample(nrow(data)),]
  
  #Create n equal folds, since the indexes are discrete, sometimes fold length differs by 1 observation
  folds <- cut(seq(1,nrow(data)),breaks=nfolds,labels=FALSE)
  
  #Perform cross validation
  mRSS= matrix(NA, nrow=50, ncol=nfolds) #Create a matrix, each column is the RSS for each test set of the folds
  for(i in 1:ncol(mRSS)){
    #Segment the data by fold 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testSet <- data[testIndexes, ]
    trainSet <- data[-testIndexes, ]
    #Use the test and train data partitions on Elastic net model
    y_train <- trainSet[,1]
    y_test <- testSet[,1]
    X_train <- scale(trainSet[,-1])
    X_test <- scale(testSet[,-1])
    
    lRSS = rep(NA, nrow(mRSS))
    for (j in 1:nrow(mRSS)){
      #Calculate RSS for each lambda value
      vBeta_new = ElasticNetMM(X_train, y_train, dEps, dAlpha, lLambda[j]) 
      lRSS[j] =sum(residuals(y_test, X_test, vBeta_new)^2) # is this fine?
    }
    mRSS[,i]=lRSS #Insert RSS values into matrix
  }
  return(mRSS)
}
#' Function to get the root mean square errors for each lambda value
#' 
#' @param mRSS matrix, residual sum of squares over k-folds for each lambda value
#' @return means list, root mean squared errors of each lambda
root_mean <- function(mRSS){
  means = sqrt(rowMeans(mRSS))
  return(means)
}

#' Function for CV that looks for optimal lambda and alpha, produce 1 plot for each alpha value
#' 
#' @param mX matrix, the predictor's data
#' @param vY vector, the vector of the outcome variable
#' @param nfolds integer, number of folds
#' @param vBeta vector, the betas
#' @param dEps double, the precision epsilon
#' @param lAlpha double, the alpha parameter
#' @param lLambda double, the lambda parameter
#' @return lRMSE_min, list of minimum root mean square error for each alpha

k_fold_plots <- function(mX, vy, nfolds, vBeta, dEps, lAlpha,lLambda){
  lRMSE_min=list()
  for (a in 1:length(lAlpha)){
    
    mRSS = k_fold_lambda(mX, vy, nfolds, vBeta, dEps, lAlpha[a],lLambda)
    
    #get root mean square errors for each lambda value
    means = root_mean(mRSS)
    
    #Get index of the min lambda
    ind = which.min(means)
    lambda_min = lLambda[ind]
    cat("Alpha is: ", lAlpha[a],". The minimum lambda is: ", lambda_min, "\n", "The minimum RMSE is: ",  min(means),"\n")
    lRMSE_min[a]=min(means)
    
    #Plot
    df <- data.frame('RMSE'=means,'Log(Lambda)'=log(lLambda))
    plot(df$Log.Lambda., df$RMSE, pch=10, main= paste(lAlpha[a]) ,ylab = "Root Mean-Squared Error", xlab= "Log Lambda",col="red") +abline(v=log(lambda_min), col="blue")
  }
  ind = which.min(lRMSE_min)
  Alpha_min = lAlpha[ind]
  cat("The alpha with minimum RMSE is: Alpha = ",Alpha_min)
  return(lRMSE_min)
}

