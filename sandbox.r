# Setup -------------------------------------------------------------------
# load packages
if(!require(pacman)){install.packages("pacman")}
p_load(tidyverse, simglm)

# set seed
set.seed(321)

# import custom functions
source("./dev/toolbox.r")

# load data
load("supermarket1996.RData")


# Synth.  Data ------------------------------------------------------------
sim_arguments <- list(
  formula = y ~ 1 + x1 + x2 + x3 + x4,
  fixed = list(
    x1 = list(var_type = 'continuous', mean = 180, sd = 50),
    x2 = list(var_type = 'continuous', mean = 75, sd = 20),
    x3 = list(var_type = 'continuous', mean = -23, sd = 4),
    x4 = list(var_type = 'continuous', mean = 1, sd = 20)),
  sample_size = 1000,
  reg_weights = c(2, 5,-0.7, 100, -23)
)

# data set
dfData <- simulate_fixed(data = NULL, sim_arguments) %>%
  generate_response(sim_arguments) %>% 
  select(-c(level1_id, random_effects, error, fixed_outcome, X.Intercept.))

# data objects
mX <- dfData %>% select(-y) %>% data.matrix()
vY <- dfData %>% select(y) %>% data.matrix()


# Real Data ---------------------------------------------------------------
vY <- supermarket1996 %>% select(GROCERY_sum) %>% as.matrix()
mX <- supermarket1996 %>% select(-c("STORE","CITY","GROCCOUP_sum","SHPINDX")) %>% as.matrix() %>% scale()

# Model -------------------------------------------------------------------

# hyperparameters
dLambda <- 0.5
dAlpha <- 0.5

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