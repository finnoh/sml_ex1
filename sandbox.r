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
  sample_size = 1000000,
  reg_weights = c(2, 5,-0.7, 100, -23)
)

# data set
dfData <- simulate_fixed(data = NULL, sim_arguments) %>%
  generate_response(sim_arguments) %>% 
  select(-c(level1_id, random_effects, error, fixed_outcome, X.Intercept.))

# data objects
mX <- dfData %>% select(-y) %>% data.matrix() %>% scale()
vY <- dfData %>% select(y) %>% data.matrix()

# Real Data ---------------------------------------------------------------
vY <- supermarket1996 %>% select(GROCERY_sum) %>% as.matrix()
mX <- supermarket1996 %>% select(-c("STORE","CITY","GROCCOUP_sum","SHPINDX", "GROCERY_sum")) %>% as.matrix() %>% scale()

# Model -------------------------------------------------------------------

# hyperparameters
dLambda <- 10
dAlpha <- 1

# loop objects
iP <- ncol(mX)
iN <- nrow(mX)
ik <- 1
dLossChange <- 0
vBeta0 <- runif(ncol(mX))
dEps <- 10e-10
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

vBeta_MM <- vBeta0

print("Beta Estimate")
print(vBeta_MM)

# Comparison with glmnet --------------------------------------------------
# uses an intercept per default
model_glm <- glmnet(x = mX, y = vY, alpha = dAlpha, lambda = dLambda, 
                    intercept = FALSE, standardize = FALSE)
vBeta_glm <- model_glm %>% coef() %>% as.matrix()
vBeta_glm <- vBeta_glm[-1, ] 

# compare estimates with each other
dfCompareBeta <- cbind(vBeta_glm, vBeta_MM)

# in toy example we get the same results after running 1e-20 = eps, however
# for more complex data sets we are not able to increase the eps further due to
# numerical issues. We suspect that the algorithm used in the glmnet package 
# (generalized linear model via penalized maximum likelihood), converges faster and
# hence delivers more precise estimates than our implementation of the elastic net
# with the MM algorithm.
