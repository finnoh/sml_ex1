# Setup -------------------------------------------------------------------
# load packages
if(!require(pacman)){install.packages("pacman")}
p_load(tidyverse, simglm, rlist, latex2exp, glmnet)

# set seed
set.seed(321)

# import custom functions
source("./dev/toolbox.r")

# load data
load("supermarket1996.RData")

# writing format objects
mytheme <- theme_bw() + theme(legend.position = "bottom")

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

# parameters
dLambda <- 10
dAlpha <- 0.5
dEps <- 10e-10

vBeta_MM <- ElasticNetMM(mX, vY, dEps, dAlpha, dLambda)

print("Beta Estimate")
print(vBeta_MM)

# Comparison with glmnet --------------------------------------------------
# uses an intercept per default
model_glm <- glmnet(x = mX, y = vY, alpha = dAlpha, lambda = dLambda, 
                    intercept = FALSE, standardize = FALSE)
vBeta_glm <- model_glm %>% coef() %>% as.matrix()
vBeta_glm <- vBeta_glm[-1, ] 
vNameCoef <- names(vBeta_glm)

dfCompareBetaTable <- CompareEstimates(vBeta_glm, vBeta_MM)

plot_coef_rmse <- dfCompareBetaTable %>% mutate(MAPE = MAPE(GLMNET, MM)) %>% 
  ggplot(aes(x = Predictor, y = MAPE)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "MAPE") +
  scale_x_discrete(breaks = vNameCoef, labels = abbreviate) +
  scale_y_continuous(labels = scales::percent) +
  mytheme +
  theme(axis.text.x = element_text(angle = 45, size = 7, vjust = 0.5))
  
# in toy example we get the same results after running 1e-20 = eps, however
# for more complex data sets we are not able to increase the eps further due to
# numerical issues. We suspect that the algorithm used in the glmnet package 
# (generalized linear model via penalized maximum likelihood), converges faster and
# hence delivers more precise estimates than our implementation of the elastic net
# with the MM algorithm.


# Development for Epsilon -------------------------------------------------
iEpsStart <- 1
iEpsEnd <- -100
iEpsStep <- -1

# epsilon steps
vEps <- 10^seq(iEpsStart, iEpsEnd, iEpsStep)
lBeta_MM <- list()
lCompare <- list()

# estiamte model for each epsilon and compare to glmnet
for (i in seq_along(vEps)) {
  vBeta_MM <- ElasticNetMM(mX, vY, vEps[i], dAlpha, dLambda)
  lCompare[[i]] <- CompareEstimates(vBeta_glm, vBeta_MM)
  lCompare[[i]]$Epsilon <- vEps[i]
}

model_glm <- glmnet(x = mX, y = vY, alpha = dAlpha, lambda = dLambda, 
                    intercept = FALSE, standardize = FALSE)
vBeta_glm <- model_glm %>% coef() %>% as.matrix()
vBeta_glm <- vBeta_glm[-1, ] 
vNameCoef <- names(vBeta_glm)

# transform to dataframe
dfBetaCompareEps <- list.stack(lCompare)
dfBetaCompareEps <- dfBetaCompareEps %>% 
  mutate(MAPE = MAPE(GLMNET, MM),
         MAE = MAE(GLMNET, MM)) %>% 
  as_tibble()

plot_MAPE_eps <- dfBetaCompareEps %>% 
  group_by(Epsilon) %>%
  summarise(MAPE = mean(MAPE),
            MAE = median(MAE)) %>%
  ggplot(aes(x = log10(Epsilon), y = MAE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous(trans = "reverse", breaks = seq(iEpsStart - 1, iEpsEnd, -10)) +
  labs(x = TeX("$log_{10} (\\epsilon)$")) +
  mytheme
plot_MAPE_eps

plot_MAPE_eps_pred <- dfBetaCompareEps %>% 
  ggplot(aes(x = log10(Epsilon), y = MAPE, group = Predictor)) +
  facet_wrap(~Predictor, scales = "free_y", ncol = 9) +
  geom_smooth() +
  scale_x_continuous(trans = "reverse", breaks = seq(iEpsStart, iEpsEnd, iEpsStep)) +
  labs(x = TeX("$log_{10} (\\epsilon)$")) +
  theme(axis.text.y = element_blank()) +
  mytheme
plot_MAPE_eps_pred
