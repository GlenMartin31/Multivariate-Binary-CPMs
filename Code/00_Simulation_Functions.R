# #######################################################################################################################

# Author of code: Glen P. Martin.

# This is code for a simulation study presented in a manuscript entitled: 
# Clinical Prediction Models to Predict the Risk of Multiple Binary Outcomes: a comparison of approaches
# Authors:
#   Glen P. Martin
#   Matthew Sperrin
#   Kym I.E. Snell
#   Iain Buchan
#   Richard Riley

# #######################################################################################################################


####-----------------------------------------------------------------------------------------
## Define a function to repeat the simulation across all iterations (for a given scenario)
####-----------------------------------------------------------------------------------------
simulation_nrun_fnc <- function(n_iter,
                                N, 
                                rho,
                                beta_1_true,
                                beta_2_true) {
  #Input: n_iter = number of iterations to repeat the simulation over,
  #       N = the number of observations to generate
  #       rho = strength of the association between the latent variables used to derive Y1 and Y2. When 
  #                   rho=0, Y1 and Y2 are independent
  #       beta_1_true & beta_2_true = vector (of length 3) specifying the true coefficients for intercept, X1 and X2
  
  ## Define an empty variable, which will be used to store the results across all iterations
  results <- NULL
  ## Repeat the simulation across iter number of iterations
  for (iter in 1:n_iter) {
    simulation_results <- simulation_singlerun_fnc(N = N, 
                                                   rho = rho,
                                                   beta_1_true = beta_1_true,
                                                   beta_2_true = beta_2_true)
    
    results <- results %>%
      bind_rows(as_tibble(simulation_results) %>% 
                  mutate("Iteration" = iter, .before = "Model"))
    
    rm(simulation_results)
  }
  ## Return results from the simulation across all iterations
  return(results)
}

####-----------------------------------------------------------------------------------------
## Define a function to run the processes within a single iteration
####-----------------------------------------------------------------------------------------
simulation_singlerun_fnc <- function(N, 
                                     rho,
                                     beta_1_true,
                                     beta_2_true) {
  #Input: N = the number of observations to generate
  #       rho = strength of the association between the latent variables used to derive Y1 and Y2. When 
  #                   rho=0, Y1 and Y2 are independent
  #       beta_1_true & beta_2_true = vector (of length 3) specifying the true coefficients for intercept, X1 and X2
  
  require(tidyverse)
  require(pROC)
  require(rjags)
  require(coda)
  require(pbivnorm)
  require(nnet)
  require(VGAM)
  require(glmnet)
  
  #Generate the development and validation data:
  CombinedData <- DataGenerating.fnc(N = N + 10000, #add 10000 additional observations for a validation set 
                                     beta_1_true = beta_1_true, 
                                     beta_2_true = beta_2_true, 
                                     rho = rho)
  #Randomly split into an IPD and validation cohort, since both are a random sample of the same underlying 
  #data generating process
  IPD <- CombinedData %>%
    sample_n(N, replace = FALSE) 
  Validation.Population <- CombinedData %>%
    filter(ID %in% IPD$ID == FALSE)
  
  ####-----------------
  ## Model Development
  ####-----------------
  
  ## Univariate Models: fit a model to each outcome independently
  Uni.m1 <- glm(Y1 ~ X1 + X2, data = IPD, family = binomial(link = "logit"))
  Uni.m2 <- glm(Y2 ~ X1 + X2, data = IPD, family = binomial(link = "logit"))
  #Predict the risk of each outcome in the validation cohort, based on the univariate models
  Univariate_PY1 <- predict(Uni.m1, newdata = Validation.Population, type = "response")
  Univariate_PY2 <- predict(Uni.m2, newdata = Validation.Population, type = "response")
  #Calculate the joint risks based on the marginal probabilities
  Validation.Population$Univariate_P11 <- Univariate_PY1 * Univariate_PY2
  Validation.Population$Univariate_P10 <- Univariate_PY1 * (1 - Univariate_PY2)
  Validation.Population$Univariate_P01 <- (1 - Univariate_PY1) * Univariate_PY2
  Validation.Population$Univariate_P00 <- (1 - Univariate_PY1) * (1 - Univariate_PY2)
  
  
  ## Stacked Regression (SR)
  #Obtain the linear predictors of each outcome from the univariate models: f_1 and f_2 in the paper
  IPD$Univariate_LP1 <- predict(Uni.m1, newdata = IPD, type = "link")
  IPD$Univariate_LP2 <- predict(Uni.m2, newdata = IPD, type = "link")
  #Use the linear predictors to fit a stacked regression model to each outcome separately
  MultivariateStackedRegression.m1 <- cv.glmnet(x = IPD %>%
                                                  select(Univariate_LP1, Univariate_LP2, X1, X2) %>%
                                                  data.matrix(),
                                                y = IPD$Y1,
                                                family = "binomial",
                                                alpha = 1)
  MultivariateStackedRegression.m2 <- cv.glmnet(x = IPD %>%
                                                  select(Univariate_LP1, Univariate_LP2, X1, X2) %>%
                                                  data.matrix(),
                                                y = IPD$Y2,
                                                family = "binomial",
                                                alpha = 1)
  #Predict the risk of each outcome in the validation cohort
  SR_PY1 <- predict(MultivariateStackedRegression.m1, 
                    newx = Validation.Population %>%
                      mutate(Univariate_LP1 = predict(Uni.m1, newdata = ., type = "link"),
                             Univariate_LP2 = predict(Uni.m2, newdata = ., type = "link")) %>%
                      select(Univariate_LP1, Univariate_LP2, X1, X2) %>%
                      data.matrix(), 
                    s = "lambda.min", 
                    type = "response")
  SR_PY2 <- predict(MultivariateStackedRegression.m2, 
                    newx = Validation.Population %>%
                      mutate(Univariate_LP1 = predict(Uni.m1, newdata = ., type = "link"),
                             Univariate_LP2 = predict(Uni.m2, newdata = ., type = "link")) %>%
                      select(Univariate_LP1, Univariate_LP2, X1, X2) %>%
                      data.matrix(), 
                    s = "lambda.min", 
                    type = "response")
  #Calculate the joint risks based on the marginal probabilities
  Validation.Population$SR_P11 <- as.vector(SR_PY1 * SR_PY2)
  Validation.Population$SR_P10 <- as.vector(SR_PY1 * (1 - SR_PY2))
  Validation.Population$SR_P01 <- as.vector((1 - SR_PY1) * SR_PY2)
  Validation.Population$SR_P00 <- as.vector((1 - SR_PY1) * (1 - SR_PY2))
  
  
  ## Probabilistic Classifier Chain (PCC)
  PCC <- Probabilistic.Classifier.Chain.fnc(DevelopmentData = IPD, 
                                            TestData = Validation.Population)
  #Extract the joint risks, which this method obtains directly
  Validation.Population$PCC_P11 <- PCC$P11
  Validation.Population$PCC_P10 <- PCC$P10
  Validation.Population$PCC_P01 <- PCC$P01
  Validation.Population$PCC_P00 <- (1 - (Validation.Population$PCC_P11 + 
                                           Validation.Population$PCC_P10 + 
                                           Validation.Population$PCC_P01))
  
  
  ## Multinomial Logistic Regression (MLR)
  MLR <- multinom(Y_Categories ~ X1 + X2, 
                  data = IPD)
  #Predict the risk of each outcome combination in the validation cohort
  MLR.predictions <- predict(MLR, newdata = Validation.Population, "probs")
  #Extract the joint risks, which this method obtains directly
  Validation.Population$MLR_P11 <- MLR.predictions[,"Y11"]
  Validation.Population$MLR_P10 <- MLR.predictions[,"Y10"]
  Validation.Population$MLR_P01 <- MLR.predictions[,"Y01"]
  Validation.Population$MLR_P00 <- MLR.predictions[,"Y00"]
  
  
  ## Multivariate Logistic Model (MLM)
  MLM <- Multivariate.Logistic.Reg.fnc(X = IPD %>% select(X1, X2), 
                                       Y1 = IPD$Y1, 
                                       Y2 = IPD$Y2)
  #Predict the risk of each outcome in the validation cohort using matrix multiplication
  MLM_PY1 <- (1 + exp(-as.numeric((Validation.Population %>%
                                     select(starts_with("X")) %>%
                                     data.matrix()) %*% MLM$beta1)))^(-1) #P(y1=1 | X,beta_1)
  MLM_PY2 <- (1 + exp(-as.numeric((Validation.Population %>%
                                     select(starts_with("X")) %>%
                                     data.matrix()) %*% MLM$beta2)))^(-1) #P(y2=1 | X,beta_2)
  #Check that the estimated value of rho12 is within the correct constraints. Note that the maximum value
  # that rho12 can take is determined by the marginal probabilities. As such, this might differ to that from
  # the IPD. As such, we here take the position that if rho12 is larger than the maximum limit within the 
  # validation set, then we assign it to be at the upper limit. This ensures that all joint probabilities lie 
  # in [0,1] and all sum to 1. The limit is calculated by the fact that P11, P10, P01 and P00 all must be >0. 
  # By the formula, rho12 can only make P10 or P01 <0 (since we subtract rho12*sqrt(...)). 
  # Working backwards from this, leads to the below if statement conditions.
  # NOTE: this is only relevant at larger value of rho 
  if (MLM$rho12 > min(min((MLM_PY1 * (1 - MLM_PY2)) / sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))),
                      min((MLM_PY2 * (1 - MLM_PY1)) / sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))))) {
    #Set at the upper limit: in most cases only very small difference between IPD estimated rho12 and this limit;
    MLM$rho12 <- min(min((MLM_PY1 * (1 - MLM_PY2)) / sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))),
                     min((MLM_PY2 * (1 - MLM_PY1)) / sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))))
  }
  #Calculate the joint risks, which this method obtains directly based on the marginal risks
  Validation.Population$MLM_P11 <- ((MLM_PY1 * MLM_PY2) + 
                                      (MLM$rho12 * sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))))
  Validation.Population$MLM_P10 <- ((MLM_PY1 * (1 - MLM_PY2)) - 
                                      (MLM$rho12 * sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))))
  Validation.Population$MLM_P01 <- ((MLM_PY2 * (1 - MLM_PY1)) - 
                                      (MLM$rho12 * sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))))
  Validation.Population$MLM_P00 <- (((1 - MLM_PY1) * (1 - MLM_PY2)) + 
                                      (MLM$rho12 * sqrt(MLM_PY1 * (1 - MLM_PY1) * MLM_PY2 * (1 - MLM_PY2))))
  
  
  ## Multivariate Probit Model (MPM)
  BayesianProbitMultivariateModel <- jags.model( file = here::here("Code", "JAGSProbitModel.txt"), 
                                                 data = list(
                                                   "n" = nrow(IPD),
                                                   "P" = 3,
                                                   "Y" = cbind(IPD$Y1, 
                                                               IPD$Y2),
                                                   "X" = cbind(1,
                                                               IPD$X1, 
                                                               IPD$X2),
                                                   'b_0' = rep(0,3),
                                                   'B_0' = diag(1, ncol = 3, nrow = 3)*0.1 #precision
                                                 ), 
                                                 inits = list("Z" = cbind(IPD$Y1, 
                                                                          IPD$Y2),
                                                              "rho" = 0),
                                                 n.chains = 1, 
                                                 n.adapt = 1000)
  #sample from the posterior distribution:
  samps <- coda.samples( BayesianProbitMultivariateModel, c('Beta','rho'), n.iter = 10000 )
  tidy.samps <- samps[[1]][5001:10000,] #set first 5000 samples as burn-in
  post.means <- colMeans(tidy.samps) #take the posterior mean
  #Predict the risk of each outcome in the validation cohort using posterior mean
  X.Beta.Y1 <- as.numeric(cbind(Validation.Population$X0, 
                                Validation.Population$X1,
                                Validation.Population$X2) %*% post.means[paste("Beta[", 1:3, ",1]", sep = "")])
  X.Beta.Y2 <- as.numeric(cbind(Validation.Population$X0, 
                                Validation.Population$X1,
                                Validation.Population$X2) %*% post.means[paste("Beta[", 1:3, ",2]", sep = "")])
  #Calculate the joint risks, which this method obtains directly based on the marginal risks and estimate of rho
  Validation.Population$MPM_P11 <- pbivnorm(x = cbind(X.Beta.Y1, X.Beta.Y2), 
                                            rho = post.means["rho"])
  Validation.Population$MPM_P10 <- pbivnorm(x = cbind(X.Beta.Y1, -X.Beta.Y2), 
                                            rho = -post.means["rho"])
  Validation.Population$MPM_P01 <- pbivnorm(x = cbind(-X.Beta.Y1, X.Beta.Y2), 
                                            rho = -post.means["rho"])
  Validation.Population$MPM_P00 <- pbivnorm(x = cbind(-X.Beta.Y1, -X.Beta.Y2),
                                            rho = post.means["rho"])
  
  
  ####-----------------
  ## Model Validation
  ####-----------------
  
  ## Extract relevant information from the validation cohort
  Predictions <- Validation.Population %>%
    select(Y1, Y2, Y_Categories,
           True_PY1, True_PY2, True_P11, True_P10, True_P01, True_P00,
           Univariate_P11, Univariate_P10, Univariate_P01, Univariate_P00,
           SR_P11, SR_P10, SR_P01, SR_P00,
           PCC_P11, PCC_P10, PCC_P01, PCC_P00,
           MLR_P11, MLR_P10, MLR_P01, MLR_P00,
           MLM_P11, MLM_P10, MLM_P01, MLM_P00,
           MPM_P11, MPM_P10, MPM_P01, MPM_P00) %>%
    mutate_at(vars(ends_with("_P11"),
                   ends_with("_P01"),
                   ends_with("_P10"),
                   ends_with("_P00")), 
              ##turn very small (practically 0 probs) to small number for entry into calibration models (in VGAM)
              ~ifelse(.<=1e-10, 1e-10, .)) %>% 
    mutate(## For each model, calculate predicted marginal risks based on the predicted joint risks:
           Univariate_Py1 = Univariate_P10 + Univariate_P11,
           Univariate_Py2 = Univariate_P01 + Univariate_P11,
           
           SR_Py1 = SR_P10 + SR_P11,
           SR_Py2 = SR_P01 + SR_P11,
           
           PCC_Py1 = PCC_P10 + PCC_P11,
           PCC_Py2 = PCC_P01 + PCC_P11,
           
           MLR_Py1 = MLR_P10 + MLR_P11,
           MLR_Py2 = MLR_P01 + MLR_P11,
           
           MLM_Py1 = MLM_P10 + MLM_P11,
           MLM_Py2 = MLM_P01 + MLM_P11,
           
           MPM_Py1 = MPM_P10 + MPM_P11,
           MPM_Py2 = MPM_P01 + MPM_P11)
  
  ##Loop across each model, and for each estimate predictive performance of the marginal and joint outcomes:
  PerformanceResults <- map_dfr(.x = list("Univariate" = Predictions$Univariate_Py1, 
                                          "SR" = Predictions$SR_Py1, 
                                          "PCC" = Predictions$PCC_Py1, 
                                          "MLR" = Predictions$MLR_Py1, 
                                          "MLM" = Predictions$MLM_Py1, 
                                          "MPM" = Predictions$MPM_Py1),
                                .f = Marginal.Performance.fnc,
                                ObservedOutcome = Predictions$Y1,
                                TrueRisk = Predictions$True_PY1,
                                OutcomeSubscript = "PY1",
                                .id = "Model") %>%
    left_join(map_dfr(.x = list("Univariate" = Predictions$Univariate_Py2, 
                                "SR" = Predictions$SR_Py2, 
                                "PCC" = Predictions$PCC_Py2, 
                                "MLR" = Predictions$MLR_Py2, 
                                "MLM" = Predictions$MLM_Py2, 
                                "MPM" = Predictions$MPM_Py2),
                      .f = Marginal.Performance.fnc,
                      ObservedOutcome = Predictions$Y2,
                      TrueRisk = Predictions$True_PY2,
                      OutcomeSubscript = "PY2",
                      .id = "Model"),
              by = "Model") %>%
    left_join(map_dfr(.x = list("Univariate" = "Univariate", 
                                "SR" = "SR", 
                                "PCC" = "PCC", 
                                "MLR" = "MLR", 
                                "MLM" = "MLM", 
                                "MPM" = "MPM"),
                      .f = Joint.Performance.fnc,
                      ObservedOutcome = Predictions$Y_Categories,
                      PredictionsDataFrame = Predictions,
                      TrueRisk_P11 = Predictions$True_P11,
                      TrueRisk_P10 = Predictions$True_P10,
                      TrueRisk_P01 = Predictions$True_P01,
                      TrueRisk_P00 = Predictions$True_P00,
                      .id = "Model"),
              by = "Model") %>%
    #Calculate relevant information to store per iteration:
    mutate("rho_observed" = cor(CombinedData$Y1, CombinedData$Y2),
           "Obs_Py1" = mean(Predictions$Y1), 
           "Obs_Py2" = mean(Predictions$Y2), 
           "Obs_P11" = sum(Predictions$Y_Categories == "Y11") / length(Predictions$Y_Categories), 
           "Obs_P10" = sum(Predictions$Y_Categories == "Y10") / length(Predictions$Y_Categories), 
           "Obs_P01" = sum(Predictions$Y_Categories == "Y01") / length(Predictions$Y_Categories)) %>%
    #Select order of variables to be consistent:
    select("Model",
           "Obs_Py1", "Obs_Py2", "Obs_P11", "Obs_P10", "Obs_P01", 
           "rho_observed",
           "CalInt_PY1", "CalSlope_PY1",
           "AUC_PY1", "MSE_PY1",
           "CalInt_PY2", "CalSlope_PY2",
           "AUC_PY2", "MSE_PY2",
           "CalInt_P11", "CalInt_P10", "CalInt_P01",
           "CalSlope_P11", "CalSlope_P10", "CalSlope_P01",
           "PDI_Overall", "PDI_P11", "PDI_P10", "PDI_P01",
           "AUC_P11", "AUC_P10", "AUC_P01",
           "MSE")
  #Return the performance results of this iteration:
  return(PerformanceResults)
}


####-----------------------------------------------------------------------------------------
## Data-generating function
####-----------------------------------------------------------------------------------------
DataGenerating.fnc <- function(N, beta_1_true, beta_2_true, rho) {
  #Input: N = the number of observations to generate
  #       beta_1_true = vector (of length 3) giving the 'true association between covariates and Y1
  #       beta_2_true = vector (of length 3) giving the 'true association between covariates and Y2
  #       rho = strength of the association between the latent variables used to derive Y1 and Y2. When rho=0, Y1 and Y2 are independent
  
  require(pbivnorm)
  
  #Define a function for the inverse logistic CDF
  Inverse.Logistic.CDF <- function(X) {
    return( log(X/(1 - X)) )
  }
  
  #Generate random vectors from standard multivariate normal distribution with inputted correlation, rho
  Z <- MASS::mvrnorm(n = N, mu = c(0,0), 
                     Sigma = toeplitz(c(1, rho)))
  Epsilon <- apply(pnorm(Z), 2, Inverse.Logistic.CDF) #apply the inverse probability transform
  
  ## Define the IPD and define the two outcomes based on the inverse probability transform
  IPD <- tibble("ID" = seq(from = 1, to = N, by = 1),
                "X0" = 1,
                "X1" = rnorm(N, 0, 1),
                "X2" = rnorm(N, 0, 1)) %>%
    mutate(LP1 = as.vector(((X0*beta_1_true[1]) + (X1*beta_1_true[2]) + (X2*beta_1_true[3]))),
           LP2 = as.vector(((X0*beta_2_true[1]) + (X1*beta_2_true[2]) + (X2*beta_2_true[3]))),
           
           Y1 = ifelse(Epsilon[,1] <= LP1, 1, 0),
           Y2 = ifelse(Epsilon[,2] <= LP2, 1, 0),
           
           Y_Categories = fct_relevel(factor(ifelse(Y1 == 0 & Y2 == 0, 
                                                    "Y00",
                                                    ifelse(Y1 == 1 & Y2 == 0,
                                                           "Y10",
                                                           ifelse(Y1 == 0 & Y2 == 1,
                                                                  "Y01",
                                                                  "Y11")))),
                                      c("Y00", "Y10", "Y01", "Y11")),
           
           True_PY1 = (1 + exp(-LP1))^(-1),
           True_PY2 = (1 + exp(-LP2))^(-1),
           True_P11 = pbivnorm(x = cbind(qnorm((1 + exp(-LP1))^(-1)), qnorm((1 + exp(-LP2))^(-1))), 
                               rho = rho),
           True_P10 = pbivnorm(x = cbind(qnorm((1 + exp(-LP1))^(-1)), -qnorm((1 + exp(-LP2))^(-1))), 
                               rho = -rho),
           True_P01 = pbivnorm(x = cbind(-qnorm((1 + exp(-LP1))^(-1)), qnorm((1 + exp(-LP2))^(-1))), 
                               rho = -rho),
           True_P00 = pbivnorm(x = cbind(-qnorm((1 + exp(-LP1))^(-1)), -qnorm((1 + exp(-LP2))^(-1))), 
                               rho = rho))
  
  # sum(IPD$Y1 == 1) / nrow(IPD); mean(IPD$True_PY1) #should match approx
  # sum(IPD$Y2 == 1) / nrow(IPD); mean(IPD$True_PY2) #should match approx
  # sum(IPD$Y_Categories == "Y11") / nrow(IPD); mean(IPD$True_P11) #should match approx
  # sum(IPD$Y_Categories == "Y10") / nrow(IPD); mean(IPD$True_P10) #should match approx
  # sum(IPD$Y_Categories == "Y01") / nrow(IPD); mean(IPD$True_P01) #should match approx
  # sum(IPD$Y_Categories == "Y00") / nrow(IPD); mean(IPD$True_P00) #should match approx
  
  return(IPD) #return the simulated IPD 
}


####-----------------------------------------------------------------------------------------
## Function to fit the Probabilistic Classifier Chain
####-----------------------------------------------------------------------------------------

Probabilistic.Classifier.Chain.fnc <- function(DevelopmentData, TestData) {
  #Inputs: DevelopmentData = data on which to derive the Probabilistic Classifier Chain CPM
  #        TestData = data on which to predict risk for each outcome
  
  ##Fit the models to the development data
  #In the first 'permutation' we start with estimating P(Y1) and then estimate P(Y2 | Y1)
  m1_Pi1 <- glm(Y1 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  m1_Pi2 <- glm(Y2 ~ X1 + X2 + Y1, data = DevelopmentData, family = binomial(link = "logit"))
  #In the second 'permutation' we start with estimating P(Y2) and then estimate P(Y1 | Y2)
  m2_Pi2 <- glm(Y2 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  m2_Pi1 <- glm(Y1 ~ X1 + X2 + Y2, data = DevelopmentData, family = binomial(link = "logit"))
  
  
  ##Now predict for each observation in the TestData
  #only extract the covaraite information (i.e. dont use the observed reponses to mimic real prediction):
  TestData <- TestData %>% select(X1, X2) 
  
  P11_m1 <- predict(m1_Pi1, newdata = TestData, type = "response") *
    predict(m1_Pi2, newdata = TestData %>% mutate(Y1 = 1), type = "response") #condition on Y1=1
  P10_m1 <- predict(m1_Pi1, newdata = TestData, type = "response") *
    (1 - predict(m1_Pi2, newdata = TestData %>% mutate(Y1 = 1), type = "response")) #condition on Y1=1
  P01_m1 <- (1 - predict(m1_Pi1, newdata = TestData, type = "response")) *
    predict(m1_Pi2, newdata = TestData %>% mutate(Y1 = 0), type = "response") #condition on Y1=0
  
  P11_m2 <- predict(m2_Pi2, newdata = TestData, type = "response") *
    predict(m2_Pi1, newdata = TestData %>% mutate(Y2 = 1), type = "response") #condition on Y2=1
  P10_m2 <- (1 - predict(m2_Pi2, newdata = TestData, type = "response")) *
    predict(m2_Pi1, newdata = TestData %>% mutate(Y2 = 0), type = "response") #condition on Y2=0
  P01_m2 <- predict(m2_Pi2, newdata = TestData, type = "response") *
    (1 - predict(m2_Pi1, newdata = TestData %>% mutate(Y2 = 1), type = "response")) #condition on Y2=1
  
  
  ## Apply the ensemble to estimate each joint risk
  P11 <- apply(cbind(P11_m1, P11_m2), 1, mean)
  P10 <- apply(cbind(P10_m1, P10_m2), 1, mean)
  P01 <- apply(cbind(P01_m1, P01_m2), 1, mean)
  
  ## Return the results
  return(list("P11" = P11,
              "P10" = P10,
              "P01" = P01))
}


####-----------------------------------------------------------------------------------------
## Numerical Optimisation of Multivariate Logistic Regression model
####-----------------------------------------------------------------------------------------
Multivariate.Logistic.Reg.fnc <- function(X, Y1, Y2) {
  
  ## Define log-likelihood function for the Multivariate Logistic Regression model
  #see Gumbel EJ. "Bivariate Logistic Distributions" (1961) for details of the likelihood
  loglike <- function(pars, X, Y1, Y2) {
    
    beta1 <- pars[1:ncol(X)]
    beta2 <- pars[(ncol(X) + 1):(2*ncol(X))]
    rho12 <- pars[(2*ncol(X)) + 1]
    
    F1 <- (1 + exp(-as.numeric(X %*% beta1)))^(-1) #P(y1=1 | X,beta_1)
    F2 <- (1 + exp(-as.numeric(X %*% beta2)))^(-1) #P(y2=1 | X,beta_2)
    S1 <- 1 - F1
    S2 <- 1 - F2
    
    P11 <- ((F1*F2) + (rho12*sqrt(F1*S1*F2*S2)))
    P10 <- ((F1*S2) - (rho12*sqrt(F1*S1*F2*S2)))
    P01 <- ((F2*S1) - (rho12*sqrt(F1*S1*F2*S2)))
    P00 <- ((S1*S2) + (rho12*sqrt(F1*S1*F2*S2))) 
    
    ll <- rep(NA, length(Y1))
    for (i in 1:length(Y1)) {
      #if any of the joint probabilities are out of range (due to constraints on rho12), then set log-likelihood as 
      # large value to ensure this combination of parameters (i.e. rho12) is not chosen as optimal
      if (P11[i] < 0 | P10[i] < 0 | P01[i] < 0 | P00[i] < 0 |
          P11[i] > 1 | P10[i] > 1 | P01[i] > 1 | P00[i] > 1) {
        ll[i] <- -Inf 
      }else {
        ll[i] <- (Y1[i]*Y2[i]*log(P11[i])) + 
          (Y1[i]*(1 - Y2[i])*log(P10[i])) + 
          ((1 - Y1[i])*Y2[i]*log(P01[i])) + 
          ((1 - Y1[i])*(1 - Y2[i])*log(P00[i]))
      }
    }
    return(-sum(ll)) #return minus log likelihood as optim minimises
  }
  
  X <- data.matrix(X) #Turn data into matrix ready for optim
  m1 <- glm(Y1 ~ X, family = binomial(link = "logit")) #obtain initial estimates of \beta_1
  m2 <- glm(Y2 ~ X, family = binomial(link = "logit")) #obtain initial estimates of \beta_2
  ## Pass relevant information to optim and use Nelder-Mead to estimate parameters
  maxlike <- optim(par = c(as.numeric(coef(m1)), 
                           as.numeric(coef(m2)), 
                           0),
                   fn = loglike,
                   method = "Nelder-Mead",
                   X = cbind(1, X), #append intercept to data
                   Y1 = Y1, Y2 = Y2,
                   control = list(maxit = 5000))
  
  ## Extract the relevant parameter estimates from optim
  beta1 <- maxlike$par[1:ncol(cbind(1, X))]
  beta2 <- maxlike$par[(ncol(cbind(1, X)) + 1):(2*ncol(cbind(1, X)))]
  rho12 <- maxlike$par[(2*ncol(cbind(1, X))) + 1]
  
  return(list("beta1" = beta1,
              "beta2" = beta2,
              "rho12" = rho12))
}


####-----------------------------------------------------------------------------------------
## Function to calculate the predictive performance of the models
####-----------------------------------------------------------------------------------------
##First, a function to calculate predictive performance of marginal (i.e. not joint) outcome risks
Marginal.Performance.fnc <- function(ObservedOutcome, PredictedRisk, TrueRisk, OutcomeSubscript) {
  
  #Note: This function should not be called directly: only from the main simulation function
  
  #Input: ObservedOutcome = a binary variable of observed outcomes
  #       PredictedRisk = a vector of predicted risks
  #       TrueRisk = data-generating risk for ObservedOutcome
  #       OutcomeSubscript = character denote marginal outcome being considered (Y1, Y2)
  
  require(pROC)
  
  LP <- log(PredictedRisk / (1 - PredictedRisk)) #map to log(prob outcome / prob of no outcome)
  
  #Estimate calibration intercept (i.e. calibration-in-the-large)
  CalInt.model <- glm(ObservedOutcome ~ offset(LP), family = binomial(link = "logit"))
  
  #Estimate calibration slope
  CalSlope.model <- glm(ObservedOutcome ~ LP, family = binomial(link = "logit"))
  
  #Discrimination
  AUC <- as.numeric(roc(response = ObservedOutcome, 
                        predictor = as.vector(PredictedRisk),
                        direction = "<",
                        levels = c(0,1))$auc)
  
  #Mean square error in predicted risks and 'true' risk
  MSE <- 1/length(ObservedOutcome) * (sum((PredictedRisk - TrueRisk)^2))
  
  ## Store performance results in a data.frame and return
  Results <- data.frame(as.numeric(coef(CalInt.model)),
                        as.numeric(coef(CalSlope.model)[2]),
                        AUC,
                        MSE)
  names(Results) <- c(paste("CalInt", OutcomeSubscript, sep = "_"),
                      paste("CalSlope", OutcomeSubscript, sep = "_"),
                      paste("AUC", OutcomeSubscript, sep = "_"),
                      paste("MSE", OutcomeSubscript, sep = "_"))
  return(Results)
}



##Second, a function to calculate predictive performance of joint outcome risks
Joint.Performance.fnc <- function(ObservedOutcome, PredictionsDataFrame,
                                  TrueRisk_P11, TrueRisk_P10, TrueRisk_P01, TrueRisk_P00, 
                                  ModelName) {
  
  #Note: This function should not be called directly: only from the main simulation function
  
  #Input: ObservedOutcome = a factor vector of observed outcomes, where each level is a 
  #               joint outcome combination (e.g. {Y1=1, Y2=1}, {Y1=1, Y2=0}, etc.). Must be
  #               factor levels of Y00, Y10, Y01, Y11 (in that order)
  #       PredictionsDataFrame = a data.frame of predicted risks. The columns for the
  #               predicted joint risks should be of the form [model name]_Pij where i/j={1,0} 
  #       TrueRisk_P11, TrueRisk_P10, TrueRisk_P01, TrueRisk_P00, = data-generating risks for each joint outcome
  #       ModelName = character variable denoting the model being tested
  
  require(VGAM)
  
  #Check factor levels (and ordering) of the ObservedOutcome
  if((all(levels(ObservedOutcome) == c("Y00", "Y10", "Y01", "Y11")) == FALSE)) {
    stop("Error: incorrect factor levels in ObservedOutcome")
  }
  
  #Extract the relevant predicted risks of the model under consideration (specified by ModelName)
  PredictionsDataFrame <- select(PredictionsDataFrame, 
                                 matches(paste(ModelName, "\\_P[0-9]+", sep = ""))) %>%
    #make sure column order matches order of factor levels for ObservedOutcome
    select(ends_with("P00"), ends_with("P10"), ends_with("P01"), ends_with("P11"))
  
  ##Calibration; follows the method and code at https://doi.org/10.1002/sim.6114
  #calibration-in-the-large
  CalInt.model <- coefficients(vgam(ObservedOutcome ~ 1, 
                                    offset = data.matrix(PredictionsDataFrame %>%
                                                           rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                                            make.names(names(PredictionsDataFrame))))) %>%
                                                           mutate(P10_Z = log(P10 / P00),
                                                                  P01_Z = log(P01 / P00),
                                                                  P11_Z = log(P11 / P00)) %>%
                                                           select(P10_Z, P01_Z, P11_Z)),
                                    family = multinomial(refLevel = "Y00")))
  CalInt_P10 <- as.numeric(CalInt.model[1])
  CalInt_P01 <- as.numeric(CalInt.model[2])
  CalInt_P11 <- as.numeric(CalInt.model[3])
  
  #calibration slope
  k <- length(levels(ObservedOutcome)) #number of outcome categories
  CalSlope.model <- coefficients(vgam(ObservedOutcome ~ P10_Z + P01_Z + P11_Z, 
                                      data = PredictionsDataFrame %>%
                                        rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                         make.names(names(PredictionsDataFrame))))) %>%
                                        mutate(P10_Z = log(P10 / P00),
                                               P01_Z = log(P01 / P00),
                                               P11_Z = log(P11 / P00)) %>%
                                        select(P10_Z, P01_Z, P11_Z),
                                      family = multinomial(refLevel = "Y00"),
                                      constraints = list("(Intercept)" = diag(1, 
                                                                              ncol = (k - 1), 
                                                                              nrow = (k - 1)),
                                                         "P10_Z" = rbind(1, 0, 0),
                                                         "P01_Z" = rbind(0, 1, 0),
                                                         "P11_Z" = rbind(0, 0, 1))))
  CalSlope_P10 <- as.numeric(CalSlope.model["P10_Z"])
  CalSlope_P01 <- as.numeric(CalSlope.model["P01_Z"])
  CalSlope_P11 <- as.numeric(CalSlope.model["P11_Z"])
  
  
  #Discrimination: the Polytomous Discrimination Index (PDI) for multivariate outcomes (random prediction has 
  #PDI of 1/k for k outcomes);
  #can transform the PDI onto the 0.5-1 scale (like AUC) using a PDI^(log(2, base=k)) transformation
  PDI <- pdi.modified(y = ObservedOutcome, 
                      d = PredictionsDataFrame, 
                      method = "prob") 
  PDI_Overall <- PDI$measure #Overall Polytomous Discrimination Index
  PDI_P11 <- PDI$table$VALUES[which(PDI$table$CATEGORIES=="Y11")] #Category-specific Polytomous Discrimination Index
  PDI_P10 <- PDI$table$VALUES[which(PDI$table$CATEGORIES=="Y10")] #Category-specific Polytomous Discrimination Index
  PDI_P01 <- PDI$table$VALUES[which(PDI$table$CATEGORIES=="Y01")] #Category-specific Polytomous Discrimination Index
  
  #Discrimination: also extract a one-versus-rest AUC:
  Y11 <- ifelse(ObservedOutcome == "Y11", 1, 0)
  Y10 <- ifelse(ObservedOutcome == "Y10", 1, 0)
  Y01 <- ifelse(ObservedOutcome == "Y01", 1, 0)
  AUC_P11 <- as.numeric(roc(response = Y11, 
                            predictor = as.vector(PredictionsDataFrame %>%
                                                    rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                                     make.names(names(PredictionsDataFrame))))) %>%
                                                    .$P11),
                            direction = "<",
                            levels = c(0,1))$auc)
  AUC_P10 <- as.numeric(roc(response = Y10, 
                            predictor = as.vector(PredictionsDataFrame %>%
                                                    rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                                     make.names(names(PredictionsDataFrame))))) %>%
                                                    .$P10),
                            direction = "<",
                            levels = c(0,1))$auc)
  AUC_P01 <- as.numeric(roc(response = Y01, 
                            predictor = as.vector(PredictionsDataFrame %>%
                                                    rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                                     make.names(names(PredictionsDataFrame))))) %>%
                                                    .$P01),
                            direction = "<",
                            levels = c(0,1))$auc)
  
  #Multivariate outcome mean square error
  MSE <- sum((as.numeric(PredictionsDataFrame %>%
                          rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                           make.names(names(PredictionsDataFrame))))) %>%
                          pull(P11)) - TrueRisk_P11)^2) +
    sum((as.numeric(PredictionsDataFrame %>%
                      rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                       make.names(names(PredictionsDataFrame))))) %>%
                      pull(P10)) - TrueRisk_P10)^2) + 
    sum((as.numeric(PredictionsDataFrame %>%
                      rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                       make.names(names(PredictionsDataFrame))))) %>%
                      pull(P01)) - TrueRisk_P01)^2) +
    sum((as.numeric(PredictionsDataFrame %>%
                      rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                       make.names(names(PredictionsDataFrame))))) %>%
                      pull(P00)) - TrueRisk_P00)^2)
  MSE <- (1/nrow(PredictionsDataFrame))*MSE
    
  
  ## Store performance results in a data.frame and return
  Results <- data.frame(CalInt_P11,
                        CalInt_P10,
                        CalInt_P01,
                        CalSlope_P11,
                        CalSlope_P10,
                        CalSlope_P01,
                        PDI_Overall,
                        PDI_P11,
                        PDI_P10,
                        PDI_P01,
                        AUC_P11,
                        AUC_P10,
                        AUC_P01,
                        MSE)
  return(Results)
}


####-----------------------------------------------------------------------------------------
## Function to calculate the Polytomous Discrimination Index (PDI) for multivariate outcomes 
# NOTE: this is adapted from the pdi() function of the mcca package, with relevant as.numeric()
#transformations added to allow calculation in large N (integer overflow errors for large N).
####-----------------------------------------------------------------------------------------
pdi.modified <- function (y, d, method = "multinom", ...) #see mcca package for input definitions
{
  require(mcca)
  y = factor(y)
  y_levels = levels(y)
  y = as.numeric(y)
  d = data.matrix(d)
  num = length(unique(y))
  pp = pm(y = y, d = d, method = method, ...)
  if (num == 3) {
    n1 = which(y == 1)
    n2 = which(y == 2)
    n3 = which(y == 3)
    pp1 = pp[n1, ]
    pp2 = pp[n2, ]
    pp3 = pp[n3, ]
    pdi1 = 0
    pdi2 = 0
    pdi3 = 0
    for (i in 1:length(n1)) {
      pdi1 = pdi1 + as.numeric(sum(pp1[i, 1] > pp2[, 1])) * 
        as.numeric(sum(pp1[i, 1] > pp3[, 1]))
    }
    for (i in 1:length(n2)) {
      pdi2 = pdi2 + as.numeric(sum(pp2[i, 2] > pp1[, 2])) * 
        as.numeric(sum(pp2[i, 2] > pp3[, 2]))
    }
    for (i in 1:length(n3)) {
      pdi3 = pdi3 + as.numeric(sum(pp3[i, 3] > pp1[, 3])) * 
        as.numeric(sum(pp3[i, 3] > pp2[, 3]))
    }
    pdi = (pdi1 + pdi2 + pdi3)/(3 * as.numeric(length(n1)) * 
                                  as.numeric(length(n2)) * 
                                  as.numeric(length(n3)))
    pdis = c(pdi1, pdi2, pdi3)/(as.numeric(length(n1)) * 
                                  as.numeric(length(n2)) * 
                                  as.numeric(length(n3)))
    df = data.frame(CATEGORIES = sapply(1:num, function(i) y_levels[i]), 
                    VALUES = pdis)
    result = list(call = match.call(), measure = pdi, table = df)
    class(result) = "mcca.pdi"
    return(result)
  }
  else if (num == 4) {
    n1 = which(y == 1)
    n2 = which(y == 2)
    n3 = which(y == 3)
    n4 = which(y == 4)
    pp1 = pp[n1, ]
    pp2 = pp[n2, ]
    pp3 = pp[n3, ]
    pp4 = pp[n4, ]
    pdi1 = 0
    pdi2 = 0
    pdi3 = 0
    pdi4 = 0
    for (i in 1:length(n1)) {
      pdi1 = pdi1 + as.numeric(sum(pp1[i, 1] > pp2[, 1])) * 
        as.numeric(sum(pp1[i,1] > pp3[, 1])) * as.numeric(sum(pp1[i, 1] > pp4[, 1]))
    }
    for (i in 1:length(n2)) {
      pdi2 = pdi2 + as.numeric(sum(pp2[i, 2] > pp1[, 2])) * 
        as.numeric(sum(pp2[i, 2] > pp3[, 2])) * as.numeric(sum(pp2[i, 2] > pp4[, 2]))
    }
    for (i in 1:length(n3)) {
      pdi3 = pdi3 + as.numeric(sum(pp3[i, 3] > pp1[, 3])) * 
        as.numeric(sum(pp3[i, 3] > pp2[, 3])) * as.numeric(sum(pp3[i, 3] > pp4[, 3]))
    }
    for (i in 1:length(n4)) {
      pdi4 = pdi4 + as.numeric(sum(pp4[i, 4] > pp1[, 4])) * 
        as.numeric(sum(pp4[i, 4] > pp2[, 4])) * as.numeric(sum(pp4[i, 4] > pp3[, 4]))
    }
    
    pdi = (pdi1 + pdi2 + pdi3 + pdi4)/(4 * as.numeric(length(n1)) * 
                                         as.numeric(length(n2)) * 
                                         as.numeric(length(n3)) * 
                                         as.numeric(length(n4)))
    pdis = c(pdi1, pdi2, pdi3, pdi4)/(as.numeric(length(n1)) * 
                                        as.numeric(length(n2)) * 
                                        as.numeric(length(n3)) * 
                                        as.numeric(length(n4)))
    df = data.frame(CATEGORIES = sapply(1:num, function(i) y_levels[i]), 
                    VALUES = pdis)
    result = list(call = match.call(), measure = pdi, table = df)
    class(result) = "mcca.pdi"
    return(result)
  }
  else if (num == 2) {
    n1 = which(y == 1)
    n2 = which(y == 2)
    pp1 = pp[n1, ]
    pp2 = pp[n2, ]
    pdi1 = 0
    pdi2 = 0
    for (i in 1:length(n1)) {
      pdi1 = pdi1 + as.numeric(sum(pp1[i, 1] > pp2[, 1])) + 
        as.numeric(sum(pp1[i, 1] == pp2[, 1]))
    }
    for (i in 1:length(n2)) {
      pdi2 = pdi2 + as.numeric(sum(pp2[i, 2] > pp1[, 2]))
    }
    pdi = (pdi1 + pdi2)/(2 * as.numeric(length(n1)) * 
                           as.numeric(length(n2)))
    pdis = c(pdi1, pdi2)/(as.numeric(length(n1)) * 
                            as.numeric(length(n2)))
    df = data.frame(CATEGORIES = sapply(1:num, function(i) y_levels[i]), 
                    VALUES = pdis)
    result = list(call = match.call(), measure = pdi, table = df)
    class(result) = "mcca.pdi"
    return(result)
  }
}
