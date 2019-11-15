####-----------------------------------------------------------------------------------------
## Define a function to call the simulation across multiple iterations
####-----------------------------------------------------------------------------------------
simulationcall.fnc <- function(iter = 100,
                               N, 
                               rho.vals = 0,
                               filename,
                               startingseed) {
  #Input: iter = number of iterations to repeat the simulation over,
  #       N = the number of observations to generate
  #       rho.vals = strength of the association between the latent variables used to derive Y1 and Y2. When 
  #                   rho=0, Y1 and Y2 are independent
  #       filename = character variable giving the name of the text file to store the results
  #       startingseed = starting seed for the random number generator
  
  library(pROC)
  library(rjags)
  library(coda)
  library(pbivnorm)
  library(nnet)
  library(VGAM)
  library(tidyverse)
  
  set.seed(startingseed)
  
  #Create text file to store all the performance results
  if (file.exists(here::here("Data", paste(filename, ".txt", sep = ""))) == TRUE) {
    #stops from accidentally overwritting a results file
    stop("filename already exists - delete or move file to continue")
  }else{
    OutputNames <- c("Iteration", 
                     "Model",
                     "Obs_Py1", "Obs_Py2", "Obs_P11", "Obs_P10", "Obs_P01", 
                     "rho_latent", 
                     "rho_Y1Y2",
                     "CalInt_PY1", "CalInt_PY2",
                     "CalSlope_PY1", "CalSlope_PY2",
                     "AUC_PY1", "AUC_PY2",
                     "MSE_PY1", "MSE_PY2",
                     "CalInt_P11", "CalInt_P10", "CalInt_P01",
                     "CalSlope_P11", "CalSlope_P10", "CalSlope_P01",
                     "AUC_P11", "AUC_P10", "AUC_P01",
                     "MSE_P11", "MSE_P10", "MSE_P01")
    write(OutputNames, 
          here::here("Data", paste(filename, ".txt", sep = "")), 
          append = FALSE, 
          sep = "|", 
          ncolumns = length(OutputNames)) 
  }
  
  
  #Set 'true' coefficients
  beta_1_true <- c(-1.0, log(2), log(1.00))
  beta_2_true <- c(-1.5, log(1.00), log(3))
  
  #Repeat the simulation over all values of rho
  for (rho in rho.vals) {
    
    #define a progress bar to monitor progress across iterations of the simulation
    pb <- winProgressBar(title = "Iteration Progress bar",
                         label = paste("Iterations 0% Completed"),
                         min = 0, max = 100, initial = 0)
    for (iteration in 1:iter) {
      
      CombinedData <- DataGenerating.fnc(N = sum(N, 10000), #add 10000 additional observations for a validation set 
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
      
      
      ## Hierarchical Related Regression (HRR)
      HRR <- Hierarchical.Related.CPM.fnc(DevelopmentData = IPD, 
                                          TestData = Validation.Population)
      #Calculate the joint risks based on the marginal probabilities
      Validation.Population$HRR_P11 <- HRR$Pi_Y1 * HRR$Pi_Y2
      Validation.Population$HRR_P10 <- HRR$Pi_Y1 * (1 - HRR$Pi_Y2)
      Validation.Population$HRR_P01 <- (1 - HRR$Pi_Y1) * HRR$Pi_Y2
      Validation.Population$HRR_P00 <- (1 - HRR$Pi_Y1) * (1 - HRR$Pi_Y2)
      
      
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
                      data = IPD %>%
                        mutate(Y_Categories = fct_relevel(factor(ifelse(Y1 == 0 & Y2 == 0, 
                                                                        "(0,0)",
                                                                        ifelse(Y1 == 1 & Y2 == 0,
                                                                               "(1,0)",
                                                                               ifelse(Y1 == 0 & Y2 == 1,
                                                                                      "(0,1)",
                                                                                      "(1,1)")))),
                                                          c("(0,0)", "(1,0)", "(0,1)", "(1,1)"))))
      #Predict the risk of each outcome combination in the validation cohort
      MLR.predictions <- predict(MLR, newdata = Validation.Population, "probs")
      #Extract the joint risks, which this method obtains directly
      Validation.Population$MLR_P11 <- MLR.predictions[,"(1,1)"]
      Validation.Population$MLR_P10 <- MLR.predictions[,"(1,0)"]
      Validation.Population$MLR_P01 <- MLR.predictions[,"(0,1)"]
      Validation.Population$MLR_P00 <- MLR.predictions[,"(0,0)"]
      
      
      ## Multivariate Logistic Model (MLM)
      MLM <- suppressWarnings(Multivariate.Logistic.Reg.fnc(X = IPD %>% select(X1, X2), 
                                                            Y1 = IPD$Y1, 
                                                            Y2 = IPD$Y2))
      #Predict the risk of each outcome in the validation cohort using matrix multiplication
      MLM_PY1 <- (1 + exp(-as.numeric((Validation.Population %>%
                                         select(starts_with("X")) %>%
                                         data.matrix()) %*% MLM$beta1)))^(-1) #P(y1=1 | X,beta_1)
      MLM_PY2 <- (1 + exp(-as.numeric((Validation.Population %>%
                                         select(starts_with("X")) %>%
                                         data.matrix()) %*% MLM$beta2)))^(-1) #P(y2=1 | X,beta_2)
      #Calculate the joint risks, which this method obtains directly based on the marginal risks
      Validation.Population$MLM_P11 <- ((MLM_PY1*MLM_PY2) + 
                                          (MLM$rho12*sqrt(MLM_PY1*(1 - MLM_PY1)*MLM_PY2*(1 - MLM_PY2))))
      Validation.Population$MLM_P10 <- ((MLM_PY1*(1 - MLM_PY2)) - 
                                          (MLM$rho12*sqrt(MLM_PY1*(1 - MLM_PY1)*MLM_PY2*(1 - MLM_PY2))))
      Validation.Population$MLM_P01 <- ((MLM_PY2*(1 - MLM_PY1)) - 
                                          (MLM$rho12*sqrt(MLM_PY1*(1 - MLM_PY1)*MLM_PY2*(1 - MLM_PY2))))
      Validation.Population$MLM_P00 <- (1 - (Validation.Population$MLM_P11 + 
                                               Validation.Population$MLM_P10 + 
                                               Validation.Population$MLM_P01))
      
      
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
                                                       'B_0' = diag(3)*0.1 #precision
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
      
      
      ## Multivariate Stacked Regression (SR)
      #Obtain the linear predictors of each outcome from the univariate models: f_1 and f_2 in the paper
      IPD$Univariate_LP1 <- predict(Uni.m1, newdata = IPD, type = "link")
      IPD$Univariate_LP2 <- predict(Uni.m2, newdata = IPD, type = "link")
      #Use the linear predictors to fit a stacked regression model to each outcome seperately
      MultivariateStackedRegression.m1 <- glm(Y1 ~ Univariate_LP1 + Univariate_LP2, 
                                              data = IPD, 
                                              family = binomial(link = "logit"))
      MultivariateStackedRegression.m2 <- glm(Y2 ~ Univariate_LP1 + Univariate_LP2, 
                                              data = IPD, 
                                              family = binomial(link = "logit"))
      #Predict the risk of each outcome in the validation cohort
      SR_PY1 <- predict(MultivariateStackedRegression.m1, 
                        newdata = Validation.Population %>%
                          mutate(Univariate_LP1 = log(Univariate_PY1/(1 - Univariate_PY1)),
                                 Univariate_LP2 = log(Univariate_PY2/(1 - Univariate_PY2))), 
                        type = "response")
      SR_PY2 <- predict(MultivariateStackedRegression.m2, 
                        newdata = Validation.Population %>%
                          mutate(Univariate_LP1 = log(Univariate_PY1/(1 - Univariate_PY1)),
                                 Univariate_LP2 = log(Univariate_PY2/(1 - Univariate_PY2))), 
                        type = "response")
      #Calculate the joint risks based on the marginal probabilities
      Validation.Population$SR_P11 <- SR_PY1 * SR_PY2
      Validation.Population$SR_P10 <- SR_PY1 * (1 - SR_PY2)
      Validation.Population$SR_P01 <- (1 - SR_PY1) * SR_PY2
      Validation.Population$SR_P00 <- (1 - SR_PY1) * (1 - SR_PY2)
      
      
      ####-----------------
      ## Model Validation
      ####-----------------
      
      ## Extract relevant information from the validation cohort
      Predictions <- Validation.Population %>%
        select(Y1, Y2,
               True_PY1, True_PY2, True_P11, True_P10, True_P01,
               Univariate_P11, Univariate_P10, Univariate_P01, Univariate_P00,
               HRR_P11, HRR_P10, HRR_P01, HRR_P00,
               PCC_P11, PCC_P10, PCC_P01, PCC_P00,
               MLR_P11, MLR_P10, MLR_P01, MLR_P00,
               MLM_P11, MLM_P10, MLM_P01, MLM_P00,
               MPM_P11, MPM_P10, MPM_P01, MPM_P00,
               SR_P11, SR_P10, SR_P01, SR_P00) %>%
        mutate_at(vars(ends_with("_P11"),
                       ends_with("_P01"),
                       ends_with("_P10"),
                       ends_with("_P00")), 
                  ##turn very small (practically 0 probs) to small number for entry into calibration models (in VGAM)
                  ~ifelse(.<=1e-10, 1e-10, .)) %>% 
        mutate(Y_Categories = fct_relevel(factor(ifelse(Y1 == 0 & Y2 == 0, 
                                                        "(0,0)",
                                                        ifelse(Y1 == 1 & Y2 == 0,
                                                               "(1,0)",
                                                               ifelse(Y1 == 0 & Y2 == 1,
                                                                      "(0,1)",
                                                                      "(1,1)")))),
                                          c("(0,0)", "(1,0)", "(0,1)", "(1,1)")), #define the observed joint outcome events
               
               ## For each model, calculate predicted marginal risks based on the predicted joint risks:
               Univariate_Py1 = Univariate_P10 + Univariate_P11,
               Univariate_Py2 = Univariate_P01 + Univariate_P11,
               
               HRR_Py1 = HRR_P10 + HRR_P11,
               HRR_Py2 = HRR_P01 + HRR_P11,
               
               PCC_Py1 = PCC_P10 + PCC_P11,
               PCC_Py2 = PCC_P01 + PCC_P11,
               
               MLR_Py1 = MLR_P10 + MLR_P11,
               MLR_Py2 = MLR_P01 + MLR_P11,
               
               MLM_Py1 = MLM_P10 + MLM_P11,
               MLM_Py2 = MLM_P01 + MLM_P11,
               
               MPM_Py1 = MPM_P10 + MPM_P11,
               MPM_Py2 = MPM_P01 + MPM_P11,
               
               SR_Py1 = SR_P10 + SR_P11,
               SR_Py2 = SR_P01 + SR_P11) %>%
        select(Y1, Y2,
               Y_Categories,
               ends_with("_Py1"), ends_with("_Py2"),
               ends_with("_P11"), ends_with("_P10"), ends_with("_P01"), ends_with("_P00"))
      
      ##Loop across each model, and for each estimate predictive performance of the marginal and joint outcomes:
      Results <- NULL
      for (ModelName in c("Univariate", "HRR", "PCC", "MLR", "MLM", "MPM", "SR")) {
        Results <- Results %>%
          bind_rows(Marginal.Performance.fnc(ObservedOutcome = Predictions$Y1, 
                                             PredictedRisk = eval(parse(text = paste("Predictions$",
                                                                                     ModelName,
                                                                                     "_Py1", sep = ""))),
                                             TrueRisk = Predictions$True_PY1,
                                             OutcomeSubscript = "PY1") %>%
                      mutate("Model" = ModelName) %>%
                      left_join(Marginal.Performance.fnc(ObservedOutcome = Predictions$Y2, 
                                                         PredictedRisk = eval(parse(text = paste("Predictions$",
                                                                                                 ModelName,
                                                                                                 "_Py2", sep = ""))),
                                                         TrueRisk = Predictions$True_PY2,
                                                         OutcomeSubscript = "PY2") %>%
                                  mutate("Model" = ModelName),
                                by = "Model") %>%
                      left_join(Joint.Performance.fnc(ObservedOutcome = Predictions$Y_Categories,
                                                      PredictionsDataFrame = Predictions,
                                                      TrueRisk_P11 = Predictions$True_P11,
                                                      TrueRisk_P10 = Predictions$True_P10,
                                                      TrueRisk_P01 = Predictions$True_P01,
                                                      ModelName = ModelName) %>%
                                  mutate("Model" = ModelName),
                                by = "Model"))
      }
      
      ####-----------------
      ## Store Performance Results
      ####-----------------
      Results <- Results %>%
        mutate("Iteration" = iteration,
               "rho_latent" = rho,
               "rho_Y1Y2" = cor(CombinedData$Y1, CombinedData$Y2),
               
               "Obs_Py1" = mean(Predictions$Y1), 
               "Obs_Py2" = mean(Predictions$Y2), 
               "Obs_P11" = sum(Predictions$Y_Categories == "(1,1)") / length(Predictions$Y_Categories), 
               "Obs_P10" = sum(Predictions$Y_Categories == "(1,0)") / length(Predictions$Y_Categories), 
               "Obs_P01" = sum(Predictions$Y_Categories == "(0,1)") / length(Predictions$Y_Categories)) %>%
        select("Iteration", 
               "Model",
               "Obs_Py1", "Obs_Py2", "Obs_P11", "Obs_P10", "Obs_P01", 
               "rho_latent", 
               "rho_Y1Y2",
               "CalInt_PY1", "CalInt_PY2",
               "CalSlope_PY1", "CalSlope_PY2",
               "AUC_PY1", "AUC_PY2",
               "MSE_PY1", "MSE_PY2",
               "CalInt_P11", "CalInt_P10", "CalInt_P01",
               "CalSlope_P11", "CalSlope_P10", "CalSlope_P01",
               "AUC_P11", "AUC_P10", "AUC_P01",
               "MSE_P11", "MSE_P10", "MSE_P01")
     
      ##Write the results to the text file before starting the next iteration
      write.table(Results,
                  here::here("Data", paste(filename, ".txt", sep = "")),
                  append = TRUE,
                  sep = "|",
                  row.names = FALSE,
                  col.names = FALSE)
      rm(Results) #clear results ready for the next iteration
      
      ##Update the progress bar
      info <- sprintf(paste("Iterations %d%% Completed"), round((iteration/(iter)*100)))
      setWinProgressBar(pb, iteration/((iter))*100, label = info)
    }
    close(pb) #close the progress bar
    
    ##Print a message to the concole to update on simulation progress
    print(paste("Simulation for rho equal to", rho, "complete", sep = " "))
  }
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
  
  ## Generate independent covariates
  X <- cbind(X0 = 1,
             X1 = rnorm(N, 0, 1), 
             X2 = rnorm(N, 0, 1))
  
  #Generate random vectors from standard multivariate normal distribution with inputted correlation, rho
  Z <- MASS::mvrnorm(n = N, mu = c(0,0), 
                     Sigma = toeplitz(c(1, rho)))
  Epsilon <- apply(pnorm(Z), 2, Inverse.Logistic.CDF) #apply the inverse probability transform
  
  ## Define the IPD and define the two outcomes based on the inverse probability transform
  IPD <- data.frame("ID" = seq(from = 1, to = N, by = 1),
                    X) %>%
    mutate(LP1 = as.vector(X %*% beta_1_true),
           LP2 = as.vector(X %*% beta_2_true),
           
           Y1 = ifelse(Epsilon[,1] <= LP1, 1, 0),
           Y2 = ifelse(Epsilon[,2] <= LP2, 1, 0),
           
           True_PY1 = (1 + exp(-LP1))^(-1),
           True_PY2 = (1 + exp(-LP2))^(-1),
           True_P11 = pbivnorm(x = cbind(qnorm((1 + exp(-LP1))^(-1)), qnorm((1 + exp(-LP2))^(-1))), 
                               rho = rho),
           True_P10 = pbivnorm(x = cbind(qnorm((1 + exp(-LP1))^(-1)), -qnorm((1 + exp(-LP2))^(-1))), 
                               rho = -rho),
           True_P01 = pbivnorm(x = cbind(-qnorm((1 + exp(-LP1))^(-1)), qnorm((1 + exp(-LP2))^(-1))), 
                               rho = -rho))
  
  # sum(IPD$Y1 == 1) / nrow(IPD); mean(IPD$True_PY1) #should match approx
  # sum(IPD$Y2 == 1) / nrow(IPD); mean(IPD$True_PY2) #should match approx
  # sum(IPD$Y1 == 1 & IPD$Y2 == 1) / nrow(IPD); mean(IPD$True_P11) #should match approx
  # sum(IPD$Y1 == 1 & IPD$Y2 == 0) / nrow(IPD); mean(IPD$True_P10) #should match approx
  # sum(IPD$Y1 == 0 & IPD$Y2 == 1) / nrow(IPD); mean(IPD$True_P01) #should match approx
  
  return(IPD) #return the simulated IPD to the main analysis function
}


####-----------------------------------------------------------------------------------------
## Function to fit the Hierarchical Related CPM
####-----------------------------------------------------------------------------------------

Hierarchical.Related.CPM.fnc <- function(DevelopmentData, TestData) {
  #Inputs: DevelopmentData = data on which to derive the Hierarchical Related CPM
  #        TestData = data on which to predict risk for each outcome (marginal), 
  
  ## In the first 'permutation' we start with estimating Y1
  m1_Pi1 <- glm(Y1 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  DevelopmentData$Pi1_m1 <- predict(m1_Pi1, 
                                    newdata = DevelopmentData, 
                                    type = "response")
  #...then use the predicted risk for Y1 to estimate Y2
  m1_Pi2 <- glm(Y2 ~ X1 + X2 + Pi1_m1, data = DevelopmentData, family = binomial(link = "logit"))
  
  
  ## In the second 'permutation' we start with estimating Y2
  m2_Pi2 <- glm(Y2 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  DevelopmentData$Pi2_m2 <- predict(m2_Pi2, 
                                    newdata = DevelopmentData, 
                                    type = "response")
  #...then use the predicted risk for Y2 to estimate Y1
  m2_Pi1 <- glm(Y1 ~ X1 + X2 + Pi2_m2, data = DevelopmentData, family = binomial(link = "logit"))
  
  
  ## Now use each permutation model to predict risk of Y1 or Y2 for each observation in the TestData
  TestData$Pi1_m1 <- predict(m1_Pi1, newdata = TestData, type = "response")
  TestData$Pi2_m1 <- predict(m1_Pi2, newdata = TestData, type = "response")
  TestData$Pi2_m2 <- predict(m2_Pi2, newdata = TestData, type = "response")
  TestData$Pi1_m2 <- predict(m2_Pi1, newdata = TestData, type = "response")
  
  ## Take an ensemble of the predicted risks from each permutation model to obtaine the
  # overall (marginal) predicted risk 
  Pi_Y1 <- apply(cbind(TestData$Pi1_m1, TestData$Pi1_m2), 1, mean)
  Pi_Y2 <- apply(cbind(TestData$Pi2_m1, TestData$Pi2_m2), 1, mean)
  
  ## Return the results
  return(list("m1_Pi1" = m1_Pi1,
              "m1_Pi2" = m1_Pi2,
              "m2_Pi1" = m2_Pi1,
              "m2_Pi2" = m2_Pi2,
              "Pi_Y1" = Pi_Y1,
              "Pi_Y2" = Pi_Y2))
}


####-----------------------------------------------------------------------------------------
## Function to fit the Probabilistic Classifier Chain
####-----------------------------------------------------------------------------------------

Probabilistic.Classifier.Chain.fnc <- function(DevelopmentData, TestData) {
  #Inputs: DevelopmentData = data on which to derive the Probabilistic Classifier Chain CPM
  #        TestData = data on which to predict risk for each outcome (marginal)
  
  ##Fit the models to the full development data to arrive at parameter estimates for future prediction
  
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
    
    ## Correct any rounding issues producing small negative numbers: turn to very small positive
    if (any(P11 < 0)) {
      P11[which(P11 < 0)] <- 10^-10
      
    }else if (any(P10 < 0)) {
      P10[which(P10 < 0)] <- 10^-10
      
    }else if (any(P01 < 0)) {
      P01[which(P01 < 0)] <- 10^-10
      
    }else if (any(P00 < 0)) {
      P00[which(P00 < 0)] <- 10^-10
      
    }
    
    ll <- (Y1*Y2*log(P11)) + (Y1*(1 - Y2)*log(P10)) + ((1 - Y1)*Y2*log(P01)) + ((1 - Y1)*(1 - Y2)*log(P00))
    
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
                   X = cbind(1, X), #apend intercept to data
                   Y1 = Y1, Y2 = Y2)
  
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
  require(pROC)
  
  LP <- log(PredictedRisk / (1 - PredictedRisk)) #convert predicted risks onto linear predictor scale
  
  #Estimate calibration intercept (i.e. calibration-in-the-large)
  CalInt.model <- glm(ObservedOutcome ~ offset(LP), family = binomial(link = "logit"))
  
  #Estimate calibration slope
  CalSlope.model <- glm(ObservedOutcome ~ LP, family = binomial(link = "logit"))
  
  #Dicrimination
  AUC <- as.numeric(roc(response = ObservedOutcome, predictor = PredictedRisk)$auc)
  
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
                                  TrueRisk_P11, TrueRisk_P10, TrueRisk_P01, ModelName) {
  
  #Note: This function should not be called directly: only from the main simulation function
  
  #Input: ObservedOutcome = a factor vector of observed outcomes, where each level is a 
  #               joint outcome category (e.g. {Y1=1, Y2=1}, {Y1=1, Y2=0}, etc.)
  #       PredictionsDataFrame = a data.frame of predicted risks. The columns for the
  #               predicted joint risks should be of the form [model name]_Pij where i/j={1,0} 
  
  require(VGAM)
  
  #Extract the relevant predicted risks
  PredictionsDataFrame <- select(PredictionsDataFrame, 
                                 matches(paste(ModelName, "\\_P[0-9]+", sep = "")))
  
  
  #Estimate calibration intercept (i.e. calibration-in-the-large)
  CalInt.model <- coefficients(vgam(ObservedOutcome ~ 1, 
                                    offset = data.matrix(PredictionsDataFrame %>%
                                                           rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                                            make.names(names(PredictionsDataFrame))))) %>%
                                                           mutate(P10_Z = log(P10 / P00),
                                                                  P01_Z = log(P01 / P00),
                                                                  P11_Z = log(P11 / P00)) %>%
                                                           select(P10_Z, P01_Z, P11_Z)),
                                    family = multinomial(refLevel = "(0,0)")))
  CalInt_P10 <- as.numeric(CalInt.model[1])
  CalInt_P01 <- as.numeric(CalInt.model[2])
  CalInt_P11 <- as.numeric(CalInt.model[3])
  
  
  #Estimate calibration slope
  k <- 4 #number of outcome categories
  CalSlope.model <- coefficients(vgam(ObservedOutcome ~ P10_Z + P01_Z + P11_Z, 
                                      data = PredictionsDataFrame %>%
                                        rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                         make.names(names(PredictionsDataFrame))))) %>%
                                        mutate(P10_Z = log(P10 / P00),
                                               P01_Z = log(P01 / P00),
                                               P11_Z = log(P11 / P00)) %>%
                                        select(P10_Z, P01_Z, P11_Z),
                                      family = multinomial(refLevel = "(0,0)"),
                                      constraints = list("(Intercept)" = diag(k - 1),
                                                         "P10_Z" = rbind(1, 0, 0),
                                                         "P01_Z" = rbind(0, 1, 0),
                                                         "P11_Z" = rbind(0, 0, 1))))
  CalSlope_P10 <- as.numeric(CalSlope.model["P10_Z"])
  CalSlope_P01 <- as.numeric(CalSlope.model["P01_Z"])
  CalSlope_P11 <- as.numeric(CalSlope.model["P11_Z"])
  
  ## NEED TO CHANGE DISCRIMINATION TO Polytomous Discrimination Index (PDI) 
  #Dicrimination: consider 1 vs. rest approach (i.e. {Y1=1, Y2=1} vs. rest, {Y1=1, Y2=0} vs. rest, etc.)
  # the Polytomous Discrimination Index (PDI) would allow a single measure, but is not posisble to compute
  # with large values of N. 
  Y11 <- ifelse(ObservedOutcome == "(1,1)", 1, 0)
  Y10 <- ifelse(ObservedOutcome == "(1,0)", 1, 0)
  Y01 <- ifelse(ObservedOutcome == "(0,1)", 1, 0)
  AUC_P11 <- as.numeric(roc(response = Y11, 
                            predictor = PredictionsDataFrame %>%
                              rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                               make.names(names(PredictionsDataFrame))))) %>%
                              .$P11)$auc)
  AUC_P10 <- as.numeric(roc(response = Y10, 
                            predictor = PredictionsDataFrame %>%
                              rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                               make.names(names(PredictionsDataFrame))))) %>%
                              .$P10)$auc)
  AUC_P01 <- as.numeric(roc(response = Y01, 
                            predictor = PredictionsDataFrame %>%
                              rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                               make.names(names(PredictionsDataFrame))))) %>%
                              .$P01)$auc)
  
  #Mean square error in predicted risks and 'true' risk
  MSE_P11 <- 1/length(Y11) * (sum((PredictionsDataFrame %>%
                                            rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                             make.names(names(PredictionsDataFrame))))) %>%
                                            .$P11 - 
                                            TrueRisk_P11)^2))
  MSE_P10 <- 1/length(Y10) * (sum((PredictionsDataFrame %>%
                                            rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                             make.names(names(PredictionsDataFrame))))) %>%
                                            .$P10 - 
                                            TrueRisk_P10)^2))
  MSE_P01 <- 1/length(Y01) * (sum((PredictionsDataFrame %>%
                                            rename_all(~(sub("([A-Z,a-z]+)_", "", 
                                                             make.names(names(PredictionsDataFrame))))) %>%
                                            .$P01 - 
                                            TrueRisk_P01)^2))
  
  ## Store performance results in a data.frame and return
  Results <- data.frame(CalInt_P11,
                        CalInt_P10,
                        CalInt_P01,
                        CalSlope_P11,
                        CalSlope_P10,
                        CalSlope_P01,
                        AUC_P11,
                        AUC_P10,
                        AUC_P01,
                        MSE_P11,
                        MSE_P10,
                        MSE_P01)
  return(Results)
}

