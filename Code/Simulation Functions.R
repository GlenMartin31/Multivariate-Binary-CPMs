####-----------------------------------------------------------------------------------------
## Define a function to call the simulation across multiple iterations
####-----------------------------------------------------------------------------------------
simulationcall.fnc <- function(iter = 100,
                               N, 
                               rho.vals = 0,
                               K = N,
                               filename,
                               startingseed) {
  #Input: iter = number of iterations to repeat the simulation over,
  #       N = the number of observations to generate
  #       rho.vals = strength of the association between the latent variables used to derive Y1 and Y2. When rho=0, Y1 and Y2 are independent
  #       K = the number of folds for k-fold cross validation for Hierarchical Related CPM. K = N gives leave-one-out CV
  #       filename = character variable giving the name of the text file to store the results
  #       startingseed = starting seed for the random number generator
  
  library(tidyverse)
  library(pROC)
  library(rjags)
  library(coda)
  library(pbivnorm)
  
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
                     "BrierScore_PY1", "BrierScore_PY2",
                     "CalInt_P11", "CalInt_P10", "CalInt_P01",
                     "CalSlope_P11", "CalSlope_P10", "CalSlope_P01",
                     "AUC_P11", "AUC_P10", "AUC_P01",
                     "BrierScore_P11", "BrierScore_P10", "BrierScore_P01")
    write(OutputNames, 
          here::here("Data", paste(filename, ".txt", sep = "")), 
          append = FALSE, 
          sep = "|", 
          ncolumns = length(OutputNames)) 
  }
  
  
  #Set 'true' coefficients
  beta_1_true <- c(-1.0, log(2), log(1.00))
  beta_2_true <- c(-1.5, log(1.00), log(3))
  
  for (rho in rho.vals) {
    
    pb <- winProgressBar(title = "Iteration Progress bar",
                         label = paste("Iterations 0% Completed"),
                         min = 0, max = 100, initial = 0)
    for (iteration in 1:iter) {
      
      CombinedData <- DataGenerating.fnc(N = sum(N, 10000), #add on 10000 additional observations to serve as a validation set 
                                         beta_1_true = beta_1_true, 
                                         beta_2_true = beta_2_true, 
                                         rho = rho)
      
      IPD <- CombinedData %>%
        sample_n(N, replace = FALSE) #IPD and validation cohorts are a random sample of the same underlying data generating mechasm.
      Validation.Population <- CombinedData %>%
        filter(ID %in% IPD$ID == FALSE)
      
      
      ####-----------------
      ## Model Development
      ####-----------------
      #Start by fitting each model independently
      Uni.m1 <- glm(Y1 ~ X1 + X2, data = IPD, family = binomial(link = "logit"))
      Uni.m2 <- glm(Y2 ~ X1 + X2, data = IPD, family = binomial(link = "logit"))
      Univariate_PY1 <- predict(Uni.m1, newdata = Validation.Population, type = "response")
      Univariate_PY2 <- predict(Uni.m2, newdata = Validation.Population, type = "response")
      Validation.Population$Univariate_P11 <- Univariate_PY1 * Univariate_PY2
      Validation.Population$Univariate_P10 <- Univariate_PY1 * (1 - Univariate_PY2)
      Validation.Population$Univariate_P01 <- (1 - Univariate_PY1) * Univariate_PY2
      
      
      #Hierarchical Related Regression
      HRR <- Hierarchical.Related.CPM.fnc(DevelopmentData = IPD, 
                                          TestData = Validation.Population,
                                          K = K)
      Validation.Population$HRR_P11 <- HRR$Pi_Y1 * HRR$Pi_Y2
      Validation.Population$HRR_P10 <- HRR$Pi_Y1 * (1 - HRR$Pi_Y2)
      Validation.Population$HRR_P01 <- (1 - HRR$Pi_Y1) * HRR$Pi_Y2
      
      
      #Probabilistic Classifier Chain
      PCC <- Probabilistic.Classifier.Chain.fnc(DevelopmentData = IPD, 
                                                TestData = Validation.Population)
      Validation.Population$PCC_P11 <- PCC$P11
      Validation.Population$PCC_P10 <- PCC$P10
      Validation.Population$PCC_P01 <- PCC$P01
      
      
      #Multivariate Logistic Model (MLM)
      MLM <- suppressWarnings(Multivariate.Logistic.Reg.fnc(X = IPD %>% select(X1, X2), 
                                                            Y1 = IPD$Y1, 
                                                            Y2 = IPD$Y2))
      MLM_PY1 <- (1 + exp(-as.numeric((Validation.Population %>%
                                         select(starts_with("X")) %>%
                                         data.matrix()) %*% MLM$beta1)))^(-1) #P(y1=1 | X,beta_1)
      MLM_PY2 <- (1 + exp(-as.numeric((Validation.Population %>%
                                         select(starts_with("X")) %>%
                                         data.matrix()) %*% MLM$beta2)))^(-1) #P(y2=1 | X,beta_2)
      Validation.Population$MLM_P11 <- ((MLM_PY1*MLM_PY2) + (MLM$rho12*sqrt(MLM_PY1*(1 - MLM_PY1)*MLM_PY2*(1 - MLM_PY2))))
      Validation.Population$MLM_P10 <- ((MLM_PY1*(1 - MLM_PY2)) - (MLM$rho12*sqrt(MLM_PY1*(1 - MLM_PY1)*MLM_PY2*(1 - MLM_PY2))))
      Validation.Population$MLM_P01 <- ((MLM_PY2*(1 - MLM_PY1)) - (MLM$rho12*sqrt(MLM_PY1*(1 - MLM_PY1)*MLM_PY2*(1 - MLM_PY2))))
      
      
      #Multivariate Probit Model (MPM)
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
      samps <- coda.samples( BayesianProbitMultivariateModel, c('Beta','rho'), n.iter = 10000 )
      tidy.samps <- samps[[1]][5001:10000,] #set first 5000 samples as burn-in
      post.means <- colMeans(tidy.samps)
      X.Beta.Y1 <- as.numeric(cbind(Validation.Population$X0, 
                                    Validation.Population$X1,
                                    Validation.Population$X2) %*% post.means[paste("Beta[", 1:3, ",1]", sep = "")])
      X.Beta.Y2 <- as.numeric(cbind(Validation.Population$X0, 
                                    Validation.Population$X1,
                                    Validation.Population$X2) %*% post.means[paste("Beta[", 1:3, ",2]", sep = "")])
      
      Validation.Population$MPM_P11 <- pbivnorm(x = cbind(X.Beta.Y1, X.Beta.Y2), 
                                                rho = post.means["rho"])
      Validation.Population$MPM_P10 <- pbivnorm(x = cbind(X.Beta.Y1, -X.Beta.Y2), 
                                                rho = -post.means["rho"])
      Validation.Population$MPM_P01 <- pbivnorm(x = cbind(-X.Beta.Y1, X.Beta.Y2), 
                                                rho = -post.means["rho"])
      
      
      #Multivariate Stacked Regression
      IPD$Univariate_LP1 <- predict(Uni.m1, newdata = IPD, type = "link")
      IPD$Univariate_LP2 <- predict(Uni.m2, newdata = IPD, type = "link")
      MultivariateStackedRegression.m1 <- glm(Y1 ~ Univariate_LP1 + Univariate_LP2, 
                                              data = IPD, 
                                              family = binomial(link = "logit"))
      MultivariateStackedRegression.m2 <- glm(Y2 ~ Univariate_LP1 + Univariate_LP2, 
                                              data = IPD, 
                                              family = binomial(link = "logit"))
      
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
      Validation.Population$SR_P11 <- SR_PY1 * SR_PY2
      Validation.Population$SR_P10 <- SR_PY1 * (1 - SR_PY2)
      Validation.Population$SR_P01 <- (1 - SR_PY1) * SR_PY2
      
      
      ####-----------------
      ## Model Validation
      ####-----------------
      
      Predictions <- Validation.Population %>%
        select(Y1, Y2,
               Univariate_P11, Univariate_P10, Univariate_P01,
               HRR_P11, HRR_P10, HRR_P01,
               PCC_P11, PCC_P10, PCC_P01,
               MLM_P11, MLM_P10, MLM_P01,
               MPM_P11, MPM_P10, MPM_P01,
               SR_P11, SR_P10, SR_P01) %>%
        mutate(Y11 = ifelse(Y1 == 1 & Y2 == 1, 1, 0),
               Y10 = ifelse(Y1 == 1 & Y2 == 0, 1, 0),
               Y01 = ifelse(Y1 == 0 & Y2 == 1, 1, 0),
               
               Univariate_Py1 = Univariate_P10 + Univariate_P11,
               Univariate_Py2 = Univariate_P01 + Univariate_P11,
               
               HRR_Py1 = HRR_P10 + HRR_P11,
               HRR_Py2 = HRR_P01 + HRR_P11,
               
               PCC_Py1 = PCC_P10 + PCC_P11,
               PCC_Py2 = PCC_P01 + PCC_P11,
               
               MLM_Py1 = MLM_P10 + MLM_P11,
               MLM_Py2 = MLM_P01 + MLM_P11,
               
               MPM_Py1 = MPM_P10 + MPM_P11,
               MPM_Py2 = MPM_P01 + MPM_P11,
               
               SR_Py1 = SR_P10 + SR_P11,
               SR_Py2 = SR_P01 + SR_P11) %>%
        select(Y1, Y2,
               Y11, 
               Y10, 
               Y01,
               ends_with("_Py1"), ends_with("_Py2"),
               ends_with("_P11"), ends_with("_P10"), ends_with("_P01"))
      
      
      Results <- NULL
      for (ModelName in c("Univariate", "HRR", "PCC", "MLM", "MPM", "SR")) {
        Results <- Results %>%
          bind_rows(Performance.fnc(ObservedOutcome = Predictions$Y1, 
                                    PredictedRisk = eval(parse(text = paste("Predictions$",
                                                                            ModelName,
                                                                            "_Py1", sep = ""))),
                                    OutcomeSubscript = "PY1") %>%
                      mutate("Model" = ModelName) %>%
                      left_join(Performance.fnc(ObservedOutcome = Predictions$Y2, 
                                                PredictedRisk = eval(parse(text = paste("Predictions$",
                                                                                        ModelName,
                                                                                        "_Py2", sep = ""))),
                                                OutcomeSubscript = "PY2") %>%
                                  mutate("Model" = ModelName),
                                by = "Model") %>%
                      left_join(Performance.fnc(ObservedOutcome = Predictions$Y11, 
                                                PredictedRisk = eval(parse(text = paste("Predictions$",
                                                                                        ModelName,
                                                                                        "_P11", sep = ""))),
                                                OutcomeSubscript = "P11") %>%
                                  mutate("Model" = ModelName),
                                by = "Model") %>%
                      left_join(Performance.fnc(ObservedOutcome = Predictions$Y10, 
                                                PredictedRisk = eval(parse(text = paste("Predictions$",
                                                                                        ModelName,
                                                                                        "_P10", sep = ""))),
                                                OutcomeSubscript = "P10") %>%
                                  mutate("Model" = ModelName),
                                by = "Model") %>%
                      left_join(Performance.fnc(ObservedOutcome = Predictions$Y01, 
                                                PredictedRisk = eval(parse(text = paste("Predictions$",
                                                                                        ModelName,
                                                                                        "_P01", sep = ""))),
                                                OutcomeSubscript = "P01") %>%
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
               "Obs_P11" = mean(Predictions$Y11), 
               "Obs_P10" = mean(Predictions$Y10), 
               "Obs_P01" = mean(Predictions$Y01)) %>%
        select("Iteration", 
               "Model",
               "Obs_Py1", "Obs_Py2", "Obs_P11", "Obs_P10", "Obs_P01", 
               "rho_latent", 
               "rho_Y1Y2",
               "CalInt_PY1", "CalInt_PY2",
               "CalSlope_PY1", "CalSlope_PY2",
               "AUC_PY1", "AUC_PY2",
               "BrierScore_PY1", "BrierScore_PY2",
               "CalInt_P11", "CalInt_P10", "CalInt_P01",
               "CalSlope_P11", "CalSlope_P10", "CalSlope_P01",
               "AUC_P11", "AUC_P10", "AUC_P01",
               "BrierScore_P11", "BrierScore_P10", "BrierScore_P01")
      
      write.table(Results,
                  here::here("Data", paste(filename, ".txt", sep = "")),
                  append = TRUE,
                  sep = "|",
                  row.names = FALSE,
                  col.names = FALSE)
      rm(Results)
      
      
      info <- sprintf(paste("Iterations %d%% Completed"), round((iteration/(iter)*100)))
      setWinProgressBar(pb, iteration/((iter))*100, label = info)
    }
    close(pb) 
    
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
  
  
  #Define a function for the inverse logistic CDF
  Inverse.Logistic.CDF <- function(X) {
    return( log(X/(1 - X)) )
  }
  
  X <- cbind(X0 = 1,
             X1 = rnorm(N, 0, 1), 
             X2 = rnorm(N, 0, 1))
  
  #Generate random vectors from standard multivariate normal distribution with inputted correlation, rho
  Z <- MASS::mvrnorm(n = N, mu = c(0,0), 
                     Sigma = toeplitz(c(1, rho)))
  Epsilon <- apply(pnorm(Z), 2, Inverse.Logistic.CDF) #aaply the inverse probability transform
  
  IPD <- data.frame("ID" = seq(from = 1, to = N, by = 1),
                    X) %>%
    mutate(LP1 = as.vector(X %*% beta_1_true),
           LP2 = as.vector(X %*% beta_2_true),
           
           Y1 = ifelse(Epsilon[,1] <= LP1, 1, 0),
           Y2 = ifelse(Epsilon[,2] <= LP2, 1, 0))
  
  return(IPD)
}


####-----------------------------------------------------------------------------------------
## Function to fit the Hierarchical Related CPM
####-----------------------------------------------------------------------------------------

Hierarchical.Related.CPM.fnc <- function(DevelopmentData, TestData, K) {
  #Inputs: DevelopmentData = data on which to derive the Hierarchical Related CPM
  #        TestData = data on which to predict risk for each outcome (marginal), 
  #        K = the number of k-fold cross validation: K = nrow(DevelopmentData) is leave-one-out CV
  
  #Perform leave-one-out cross validation to obtain estimates of pi1 and pi2 from each 
  # 'initial' model in each 'permutation' stage
  DevelopmentData$Pi1_m1 <- NA
  DevelopmentData$Pi2_m2 <- NA
  
  folds <- cut(seq(1,nrow(DevelopmentData)), breaks = K, labels = FALSE)
  cv.pb <- txtProgressBar()
  print("Cross Validation Progress for Hierarchical Related CPM")
  for (i in 1:K) {
    testIndexes <- which(folds == i, arr.ind = TRUE)
    
    CV.Data <- DevelopmentData %>% slice(-testIndexes)
    
    m1_Pi1 <- glm(Y1 ~ X1 + X2, data = CV.Data, family = binomial(link = "logit"))
    DevelopmentData$Pi1_m1[testIndexes] <- predict(m1_Pi1, 
                                                   newdata = DevelopmentData %>% 
                                                     slice(testIndexes), 
                                                   type = "response")
    
    
    m2_Pi2 <- glm(Y2 ~ X1 + X2, data = CV.Data, family = binomial(link = "logit"))
    DevelopmentData$Pi2_m2[testIndexes] <- predict(m2_Pi2, 
                                                   newdata = DevelopmentData %>%
                                                     slice(testIndexes), 
                                                   type = "response")
    
    setTxtProgressBar(cv.pb, i/K)
  }; close(cv.pb)
  ##Now fit the models to the full development data to arrive at parameter estimates for future prediction
  #In the first 'permutation' we start with estimating Y1
  m1_Pi1 <- glm(Y1 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  m1_Pi2 <- glm(Y2 ~ X1 + X2 + Pi1_m1, data = DevelopmentData, family = binomial(link = "logit"))
  #In the second 'permutation' we start with estimating Y2
  m2_Pi2 <- glm(Y2 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  m2_Pi1 <- glm(Y1 ~ X1 + X2 + Pi2_m2, data = DevelopmentData, family = binomial(link = "logit"))
  
  
  #Now predict for each observation in the TestData
  TestData$Pi1_m1 <- predict(m1_Pi1, newdata = TestData, type = "response")
  TestData$Pi2_m1 <- predict(m1_Pi2, newdata = TestData, type = "response")
  TestData$Pi2_m2 <- predict(m2_Pi2, newdata = TestData, type = "response")
  TestData$Pi1_m2 <- predict(m2_Pi1, newdata = TestData, type = "response")
  
  Pi_Y1 <- apply(cbind(TestData$Pi1_m1, TestData$Pi1_m2), 1, mean)
  Pi_Y2 <- apply(cbind(TestData$Pi2_m1, TestData$Pi2_m2), 1, mean)
  
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
  #Inputs: DevelopmentData = data on which to derive the Hierarchical Related CPM
  #        TestData = data on which to predict risk for each outcome (marginal)
  
  ##Fit the models to the full development data to arrive at parameter estimates for future prediction
  
  #In the first 'permutation' we start with estimating P(Y1) and then estimate P(Y2 | Y1)
  m1_Pi1 <- glm(Y1 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  m1_Pi2 <- glm(Y2 ~ X1 + X2 + Y1, data = DevelopmentData, family = binomial(link = "logit"))
  #In the second 'permutation' we start with estimating P(Y2) and then estimate P(Y1 | Y2)
  m2_Pi2 <- glm(Y2 ~ X1 + X2, data = DevelopmentData, family = binomial(link = "logit"))
  m2_Pi1 <- glm(Y1 ~ X1 + X2 + Y2, data = DevelopmentData, family = binomial(link = "logit"))
  
  
  ##Now predict for each observation in the TestData
  
  #only extract the covaraite information (i.e. dont use the observed reponses to mimic prediction time):
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
  
  
  P11 <- apply(cbind(P11_m1, P11_m2), 1, mean)
  P10 <- apply(cbind(P10_m1, P10_m2), 1, mean)
  P01 <- apply(cbind(P01_m1, P01_m2), 1, mean)
  
  return(list("P11" = P11,
              "P10" = P10,
              "P01" = P01))
}



####-----------------------------------------------------------------------------------------
## Numerical Optimisation of Multivariate Logistic Regression model
####-----------------------------------------------------------------------------------------
Multivariate.Logistic.Reg.fnc <- function(X, Y1, Y2) {
  
  #Define log-likelihood function for the Multivariate Logistic Regression model
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
  
  X <- data.matrix(X)
  
  m1 <- glm(Y1 ~ X, family = binomial(link = "logit"))
  m2 <- glm(Y2 ~ X, family = binomial(link = "logit"))
  
  
  maxlike <- optim(par = c(as.numeric(coef(m1)), 
                           as.numeric(coef(m2)), 
                           0),
                   fn = loglike,
                   method = "Nelder-Mead",
                   X = cbind(1, X), #apend intercept to data
                   Y1 = Y1, Y2 = Y2)
  
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
Performance.fnc <- function(ObservedOutcome, PredictedRisk, OutcomeSubscript) {
  require(pROC)
  
  LP <- log(PredictedRisk / (1 - PredictedRisk))
  
  CalInt.model <- glm(ObservedOutcome ~ offset(LP), family = binomial(link = "logit"))
  CalSlope.model <- glm(ObservedOutcome ~ LP, family = binomial(link = "logit"))
  AUC <- as.numeric(roc(response = ObservedOutcome, predictor = PredictedRisk)$auc)
  
  BrierScore <- 1/length(ObservedOutcome) * (sum((PredictedRisk - ObservedOutcome)^2))
  
  Results <- data.frame(as.numeric(coef(CalInt.model)),
                        as.numeric(coef(CalSlope.model)[2]),
                        AUC,
                        BrierScore)
  names(Results) <- c(paste("CalInt", OutcomeSubscript, sep = "_"),
                      paste("CalSlope", OutcomeSubscript, sep = "_"),
                      paste("AUC", OutcomeSubscript, sep = "_"),
                      paste("BrierScore", OutcomeSubscript, sep = "_"))
  return(Results)
}