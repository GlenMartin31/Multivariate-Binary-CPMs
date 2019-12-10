####-----------------------------------------------------------------------------------------
## Run the Simulation
####-----------------------------------------------------------------------------------------

#First run all the functions located in here::here("Code", "Simulation", "Simulation Functions.R")

simulationcall.fnc(iter = 100,
                   N = 5000, 
                   rho.vals = c(0.95, 0.75, 0.5, 0.25, 0),
                   beta_1_true = c(-1.0, log(2), log(1.00)),
                   beta_2_true = c(-1.5, log(1.00), log(3)),
                   filename = "SimulationResults",
                   startingseed = 4651)


simulationcall.fnc(iter = 100,
                   N = 5000, 
                   rho.vals = c(0.95, 0.75, 0.5, 0.25, 0),
                   beta_1_true = c(-3.0, log(2), log(1.00)),
                   beta_2_true = c(-3.5, log(1.00), log(3)),
                   filename = "SimulationResults_LowerMarginalProportion",
                   startingseed = 9473)
