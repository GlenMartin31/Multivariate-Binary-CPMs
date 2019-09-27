####-----------------------------------------------------------------------------------------
## Run the Simulation
####-----------------------------------------------------------------------------------------

#First run all the functions located in here::here("Code", "Simulation", "Simulation Functions.R")

simulationcall.fnc(iter = 100,
                   N = 5000, 
                   rho.vals = c(0, 0.2, 0.5, 0.75, 0.95),
                   filename = "SimulationResults",
                   startingseed = 4651)
