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
## This script runs the simulations across all scenarios: runs in parallel using furrr
####-----------------------------------------------------------------------------------------

#Load the simulation functions
source(here::here("Code", "00_Simulation_Functions.R"))

library(tidyverse)
library(furrr)

# Define a dataset that includes all combinations of simulation parameters (i.e. simulation cases)
sims_parameters <- crossing(
  N = 5000,
  rho = c(0.95, 0.75, 0.5, 0.25, 0)
)

# number of repeats per scenario
n_rep <- 100

#Set up parallel process
plan(multiprocess)

#Run the main simulation study:
set.seed(46513)
sims_main <- sims_parameters %>%
  mutate(results = future_pmap(.l = list(N = N,
                                         rho = rho),
                               .f = simulation_nrun_fnc,
                               n_iter = n_rep,
                               beta_1_true = c(-1.0, log(2), log(1.00)),
                               beta_2_true = c(-1.5, log(1.00), log(3)),
                               .progress = TRUE)
         )
write_rds(sims_main, path = here::here("Data", "sims_main.RDS"))

#Run a sensitivity analysis with a lower (marginal) outcome proportion for Y1 and Y2:
set.seed(94735)
sims_sensanalysis <- sims_parameters %>%
  mutate(results = future_pmap(.l = list(N = N,
                                         rho = rho),
                               .f = simulation_nrun_fnc,
                               n_iter = n_rep,
                               
                               beta_1_true = c(-3.0, log(2), log(1.00)),
                               beta_2_true = c(-3.5, log(1.00), log(3)),
                               .progress = TRUE)
  )

write_rds(sims_sensanalysis, path = here::here("Data", "sims_sensanalysis.RDS"))
