####-----------------------------------------------------------------------------------------
## Run the Simulation
####-----------------------------------------------------------------------------------------

#First run all the functions located in here::here("Code", "Simulation", "Simulation Functions.R")

simulationcall.fnc(iter = 100,
                   N = 5000, 
                   rho.vals = c(0, 0.2, 0.5, 0.75, 0.95),
                   K = 100,
                   filename = "SimulationResults",
                   startingseed = 4651)



library(tidyverse)
Results <- read.table(here::here("Data", "SimulationResults.txt"), sep = "|", header = TRUE)

#Summarise the relationship between rho_latent and cor(Y1, Y2), as well as the
#observed outcome rates
View(Results %>%
       group_by(rho_latent) %>%
       summarise(mean(rho_Y1Y2),
                 mean(Obs_Py1),
                 mean(Obs_Py2),
                 mean(Obs_P11),
                 mean(Obs_P10),
                 mean(Obs_P01)))

mean(Results$Obs_Py1)
mean(Results$Obs_Py2)



#Now plot the performance results across all simulation scenarios
SummaryResults <- Results %>%
  select(Iteration, Model, rho_latent, 
         starts_with("CalInt_"),
         starts_with("CalSlope_"),
         starts_with("AUC_"), 
         starts_with("BrierScore_")) %>%
  group_by(Model, rho_latent) %>%
  summarise_all(list("mean" = mean, "sd" = sd)) %>%
  select(-starts_with("Iteration"))


#Calibration Intercept plot
SummaryResults %>% 
  select(Model, rho_latent, starts_with("CalInt_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalInt_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed") +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  xlab("Model") + ylab("Calibration Intercept")


#Calibration Slope plot
SummaryResults %>% 
  select(Model, rho_latent, starts_with("CalSlope_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalSlope_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed") +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  xlab("Model") + ylab("Calibration Slope")


#AUC plot
SummaryResults %>% 
  select(Model, rho_latent, starts_with("AUC_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "AUC_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed") +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  xlab("Model") + ylab("AUC")



#Brier Score plot
SummaryResults %>% 
  select(Model, rho_latent, starts_with("BrierScore_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "BrierScore_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed") +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  xlab("Model") + ylab("Brier Score")
