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
library(tidyverse)

##Load the simulation results
sims_main <- read_rds(here::here("Data", "sims_main.RDS"))
sims_sensanalysis <- read_rds(here::here("Data", "sims_sensanalysis.RDS"))

###-------------------------------------------------------------------------------------
## Combine Simulation scenarios
###-------------------------------------------------------------------------------------
sims_all <- sims_main %>%
  mutate("Simulation" = "main") %>%
  bind_rows(sims_sensanalysis %>%
              mutate("Simulation" = "Sensitivity")) %>%
  mutate("Simulation" = factor(Simulation, levels = c("main", "Sensitivity")),
         
         results = map(results, as_tibble)) %>%
  unnest(cols = c(results)) 


###-------------------------------------------------------------------------------------
## Summarise the relationship between rho_observed and cor(Y1, Y2), as well as the
# observed outcome proportion
###-------------------------------------------------------------------------------------
View(sims_all %>%
       group_by(rho, Simulation) %>%
       summarise(round(mean(rho_observed), 3),
                 round(mean(Obs_P11), 3),
                 round(mean(Obs_P10), 3),
                 round(mean(Obs_P01), 3),
                 .groups = "drop"))

mean(Results$Obs_Py1)
mean(Results$Obs_Py2)


###-------------------------------------------------------------------------------------
## Summarise results in each simulation scenario across iterations
###-------------------------------------------------------------------------------------
SummaryResults <- sims_all %>%
  select(Simulation, Iteration, Model, rho, N,
         starts_with("CalInt_"),
         starts_with("CalSlope_"),
         starts_with("AUC_"), 
         starts_with("PDI_"),
         "MSE_PY1",
         "MSE_PY2",
         "MSE") %>%
  rename("MSE_Joint" = "MSE") %>%
  group_by(Model, Simulation, rho, N) %>%
  summarise_all(list("mean" = mean, "sd" = sd)) %>%
  ungroup() %>%
  select(-starts_with("Iteration")) 

###-------------------------------------------------------------------------------------
## Calibration Intercept plots
###-------------------------------------------------------------------------------------

#Joint Risks: main simulation
SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("CalInt_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalInt_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration-in-the-large")

# ggsave(here::here("Outputs", "Manuscript", "Fig1_CalibrationInt.tiff"), 
#        height = 7, width = 8, units = "in", dpi = 300)


#Joint Risks: sensitivity simulation
SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>%
  select(Model, rho, starts_with("CalInt_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalInt_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration-in-the-large")



#Marginal Risks: main simulation
SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("CalInt_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalInt_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration-in-the-large")


#Marginal Risks: sensitivity simulation
SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>%
  select(Model, rho, starts_with("CalInt_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalInt_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration-in-the-large")



###-------------------------------------------------------------------------------------
## Calibration Slope plots
###-------------------------------------------------------------------------------------

#Joint Risks: main simulation
SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("CalSlope_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalSlope_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration Slope")

# ggsave(here::here("Outputs", "Manuscript", "Fig2_CalibrationSlope.tiff"), 
#        height = 7, width = 8, units = "in", dpi = 300)

#Joint Risks: sensitivity simulation
SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>%
  select(Model, rho, starts_with("CalSlope_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalSlope_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration Slope")


#Marginal Risks: main simulation
SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("CalSlope_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalSlope_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration Slope")


#Marginal Risks: sensitivity simulation
SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>%
  select(Model, rho, starts_with("CalSlope_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalSlope_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration Slope")

###-------------------------------------------------------------------------------------
## PDI plots
###-------------------------------------------------------------------------------------

#Main simulation
SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("PDI_Overall_")) %>%
  rename("PDI_mean" = "PDI_Overall_mean",
         "PDI_sd" = "PDI_Overall_sd") %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate","Quantity", "PDI_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~rho, scales = "fixed", nrow = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = (1/4), slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Overall PDI") +
  ylim(0.25,0.75)

# ggsave(here::here("Outputs", "Manuscript", "Fig3_PDI.tiff"), 
#        height = 7, width = 8, units = "in", dpi = 300)

SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("PDI_")) %>%
  select(-starts_with("PDI_Overall_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "PDI_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = (1/4), slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Outcome-Specific PDI") +
  ylim(0.25,1)



#Sensitivity simulation
SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>%
  select(Model, rho, starts_with("PDI_Overall_")) %>%
  rename("PDI_mean" = "PDI_Overall_mean",
         "PDI_sd" = "PDI_Overall_sd") %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate","Quantity", "PDI_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~rho, scales = "fixed", nrow = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = (1/4), slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Overall PDI") +
  ylim(0.25,0.75)

###-------------------------------------------------------------------------------------
## AUC plots for marginal risk estimates
###-------------------------------------------------------------------------------------

SummaryResults %>% 
  filter(Simulation == "main") %>% 
  select(Model, rho, starts_with("AUC_")) %>%
  gather(key = "Estimate", value, -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "AUC_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("AUC")


SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>% 
  select(Model, rho, starts_with("AUC_")) %>%
  gather(key = "Estimate", value, -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "AUC_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("AUC")


###-------------------------------------------------------------------------------------
## MSE plots
###-------------------------------------------------------------------------------------

#Joint risks: Main simulation
SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("MSE_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "MSE_([A-Za-z0-9]+)_([a-z]+)") %>%
  filter(Outcome == "Joint") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~rho, scales = "fixed", nrow = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = (1/4), slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("MSE")

# ggsave(here::here("Outputs", "Manuscript", "Fig4_MSE.tiff"), 
#        height = 7, width = 8, units = "in", dpi = 300)


#Joint risks: Sensitivity simulation
SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>%
  select(Model, rho, starts_with("MSE_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "MSE_([A-Za-z0-9]+)_([a-z]+)") %>%
  filter(Outcome == "Joint") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~rho, scales = "fixed", nrow = 1) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = (1/4), slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("MSE")


#Marginal risks: Main simulation
SummaryResults %>% 
  filter(Simulation == "main") %>%
  select(Model, rho, starts_with("MSE_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "MSE_([A-Za-z0-9]+)_([a-z]+)") %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = (1/4), slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("MSE")


#Marginal risks: Sensitivity simulation
SummaryResults %>% 
  filter(Simulation == "Sensitivity") %>%
  select(Model, rho, starts_with("MSE_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho) %>%
  extract("Estimate", c("Outcome", "Quantity"), "MSE_([A-Za-z0-9]+)_([a-z]+)") %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (qnorm(0.975)*sd)),
         upper = (mean + (qnorm(0.975)*sd))) %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "SR", "PCC", "MLR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = (1/4), slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("MSE")




