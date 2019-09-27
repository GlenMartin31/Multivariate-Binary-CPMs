library(tidyverse)
Results <- read.table(here::here("Data", "SimulationResults.txt"), sep = "|", header = TRUE)

#Summarise the relationship between rho_latent and cor(Y1, Y2), as well as the
#observed outcome rates
View(Results %>%
       group_by(rho_latent) %>%
       summarise(round(mean(rho_Y1Y2), 3),
                 round(mean(Obs_P11), 3),
                 round(mean(Obs_P10), 3),
                 round(mean(Obs_P01), 3)))

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


####-----------------------------------------------------------------------------------------
## Calibration Intercept plots
####-----------------------------------------------------------------------------------------

#Joint Risks
SummaryResults %>% 
  select(Model, rho_latent, starts_with("CalInt_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalInt_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration Intercept")

ggsave(here::here("Outputs", "Manuscript", "Fig1_CalibrationInt.tiff"), 
       height = 7, width = 8, units = "in", dpi = 300)


#Marginal Risks
SummaryResults %>% 
  select(Model, rho_latent, starts_with("CalInt_")) %>%
  gather(key = "Estimate", value = "value", -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalInt_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Model") + ylab("Calibration Intercept")



####-----------------------------------------------------------------------------------------
## Calibration Slope plots
####-----------------------------------------------------------------------------------------

#Joint Risks
SummaryResults %>% 
  select(Model, rho_latent, starts_with("CalSlope_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalSlope_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("Calibration Slope")

ggsave(here::here("Outputs", "Manuscript", "Fig2_CalibrationSlope.tiff"), 
       height = 7, width = 8, units = "in", dpi = 300)


#Marginal Risks
SummaryResults %>% 
  select(Model, rho_latent, starts_with("CalSlope_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "CalSlope_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed") +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("Calibration Slope")


####-----------------------------------------------------------------------------------------
## AUC plots
####-----------------------------------------------------------------------------------------

#Joint Risk
SummaryResults %>% 
  select(Model, rho_latent, starts_with("AUC_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "AUC_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("AUC")

ggsave(here::here("Outputs", "Manuscript", "Fig3_AUC.tiff"), 
       height = 7, width = 8, units = "in", dpi = 300)


#Marginal Risk
SummaryResults %>% 
  select(Model, rho_latent, starts_with("AUC_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "AUC_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("AUC")


####-----------------------------------------------------------------------------------------
## Brier Score plots
####-----------------------------------------------------------------------------------------

#Joint Risk
SummaryResults %>% 
  select(Model, rho_latent, starts_with("BrierScore_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "BrierScore_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "P01" |
           Outcome == "P10" |
           Outcome == "P11") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 3) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("Brier Score")

ggsave(here::here("Outputs", "Manuscript", "Fig4_BrierScore.tiff"), 
       height = 7, width = 8, units = "in", dpi = 300)

#Marginal Risk
SummaryResults %>% 
  select(Model, rho_latent, starts_with("BrierScore_")) %>%
  gather(key = "Estimate", value, -Model, -rho_latent) %>%
  extract("Estimate", c("Outcome", "Quantity"), "BrierScore_([A-Z0-9]+)_([a-z]+)") %>%
  spread(Quantity, value) %>%
  mutate(lower = (mean - (1.96*sd)),
         upper = (mean + (1.96*sd))) %>%
  filter(Outcome == "PY1" |
           Outcome == "PY2") %>%
  ungroup() %>%
  mutate(Model = fct_relevel(Model,
                             c("Univariate", "HRR", "PCC", "MLR", "SR", "MLM", "MPM"))) %>%
  ggplot(aes(x = Model, y = mean)) +
  facet_wrap(~Outcome + rho_latent, scales = "fixed", nrow = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Model") + ylab("Brier Score")

