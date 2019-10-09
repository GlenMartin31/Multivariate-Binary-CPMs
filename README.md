# Multivariate Logistic CPMs
Repo for the paper entitled "Developing Clinical Prediction Models for Multiple Binary Outcomes: a comparison of multivariate approaches", which is in review stage. The paper explores statistical methodology to allow the development of clinical prediction models (CPMs) for two binary outcomes simultaneously, rather than modelling as seperate outcomes. 

The repo contains the coding scripts and results from the simulation study described in the paper as follows:
## Code sub-folder
This contains the R scripts used in the simulation study; namely the R scripts entitled "Simulation Functions", "Simulation Call" and "Simulation Analysis", which should be run in this order. Additionally, this folder also contains the SQL script used to extract the cohort from the MIMIC-III database, which was used as part of the empirical study described in the paper.

## Data sub-folder
This contains a text file of the results of the simulation, as outlined and described in the paper. This text file can be loaded into R and analysed as per the "Simulation Analysis.R" file, which is located in the "Code" sub-folder
