###
### Alison C. Ketz 08/31/2021
###
### Transmission following Heisey et al (2010) rcode
### for fitting FOI model given the cross-sectional harvest data
###

###########################################################
### Preliminaries
###########################################################

rm(list = ls())

setwd("~/Documents/integrate_s_foi/s_foi_v5/joint_survival_foi")

library(viridis)
library(RColorBrewer)
library(Hmisc)
library(lubridate)
library(readxl)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(xtable)
library(nimble)
library(tidyverse)
library(dplyr)
library(lattice)
library(foreign)
library(rgdal)
library(rgeos)
library(zoo)
library(spdep)
library(parallel)
library(doParallel)
library(coda)
library(INLA)
library(sf)
library(terra)
library(splines)
library(MetBrewer)
library(ggforce)

###########################################################
### Source summary function for posteriors
###########################################################

source("summarize.R")

###########################################################
# Load Data
###########################################################

source("02_load_clean_data_foi.R")
source("03_load_clean_data_survival.R")

###############################################################
### Load/clean Spatially referenced data
###############################################################

source("04_homerange_data_survival.R")

###########################################################
### Format data for fitting age-period survival models
###########################################################

source("05_format_data_survival.R")

##########################################################
### Setup collar data for FOI + Survival
##########################################################

source("06_format_data_combine_foi_surv.R")

###########################################################
### Setup consts etc for running the model
###########################################################

source("07_prelim_survival.R")
source("08_prelim_foi.R")
source("09_prelim_collar_foi.R")

###########################################################
### Run model
###########################################################

# source("10_distributions_check.R")
source("10_distributions_fast.R")

###########################################################
### Run model
###########################################################

source("11_run_model.R")

###########################################################
### Post processing
###########################################################

source("12_post_process.R")
