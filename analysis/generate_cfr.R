# 05_generate_cfr.R
#------------------------------------------------------------------------------#
# This script generates MCMA's Covid-19 time varying case fatality rates by    #
# age group.                                                                   #
#                                                                              #
# Author:                                                                      #
#     Valeria Gracia Olvera, <valeria.gracia@cide.edu>                         #
#                                                                              #
#------------------------------------------------------------------------------#
rm(list = ls())

# Load libraries and functions -------------------------------------------

# Libraries
library(survival)   # To compute mortality risk
library(splines)    # To compute mortality risk
library(tidyverse)
library(data.table)

# Functions
source("R/01_cumulative_mortality_risk.R")


# Relative risk of death  -------------------------------------------------

# Relative risk of death: from SSA data
l_rr_death <- list(SQ = 1,
                   SA = 0.06524525)

# State-specific CFR ------------------------------------------------------

for(n_proj_type in c("SQ", "SA")){ # n_proj_type <- "SQ"
  
  # Load state specific data
  n_time_stamp <- "2020-12-13"
  load(paste0("data/MCMA_calibration_data_",n_proj_type,"_",n_time_stamp,".RData"))
  
  # Get Mortality Risk (Author: Jorge Roa)
  # df_MortRisk_hat_raw <- get_cum_mort_risk(save_data = TRUE)
  
  df_MortRisk_hat_raw <- data.table::fread(paste0("data/MortalityRisk_MCMA_",n_time_stamp,".csv"))
  
  # Plot
  # plot_mort_risk(df_MortRisk_hat_raw, n_MR_day = 30, save_plot = FALSE)
  
  # Load population size from SC-COSMO
  load(paste0("data/df_pop_size_",n_proj_type,".rda"))
  
  # Generate matrix of estimated time-varying CFR
    n_death_delay <- 9 # based on mean time-to-death: mean(ssa_covid_surv_pos$time_dx[!is.na(ssa_covid_surv_pos$date_dead)])
  
  # Proportion of adults: from model
  df_MortRisk_hat <- df_MortRisk_hat_raw %>%
    filter(date <= n_date_end_calib &
             date >= n_date_ini &
             MR == "Mortality Risk Day 30" &
             type == "Estimated") %>%
    dplyr::select(-lb, -ub, -type) %>%
    mutate(p_death_adul = value/(df_pop_size$prop_adults +
                                   ((1-df_pop_size$prop_adults)*l_rr_death[[n_proj_type]]
                                   )),
           p_death_kids = p_death_adul*l_rr_death[[n_proj_type]],
           p_death = df_pop_size$prop_adults*p_death_adul +
             (1 - df_pop_size$prop_adults)*p_death_kids)    # Must be equal to value
  
  df_MortRisk_hat <- df_MortRisk_hat[order(df_MortRisk_hat$date),]
  
  # Manual calibration: from 11/01 forward implements a 90% reduction on MR
  n_date_last_cfr <- "2020-11-01"
  n_day_last_cfr <- df_MortRisk_hat$day_dx[which(df_MortRisk_hat$date == n_date_last_cfr)]
  red_serie <- seq(1,0.9,
                   length.out = last(df_MortRisk_hat$day_dx) + 1 - n_day_last_cfr)
  
  # Generate matrix of estimated time-varying CFR
  if(n_proj_type == "SA"){
    n_mult_deaths <- 0.9   # multiplier for the series
  }else{
    n_mult_deaths <- 1.05  # multiplier for the series
  }
  
  m_cfr <- cbind(matrix(rep(0, n_lag_inf*8), ncol = n_lag_inf, nrow = 8),
                 rbind(matrix(rep(df_MortRisk_hat$p_death_kids[1:n_day_last_cfr], each = 2),
                              nrow = 2, byrow=F),
                       matrix(rep(df_MortRisk_hat$p_death_adul[1:n_day_last_cfr], each = 6),
                              nrow = 6, byrow=F)),
                 rbind(matrix(rep(df_MortRisk_hat$p_death_kids[which(df_MortRisk_hat$date == n_date_last_cfr)]*red_serie, each = 2),
                              nrow = 2, byrow = F),
                       matrix(rep(df_MortRisk_hat$p_death_adul[which(df_MortRisk_hat$date == n_date_last_cfr)]*red_serie, each = 6),
                              nrow = 6, byrow = F)))*n_mult_deaths
  
  save(m_cfr,
       n_death_delay,
       file = paste0("data/m_cfr_",n_proj_type,".RData"))
  
}

