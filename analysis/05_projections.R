#05_projections.R
#-------------------------------------------------------------------------# 
# This script generates epidemiological projections from the the Mexican  #
# Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO).                  #
#                                                                         # 
# Authors:                                                                #
#     - Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>         #
#     - Andrea Luviano, MD, MPH                                           #
#     - Jeremy D Goldhaber-Fiebert, PhD                                   #
#     - Valeria Gracia Olvera, MsC, <valeria.gracia@cide.edu>             #
#-------------------------------------------------------------------------#

rm(list = ls()) # to clean the workspace
devtools::load_all(".")

# Load packages and functions --------------------------------------------

# Packages
library(tidyverse)
library(lubridate)
library(dampack)
library(broom)
library(sccosmomcma)
library(doParallel)

# Functions
source("R/03_target_functions.R")
source("R/03_calibration_functions.R")
source("R/04_validation_functions.R")
source("R/05_projections_functions.R")


# Start projections -------------------------------------------------------

# Time stamp
n_time_stamp <- as.Date("2020-12-13")

for(n_proj_type in c("SQ", "SA")){ # n_proj_type = "SQ"
  
  # Load data and inputs ----------------------------------------------------
  
  # State-specific data
  load(paste0("data/MCMA_calibration_data_",n_proj_type,"_",n_time_stamp,".RData"))
  
  # Load hospitalization data and calibrated parameters
  load(paste0("data/hosp_prop_data_",n_proj_type,"_",n_time_stamp,".RData"))
  load(paste0("output/03_map_output_hosp_NM_MCMA_",n_proj_type,"_",n_time_stamp,".Rdata"))
  
  # Load state specific time varying CFR
  load(paste0("data/m_cfr_",n_proj_type,".RData"))
  
  # Calibrated parameters
  load(paste0("./output/03_map_output_IMIS_", abbrev_state,"_",n_proj_type,"_",n_time_stamp,".RData"))
  load(paste0("output/03_summary_posterior_IMIS_",abbrev_state,"_",n_proj_type,"_",n_time_stamp,".RData"))
  
  ### Posterior distribution of calibrated parameters
  v_calib_post <- df_calib_post_map_IMIS[, "map"]
  names(v_calib_post) <- rownames(df_calib_post_map_IMIS)
  
  
  # Projection parameters ---------------------------------------------------
  
  # Number of age groups
  n_ages <- l_params_all$n_ages
  
  # First day of calibration
  n_date_ini <- first(l_targets$cases$Date)
  
  # Date end calibration
  n_date_end_calib <- last(l_targets$cases$Date)
  
  # Number of days to project for selected state
  n_t_calib <- as.numeric(n_date_end_calib - n_date_ini)
  
  # Vector with NPI times
  v_n_date_NPI   <- c(n_date_NPI, n_date_NPI_2, n_date_NPI_3, n_date_NPI_4, n_date_NPI_5)
  v_n_date0_NPI  <- as.numeric(v_n_date_NPI - n_date_ini) 
  v_n_offset_NPI <- c(0, diff(v_n_date0_NPI))
  
  # Select state to generate projections
  v_states_project <- state_i
  
  # Number of days to project for selected state
  n_date_end_project <- as.Date("2021-03-07")
  n_t_project <- as.numeric(n_date_end_project - n_date_ini) #+ n_lag_inf
  
  
  # Run projections ---------------------------------------------------------
  
  # Set length of last value to be projected
  n_t_proj_calib <- n_t_project - n_t_calib
  
  # Time varying CFR
  m_cfr_proj <- cbind(m_cfr,
                      matrix(rep(m_cfr[,dim(m_cfr)[2]], n_t_proj_calib), nrow = n_ages))
  
  # Hospitalization proportions
  m_p_nonicu_hosp_proj <- cbind(matrix(rep(0,n_lag_inf*n_ages), nrow = n_ages), 
                                m_p_nonicu_hosp,
                                matrix(rep(m_p_nonicu_hosp[,dim(m_p_nonicu_hosp)[2]], n_t_proj_calib), nrow = n_ages))
  
  m_p_icu_hosp_proj <- cbind(matrix(rep(0,n_lag_inf*n_ages), nrow = n_ages),
                             m_p_icu_hosp,
                             matrix(rep(m_p_icu_hosp[,dim(m_p_icu_hosp)[2]], n_t_proj_calib), nrow = n_ages))
  
  m_p_tot_hosp_proj <- cbind(matrix(rep(0,n_lag_inf*n_ages), nrow = n_ages),
                             m_p_tot_hosp,
                             matrix(rep(m_p_tot_hosp[,dim(m_p_tot_hosp)[2]], n_t_proj_calib), nrow = n_ages))
  
  ## Run probabilistic projections ------------------------------------------
  
  l_out_mex_total_prob <- project_interventions_probabilistic(m_calib_post  = m_calib_post,
                                                              n_date_ini    = n_date_ini, 
                                                              v_n_date0_NPI = v_n_date0_NPI,
                                                              n_t_calib     = n_t_calib,
                                                              n_t_project   = n_t_project, 
                                                              n_lag_inf     = n_lag_inf)
  
  
  df_out_mex_total_prob <- l_out_mex_total_prob$df_summ
  df_out_mex_total_prob_all <- l_out_mex_total_prob$df_all
  
  # Save all probabilistic projections
  df_out_mex_total_prob_all <- df_out_mex_total_prob_all %>%
    mutate(proj_type = n_proj_type)
  
  save(df_out_mex_total_prob_all,
       file = paste0("output/05_projections_probabilistic_all_",n_proj_type,"_",abbrev_state,"_", n_time_stamp, ".RData"))
  
  rm(df_out_mex_total_prob_all)
  
  # Wrangle projections
  df_out_mex_total_prob <- df_out_mex_total_prob %>%
    mutate(state = state_i,
           s_code = abbrev_state,
           time_stamp = n_time_stamp) %>%
    ungroup()
  
  df_out_mex_total_prob$Intervention <- as.factor(df_out_mex_total_prob$Intervention)
  df_out_mex_total_prob$date_NPI     <-  n_date_NPI
  
  # Add state name in spanish
  df_out_mex_total_prob$state_esp <-  "ZMVM"
  
  df_out_mex_total_prob <- df_out_mex_total_prob %>%
    mutate(proj_type = n_proj_type)
  
  
  save(df_out_mex_total_prob, 
       file = paste0("output/05_projections_probabilistic_",n_proj_type,"_",abbrev_state,"_", n_time_stamp, ".RData"))
  
  rm(df_out_mex_total_prob)
  
}



