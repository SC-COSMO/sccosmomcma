#03_hosp_calibration.R
#------------------------------------------------------------------------------#
# This script calibrates length of stay mean and standard deviation to hospi-  #
# talization prevalence in Mexico City Metropolitan Area.                      #
#                                                                              #
# Author:                                                                      #
#     Valeria Gracia Olvera, MSc, valeria.gracia@cide.edu                      #
#                                                                              #
#------------------------------------------------------------------------------#

rm(list=ls())


# Load libraries and functions --------------------------------------------

# Libraries
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(dampack)
library(epitools)
library(readxl)
library(sccosmomcma)
library(tidyr)
library(lubridate)
library(betareg)

library(foreach)
library(doParallel)

# Functions
source("R/03_hosp_target_functions.R")
source("R/03_hosp_calibration_functions.R")
source("R/03_hosp_prop_functions.R")

# Relative risk of hospitalizations ---------------------------------------

# Relative risk of hospitalization: from SSA data
l_rr_hosp <- list(SQ = 1,
                  SA = 0.3625)

# Relative risk of needing ventilator when hospitalized: from SSA data
l_rr_hosp_vent <- list(SQ = 1,
                       SA = 0.37)

# Relative risk of not needing ventilator when hospitalized: from SSA data
l_rr_hosp_novent <- list(SQ = 1,
                         SA = 0.63)

# Start hosp calibration for each type -----------------------------------

for(n_proj_type in c("SQ", "SA")){ # n_proj_type = "SQ"
  
  # Calibration tye and log file functions
  GLOBAL_CALIB_TYPE      = "NM"
  GLOBAL_LOGGING_ENABLED = TRUE
  
  # Save data
  save_data_hosp <- TRUE
  
  get_log_file_name <- function() {
    paste0("temp/log_hosp_calib_mx_",GLOBAL_CALIB_TYPE,"_",n_proj_type,".txt")
  }
  
  init_log_file <- function(log_flag = TRUE) {
    if (log_flag == TRUE) {
      log_file_name <- get_log_file_name()
      if (file.exists(log_file_name)) {
        #Delete file if it exists
        file.remove(log_file_name)
      }
      writeLines(c(""), log_file_name)
    }
  }
  
  write_log_file <- function(msg, log_flag = TRUE) {
    if (log_flag == TRUE) {
      cat(msg, file = get_log_file_name(), append = TRUE)
    } else {
      print(msg)
    }
  }
  
  # Initialize log file
  init_log_file(log_flag = GLOBAL_LOGGING_ENABLED)
  
  # Load data and functions -------------------------------------------------
  
  # State-specific data
  load(paste0("data/MCMA_calibration_data_",n_proj_type,"_",n_time_stamp,".RData"))
  
  # Load state specific time varying CFR
  load(paste0("data/m_cfr_",n_proj_type,".RData"))
  
  # Load population size from SC-COSMO
  load(paste0("data/df_pop_size_",n_proj_type,".rda"))
  
  # Load m_params_calib and v_map
  load(paste0("output/03_map_output_IMIS_MCMA_",n_proj_type,"_","2020-12-13.RData"))
  v_map <- colMeans(m_calib_post)

  # Generate targets --------------------------------------------------------
  
  l_hosp_targets <- gen_hosp_targets(n_time_stamp = "2020-12-21",
                                     n_date_last  = "2020-12-07"#n_date_last
                                     )
  
  # # Plot targets
  # plot_hosp_targets(l_hosp_targets, save_plot = F)
  
  ## List of number of observations for each target
  l_dates_hosp_targets <- list(hosp   = c(first(l_hosp_targets$hosp$Date),
                                          last(l_hosp_targets$hosp$Date)),
                               vent   = c(first(l_hosp_targets$vent$Date),
                                          last(l_hosp_targets$vent$Date)),
                               novent = c(first(l_hosp_targets$novent$Date),
                                          last(l_hosp_targets$novent$Date))
                               
  )
  
  
  # Set calibration parameters ----------------------------------------------
  
  # Parameters to calibrate
  v_params_calib <- c(m_r_exit_tot    = 0.32,
                      m_r_exit_nonicu = 0.32,
                      m_r_exit_icu    = 0.13,
                      m_sigma_tot     = 2,
                      m_sigma_nonicu  = 2,
                      m_sigma_icu     = 1.5)
  
  v_param_names <- c("m_r_exit_tot",
                     "m_r_exit_nonicu",
                     "m_r_exit_icu",
                     "m_sigma_tot",
                     "m_sigma_nonicu",
                     "m_sigma_icu")
  
  # First day of calibration
  n_date_ini_hosp <- first(as.Date(l_hosp_targets$hosp$Date))
  
  # Last day of calibration
  n_date_end_hosp <- last(as.Date(l_hosp_targets$hosp$Date))
  
  # Total duration of calibration
  n_t_hosp_calib <- as.numeric(n_date_end_hosp - n_date_ini_hosp)
  
  # Define Status quo interventions
  i1 <- make_intervention(intervention_type = "StatusQuo",
                          time_start = 0,
                          time_stop  = n_date0_NPI + n_lag_inf)
  i2 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = n_date0_NPI + n_lag_inf,
                          time_stop  = n_date0_NPI + n_lag_inf + n_offset_NPI_2,
                          intervention_factor = v_map["r_soc_dist_factor"],
                          intervention_change_rate = 0.5)
  i3 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = n_date0_NPI + n_lag_inf + n_offset_NPI_2,
                          time_stop  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3,
                          intervention_factor = v_map["r_soc_dist_factor_2"],
                          intervention_change_rate = 0.5)
  i4 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3,
                          time_stop   = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + n_offset_NPI_4,
                          intervention_factor = v_map["r_soc_dist_factor_3"],
                          intervention_change_rate = 0.5)
  i5 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + n_offset_NPI_4,
                          time_stop   = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + n_offset_NPI_4 + n_offset_NPI_5,
                          intervention_factor = v_map["r_soc_dist_factor_4"],
                          intervention_change_rate = 0.5)
  i6 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + n_offset_NPI_4 + n_offset_NPI_5,
                          time_stop   = n_t_calib + n_lag_inf + 1,
                          intervention_factor = v_map["r_soc_dist_factor_5"],
                          intervention_change_rate = 0.5)
  
  ### Add interventions
  l_interventions <- add_intervention(interventions = NULL, intervention = i1)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
  
  ### Load parameters
  l_params_init <- sccosmomcma::load_params_init(
    n_t              = n_t_calib + n_lag_inf,  # Number of days
    ctry             = "Mexico",
    ste              = state_i,
    cty              = state_i, 
    v_reduced_sus    = v_reduced_sus, 
    r_beta           = v_map["r_beta"], 
    l_nu_exp2_dx     = add_period(l_period_def = NULL,
                                  l_period_add = make_period(
                                    functional_form = "general logit",
                                    time_start = 0,
                                    time_stop = n_t_calib + n_lag_inf,
                                    val_start = as.numeric(v_map["r_nu_exp2_dx_lb"]),
                                    val_end   = as.numeric(v_map["r_nu_exp2_dx_ub"]),
                                    v_logit_change_rate = as.numeric(v_map["r_nu_exp2_dx_rate"]),
                                    v_logit_change_mid  = as.numeric(v_map["n_nu_exp2_dx_mid"]))),
    l_nu_inf2_dx     = add_period(l_period_def = NULL,
                                  l_period_add = make_period(
                                    functional_form = "general logit",
                                    time_start = 0,
                                    time_stop = n_t_calib + n_lag_inf,
                                    val_start = as.numeric(v_map["r_nu_exp2_dx_lb"]),
                                    val_end   = as.numeric(v_map["r_nu_exp2_dx_ub"]),
                                    v_logit_change_rate = as.numeric(v_map["r_nu_exp2_dx_rate"]),
                                    v_logit_change_mid  = as.numeric(v_map["n_nu_exp2_dx_mid"]))),
    v_inf_init_ages  = v_inf_init_ages,
    l_contact_info   = l_contact_matrices,
    l_interventions  = l_interventions,
    n_hhsize         = n_hhsize,
    r_tau            = v_map["r_tau"],
    r_omega          = 0, #1/200
    l_cfr            = get_non_const_multiage_list(v_time_stop = 1:(n_t_calib+n_lag_inf), m_ageval = m_cfr)
    # m_r_exit_tot    = 0.32,
    # m_r_exit_icu    = 0.13,
    # m_r_exit_nonicu = 0.32,
    # m_sigma_tot     = 2,
    # m_sigma_nonicu  = 2,
    # m_sigma_icu     = 1.5
  ) 
  
  # Load all parameters values
  l_params_all <- sccosmomcma::load_all_params(l_params_init = l_params_init)
  
  
  # Hospitalization proportions ---------------------------------------------
  
  # l_hosp_prop <- get_MCMA_prop_hosp(n_date_ini  = n_date_ini_hosp,
  #                                   n_date_last = n_date_last)
  
  load("data/l_hosp_prop_2020-12-13.RData")
  
  
  ## Total hospitalizations -------------------------------------------------
  
  df_hosp_prop <- l_hosp_prop$ssa
  
  # Beta regression
  df_p_hosp <- df_hosp_prop %>%
    filter(p_hosp > 0) %>%
    select(date, date0, p_hosp)
  
  fit_hosp_linear <- betareg(p_hosp ~ date0, data = df_p_hosp)
  fit_hosp_logit  <- betareg(p_hosp ~ splines::ns(date0, 4), data = df_p_hosp)
  fit_hosp_loglog <- betareg(p_hosp ~ splines::ns(date0, 4), data = df_p_hosp, link = "loglog")
  df_AIC <- AIC(fit_hosp_linear, fit_hosp_logit, fit_hosp_loglog)
  
  best_model <- rownames(df_AIC[which(df_AIC$AIC == min(df_AIC$AIC)),])
  
  if(best_model == "fit_hosp_linear"){
    ## Linear model
    df_p_hosp$hosp_hat  <- predict(fit_hosp_linear, type = "response")
    df_p_hosp$hosp_hat_lb  <- predict(fit_hosp_linear, type = "quantile", at = 0.025)
    df_p_hosp$hosp_hat_ub  <- predict(fit_hosp_linear, type = "quantile", at = 0.975)
  } else if(best_model == "fit_hosp_logit"){
    ## Logit model
    df_p_hosp$hosp_hat  <- predict(fit_hosp_logit, type = "response")
    df_p_hosp$hosp_hat_lb  <- predict(fit_hosp_logit, type = "quantile", at = 0.025)
    df_p_hosp$hosp_hat_ub  <- predict(fit_hosp_logit, type = "quantile", at = 0.975)
  } else if(best_model == "fit_hosp_loglog"){
    ## Log-log model
    df_p_hosp$hosp_hat  <- predict(fit_hosp_loglog, type = "response")
    df_p_hosp$hosp_hat_lb  <- predict(fit_hosp_loglog, type = "quantile", at = 0.025)
    df_p_hosp$hosp_hat_ub  <- predict(fit_hosp_loglog, type = "quantile", at = 0.975)
  }
  # Data.frame to plot observed vd estimated
  df_plot_phosp <- rbind.data.frame(data.frame(Outcome = "Total hospitalizations",
                                               state   = "MCMA",
                                               date  = df_hosp_prop$date,
                                               date0 = df_hosp_prop$date0,
                                               prop  = df_hosp_prop$p_hosp,
                                               type  = "Observed",
                                               lb    = NA,
                                               ub    = NA),
                                    data.frame(Outcome = "Total hospitalizations",
                                               state   = "MCMA",
                                               date  = df_p_hosp$date,
                                               date0 = df_p_hosp$date0,
                                               prop  = df_p_hosp$hosp_hat,
                                               type  = "Estimated",
                                               lb    = df_p_hosp$hosp_hat_lb,
                                               ub    = df_p_hosp$hosp_hat_ub))
  
  
  # Plot hosp proportions
  # plot_hosp_prop(df_plot_phosp, save_plot = FALSE)
  
  # Save data.frame
  # df_p_hosp_adjusted <- subset(df_plot_phosp, date>="2020-03-24")
  # save(df_p_hosp_adjusted, file = "data/df_p_hosp_adjusted.Rdata")
  
  # Generate matrix of estimated time-varying p_hosp with respective
  # reduction on hospitalizations
  
  df_p_hosp_mod <- df_p_hosp %>%
    complete(date = seq.Date(from = as.Date(n_date_ini),
                             to   = as.Date(n_date_end_hosp),
                             by   = "day"),
             fill = list(date0 = NA,
                         p_hosp = 0,
                         p_hosp_adul = 0,
                         hosp_hat = 0,
                         hosp_hat_lb = 0,
                         hosp_hat_ub = 0
             )) %>%
    mutate(p_hosp_adul = hosp_hat/(df_pop_size$prop_adults + 
                                     ((1-df_pop_size$prop_adults)*l_rr_hosp[[n_proj_type]])),
           p_hosp_kids = p_hosp_adul*l_rr_hosp[[n_proj_type]],
           p_hosp_tot = df_pop_size$prop_adults*p_hosp_adul +
             (1-df_pop_size$prop_adults)*p_hosp_kids,   # Must be equal to hosp_hat
           date0 = as.numeric(date - date[1])) %>%
    filter(date <= n_date_end_calib)
  
  # Order data.frame
  df_p_hosp_mod <- df_p_hosp_mod[order(df_p_hosp_mod$date),]
  
  # Matrix of hospitalization proportions
  m_p_tot_hosp <- cbind(#matrix(rep(0,n_lag_inf*n_ages), nrow = 8, byrow = F),
                        #matrix(rep(0,(n_t_calib - n_t_hosp_calib)*n_ages), nrow = 8, byrow = F),
                        rbind(matrix(rep(df_p_hosp_mod$p_hosp_kids,each = 2), nrow = 2, byrow=F),
                              matrix(rep(df_p_hosp_mod$p_hosp_adul,each = 6), nrow = 6, byrow=F)))
  
  
  ## Hospitalizations without ventilator ---------------------------------------
  
  df_novent_hosp_prop <- l_hosp_prop$adip
  
  df_novent_hosp_prop$date <- as.Date(df_novent_hosp_prop$date)
  
  # Beta regression
  df_p_novent_hosp <- df_novent_hosp_prop %>%
    filter(p_novent_hosp > 0) %>%
    select(date, date0, p_novent_hosp)
  
  fit_hosp_novent_linear <- betareg(p_novent_hosp ~ date0, data = df_p_novent_hosp)
  fit_hosp_novent_logit  <- betareg(p_novent_hosp ~ splines::ns(date0, 4), data = df_p_novent_hosp)
  fit_hosp_novent_loglog <- betareg(p_novent_hosp ~ splines::ns(date0, 4), data = df_p_novent_hosp, link = "loglog")
  
  df_AIC <- AIC(fit_hosp_novent_linear, fit_hosp_novent_logit, fit_hosp_novent_loglog)
  
  best_model <- rownames(df_AIC[which(df_AIC$AIC == min(df_AIC$AIC)),])
  
  if(best_model == "fit_hosp_novent_linear"){
    ## Linear model
    df_p_novent_hosp$novent_hat  <- predict(fit_hosp_novent_linear, type = "response")
    df_p_novent_hosp$novent_hat_lb  <- predict(fit_hosp_novent_linear, type = "quantile", at = 0.025)
    df_p_novent_hosp$novent_hat_ub  <- predict(fit_hosp_novent_linear, type = "quantile", at = 0.975)
  } else if(best_model == "fit_hosp_novent_logit"){
    ## Logit model
    df_p_novent_hosp$novent_hat  <- predict(fit_hosp_novent_logit, type = "response")
    df_p_novent_hosp$novent_hat_lb  <- predict(fit_hosp_novent_logit, type = "quantile", at = 0.025)
    df_p_novent_hosp$novent_hat_ub  <- predict(fit_hosp_novent_logit, type = "quantile", at = 0.975)
  } else if(best_model == "fit_hosp_novent_loglog"){
    ## Log-log model
    df_p_novent_hosp$novent_hat  <- predict(fit_hosp_novent_loglog, type = "response")
    df_p_novent_hosp$novent_hat_lb  <- predict(fit_hosp_novent_loglog, type = "quantile", at = 0.025)
    df_p_novent_hosp$novent_hat_ub  <- predict(fit_hosp_novent_loglog, type = "quantile", at = 0.975)
  }
  
  # Data.frame to plot observed vs estimated
  df_plot_novent_phosp <- rbind.data.frame(data.frame(Outcome = "Hospitalizations without ventilator",
                                                      state   = "MCMA",
                                                      date   = df_novent_hosp_prop$date,
                                                      date0  = df_novent_hosp_prop$date0,
                                                      prop   = df_novent_hosp_prop$p_novent_hosp,
                                                      type   = "Observed",
                                                      lb     = NA,
                                                      ub     = NA),
                                           data.frame(Outcome = "Hospitalizations without ventilator",
                                                      state   = "MCMA",
                                                      date   = df_p_novent_hosp$date,
                                                      date0  = df_p_novent_hosp$date0,
                                                      prop   = df_p_novent_hosp$novent_hat,
                                                      type   = "Estimated",
                                                      lb     = df_p_novent_hosp$novent_hat_lb,
                                                      ub     = df_p_novent_hosp$novent_hat_ub))
  
  # Plot hosp proportions
  # plot_hosp_prop(df_plot_novent_phosp, save_plot = FALSE)
  
  # Save data.frame
  # df_p_novent_adjusted <- df_plot_novent_phosp
  # save(df_p_novent_adjusted, file = "data/df_p_novent_adjusted.Rdata")
  
  # Generate matrix of estimated time-varying p_hosp with respective
  # reduction on hospitalizations
  
  df_p_novent_hosp_mod <- df_p_novent_hosp %>%
    complete(date = seq.Date(from = as.Date(n_date_ini),
                             to   = as.Date(n_date_end_hosp),
                             by   = "day"),
             fill = list(date0 = NA,
                         p_novent_hosp = 0,
                         p_hosp_novent_adul = 0,
                         novent_hat = 0,
                         novent_hat_lb = 0,
                         novent_hat_ub = 0
             )) %>%
    mutate(p_hosp_novent_adul = novent_hat/(df_pop_size$prop_adults + 
                                              ((1-df_pop_size$prop_adults)*l_rr_hosp_novent[[n_proj_type]])),
           p_hosp_novent_kids = p_hosp_novent_adul*l_rr_hosp_novent[[n_proj_type]],
           p_hosp_novent = df_pop_size$prop_adults*p_hosp_novent_adul +
             (1 - df_pop_size$prop_adults)*p_hosp_novent_kids,
           date0 = as.numeric(date - date[1]))
  
  # Order data.frame
  df_p_novent_hosp_mod <- df_p_novent_hosp_mod[order(df_p_novent_hosp_mod$date),]
  
  # Matrix
  m_p_nonicu_hosp <-  cbind(#matrix(rep(0,n_lag_inf*n_ages), nrow = 8, byrow = F),
                           #matrix(rep(0,(n_t_calib - n_t_hosp_calib)*n_ages), nrow = 8, byrow = F),
                           rbind(matrix(rep(df_p_novent_hosp_mod$p_hosp_novent_kids,each = 2), nrow = 2, byrow=F),
                                 matrix(rep(df_p_novent_hosp_mod$p_hosp_novent_adul,each = 6), nrow = 6, byrow=F)))
  
  
  ## Hospitalizations with ventilator ---------------------------------------
  
  df_vent_hosp_prop <- l_hosp_prop$adip
  
  df_vent_hosp_prop$date <- as.Date(df_vent_hosp_prop$date)
  
  # Beta regression
  df_p_vent_hosp <- df_vent_hosp_prop %>%
    filter(p_vent_hosp > 0) %>%
    select(date, date0, p_vent_hosp)
  
  fit_hosp_vent_linear <- betareg(p_vent_hosp ~ date0, data = df_p_vent_hosp)
  fit_hosp_vent_logit  <- betareg(p_vent_hosp ~ splines::ns(date0, 4), data = df_p_vent_hosp)
  fit_hosp_vent_loglog <- betareg(p_vent_hosp ~ splines::ns(date0, 4), data = df_p_vent_hosp, link = "loglog")
  
  df_AIC <- AIC(fit_hosp_vent_linear, fit_hosp_vent_logit, fit_hosp_vent_loglog)
  
  best_model <- rownames(df_AIC[which(df_AIC$AIC == min(df_AIC$AIC)),])
  
  if(best_model == "fit_hosp_vent_linear"){
    ## Linear model
    df_p_vent_hosp$vent_hat  <- predict(fit_hosp_vent_linear, type = "response")
    df_p_vent_hosp$vent_hat_lb  <- predict(fit_hosp_vent_linear, type = "quantile", at = 0.025)
    df_p_vent_hosp$vent_hat_ub  <- predict(fit_hosp_vent_linear, type = "quantile", at = 0.975)
  } else if(best_model == "fit_hosp_vent_logit"){
    ## Logit model
    df_p_vent_hosp$vent_hat  <- predict(fit_hosp_vent_logit, type = "response")
    df_p_vent_hosp$vent_hat_lb  <- predict(fit_hosp_vent_logit, type = "quantile", at = 0.025)
    df_p_vent_hosp$vent_hat_ub  <- predict(fit_hosp_vent_logit, type = "quantile", at = 0.975)
  } else if(best_model == "fit_hosp_vent_loglog"){
    ## Log-log model
    df_p_vent_hosp$vent_hat  <- predict(fit_hosp_vent_loglog, type = "response")
    df_p_vent_hosp$vent_hat_lb  <- predict(fit_hosp_vent_loglog, type = "quantile", at = 0.025)
    df_p_vent_hosp$vent_hat_ub  <- predict(fit_hosp_vent_loglog, type = "quantile", at = 0.975)
  }
  
  # Data.frame to plot observed vs estimated
  df_plot_vent_phosp <- rbind.data.frame(data.frame(Outcome = "Hospitalizations with ventilator",
                                                    state   = "MCMA",
                                                    date   = df_vent_hosp_prop$date,
                                                    date0  = df_vent_hosp_prop$date0,
                                                    prop   = df_vent_hosp_prop$p_vent_hosp,
                                                    type   = "Observed",
                                                    lb     = NA,
                                                    ub     = NA),
                                         data.frame(Outcome = "Hospitalizations with ventilator",
                                                    state   = "MCMA",
                                                    date   = df_p_vent_hosp$date,
                                                    date0  = df_p_vent_hosp$date0,
                                                    prop   = df_p_vent_hosp$vent_hat,
                                                    type   = "Estimated",
                                                    lb     = df_p_vent_hosp$vent_hat_lb,
                                                    ub     = df_p_vent_hosp$vent_hat_ub))
  
  # Plot hosp proportions
  # plot_hosp_prop(df_plot_vent_phosp, save_plot = FALSE)
  
  # Save data.frame
  # df_p_vent_adjusted <- df_plot_vent_phosp
  # save(df_p_vent_adjusted, file = "data/df_p_vent_adjusted.Rdata")
  
  
  # Generate matrix of estimated time-varying p_hosp with respective
  # reduction on hospitalizations
  
  df_p_vent_hosp_mod <- df_p_vent_hosp %>%
    complete(date = seq.Date(from = as.Date(n_date_ini),
                             to   = as.Date(n_date_end_hosp),
                             by   = "day"),
             fill = list(date0 = NA,
                         p_vent_hosp = 0,
                         p_hosp_vent_adul = 0,
                         vent_hat = 0,
                         vent_hat_lb = 0,
                         vent_hat_ub = 0
             )) %>%
    mutate(p_hosp_vent_adul = vent_hat/(df_pop_size$prop_adults + 
                                          ((1-df_pop_size$prop_adults)*l_rr_hosp_vent[[n_proj_type]])),
           p_hosp_vent_kids = p_hosp_vent_adul*l_rr_hosp_vent[[n_proj_type]],
           p_hosp_vent = df_pop_size$prop_adults*p_hosp_vent_adul +
             (1 - df_pop_size$prop_adults)*p_hosp_vent_kids,
           date0 = as.numeric(date - date[1]))
  
  # Order data.frame
  df_p_vent_hosp_mod <- df_p_vent_hosp_mod[order(df_p_vent_hosp_mod$date),]
  
  # Matrix
  m_p_icu_hosp <-  cbind(#matrix(rep(0,n_lag_inf*n_ages), nrow = 8, byrow = F),
                         #matrix(rep(0,(n_t_calib - n_t_hosp_calib)*n_ages), nrow = 8, byrow = F),
                         rbind(matrix(rep(df_p_vent_hosp_mod$p_hosp_vent_kids,each = 2), nrow = 2, byrow=F),
                               matrix(rep(df_p_vent_hosp_mod$p_hosp_vent_adul,each = 6), nrow = 6, byrow=F)))
  
  
  # Save hospitalization data ----------------------------------------------
  
  if(save_data_hosp){
    save(m_p_tot_hosp,
         m_p_nonicu_hosp,
         m_p_icu_hosp,
         n_date_ini_hosp,
         n_date_end_hosp,
         n_t_hosp_calib,
         l_hosp_targets,
         file = paste0("data/hosp_prop_data_",n_proj_type,"_",n_time_stamp,".RData"))
  }
  
  
  # Test calibration functions ----------------------------------------------
  
  # source("R/03_hosp_calibration_functions.R")
  # 
  # l_res <- hosp_calibration_out(v_params_calib       = v_params_calib,
  #                               l_params_all         = l_params_all,
  #                               n_lag_inf            = 14,
  #                               n_lag_conf           = 0,
  #                               l_dates_hosp_targets = l_dates_hosp_targets)
  # 
  # hosp_log_lik(v_params             = v_params_calib,
  #              l_params_all         = l_params_all,
  #              n_lag_inf            = 14,
  #              n_lag_conf           = 0,
  #              l_dates_hosp_targets = l_dates_hosp_targets)
  # 
  # hosp_log_lik_opt(v_params             = v_params_calib,
  #                  l_params_all         = l_params_all,
  #                  n_lag_inf            = 14,
  #                  n_lag_conf           = 0,
  #                  l_dates_hosp_targets = l_dates_hosp_targets)
  # 
  # hosp_log_lik_par(v_params             = v_params_calib,
  #                  l_params_all         = l_params_all,
  #                  n_lag_inf            = 14,
  #                  n_lag_conf           = 0,
  #                  l_dates_hosp_targets = l_dates_hosp_targets,
  #                  log_lik_offset       = 0)
  # 
  # hosp_log_lik(v_params             = hosp_sample_prior(6),
  #              l_params_all         = l_params_all,
  #              n_lag_inf            = 14,
  #              n_lag_conf           = 0,
  #              l_dates_hosp_targets = l_dates_hosp_targets)
  # 
  # hosp_log_lik_opt(v_params             = hosp_sample_prior(6),
  #                  l_params_all         = l_params_all,
  #                  n_lag_inf            = 14,
  #                  n_lag_conf           = 0,
  #                  l_dates_hosp_targets = l_dates_hosp_targets)
  # 
  # hosp_log_prior_opt(v_params = hosp_sample_prior(6))
  # 
  # hosp_log_lik_par(v_params             = hosp_sample_prior(6),
  #                  l_params_all         = l_params_all,
  #                  n_lag_inf            = 14,
  #                  n_lag_conf           = 0,
  #                  l_dates_hosp_targets = l_dates_hosp_targets,
  #                  log_lik_offset       = 0)
  # 
  # hosp_likelihood(v_params_calib)
  # 
  # hosp_likelihood(hosp_sample_prior(5))
  
  # hosp_log_post_opt(v_params             = hosp_sample_prior(6),
  #                   l_params_all         = l_params_all,
  #                   n_lag_inf            = 14,
  #                   n_lag_conf           = 0,
  #                   l_dates_hosp_targets = l_dates_hosp_targets,
  #                   log_lik_offset       = 0)
  
  
  # Find a maximum-a-posteriori: NM parallel --------------------------------
  
  ## Set seed
  set.seed(20200522)
  
  # Number and parameters set to calibrate
  n_param_nm <- 200
  m_params_calib <- hosp_sample_prior(n_samp = n_param_nm)
  
  # Parallel set up
  no_cores <- detectCores() - 2   # detect number of cores
  cl <- makeCluster(no_cores)       # initialize cluster object
  registerDoParallel(cl)
  opts <- list(attachExportEnv = TRUE)
  
  # Start paralleled iterations
  l_out_map_multi <- foreach(parallel_i = 1:n_param_nm, .combine = rbind, .export = ls(globalenv()),
                             .packages=c("sccosmomcma",
                                         "ggplot2",
                                         "tidyverse",
                                         "dplyr",
                                         "lubridate",
                                         "epitools"),
                             .options.snow = opts) %dopar% {
                               rerun <-  FALSE
                               write_log_file(msg=paste(Sys.time(),": +++ INITIATING iteration:",parallel_i,"\n") ,
                                              log_flag=GLOBAL_LOGGING_ENABLED)
                               
                               #t0 <- Sys.time()
                               l_out_map_par <- optim(par                  = m_params_calib[parallel_i, ], 
                                                      fn                   = hosp_log_post_opt, 
                                                      l_params_all         = l_params_all, 
                                                      n_lag_inf            = n_lag_inf, 
                                                      n_lag_conf           = n_lag_conf, 
                                                      l_dates_hosp_targets = l_dates_hosp_targets,
                                                      control              = list(fnscale = -1, 
                                                                                  maxit = 500, 
                                                                                  trace = 3, 
                                                                                  abstol  = 1e-1,
                                                                                  reltol  = 1e-4), 
                                                      hessian              = FALSE)
                               
                               write_log_file(msg=paste(Sys.time(),":  ++ AFTER OPTIMIZATON:",parallel_i,"\n"),
                                              log_flag=GLOBAL_LOGGING_ENABLED)
                               
                               if(l_out_map_par$convergence!=0){
                                 rerun <- TRUE
                                 write_log_file(msg=paste(Sys.time(),":  ++ FIRST OPTIMIZATON DID NOT CONVERGE:",parallel_i,l_out_map_par$convergence,
                                                          cat(msg=l_out_map_par$par, collapse=" "),"\n") ,log_flag=GLOBAL_LOGGING_ENABLED)
                                 
                                 write_log_file(msg=paste(Sys.time(),":  ++ BEFORE 2nd OPTIMIZATON:",parallel_i,"\n"),
                                                log_flag=GLOBAL_LOGGING_ENABLED)
                                 # print(paste(state_i, "did not converge first 500"))
                                 
                                 l_out_map_par <- optim(par                  = m_params_calib[parallel_i, ], 
                                                        fn                   = hosp_log_post_opt, 
                                                        l_params_all         = l_params_all, 
                                                        n_lag_inf            = n_lag_inf, 
                                                        n_lag_conf           = n_lag_conf, 
                                                        l_dates_hosp_targets = l_dates_hosp_targets,
                                                        control              = list(fnscale = -1, 
                                                                                    maxit = 500, 
                                                                                    trace = 3, 
                                                                                    abstol  = 1e-1,
                                                                                    reltol  = 1e-4), 
                                                        hessian              = FALSE)
                                 
                                 write_log_file(msg=paste(Sys.time(),":  ++ AFTER 2nd OPTIMIZATON:",parallel_i,"\n"),
                                                log_flag=GLOBAL_LOGGING_ENABLED)
                                 
                               }
                               
                               
                               #t1 <- Sys.time()
                               df_temp <- data.frame(`Time stamp` = n_time_stamp,
                                                     LastDay      = n_date_last,
                                                     state        = state_i,
                                                     Run          = parallel_i,
                                                     Rerun        = rerun,
                                                     Value        = l_out_map_par$value,
                                                     Convergence  = l_out_map_par$convergence,
                                                     Counts       = l_out_map_par$counts[1],
                                                     Parameter    = v_param_names,
                                                     map          = l_out_map_par$par,
                                                     lb           = NA,
                                                     ub           = NA,
                                                     se           = NA)
                               #WriteXLS(df_temp, 
                               #          paste0("temp/df_results_iter_",parallel_i,"_",abbrev_state,"_",n_time_stamp,".xls"), 
                               #         row.names=TRUE, perl = "C:/Strawberry/perl/bin/perl.exe")  
                               
                               write_log_file(msg=paste(Sys.time(),": +++ FINISHING iteration:",parallel_i,"\n") ,
                                              log_flag=GLOBAL_LOGGING_ENABLED)
                               df_temp
                             }
  
  stopCluster(cl)
  
  max_value <- max(l_out_map_multi$Value)
  v_hosp_map <- l_out_map_multi[l_out_map_multi$Value == max_value, ]$map
  names(v_hosp_map) <- v_param_names
  
  # Matrix with parameters
  df_calib_post_map_NMpar <- data.frame(`Time stamp` = n_time_stamp,
                                        LastDay      = n_date_last,
                                        state        = state_i,
                                        Method       = "NM parallel",
                                        value        = max_value,
                                        map          = v_hosp_map,
                                        lb           = NA,
                                        ub           = NA,
                                        se           = NA)
  rownames(df_calib_post_map_NMpar) <- v_param_names
  
  
  # Store posterior and MAP from NM_par calibration
  save(v_hosp_map,
       l_out_map_multi,
       df_calib_post_map_NMpar,
       n_date_ini_hosp,
       n_date_end_hosp,
       file = paste0("output/03_map_output_hosp_NM_",
                     abbrev_state,"_",n_proj_type,"_",n_time_stamp,".RData"))
  
  
  
  # Validation --------------------------------------------------------------
  
  # Load parameters
  l_params_init <- sccosmomcma::load_params_init(
    n_t              = n_t_calib + n_lag_inf,  # Number of days
    ctry             = "Mexico",
    ste              = state_i,
    cty              = state_i, 
    v_reduced_sus    = v_reduced_sus,
    r_beta           = v_map["r_beta"], #m_calib_post[1, "r_beta"]
    l_nu_exp2_dx     = add_period(l_period_def = NULL,
                                  l_period_add = make_period(
                                    functional_form = "general logit",
                                    time_start = 0,
                                    time_stop = n_t_calib + n_lag_inf,
                                    val_start = as.numeric(v_map["r_nu_exp2_dx_lb"]),
                                    val_end   = as.numeric(v_map["r_nu_exp2_dx_ub"]),
                                    v_logit_change_rate = as.numeric(v_map["r_nu_exp2_dx_rate"]),
                                    v_logit_change_mid  = as.numeric(v_map["n_nu_exp2_dx_mid"]))),
    l_nu_inf2_dx     = add_period(l_period_def = NULL,
                                  l_period_add = make_period(
                                    functional_form = "general logit",
                                    time_start = 0,
                                    time_stop = n_t_calib + n_lag_inf,
                                    val_start = as.numeric(v_map["r_nu_exp2_dx_lb"]),
                                    val_end   = as.numeric(v_map["r_nu_exp2_dx_ub"]),
                                    v_logit_change_rate = as.numeric(v_map["r_nu_exp2_dx_rate"]),
                                    v_logit_change_mid  = as.numeric(v_map["n_nu_exp2_dx_mid"]))),
    v_inf_init_ages  = v_inf_init_ages,
    l_contact_info   = l_contact_matrices,
    l_interventions  = l_interventions,
    n_hhsize         = n_hhsize,
    r_tau            = v_map["r_tau"],
    r_omega          = 0, #1/200
    l_cfr            = get_non_const_multiage_list(v_time_stop = 1:(n_t_calib+n_lag_inf), m_ageval = m_cfr),
    m_r_exit_tot     = v_hosp_map["m_r_exit_tot"],
    m_r_exit_icu     = v_hosp_map["m_r_exit_icu"],
    m_r_exit_nonicu  = v_hosp_map["m_r_exit_nonicu"],
    m_sigma_tot      = v_hosp_map["m_sigma_tot"],
    m_sigma_nonicu   = v_hosp_map["m_sigma_nonicu"],
    m_sigma_icu      = v_hosp_map["m_sigma_icu"]
  ) 
  
  # Load all parameter values
  l_params_all <- sccosmomcma::load_all_params(l_params_init = l_params_init)
  #l_out_cosmo <- cosmo(l_params_all = l_params_all)
  
  
  # Run model-predicted outputs
  df_out_cosmo_map_bc <- hosp_calibration_out(v_params_calib       = v_hosp_map, 
                                              l_params_all         = l_params_all,
                                              n_lag_inf            = n_lag_inf, 
                                              n_lag_conf           = 0,
                                              l_dates_hosp_targets = l_dates_hosp_targets)
  
  df_out_map <- data.frame(type = "Model",
                           series = c(rep(unique(l_hosp_targets$df_all_targets$series)[3], length(df_out_cosmo_map_bc$H_prev$value)),
                                      rep(unique(l_hosp_targets$df_all_targets$series)[2], length(df_out_cosmo_map_bc$H_novent_prev$value)),
                                      rep(unique(l_hosp_targets$df_all_targets$series)[1], length(df_out_cosmo_map_bc$H_vent_prev$value))),
                           Date = c(seq.Date(from = l_dates_hosp_targets$hosp[1], 
                                             to = l_dates_hosp_targets$hosp[2], by = "days"),
                                    seq.Date(from = l_dates_hosp_targets$novent[1], 
                                             to = l_dates_hosp_targets$novent[2], by = "days"),
                                    seq.Date(from = l_dates_hosp_targets$vent[1], 
                                             to = l_dates_hosp_targets$vent[2], by = "days")),
                           Date0 = c(l_hosp_targets$df_all_targets$Date0[l_hosp_targets$df_all_targets$Target=="Total hospitalizations"],
                                     l_hosp_targets$df_all_targets$Date0[l_hosp_targets$df_all_targets$Target=="Hospitalizations without ventilator"],
                                     l_hosp_targets$df_all_targets$Date0[l_hosp_targets$df_all_targets$Target=="Hospitalizations with ventilator"]),
                           value = c(df_out_cosmo_map_bc$H_prev$value, 
                                     df_out_cosmo_map_bc$H_novent_prev$value,
                                     df_out_cosmo_map_bc$H_vent_prev$value),
                           lb = NA,
                           ub = NA 
  )
  
  df_targets_out <- dplyr::bind_rows(l_hosp_targets$hosp[, c("type", "series", "Date", "Date0", "value", "lb", "ub")],
                                     l_hosp_targets$vent[, c("type", "series", "Date", "Date0", "value", "lb", "ub")],
                                     l_hosp_targets$novent[, c("type", "series", "Date", "Date0", "value", "lb", "ub")],
                                     df_out_map)
  
  # Plot target vs projections on calendar days since first confirmed case 
  ggplot(subset(df_targets_out, (type == "Model" & Date <= "2020-11-30") | 
                  (type == "Target" & Date <= "2020-11-30") ),
         aes(x = Date, y = value, 
             ymin = lb, ymax = ub,
             color = type, shape = type)) +
    facet_wrap(~series, scales = "free_y") + 
    geom_point(size = 3) + # for target: shape = 8
    geom_errorbar() +
    scale_shape_manual(values = c(8, 1)) +
    scale_color_manual(values = c("red", "black")) +
    scale_x_date(breaks = number_ticks(8)) +
    scale_y_continuous(breaks = number_ticks(6)) +
    labs(title = paste("Hospital occupancy of COVID-19 in",state_i,"Mexico"),
         x = "Days since first case", y = "Hospital occupancy") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 60, hjust = 1, size = 10))
  
  # Save plot
  ggsave(paste0("figs/04_validation_targets_vs_model_hosp_",GLOBAL_CALIB_TYPE,"_",abbrev_state,"_",
                n_proj_type,"_", n_time_stamp,".pdf"), 
         width = 12, height = 8)
  
}
