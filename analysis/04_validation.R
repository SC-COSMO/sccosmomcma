# 04_validation.R
#------------------------------------------------------------------------------# 
# This script conducts an internal validation of the Mexican Stanford-CIDE     #
# COronavirus Simulation MOdel (SC-COSMO) by comparing the model-predicted     #
# outputs evaluated at the calibrated parameters vs the calibration targets.   #
#                                                                              # 
# Authors:                                                                     #
#     - Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>              #
#     - Andrea Luviano, MD, MPH                                                #
#     - Jeremy D Goldhaber-Fiebert, PhD                                        #                                                       #
#     - Valeria Gracia Olvera, MsC, <valeria.gracia@cide.edu>                  #
#------------------------------------------------------------------------------# 

rm(list = ls()) # to clean the workspace


# Load packages, functions and inputs -------------------------------------

# Packages
library(sccosmomcma)
library(deSolve)
library(tidyverse)
library(lubridate)
library(doParallel)
library(doParallel)  # Parallel setup 
library(foreach)     # Parallel setup 

# Functions
source("R/03_target_functions.R")
source("R/03_calibration_functions.R")
source("R/04_validation_functions.R")


# Start validation --------------------------------------------------------

# Number of age groups
n_ages <- 8

# Reduced susceptibility
reduced_sus <- 0.25 # based on : https://www.nature.com/articles/s41591-020-0962-9 abstract and figure 1 panel B
l_reduced_sus <- list(SQ = rep(1, n_ages), # SQ
                      SA = c(rep(reduced_sus, 2), rep(1, (n_ages-2)))) # SA
                           

for(n_proj_type in c("SQ", "SA")){ # n_proj_type = "SQ"
        
        # Load data ---------------------------------------------------------------

        # Load state specific data
        n_time_stamp <- "2020-12-13"
        load(paste0("data/MCMA_calibration_data_",n_proj_type,"_",n_time_stamp,".RData"))
        
        # Load state specific time varying CFR
        load(paste0("data/m_cfr_",n_proj_type,".RData"))

        # Load calibrated parameters
        files <- list.files(paste0("./output/"), 
                            pattern = paste0("03_map_output_IMIS_",abbrev_state,"_",n_proj_type,"_",n_time_stamp))
        
        load(file=paste0("./output/",files))
        
        v_map <- df_calib_post_map_IMIS$map
        names(v_map) <- v_param_names
        
        ## Load all parameters and initialize model ----------------------------
        
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
                                time_stop   = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 +
                                        n_offset_NPI_4,
                                intervention_factor = v_map["r_soc_dist_factor_3"],
                                intervention_change_rate = 0.5)
        i5 <- make_intervention(intervention_type = "SocialDistancing",
                                time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + 
                                        n_offset_NPI_4,
                                time_stop   = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + 
                                        n_offset_NPI_4 + n_offset_NPI_5,
                                intervention_factor = v_map["r_soc_dist_factor_4"],
                                intervention_change_rate = 0.5)
        i6 <- make_intervention(intervention_type = "SocialDistancing",
                                time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 +
                                        n_offset_NPI_4 + n_offset_NPI_5,
                                time_stop   = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + 
                                        n_offset_NPI_4 + n_offset_NPI_5 + n_offset_NPI_6,
                                intervention_factor = v_map["r_soc_dist_factor_5"],
                                intervention_change_rate = 0.5)
        i7 <- make_intervention(intervention_type = "SocialDistancing",
                                time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + 
                                        n_offset_NPI_4 + n_offset_NPI_5 + n_offset_NPI_6,
                                time_stop   = n_t_calib + n_lag_inf + 1,
                                intervention_factor = v_map["r_soc_dist_factor_6"],
                                intervention_change_rate = 0.5)
        # i8 <- make_intervention(intervention_type = "SocialDistancing",
        #                         time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + 
        #                                 n_offset_NPI_4 + n_offset_NPI_5 + n_offset_NPI_6 + n_offset_NPI_7 + n_offset_NPI_8,
        #                         time_stop   = n_t_calib + n_lag_inf + 1,
        #                         intervention_factor = v_map["r_soc_dist_factor_7"],
        #                         intervention_change_rate = 0.5)
        # i9 <- make_intervention(intervention_type = "SocialDistancing",
        #                         time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + 
        #                                 n_offset_NPI_4 + n_offset_NPI_5 + n_offset_NPI_6 + n_offset_NPI_7 + n_offset_NPI_8,
        #                         time_stop   = n_t_calib + n_lag_inf + 1,
        #                         intervention_factor = v_map["r_soc_dist_factor_8"],
        #                         intervention_change_rate = 0.5)
        
        ### Add interventions
        l_interventions <- add_intervention(interventions = NULL, intervention = i1)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i7)
        # l_interventions <- add_intervention(interventions = l_interventions, intervention = i8)
        # l_interventions <- add_intervention(interventions = l_interventions, intervention = i9)
        
        ### Time varying parameters
        l_params_init <- sccosmomcma::load_params_init(n_t = n_t_calib + n_lag_inf, # Number of days
                                                       ctry = "Mexico",
                                                       ste  = state_i,
                                                       cty  = state_i,
                                                       v_reduced_sus = v_reduced_sus,
                                                       r_beta = v_map["r_beta"],
                                                       l_nu_exp2_dx = add_period(l_period_def = NULL,
                                                                                 l_period_add = make_period(
                                                                                         functional_form = "general logit",
                                                                                         time_start = 0,
                                                                                         time_stop = n_t_calib + n_lag_inf,
                                                                                         val_start = as.numeric(v_map["r_nu_exp2_dx_lb"]),
                                                                                         val_end   = as.numeric(v_map["r_nu_exp2_dx_ub"]),
                                                                                         v_logit_change_rate = as.numeric(v_map["r_nu_exp2_dx_rate"]),
                                                                                         v_logit_change_mid  = as.numeric(v_map["n_nu_exp2_dx_mid"]))),
                                                       l_nu_inf2_dx = add_period(l_period_def = NULL,
                                                                                 l_period_add = make_period(
                                                                                         functional_form = "general logit",
                                                                                         time_start = 0,
                                                                                         time_stop = n_t_calib + n_lag_inf,
                                                                                         val_start = as.numeric(v_map["r_nu_exp2_dx_lb"]),
                                                                                         val_end   = as.numeric(v_map["r_nu_exp2_dx_ub"]),
                                                                                         v_logit_change_rate = as.numeric(v_map["r_nu_exp2_dx_rate"]),
                                                                                         v_logit_change_mid  = as.numeric(v_map["n_nu_exp2_dx_mid"]))),
                                                       v_inf_init_ages  = v_inf_init_ages,
                                                       l_contact_info  = l_contact_matrices,
                                                       l_interventions = l_interventions,
                                                       n_hhsize = n_hhsize,
                                                       r_tau   = v_map["r_tau"],
                                                       r_omega = 0, 
                                                       l_cfr   = get_non_const_multiage_list(v_time_stop = 1:(n_t_calib+n_lag_inf), m_ageval = m_cfr)
        )
        
        ## Load all parameter values
        l_params_all <- sccosmomcma::load_all_params(l_params_init = l_params_init)
        
        ## Targets
        df_valid_targets <- rbind.data.frame(l_targets$cases,
                                             l_targets$cases_inc,
                                             l_targets$deaths,
                                             l_targets$deaths_inc)
        
        
        ## Compute model-predicted outputs -------------------------------------
        
        l_out_post_summ <- validate_out(m_calib_post    = m_calib_post, 
                                        l_params_all    = l_params_all, 
                                        n_date_ini      = n_date_ini, 
                                        n_t_calib       = n_t_calib, 
                                        n_lag_inf       = n_lag_inf,
                                        l_dates_targets = l_dates_targets)
        
        # Save model-predicted outputs
        save(l_fit_imis,
             m_calib_post,
             v_calib_post_mode,
             df_samp_prior_post,
             df_calib_post_map_IMIS,
             gg_post_imis,
             m_cfr,
             l_out_post_summ,
             file = paste0("output/04_validation_IMIS_",abbrev_state,"_",n_proj_type,"_",n_time_stamp,".RData"))
        
        
        ### Plot model-predicted outputs vs. targets ---------------------------
        source("R/04_validation_functions.R")
        gg_validate_imis <- plot_model_out_vs_targets(df_all_targets = df_valid_targets,
                                                      df_model_out   = l_out_post_summ$df_out_all_summ)
        
        # Save plot
        ggsave(filename = paste0("figs/04_validation_targets_vs_model_IMIS_",abbrev_state,"_",n_proj_type,"_", n_time_stamp,".pdf"),
               plot     = gg_validate_imis$gg_model_out_vs_targets, 
               width    = 10, 
               height   = 8)
}

