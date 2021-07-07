# pop_size.R
#------------------------------------------------------------------------------#
# This script generates MCMA population size and proportion by age groups      #
# from the SC-COSMO model.                                                     #
#                                                                              #
# Authors:                                                                     #
#      - Valeria Gracia Olvera, MsC, <valeria.gracia@cide.edu>                 #
#                                                                              #
#------------------------------------------------------------------------------#

rm(list = ls())

# Libraries
library(tidyverse)
library(dplyr)

# Time stamp
n_time_stamp <- "2020-12-13"

for(n_proj_type in c("SQ", "SA")){
  # Load state specific data
  load(paste0("data/MCMA_calibration_data_",n_proj_type,"_",n_time_stamp,".RData"))
  
  
  # Load calibrated parameters
  
  files <- list.files(paste0("./output/"), 
                      pattern = paste0("03_map_output_IMIS_",abbrev_state,"_",n_proj_type,"_",n_time_stamp))
  
  load(file=paste0("./output/",files))
  
  v_map <- df_calib_post_map_IMIS$map
  names(v_map) <- v_param_names
  
  
  # Load all parameters and initialize model 
  
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
  
  l_out_cosmo <- sccosmomcma::cosmo(l_params_all = l_params_all)
  
  # Age groups names
  v_init_age_grps <- c(0, 5, 15, 25, 45, 55, 65, 70)
  v_names_ages <- ordered(c(paste(v_init_age_grps[-length(v_init_age_grps)], 
                                  (v_init_age_grps[-1]-1), sep = "-"), 
                            paste0(v_init_age_grps[length(v_init_age_grps)], "+")),
                          c(paste(v_init_age_grps[-length(v_init_age_grps)], 
                                  (v_init_age_grps[-1]-1), sep = "-"), 
                            paste0(v_init_age_grps[length(v_init_age_grps)], "+")))
  v_names_age_groups <- paste(v_names_ages, "years")
  
  ### Population
  df_pop_size_raw <- calc_popsize_totals(l_out_cosmo)
  names(df_pop_size_raw)[2:9] <- v_names_age_groups
  
  df_pop_size <- df_pop_size_raw[-c(1:n_lag_inf),] %>%
    mutate(state = state_i,
           date = seq.Date(from = as.Date(n_date_ini),
                           to = as.Date(n_date_end_calib),
                           by = "day"),
           prop_kids = (`0-4 years` + `5-14 years`)/All,
           prop_adults = 1 - prop_kids)
  
  # df_pop_size <- df_pop_size_raw[-c(1:n_lag_inf),] %>%
  #   mutate(state = state_i,
  #          date = seq.Date(from = as.Date(n_date_ini),
  #                          to = as.Date(n_date_end_calib),
  #                          by = "day")) %>%
  #   pivot_longer(names_to = "age_groups", cols = 2:9, values_to = "age_pop") %>%
  #   mutate(proj_type = n_proj_type,
  #          age_prop = age_pop/All) %>%
  #   rename(tot_pop = All) %>%
  #   select(state, proj_type, age_groups, date, time, tot_pop, age_pop, age_prop)
  
  save(df_pop_size, file = paste0("data/df_pop_size_",n_proj_type,".rda"))
}
