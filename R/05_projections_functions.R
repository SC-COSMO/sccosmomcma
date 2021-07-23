#' Generates interventions input list for performing projections
#' with the SC-COSMO model
#'
#' \code{get_projection_scenarios} generates interventions input 
#' list for performing projections with the SC-COSMO model 
#' for selected state.
#' @param n_t Simulation total time.
#' @param v_soc_dist_factor Vector with social distancing multipliers for the 
#' different mobility segments calibrated to.
#' @param v_mean_soc_dist_factor Vector with mean social distancing multipliers 
#' for different mobility segments calibrated to.
#' @param v_n_date0_NPI Vector with the time steps (0 = \code{date_init}) at 
#' which effect of NPI changed in the calibration period.
#' @param date_proj0 the time step (0 = \code{date_init}) where projection starts
#' (calibration period is over).
#' @return 
#' A list of named (scenarios) of intervention lists formatted for
#' input into the SC-COSMO model.
#' @export
get_projection_scenarios <- function(n_t , 
                                     v_soc_dist_factor, 
                                     v_mean_soc_dist_factor,
                                     v_n_date0_NPI,
                                     date_proj0) {
  
  v_ind_red_factor <- 1 - v_mean_soc_dist_factor
  
  l_projection_scenarios <- list()
  
# No Intervention ---------------------------------------------------------

  # ## Begin from calibration start 
  # i1 <- make_intervention(intervention_type = "StatusQuo",
  #                         time_start = 0,
  #                         time_stop = v_n_date0_NPI[1] + n_lag_inf)
  # 
  # i2 <- make_intervention(intervention_type = "SocialDistancing",
  #                         time_start = v_n_date0_NPI[1] + n_lag_inf,
  #                         time_stop  = v_n_date0_NPI[2] + n_lag_inf,
  #                         intervention_factor = 1,
  #                         intervention_change_rate = 0.5,
  #                         resume_school = FALSE)
  # 
  # i3 <- make_intervention(intervention_type = "SocialDistancing",
  #                         time_start = v_n_date0_NPI[2] + n_lag_inf,
  #                         time_stop  = v_n_date0_NPI[3] + n_lag_inf,
  #                         intervention_factor = 1,
  #                         intervention_change_rate = 0.5,
  #                         resume_school = FALSE)
  # 
  # i4 <- make_intervention(intervention_type = "SocialDistancing",
  #                         time_start = v_n_date0_NPI[3] + n_lag_inf,
  #                         time_stop  = v_n_date0_NPI[4] + n_lag_inf,
  #                         intervention_factor = 1,
  #                         intervention_change_rate = 0.5,
  #                         resume_school = FALSE)
  # 
  # i5 <- make_intervention(intervention_type = "SocialDistancing",
  #                         time_start = v_n_date0_NPI[4] + n_lag_inf,
  #                         time_stop  = v_n_date0_NPI[5] + n_lag_inf,
  #                         intervention_factor = 1,
  #                         intervention_change_rate = 0.5,
  #                         resume_school = FALSE)
  # 
  # i6 <- make_intervention(intervention_type = "SocialDistancing",
  #                         time_start = v_n_date0_NPI[5] + n_lag_inf,
  #                         time_stop  = date_proj0 + n_lag_inf,
  #                         intervention_factor = 1,
  #                         intervention_change_rate = 0.5,
  #                         resume_school = FALSE)
  # 
  # ## Start projection
  # i7 <- make_intervention(intervention_type = "SocialDistancing",
  #                         time_start = date_proj0 + n_lag_inf,
  #                         time_stop = n_t + 1,
  #                         intervention_factor = 1,
  #                         intervention_change_rate = 0.5,
  #                         resume_school = FALSE)
  # 
  # l_interventions <- add_intervention(interventions = NULL,            intervention = i1)
  # l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
  # l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
  # l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
  # l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
  # l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
  # l_interventions <- add_intervention(interventions = l_interventions, intervention = i7)
  # 
  # name_int <- "No NPIs implemented"
  # 
  # l_projection_scenarios[[name_int]] <- l_interventions
  # 

# Base Case ---------------------------------------------------------------

  ## Begin from calibration start 
  i1 <- make_intervention(intervention_type = "StatusQuo",
                          time_start = 0,
                          time_stop = v_n_date0_NPI[1] + n_lag_inf)
  
  i2 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = v_n_date0_NPI[1] + n_lag_inf, # date_sd0
                          time_stop  = v_n_date0_NPI[2] + n_lag_inf,
                          intervention_factor = v_soc_dist_factor[1],
                          intervention_change_rate = 0.5,
                          resume_school = FALSE)
  
  i3 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = v_n_date0_NPI[2] + n_lag_inf,
                          time_stop  = v_n_date0_NPI[3] + n_lag_inf,
                          intervention_factor = v_soc_dist_factor[2],
                          intervention_change_rate = 0.5,
                          resume_school = FALSE)
  
  i4 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = v_n_date0_NPI[3] + n_lag_inf,
                          time_stop  = v_n_date0_NPI[4] + n_lag_inf,
                          intervention_factor = v_soc_dist_factor[3],
                          intervention_change_rate = 0.5,
                          resume_school = FALSE)
  
  i5 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = v_n_date0_NPI[4] + n_lag_inf,
                          time_stop  = v_n_date0_NPI[5] + n_lag_inf,
                          intervention_factor = v_soc_dist_factor[4],
                          intervention_change_rate = 0.5,
                          resume_school = FALSE)
  
  i6 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = v_n_date0_NPI[5] + n_lag_inf,
                          time_stop  = date_proj0 + n_lag_inf,
                          intervention_factor = v_soc_dist_factor[5],
                          intervention_change_rate = 0.5,
                          resume_school = FALSE)
  ## Start projection
  i7 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = date_proj0 + n_lag_inf,
                          time_stop = n_t + 1,
                          intervention_factor = v_soc_dist_factor[5],
                          intervention_change_rate = 0.5,
                          resume_school = FALSE)
  
  l_interventions <- add_intervention(interventions = NULL,            intervention = i1)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i7)
  
  name_int <- "Social distancing: status quo; Schooling: not in-person; Holiday bump: no"
  
  l_projection_scenarios[[name_int]] <- l_interventions
  

# Interventions -----------------------------------------------------------

  # Create vector with effective contact rates
  v_ind_red_factor_proj <- c(BaseCase   = as.numeric(v_soc_dist_factor[5]),        # continue with estimated SD
                             IncreaseSD = as.numeric(min(v_soc_dist_factor))       # increase SD by applying the min SD observed
  )
  
  school_factor <-  as.numeric(v_soc_dist_factor[5])

## Scenario: STATUS QUO --------------------------------------------------
  
  for(int_factor in v_ind_red_factor_proj){ # int_factor = v_soc_dist_factor[5]
    name_int_red_factor <- names(v_ind_red_factor_proj)[which(v_ind_red_factor_proj == int_factor)]
    
    for(resume_school in c(T,F)){ # resume_school = T
      if(resume_school == F & name_int_red_factor == "BaseCase"){
        next
      }
      
      ### Intervention: resume_school and/or social distancing
      n_date_NPI_proj  <- "2021-01-10"
      n_date0_NPI_proj <- as.numeric(difftime(as.Date(n_date_NPI_proj),
                                              n_date_ini, units = "days"))
      
      
      i1 <- make_intervention(intervention_type = "StatusQuo",
                              time_start = 0,
                              time_stop = v_n_date0_NPI[1] + n_lag_inf)
      
      i2 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[1] + n_lag_inf, # date_sd0
                              time_stop  = v_n_date0_NPI[2] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[1],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i3 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[2] + n_lag_inf,
                              time_stop  = v_n_date0_NPI[3] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[2],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i4 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[3] + n_lag_inf,
                              time_stop  = v_n_date0_NPI[4] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[3],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i5 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[4] + n_lag_inf,
                              time_stop  = v_n_date0_NPI[5] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[4],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i6 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[5] + n_lag_inf,
                              time_stop  = date_proj0 + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[5],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      ## Start projection 
      # Continue with social distancing until 2021-01-11
      i7 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = date_proj0 + n_lag_inf,
                              time_stop = n_date0_NPI_proj + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[5],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      if(resume_school){

        # Add interventions
        i8 <- make_intervention(intervention_type = "SocialDistancing",
                                time_start = n_date0_NPI_proj + n_lag_inf,
                                time_stop = n_t + 1,
                                intervention_factor = int_factor,
                                intervention_change_rate = 0.5,
                                resume_school = resume_school,
                                school_intervention_factor = school_factor)
        
        l_interventions <- add_intervention(interventions = NULL,            intervention = i1)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i7)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i8)
        
        if(name_int_red_factor == "IncreaseSD"){
          name_int <- "Social distancing: stricter; Schooling: in-person; Holiday bump: no"
        }
        if(name_int_red_factor == "BaseCase"){
          name_int <- "Social distancing: status quo; Schooling: in-person; Holiday bump: no"
        }
        l_projection_scenarios[[name_int]] <- l_interventions
        
      }else{
        
        # Add interventions: school and/or work at 50% and 75%
        i8 <- make_intervention(intervention_type = "SocialDistancing",
                                time_start = n_date0_NPI_proj + n_lag_inf,
                                time_stop = n_t + 1,
                                intervention_factor = int_factor,
                                intervention_change_rate = 0.5,
                                resume_school = resume_school)
        
        l_interventions <- add_intervention(interventions = NULL,            intervention = i1)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i7)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i8)
        
        if(name_int_red_factor == "IncreaseSD"){
          name_int <- "Social distancing: stricter; Schooling: not in-person; Holiday bump: no"
        }
        if(name_int_red_factor == "BaseCase"){
          name_int <- "Social distancing: status quo; Schooling: not in-person; Holiday bump: no"
        }
        
        l_projection_scenarios[[name_int]] <- l_interventions
      }
    }
  }

## Scenario: increase in contacts - HOLIDAYS ------------------------------
  
  for(int_factor in v_ind_red_factor_proj){
    name_int_red_factor <- names(v_ind_red_factor_proj)[which(v_ind_red_factor_proj == int_factor)]
    for(resume_school in c(T,F)){
      
      ### Intervention: resume_school and/or social distancing
      n_date_NPI_proj  <- "2021-01-10"
      n_date0_NPI_proj <- as.numeric(difftime(as.Date(n_date_NPI_proj),
                                              n_date_ini, units = "days"))
      
      ### Intervention: increase effective contacts on holidays
      sd_holidays <- min(max(v_soc_dist_factor[5]+0.30, 0),1)
      
      n_date_holidays_init  <- "2020-12-24"
      n_date0_holidays_init <- as.numeric(difftime(as.Date(n_date_holidays_init),
                                                   n_date_ini, units = "days"))
      
      n_date_holidays_end  <- "2021-01-06"
      n_date0_holidays_end <- as.numeric(difftime(as.Date(n_date_holidays_end),
                                                  n_date_ini, units = "days"))
      
      
      i1 <- make_intervention(intervention_type = "StatusQuo",
                              time_start = 0,
                              time_stop = v_n_date0_NPI[1] + n_lag_inf)
      
      i2 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[1] + n_lag_inf, # date_sd0
                              time_stop  = v_n_date0_NPI[2] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[1],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i3 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[2] + n_lag_inf,
                              time_stop  = v_n_date0_NPI[3] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[2],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i4 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[3] + n_lag_inf,
                              time_stop  = v_n_date0_NPI[4] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[3],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i5 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[4] + n_lag_inf,
                              time_stop  = v_n_date0_NPI[5] + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[4],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      i6 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = v_n_date0_NPI[5] + n_lag_inf,
                              time_stop  = date_proj0 + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[5],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      ## Start projection 
      # Continue with social distancing until 2020-12-24
      i7 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = date_proj0 + n_lag_inf,
                              time_stop = n_date0_holidays_init + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[5],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      # Decrease SD on holidays from 2020-12-24 to 2021-01-06
      i8 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = n_date0_holidays_init + n_lag_inf,
                              time_stop = n_date0_holidays_end + n_lag_inf,
                              intervention_factor = sd_holidays,
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      # Return to estimated SD at 12/07 from 2021-01-06 to 2021-01-10
      i9 <- make_intervention(intervention_type = "SocialDistancing",
                              time_start = n_date0_holidays_end + n_lag_inf,
                              time_stop = n_date0_NPI_proj + n_lag_inf,
                              intervention_factor = v_soc_dist_factor[5],
                              intervention_change_rate = 0.5,
                              resume_school = FALSE)
      
      if(resume_school){
        
        # Add interventions: school and/or work at 50% and 75%
        i10 <- make_intervention(intervention_type = "SocialDistancing",
                                 time_start = n_date0_NPI_proj + n_lag_inf,
                                 time_stop = n_t + 1,
                                 intervention_factor = int_factor,
                                 intervention_change_rate = 0.5,
                                 resume_school = resume_school,
                                 school_intervention_factor = school_factor)
        
        
        l_interventions <- add_intervention(interventions = NULL,            intervention = i1)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i7)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i8)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i9)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i10)
        
        if(name_int_red_factor == "IncreaseSD"){
          name_int <- "Social distancing: stricter; Schooling: in-person; Holiday bump: yes"
        }
        if(name_int_red_factor == "BaseCase"){
          name_int <- "Social distancing: status quo; Schooling: in-person; Holiday bump: yes"
        }
        l_projection_scenarios[[name_int]] <- l_interventions
    
      }else{
        
        # Add interventions: school and/or work at 50% and 75%
        i10 <- make_intervention(intervention_type = "SocialDistancing",
                                 time_start = n_date0_NPI_proj + n_lag_inf,
                                 time_stop = n_t + 1,
                                 intervention_factor = int_factor,
                                 intervention_change_rate = 0.5,
                                 resume_school = resume_school)
        
        l_interventions <- add_intervention(interventions = NULL,            intervention = i1)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i7)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i8)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i9)
        l_interventions <- add_intervention(interventions = l_interventions, intervention = i10)
        
        if(name_int_red_factor == "IncreaseSD"){
          name_int <- "Social distancing: stricter; Schooling: not in-person; Holiday bump: yes"
        }
        if(name_int_red_factor == "BaseCase"){
          name_int <- "Social distancing: status quo; Schooling: not in-person; Holiday bump: yes"
        }
        l_projection_scenarios[[name_int]] <- l_interventions
      }
    }
  }
  
  return(l_projection_scenarios)
}

#' Calculate epidemiological outcomes projections from SC-COSMO
#'
#' \code{project_epi_out} projects epidemiological outcomes from the 
#' SC-COSMO model for selected states in Mexico.
#' 
#' @param v_states_project Vector specifying the name of the states to be
#' projected.
#' @param l_params_all List of all SC-COSMO model parameters.
#' @param n_date_ini Initial calendar date of the simulation.
#' @param n_lag_inf Lag in time series of infectious individuals.
#' @return 
#' A list with epidemiological outcomes.
#' @export
project_epi_out <- function(v_states_project,
                            l_params_all,
                            n_date_ini,
                            n_lag_inf = NULL){ # User defined
  
  df_DXCumtot      <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_DXInctot      <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_DXtot         <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_DXDcov_tot    <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_DXDIncCov_tot <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_Hosp_tot      <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_NonICU_tot    <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_ICU_tot       <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_InfCumtot     <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_InfCumprop    <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_InfInctot     <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_Inftot        <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_ExpInf_tot    <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_Dcov_tot      <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_Rec_tot       <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_Rec_prop      <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_Rt_tot        <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_CDR_tot       <- data.frame(county = NULL, Outcome = NULL, time = NULL, value = NULL)
  df_DXCumages     <- c()
  df_DXIncages     <- c()
  df_InfCumages    <- c()
  
  # print(paste0(l_interventions[[4]]$intervention_factor,
  #              "; r_beta = ", round(l_params_all$r_beta, 3), 
  #              "; r_tau = ", round(l_params_all$r_tau, 3), 
  #              # "; r_nu_exp2_dx_lb = ", round(l_params_all$r_nu_exp2_dx_lb, 3), 
  #              # "; r_nu_exp2_dx_ub = ", round(l_params_all$r_nu_exp2_dx_ub, 3),
  #              # "; r_nu_exp2_dx_rate = ", round(l_params_all$r_nu_exp2_dx_rate, 3), 
  #              # "; n_nu_exp2_dx_mid = ", round(l_params_all$n_nu_exp2_dx_mid, 3),
  #              "; n_date_ini = ", n_date_ini))
  # 
  ### Run SC-COSMO model with updated calibrated parameters
  l_out_cosmo <- sccosmomcma::cosmo(l_params_all = l_params_all)
  
  ### Population
  df_popize <- calc_popsize_totals(l_out_cosmo)
  
  ####### Epidemiological Output ###########################################
  v_dates  <- n_date_ini + 0:n_t_project
  v_dates0 <- 0:sum(n_t_project)
  
  #### Cumulative infections total  ####
  df_infcum_ages <- calc_infcum_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_InfCumtot <- bind_rows(df_InfCumtot,
                              data.frame(county  =  v_states_project,
                                         Outcome = "Cumulative infections",
                                         dates   = v_dates,
                                         dates0  = v_dates0,
                                         time    = df_infcum_ages$time,
                                         value   = df_infcum_ages[, "All"])) # (l_params_all$n_ages + 2)
    df_InfCumages <- bind_rows(df_InfCumages,
                               data.frame(county =  v_states_project,
                                          Outcome = "Cumulative infections",
                                          dates = v_dates,
                                          dates0 = v_dates0,
                                          df_infcum_ages,
                                          check.names = FALSE))
  } else {
    df_InfCumtot <- bind_rows(df_InfCumtot,
                              data.frame(county =  v_states_project,
                                         Outcome = "Cumulative infections",
                                         dates = v_dates,
                                         dates0 = v_dates0,
                                         time = df_infcum_ages$time[-c(1:n_lag_inf)],
                                         value = df_infcum_ages[-c(1:n_lag_inf), "All"])) # (l_params_all$n_ages + 2)
    df_InfCumages <- bind_rows(df_InfCumages,
                               data.frame(county =  v_states_project,
                                          Outcome = "Cumulative infections",
                                          dates = v_dates,
                                          dates0 = v_dates0,
                                          df_infcum_ages[-c(1:n_lag_inf), ],
                                          check.names = FALSE))
  }
  
  #### Cumulative infections proportion  ####
  if(is.null(n_lag_inf)){
    df_InfCumprop <- bind_rows(df_InfCumprop,
                               data.frame(county  =  v_states_project,
                                          Outcome = "Cumulative infections proportion",
                                          dates   = v_dates,
                                          dates0  = v_dates0,
                                          time    = df_infcum_ages$time,
                                          value   = df_infcum_ages[, "All"]/df_popize[, ncol(df_popize)])) # (l_params_all$n_ages + 2)
  } else {
    df_InfCumprop <- bind_rows(df_InfCumprop,
                               data.frame(county =  v_states_project,
                                          Outcome = "Cumulative infections proportion",
                                          dates = v_dates,
                                          dates0 = v_dates0,
                                          time = df_infcum_ages$time[-c(1:n_lag_inf)],
                                          value = df_infcum_ages[-c(1:n_lag_inf), "All"]/df_popize[-c(1:n_lag_inf), ncol(df_popize)])) # (l_params_all$n_ages + 2)
  }
  
  #### Incident infections total ####
  df_infinc_ages <- calc_infinc_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_InfInctot <- bind_rows(df_InfInctot,
                              data.frame(county =  v_states_project,
                                         Outcome = "Incident infections",
                                         dates = v_dates,
                                         dates0 = v_dates0,
                                         time = df_infinc_ages$time,
                                         value = df_infinc_ages[, "All"])) # (l_params_all$n_ages + 2)
  } else {
    df_InfInctot <- bind_rows(df_InfInctot,
                              data.frame(county =  v_states_project,
                                         Outcome = "Incident infections",
                                         dates = v_dates,
                                         dates0 = v_dates0,
                                         time = df_infinc_ages$time[-c(1:n_lag_inf)],
                                         value = df_infinc_ages[-c(1:n_lag_inf), "All"])) # (l_params_all$n_ages + 2)
  }
  
  #### Prevalence: Total COVID Infections  ####
  df_inftot_ages <- calc_inf_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_Inftot <- bind_rows(df_Inftot,
                           data.frame(county =  v_states_project,
                                      Outcome = "Prevalent infections",
                                      time = df_inftot_ages$time,
                                      dates = v_dates,
                                      dates0 = v_dates0,
                                      value = df_inftot_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_Inftot <- bind_rows(df_Inftot,
                           data.frame(county =  v_states_project,
                                      Outcome = "Prevalent infections",
                                      dates = v_dates,
                                      dates0 = v_dates0,
                                      time = df_inftot_ages$time[-c(1:n_lag_inf)],
                                      value = df_inftot_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }
  
  #### Prevalence: Total COVID Infections (Es and Is)  ####
  df_expinftot_ages <- calc_expinf_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_ExpInf_tot <- bind_rows(df_ExpInf_tot,
                               data.frame(county =  v_states_project,
                                          Outcome = "Infections (Es and Is)",
                                          time = df_expinftot_ages$time,
                                          value = df_expinftot_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_ExpInf_tot <- bind_rows(df_ExpInf_tot,
                               data.frame(county =  v_states_project,
                                          Outcome = "Infections (Es and Is)",
                                          time = df_expinftot_ages$time[-c(1:n_lag_inf)],
                                          value = df_expinftot_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }
  
  #### Cumulative detected cases total ####
  df_dxcum_ages <- calc_dxcum_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_DXCumtot <- bind_rows(df_DXCumtot,
                             data.frame(county =  v_states_project,
                                        Outcome = "Cumulative detected cases",
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        time = df_dxcum_ages$time,
                                        value = df_dxcum_ages[, (l_params_all$n_ages + 2)]))
    df_DXCumages <- bind_rows(df_DXCumages,
                              data.frame(county =  v_states_project,
                                         Outcome = "Cumulative detected cases",
                                         dates = v_dates,
                                         dates0 = v_dates0,
                                         df_dxcum_ages,
                                         check.names = FALSE))
  } else {
    df_DXCumtot <- bind_rows(df_DXCumtot,
                             data.frame(county =  v_states_project,
                                        Outcome = "Cumulative detected cases",
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        time = df_dxcum_ages$time[-c(1:n_lag_inf)],
                                        value = df_dxcum_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
    df_DXCumages <- bind_rows(df_DXCumages,
                              data.frame(county =  v_states_project,
                                         Outcome = "Cumulative detected cases",
                                         dates = v_dates,
                                         dates0 = v_dates0,
                                         df_dxcum_ages[-c(1:n_lag_inf), ],
                                         check.names = FALSE))
  }
  
  #### Incident detected total cases  ####
  df_dxinc_ages <- calc_dxinc_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_DXInctot <- bind_rows(df_DXInctot,
                             data.frame(county =  v_states_project,
                                        Outcome = "Incident detected cases",
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        time = df_dxinc_ages$time,
                                        value = df_dxinc_ages[, (l_params_all$n_ages + 2)]))
    df_DXIncages <- bind_rows(df_DXIncages,
                              data.frame(county =  v_states_project,
                                         Outcome = "Incident detected cases",
                                         dates = v_dates,
                                         dates0 = v_dates0,
                                         df_dxinc_ages,
                                         check.names = FALSE))
  } else {
    df_DXInctot <- bind_rows(df_DXInctot,
                             data.frame(county =  v_states_project,
                                        Outcome = "Incident detected cases",
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        time = df_dxinc_ages$time[-c(1:n_lag_inf)],
                                        value = df_dxinc_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
    df_DXIncages <- bind_rows(df_DXIncages,
                              data.frame(county =  v_states_project,
                                         Outcome = "Incident detected cases",
                                         dates = v_dates,
                                         dates0 = v_dates0,
                                         df_dxinc_ages[-c(1:n_lag_inf), ],
                                         check.names = FALSE))
  }
  
  #### Prevalent detected cases total ####
  df_dx_ages <- calc_dx_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_DXtot <- bind_rows(df_DXtot,
                          data.frame(county =  v_states_project,
                                     Outcome = "Detected cases",
                                     dates = v_dates,
                                     dates0 = v_dates0,
                                     time = df_dx_ages$time,
                                     value = df_dx_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_DXtot <- bind_rows(df_DXtot,
                          data.frame(county =  v_states_project,
                                     Outcome = "Detected cases",
                                     dates = v_dates,
                                     dates0 = v_dates0,
                                     time = df_dx_ages$time[-c(1:n_lag_inf)],
                                     value = df_dx_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }
  
  #### Prevalence Recovered total ####
  df_Rec_ages <- calc_rec_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_Rec_tot <- bind_rows(df_Rec_tot,
                            data.frame(county  =  v_states_project,
                                       Outcome = "Recovered prevalence",
                                       dates   = v_dates,
                                       dates0  = v_dates0,
                                       time    = df_Rec_ages$time,
                                       value   = df_Rec_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_Rec_tot <- bind_rows(df_Rec_tot,
                            data.frame(county  =  v_states_project,
                                       Outcome = "Recovered prevalence",
                                       dates   = v_dates,
                                       dates0  = v_dates0,
                                       time    = df_Rec_ages$time[-c(1:n_lag_inf)],
                                       value   = df_Rec_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
    
  }
  
  #### Prevalence Recovered proportion ####
  if(is.null(n_lag_inf)){
    df_Rec_prop <- bind_rows(df_Rec_prop,
                             data.frame(county  =  v_states_project,
                                        Outcome = "Recovered prevalence proportion",
                                        dates   = v_dates,
                                        dates0  = v_dates0,
                                        time    = df_Rec_ages$time,
                                        value   = df_Rec_ages[, (l_params_all$n_ages + 2)]/df_popize[, ncol(df_popize)]))
  } else {
    df_Rec_prop <- bind_rows(df_Rec_prop,
                             data.frame(county  =  v_states_project,
                                        Outcome = "Recovered prevalence proportion",
                                        dates   = v_dates,
                                        dates0  = v_dates0,
                                        time    = df_Rec_ages$time[-c(1:n_lag_inf)],
                                        value   = df_Rec_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]/df_popize[-c(1:n_lag_inf), ncol(df_popize)]))
    
  }
  
  #### Case Detection Ratio  ####
  df_CDRdenom <- df_ExpInf_tot %>%
    rename(denom = value)
  
  df_CDR_tot <- merge((df_DXtot %>% select(-Outcome)), (df_CDRdenom %>% select(-Outcome))) %>%
    mutate(value_orig = value) %>%
    mutate(value = value_orig/denom) %>%
    mutate(Outcome = "CDR proportion") %>%
    select(-value_orig, -denom)
  
  df_CDR_tot <- df_CDR_tot[, c("county", "Outcome", "dates", "dates0", "time", "value")]
  #relocate(county, Outcome, dates, dates0, time, value)
  
  #### Incident COVID19 deaths infections ####
  df_DXDIncCov_ages <- calc_incdeathsdx_totals(l_out_cosmo)
  
  if(is.null(n_lag_inf)){
    df_DXDIncCov_tot <- bind_rows(df_DXDIncCov_tot,
                                  data.frame(county =  v_states_project,
                                             Outcome = "Incident COVID19 deaths infections",
                                             time = df_DXDIncCov_ages$time,
                                             dates = v_dates,
                                             dates0 = v_dates0,
                                             value = df_DXDIncCov_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_DXDIncCov_tot <- bind_rows(df_DXDIncCov_tot,
                                  data.frame(county =  v_states_project,
                                             Outcome = "Incident COVID19 deaths infections",
                                             dates = v_dates,
                                             dates0 = v_dates0,
                                             time = df_DXDIncCov_ages$time[-c(1:n_lag_inf)],
                                             value = df_DXDIncCov_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }
  
  ### Apply delay in deaths
  df_DXDIncCov_tot <- df_DXDIncCov_tot %>%
    mutate(dates = dates + n_death_delay) %>%
    filter(dates <= n_date_end_project & dates >= l_dates_targets$deaths[1]) %>%
    complete(dates = seq.Date(from = as.Date(l_dates_targets$deaths[1]),
                              to   = as.Date(n_date_end_project),
                              by   = "day"),
             fill = list(county = v_states_project,
                         Outcome    = "Incident COVID19 deaths infections",
                         value = df_DXDIncCov_tot$value[1]
             )) %>%
    mutate(dates0 = as.numeric(dates - dates[1]),
           time = n_lag_inf:(as.Date(n_date_end_project) - 
                               as.Date(l_dates_targets$deaths[1]) + n_lag_inf))
  
  #### Cumulative COVID Deaths  ####
  df_DCov_ages <- calc_deaths_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_Dcov_tot <- bind_rows(df_Dcov_tot,
                             data.frame(county =  v_states_project,
                                        Outcome = "COVID deaths",
                                        time = df_DCov_ages$time,
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        value = df_DCov_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_Dcov_tot <- bind_rows(df_Dcov_tot,
                             data.frame(county =  v_states_project,
                                        Outcome = "COVID deaths",
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        time = df_DCov_ages$time[-c(1:n_lag_inf)],
                                        value = df_DCov_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }
  
  ### Apply delay in deaths
  df_Dcov_tot <- df_Dcov_tot %>%
    mutate(dates = dates + n_death_delay) %>%
    filter(dates <= n_date_end_project & dates >= l_dates_targets$deaths[1]) %>%
    complete(dates = seq.Date(from = as.Date(l_dates_targets$deaths[1]),
                              to   = as.Date(n_date_end_project),
                              by   = "day"),
             fill = list(county = v_states_project,
                         Outcome    = "COVID deaths",
                         value = df_Dcov_tot$value[1]
             )) %>%
    mutate(dates0 = as.numeric(dates - dates[1]),
           time = n_lag_inf:(as.Date(n_date_end_project) - 
                               as.Date(l_dates_targets$deaths[1]) + n_lag_inf))
  
  ###################### Effective Reproduction Number Rt ######################
  rt_start <- 1
  system.time(df_Rt_raw <- calc_reproduction_number_wt(l_out_cosmo,
                                                       v_time = rt_start:(l_params_all$n_t),
                                                       nsim_chosen = 100))
  
  if(is.null(n_lag_inf)){
    df_Rt_tot <- bind_rows(df_Rt_tot,
                           data.frame(county  = v_states_project,
                                      Outcome = "R effective",
                                      dates   = v_dates[-c(1:(rt_start + 1))],
                                      dates0  = v_dates0[-c(1:(rt_start + 1))],
                                      time    = df_Rt_raw$time,
                                      value   = df_Rt_raw$Rt))
  } else {
    df_Rt_tot <- bind_rows(df_Rt_tot,
                           data.frame(county  = v_states_project,
                                      Outcome = "R effective",
                                      dates   = v_dates[-c(1:(rt_start + 1))],
                                      dates0  = v_dates0[-c(1:(rt_start + 1))],
                                      time    = df_Rt_raw$time[-c(1:(n_lag_inf))],
                                      value   = df_Rt_raw$Rt[-c(1:(n_lag_inf))]))
    
  }
  
  # print(paste0("Rt calculated in ", 
  #              round(Rt_time[3]/60, 2), " minutes"))
  
  ############################# Hospitalizations ##############################
  #### PREPPING HOSPITAL CALCS
  l_hosp <- prep_dx_hospitalizations(l_out_cosmo, use_prevalence = FALSE)
  
  #### All Hospitalization Prevalence ####
  df_Hosp_ages <- calc_dx_hosp(l_hosp, l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_Hosp_tot <- bind_rows(df_Hosp_tot,
                             data.frame(county =  v_states_project, 
                                        Outcome = "Total hospitalizations",
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        time = df_Hosp_ages$time,
                                        value = df_Hosp_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_Hosp_tot <- bind_rows(df_Hosp_tot,
                             data.frame(county =  v_states_project, 
                                        Outcome = "Total hospitalizations",
                                        dates = v_dates,
                                        dates0 = v_dates0,
                                        time = df_Hosp_ages$time[-c(1:n_lag_inf)],
                                        value = df_Hosp_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }    
  
  #### Hospitalizations without ventilator ####
  df_NonICU_ages <- calc_dx_hosp_nonicu(l_hosp, l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_NonICU_tot <- bind_rows(df_NonICU_tot,
                               data.frame(county =  v_states_project, 
                                          Outcome = "Hospitalizations without ventilator", 
                                          dates = v_dates,
                                          dates0 = v_dates0,
                                          time = df_NonICU_ages$time,
                                          value = df_NonICU_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_NonICU_tot <- bind_rows(df_NonICU_tot,
                               data.frame(county =  v_states_project, 
                                          Outcome = "Hospitalizations without ventilator", 
                                          dates = v_dates,
                                          dates0 = v_dates0,
                                          time = df_NonICU_ages$time[-c(1:n_lag_inf)],
                                          value = df_NonICU_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }
  
  #### ICU Prevalence ####
  df_ICU_ages <- calc_dx_hosp_icu(l_hosp, l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_ICU_tot <- bind_rows(df_ICU_tot,
                            data.frame(county =  v_states_project, 
                                       Outcome = "Hospitalizations with ventilator", 
                                       dates = v_dates,
                                       dates0 = v_dates0,
                                       time = df_ICU_ages$time,
                                       value = df_ICU_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_ICU_tot <- bind_rows(df_ICU_tot,
                            data.frame(county =  v_states_project, 
                                       Outcome = "Hospitalizations with ventilator", 
                                       dates = v_dates,
                                       dates0 = v_dates0,
                                       time = df_ICU_ages$time[-c(1:n_lag_inf)],
                                       value = df_ICU_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }     
  
  #### DXCOVID19 deaths ####
  df_DXDCov_ages <- calc_deathsdx_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_DXDcov_tot <- bind_rows(df_DXDcov_tot,
                               data.frame(county =  v_states_project, 
                                          Outcome = "Detected COVID deaths",
                                          dates = v_dates,
                                          dates0 = v_dates0,
                                          time = df_DXDCov_ages$time,
                                          value = df_DXDCov_ages[, (l_params_all$n_ages + 2)]))
  } else {
    df_DXDcov_tot <- bind_rows(df_DXDcov_tot,
                               data.frame(county =  v_states_project, 
                                          Outcome = "Detected COVID deaths", 
                                          dates = v_dates,
                                          dates0 = v_dates0,
                                          time = df_DXDCov_ages$time[-c(1:n_lag_inf)],
                                          value = df_DXDCov_ages[-c(1:n_lag_inf), (l_params_all$n_ages + 2)]))
  }
  
  ### Apply delay in deaths
  df_DXDcov_tot <- df_DXDcov_tot %>%
    mutate(dates = dates + n_death_delay) %>%
    filter(dates <= n_date_end_project & dates >= l_dates_targets$deaths[1]) %>%
    complete(dates = seq.Date(from = as.Date(l_dates_targets$deaths[1]),
                              to   = as.Date(n_date_end_project),
                              by   = "day"),
             fill = list(county = v_states_project,
                         Outcome    = "Detected COVID deaths",
                         value = df_DXDcov_tot$value[1]
             )) %>%
    mutate(dates0 = as.numeric(dates - dates[1]),
           time = n_lag_inf:(as.Date(n_date_end_project) - 
                               as.Date(l_dates_targets$deaths[1]) + n_lag_inf))
  
  ### Combine outputs and generate dates for each county
  df_out_mex_tot <- bind_rows(df_InfCumtot,
                              df_InfCumprop,
                              df_InfInctot,                              
                              df_Inftot,
                              df_Dcov_tot,
                              df_DXCumtot,
                              df_DXInctot,
                              df_DXtot,
                              df_Rec_tot,
                              df_Rec_prop,
                              df_Rt_tot,
                              df_Hosp_tot,
                              df_NonICU_tot,
                              df_ICU_tot,
                              df_DXDcov_tot,
                              df_CDR_tot,
                              df_DXDIncCov_tot) 
  
  df_out_mex_ages <- bind_rows(df_InfCumages,
                               df_DXCumages, 
                               df_DXIncages)
  
  ### Return data.frame
  return(list(Total = df_out_mex_tot,
              Ages  = df_out_mex_ages))
}

acomb <- function(...) abind::abind(..., along=3)

#' Probabilistic projections of interventions
#'
#' \code{project_interventions_probabilistic} produces probabilistic projections 
#' of interventions
#'
#' @param m_calib_post Matrix with calibrated parameters from posterior 
#' distribution.
#' @param n_date_ini Initial date of calibration.
#' @param v_n_date0_NPI Vector with the time steps (\code{0 = date_init}) at 
#' which effect of NPI changed in the calibration period.
#' @param n_t_calib Number of calibration days.
#' @param n_t_project Number of projection days.
#' @param n_lag_inf Lag in time series of infectious individuals.
#' @return
#' A list with probabilistic projections of interventions.
#' @export
project_interventions_probabilistic <- function(m_calib_post,
                                                n_date_ini,
                                                v_n_date0_NPI,
                                                n_t_calib,
                                                n_t_project,
                                                n_lag_inf){
  if(is.null(dim(m_calib_post))) { # If vector, change to matrix
    m_calib_post <- t(m_calib_post) 
  }
  
  ### Number of posterior samples
  n_samp <- nrow(m_calib_post)
  
  #### Compute model-predicted outputs for al interventions for each sample of posterior distribution ####
  v_mean_soc_dist_factor <- colMeans(m_calib_post[, c("r_soc_dist_factor", 
                                                      "r_soc_dist_factor_2", 
                                                      "r_soc_dist_factor_3",
                                                      "r_soc_dist_factor_4",
                                                      "r_soc_dist_factor_5"), drop = FALSE])
  
  ### Get OS
  os <- get_os()
  print(paste0("Parallelized projections on ", os))
  
  ### Get cores
  no_cores <- 50
  
  ### Evaluate model at each posterior sample and store resultsl
  if(os == "macosx" | os == "linux"){
    cl <- makeForkCluster(no_cores) 
    registerDoParallel(cl)
    
    time_foreach <- system.time(
      df_out_projection_post_all <- foreach(i = 1:n_samp, .combine = c) %dopar% { # i = 1
        # ### Progress bar
        # if(!exists("pb")) pb <- tcltk::tkProgressBar(title = "Parallel task for Target coverage",
        #                                              min = 1, max = n_samp)
        # info <- sprintf("%s%% done", round(i/n_samp*100))
        # tcltk::setTkProgressBar(pb, i, label = sprintf("Progress of simulations (%s)", info))
        
        # write_log_file(msg=paste(Sys.time(),":  INITIATING iteration:",i,"\n"),log_flag=GLOBAL_LOGGING_ENABLED)
        
        ### Call projection scenarios
        l_interventions_scenarios <- get_projection_scenarios(n_t = n_t_project + n_lag_inf, 
                                                              v_soc_dist_factor = m_calib_post[i, 
                                                                                               c("r_soc_dist_factor", 
                                                                                                 "r_soc_dist_factor_2",
                                                                                                 "r_soc_dist_factor_3",
                                                                                                 "r_soc_dist_factor_4",
                                                                                                 "r_soc_dist_factor_5")],
                                                              v_mean_soc_dist_factor = v_mean_soc_dist_factor,
                                                              v_n_date0_NPI = v_n_date0_NPI, 
                                                              date_proj0    = n_t_calib + n_lag_inf)
        
        df_out_mex_total  <- c()
        df_out_mex_ages   <- c()
        ### Iterate over projection scenarios
        for(scenario_name in names(l_interventions_scenarios)) { # scenario_name <- names(l_interventions_scenarios)[3]
          print(paste0("Running scenario ", scenario_name))
          l_interventions <- l_interventions_scenarios[[scenario_name]]
          ### Initialize parameters
          l_params_init <- sccosmomcma::load_params_init(n_t = n_t_project + n_lag_inf,  # Number of days
                                                     ctry = "Mexico",
                                                     ste  = v_states_project,
                                                     cty  = v_states_project, 
                                                     r_beta = m_calib_post[i,"r_beta"],
                                                     l_nu_exp2_dx = add_period(l_period_def = NULL,
                                                                               l_period_add = make_period(
                                                                                 functional_form = "general logit",
                                                                                 time_start = 0,
                                                                                 time_stop = n_t_project + n_lag_inf,
                                                                                 val_start = as.numeric(m_calib_post[i,"r_nu_exp2_dx_lb"]),
                                                                                 val_end   = as.numeric(m_calib_post[i,"r_nu_exp2_dx_ub"]),
                                                                                 v_logit_change_rate = as.numeric(m_calib_post[i,"r_nu_exp2_dx_rate"]),
                                                                                 v_logit_change_mid  = as.numeric(m_calib_post[i,"n_nu_exp2_dx_mid"]))),
                                                     l_nu_inf2_dx = add_period(l_period_def = NULL,
                                                                               l_period_add = make_period(
                                                                                 functional_form = "general logit",
                                                                                 time_start = 0,
                                                                                 time_stop = n_t_project + n_lag_inf,
                                                                                 val_start = as.numeric(m_calib_post[i,"r_nu_exp2_dx_lb"]),
                                                                                 val_end   = as.numeric(m_calib_post[i,"r_nu_exp2_dx_ub"]),
                                                                                 v_logit_change_rate = as.numeric(m_calib_post[i,"r_nu_exp2_dx_rate"]),
                                                                                 v_logit_change_mid  = as.numeric(m_calib_post[i,"n_nu_exp2_dx_mid"]))),
                                                     v_inf_init_ages  = v_inf_init_ages,
                                                     l_contact_info  = l_contact_matrices,
                                                     l_interventions = l_interventions,
                                                     n_hhsize = n_hhsize,
                                                     r_tau   = m_calib_post[i,"r_tau"],
                                                     r_omega = 0, #1/200
                                                     l_cfr = get_non_const_multiage_list(v_time_stop = 1:(n_t_project+n_lag_inf), m_ageval = m_cfr_proj),
                                                     m_r_exit_tot    = v_hosp_map["m_r_exit_tot"],
                                                     m_r_exit_icu    = v_hosp_map["m_r_exit_icu"],
                                                     m_r_exit_nonicu = v_hosp_map["m_r_exit_nonicu"],
                                                     m_sigma_tot     = v_hosp_map["m_sigma_tot"],
                                                     m_sigma_nonicu  = v_hosp_map["m_sigma_nonicu"],
                                                     m_sigma_icu     = v_hosp_map["m_sigma_icu"]
          )
          ## Load all parameter values
          l_params_all <- sccosmomcma::load_all_params(l_params_init = l_params_init)
          
          df_out_scenario <- project_epi_out(v_states_project = v_states_project,
                                             l_params_all     = l_params_all,
                                             n_date_ini       = n_date_ini,
                                             n_lag_inf        = n_lag_inf)
          
          ### Store interventions for total population and by age groups
          df_out_scenario_tot  <- df_out_scenario$Total
          df_out_scenario_ages <- df_out_scenario$Ages
          
          ### Add intervention names
          df_out_scenario_tot$Intervention  <- scenario_name
          df_out_scenario_ages$Intervention <- scenario_name
          
          ### Add intervention and base-case type
          df_out_scenario_tot$intervention_type = ""
          df_out_scenario_tot$BaseCase_type = ""
          
          if(scenario_name == "No NPIs implemented"){
            df_out_scenario_tot$intervention_type <- "NoNPI"
            df_out_scenario_tot$BaseCase_type <- NA
            
          }else if(scenario_name == "Social distancing: status quo; Schooling: not in-person; Holiday bump: no"){
            df_out_scenario_tot$intervention_type <- "BaseCase"
            df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
            
          }else if(scenario_name == "Social distancing: status quo; Schooling: not in-person; Holiday bump: yes"){
            df_out_scenario_tot$intervention_type <- "BaseCase"
            df_out_scenario_tot$BaseCase_type <-  "Holidays"
            
          }else if(scenario_name == "Social distancing: status quo; Schooling: in-person; Holiday bump: no"){
            df_out_scenario_tot$intervention_type <- "SchoolSD"
            df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
            
          }else if(scenario_name == "Social distancing: status quo; Schooling: in-person; Holiday bump: yes"){
            df_out_scenario_tot$intervention_type <- "SchoolSD"
            df_out_scenario_tot$BaseCase_type <-  "Holidays"
            
          }else if(scenario_name == "Social distancing: stricter; Schooling: not in-person; Holiday bump: no"){
            df_out_scenario_tot$intervention_type <- "IncreaseSD"
            df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
            
          }else if(scenario_name == "Social distancing: stricter; Schooling: not in-person; Holiday bump: yes"){
            df_out_scenario_tot$intervention_type <- "IncreaseSD"
            df_out_scenario_tot$BaseCase_type <-  "Holidays"
            
          }else if(scenario_name == "Social distancing: stricter; Schooling: in-person; Holiday bump: no"){
            df_out_scenario_tot$intervention_type <- "IncreaseSDSchoolSD"
            df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
            
          }else if(scenario_name == "Social distancing: stricter; Schooling: in-person; Holiday bump: yes"){
            df_out_scenario_tot$intervention_type <- "IncreaseSDSchoolSD"
            df_out_scenario_tot$BaseCase_type <-  "Holidays"
            
          }
          
          # Combine scenarios
          df_out_scenario_tot <- df_out_scenario_tot %>%
            mutate(type = "Model-predicted",
                   simulation = i)
          df_out_scenario_ages <- df_out_scenario_ages %>%
            mutate(type = "Model-predicted",
                   simulation = i)
          
          df_out_mex_total  <- bind_rows(df_out_mex_total, 
                                         df_out_scenario_tot)
          # df_out_mex_ages <- bind_rows(df_out_mex_ages, 
          #                              df_out_scenario_ages)
          
        }

        # Return data.frame
        df_out_mex_total
      }
    )
    
  }else if(os == "windows"){
    cl <- makeCluster(no_cores)       # initialize cluster object
    registerDoParallel(cl)
    opts <- list(attachExportEnv = TRUE)
    
    time_foreach <- system.time(
      df_out_projection_post_all <- foreach(i = 1:n_samp, .combine = rbind, .export = ls(globalenv()), # i = 1
                                            .packages=c("sccosmomcma",
                                                        "tidyverse",
                                                        "dplyr",
                                                        "lubridate",
                                                        "dampack",
                                                        "epitools"),
                                            .options.snow = opts) %dopar% { # i = 1
                                              # ### Progress bar
                                              # if(!exists("pb")) pb <- tcltk::tkProgressBar(title = "Parallel task for Target coverage",
                                              #                                              min = 1, max = n_samp)
                                              # info <- sprintf("%s%% done", round(i/n_samp*100))
                                              # tcltk::setTkProgressBar(pb, i, label = sprintf("Progress of simulations (%s)", info))
                                              
                                              # write_log_file(msg=paste(Sys.time(),":  INITIATING iteration:",i,"\n"),log_flag=GLOBAL_LOGGING_ENABLED)
                                              
                                              ### Call projection scenarios
                                              l_interventions_scenarios <- get_projection_scenarios(n_t = n_t_project + n_lag_inf, 
                                                                                                    v_soc_dist_factor = m_calib_post[i, 
                                                                                                                                     c("r_soc_dist_factor", 
                                                                                                                                       "r_soc_dist_factor_2",
                                                                                                                                       "r_soc_dist_factor_3",
                                                                                                                                       "r_soc_dist_factor_4",
                                                                                                                                       "r_soc_dist_factor_5")],
                                                                                                    v_mean_soc_dist_factor = v_mean_soc_dist_factor,
                                                                                                    v_n_date0_NPI = v_n_date0_NPI, 
                                                                                                    date_proj0    = n_t_calib + n_lag_inf)
                                              
                                              df_out_mex_total  <- c()
                                              df_out_mex_ages   <- c()
                                              ### Iterate over projection scenarios
                                              for(scenario_name in names(l_interventions_scenarios)) { # scenario_name <- names(l_interventions_scenarios)[1]
                                                print(paste0("Running scenario ", scenario_name))
                                                l_interventions <- l_interventions_scenarios[[scenario_name]]
                                                ### Initialize parameters
                                                l_params_init <- sccosmomcma::load_params_init(n_t = n_t_project + n_lag_inf,  # Number of days
                                                                                           ctry = "Mexico",
                                                                                           ste  = v_states_project,
                                                                                           cty  = v_states_project, 
                                                                                           v_reduced_sus = v_reduced_sus,
                                                                                           r_beta = m_calib_post[i,"r_beta"],
                                                                                           l_nu_exp2_dx = add_period(l_period_def = NULL,
                                                                                                                     l_period_add = make_period(
                                                                                                                       functional_form = "general logit",
                                                                                                                       time_start = 0,
                                                                                                                       time_stop = n_t_project + n_lag_inf,
                                                                                                                       val_start = as.numeric(m_calib_post[i,"r_nu_exp2_dx_lb"]),
                                                                                                                       val_end   = as.numeric(m_calib_post[i,"r_nu_exp2_dx_ub"]),
                                                                                                                       v_logit_change_rate = as.numeric(m_calib_post[i,"r_nu_exp2_dx_rate"]),
                                                                                                                       v_logit_change_mid  = as.numeric(m_calib_post[i,"n_nu_exp2_dx_mid"]))),
                                                                                           l_nu_inf2_dx = add_period(l_period_def = NULL,
                                                                                                                     l_period_add = make_period(
                                                                                                                       functional_form = "general logit",
                                                                                                                       time_start = 0,
                                                                                                                       time_stop = n_t_project + n_lag_inf,
                                                                                                                       val_start = as.numeric(m_calib_post[i,"r_nu_exp2_dx_lb"]),
                                                                                                                       val_end   = as.numeric(m_calib_post[i,"r_nu_exp2_dx_ub"]),
                                                                                                                       v_logit_change_rate = as.numeric(m_calib_post[i,"r_nu_exp2_dx_rate"]),
                                                                                                                       v_logit_change_mid  = as.numeric(m_calib_post[i,"n_nu_exp2_dx_mid"]))),
                                                                                           v_inf_init_ages  = v_inf_init_ages,
                                                                                           l_contact_info  = l_contact_matrices,
                                                                                           l_interventions = l_interventions,
                                                                                           n_hhsize = n_hhsize,
                                                                                           r_tau   = m_calib_post[i,"r_tau"],
                                                                                           r_omega = 0, #1/200
                                                                                           l_cfr = get_non_const_multiage_list(v_time_stop = 1:(n_t_project+n_lag_inf), m_ageval = m_cfr_proj),
                                                                                           m_r_exit_tot    = v_hosp_map["m_r_exit_tot"],
                                                                                           m_r_exit_icu    = v_hosp_map["m_r_exit_icu"],
                                                                                           m_r_exit_nonicu = v_hosp_map["m_r_exit_nonicu"],
                                                                                           m_sigma_tot     = v_hosp_map["m_sigma_tot"],
                                                                                           m_sigma_nonicu  = v_hosp_map["m_sigma_nonicu"],
                                                                                           m_sigma_icu     = v_hosp_map["m_sigma_icu"]
                                                )
                                                ## Load all parameter values
                                                l_params_all <- load_all_params(l_params_init = l_params_init)
                                                
                                                df_out_scenario <- project_epi_out(v_states_project = v_states_project,
                                                                                   l_params_all     = l_params_all,
                                                                                   n_date_ini       = n_date_ini,
                                                                                   n_lag_inf        = n_lag_inf)
                                                
                                                ### Store interventions for total population and by age groups
                                                df_out_scenario_tot  <- df_out_scenario$Total
                                                df_out_scenario_ages <- df_out_scenario$Ages
                                                
                                                ### Add intervention names
                                                df_out_scenario_tot$Intervention  <- scenario_name
                                                df_out_scenario_ages$Intervention <- scenario_name
                                                
                                                ### Add intervention and base-case type
                                                df_out_scenario_tot$intervention_type = ""
                                                df_out_scenario_tot$BaseCase_type = ""
                                                
                                                if(scenario_name == "No NPIs implemented"){
                                                  df_out_scenario_tot$intervention_type <- "NoNPI"
                                                  df_out_scenario_tot$BaseCase_type <- NA
                                                  
                                                }else if(scenario_name == "Social distancing: status quo; Schooling: not in-person; Holiday bump: no"){
                                                  df_out_scenario_tot$intervention_type <- "BaseCase"
                                                  df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
                                                  
                                                }else if(scenario_name == "Social distancing: status quo; Schooling: not in-person; Holiday bump: yes"){
                                                  df_out_scenario_tot$intervention_type <- "BaseCase"
                                                  df_out_scenario_tot$BaseCase_type <-  "Holidays"
                                                  
                                                }else if(scenario_name == "Social distancing: status quo; Schooling: in-person; Holiday bump: no"){
                                                  df_out_scenario_tot$intervention_type <- "SchoolSD"
                                                  df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
                                                  
                                                }else if(scenario_name == "Social distancing: status quo; Schooling: in-person; Holiday bump: yes"){
                                                  df_out_scenario_tot$intervention_type <- "SchoolSD"
                                                  df_out_scenario_tot$BaseCase_type <-  "Holidays"
                                                  
                                                }else if(scenario_name == "Social distancing: stricter; Schooling: not in-person; Holiday bump: no"){
                                                  df_out_scenario_tot$intervention_type <- "IncreaseSD"
                                                  df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
                                                  
                                                }else if(scenario_name == "Social distancing: stricter; Schooling: not in-person; Holiday bump: yes"){
                                                  df_out_scenario_tot$intervention_type <- "IncreaseSD"
                                                  df_out_scenario_tot$BaseCase_type <-  "Holidays"
                                                  
                                                }else if(scenario_name == "Social distancing: stricter; Schooling: in-person; Holiday bump: no"){
                                                  df_out_scenario_tot$intervention_type <- "IncreaseSDSchoolSD"
                                                  df_out_scenario_tot$BaseCase_type <-  "StatusQuo"
                                                  
                                                }else if(scenario_name == "Social distancing: stricter; Schooling: in-person; Holiday bump: yes"){
                                                  df_out_scenario_tot$intervention_type <- "IncreaseSDSchoolSD"
                                                  df_out_scenario_tot$BaseCase_type <-  "Holidays"
                                                  
                                                }
                                                
                                                ######## Combine scenarios #######
                                                df_out_scenario_tot <- df_out_scenario_tot %>%
                                                  mutate(type = "Model-predicted",
                                                         simulation = i)
                                                df_out_scenario_ages <- df_out_scenario_ages %>%
                                                  mutate(type = "Model-predicted",
                                                         simulation = i)
                                                
                                                df_out_mex_total  <- bind_rows(df_out_mex_total, 
                                                                               df_out_scenario_tot)
                                                # df_out_mex_ages <- bind_rows(df_out_mex_ages, 
                                                #                              df_out_scenario_ages)
                                                
                                              }

                                              # Return data.frame
                                              df_out_mex_total
                                            }
    )
  }
  stopCluster(cl)
  print(paste0("Model evaluated ", scales::comma(n_samp), " times in ", 
               round(time_foreach[3]/60, 2), " minutes"))
  
  df_out_projection_post_all_summ <- df_out_projection_post_all %>% 
    group_by(type, Intervention, intervention_type, BaseCase_type, Outcome, dates) %>%
    summarise(mean = mean(value),
              median = quantile(value, probs = 0.5, names = FALSE),
              sd = sd(value),
              lb = quantile(value, probs = 0.025, names = FALSE),
              ub = quantile(value, probs = 0.975, names = FALSE)
              ) 
  
  colnames(df_out_projection_post_all_summ)[colnames(df_out_projection_post_all_summ)=="mean"] <- "value"
  return(list(df_all  = df_out_projection_post_all,
              df_summ = df_out_projection_post_all_summ))
}
