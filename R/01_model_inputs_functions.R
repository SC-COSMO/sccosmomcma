#' Base-case inital parameter set
#'
#' \code{load_params_init} generates the initial values of the SC-COSMO parameters
#' 
#' @param n_t Time horizon in day.
#' @param time_step Model evaluations per day.
#' @param comp Logical. Should the compiled version of SC-COSMO be used?
#' @param v_init_age_grps Vector with initial ages of age groups.
#' @param v_inf_init_ages Vector with number of individuals in first .
#' infectious class for each age group. Default will assign 1 for each age 
#' group.
#' @param ctry Country.
#' @param ste State.
#' @param cty County or municipality.
#' @param l_contact_info A list with overall and setting specific contact 
#' matrices, typically obtained from `sccosmoData`.
#' @param r_birth Daily birth rate (crude).
#' @param r_beta Transmission parameters.
#' @param v_sigma Age-specific daily rate of progression of exposed individuals 
#' (Latent period)
#' @param v_gamma Age-specific daily rate of recovery of infectious individuals 
#' (Infectiousness period)
#' @param v_omega Age-specific daily rate of waning immunity of recovered 
#' individuals
#' @param l_nu_exp1_dx list specifying time-varying case detection when exposed l=1
#' @param l_nu_exp2_dx list specifying time-varying case detection when exposed l=2
#' @param l_nu_inf1_dx list specifying time-varying case detection when infectious l=1
#' @param l_nu_inf2_dx list specifying time-varying case detection when infectious l=2
#' @param l_phi_exp1_dx list specifying time-varying screen detection when exposed l=1
#' @param l_phi_exp2_dx list specifying time-varying screen detection when exposed l=2
#' @param l_phi_inf1_dx list specifying time-varying screen detection when infectious l=1
#' @param l_phi_inf2_dx list specifying time-varying screen detection when infectious l=2
#' @param v_cfr Age-specific case fatality rate. Default from \code{data-raw} 
#' folder.
#' @param v_ifr Age-specific infection fatality rate. Default from \code{data-raw} 
#' folder.
#' @param l_interventions a list of intervention "objects" defining type, timing, and
#' effects of interventions.
#' @param v_idx_scale_factor Age-specific infectiousness reduction factor on detected cases
#' @param v_alpha_dx Age-specific reduction in mortality on detected infectious vs undetected infectious
#' @param l_l_hospitalizations list of lists specifying age-specific,
#'                             time-varying probability of hospitalizations
#' @param v_p_hosp Age-specific Probability of hospitalization for detected (prevalent or 
#' incident) cases     
#' @param v_p_s_hosp Age-specific Proportion of hospitalizations that are 
#' severe based on NEJM definition
#' @param v_p_icu_s_hosp Age-specific Probability of going to ICU given that 
#' you are severe hospitalization 
#' @param v_p_icu_ns_hosp Age-specific Probability of going to ICU given that 
#' you are non-severe hospitalization
#' @param m_r_exit_hns Age-specific rate of non-severe hospital exit
#' @param m_r_exit_hs Age-specific rate of severe hospital exit
#' @param m_r_exit_icu Age-specific rate of ICU hospital exit
#' @param m_sigma_hns Age-specific sigma (gamma distrib) non-severe hosp dur
#' @param m_sigma_hs  Age-specific sigma (gamma distrib) severe hosp dur
#' @param m_sigma_icu Age-specific sigma (gamma distrib) ICU hosp dur
#' @param n_hhsize Household size (integer)
#' @param r_tau Transmission rate in the household model
#' @param r_sigma Daily rate of progression of exposed individuals 
#' (Latent period). Default to the first entry of the age-specific v_sigma
#' @param r_gamma Daily rate of recovery of infectious individuals 
#' (Infectiousness period). Default to the first entry of the age-specific 
#' v_gamma
#' @param r_omega Waning rate for the household model (default to the first 
#' entry of the age-specific v_omega)
#' Default will assign a reduction factor of 1 (i.e., no social distancing).
#' @return 
#' List of all parameters 
#' @export
load_params_init <- function(
  n_t       = 30,
  time_step = 1,
  comp      = TRUE,
  v_init_age_grps = c(0, 5, 15, 25, 45, 55, 65, 70),
  v_inf_init_ages = NULL,
  ctry = "Mexico", 
  ste  = "Mexico City", 
  cty  = "Mexico City",
  l_contact_info = l_contact_matrices,
  # v_beta   = rep(0.0536, n_ages), #rep(0.036, n_ages)
  r_birth  = 0, # daily birth rate (crude)
  r_beta   = 0.0536, #rep(0.036, n_ages) # transmission rate
  v_sigma  = rep(1/3, length(v_init_age_grps)), # infectious (Latent period)
  v_gamma  = rep(1/3.1, length(v_init_age_grps)), # recover (Infectiousness)
  v_omega  = rep(0, length(v_init_age_grps)), # immunity (Waning); for the household model, we assume that the first entry corresponds to the household omega
  l_nu_exp1_dx = add_period(l_period_def = NULL, 
                                l_period_add = make_period(
                                        functional_form = "constant",
                                        time_start = 0,
                                        time_stop = n_t,
                                        val_start = 0,
                                        val_end   = 0)),
  l_nu_exp2_dx = add_period(l_period_def = NULL, 
                            l_period_add = make_period(
                                    functional_form = "constant",
                                    time_start = 0,
                                    time_stop = n_t,
                                    val_start = 1/12,
                                    val_end   = 1/12)),
  l_nu_inf1_dx = add_period(l_period_def = NULL, 
                            l_period_add = make_period(
                              functional_form = "constant",
                              time_start = 0,
                              time_stop = n_t,
                              val_start = 0,
                              val_end   = 0)),
  l_nu_inf2_dx = add_period(l_period_def = NULL, 
                            l_period_add = make_period(
                              functional_form = "constant",
                              time_start = 0,
                              time_stop = n_t,
                              val_start = 1/12,
                              val_end   = 1/12)),
  l_phi_exp1_dx = add_period(l_period_def = NULL, 
                             l_period_add = make_period(
                              functional_form = "constant",
                              time_start = 0,
                              time_stop = n_t,
                              val_start = 0,
                              val_end   = 0)),
  l_phi_exp2_dx = add_period(l_period_def = NULL, 
                             l_period_add = make_period(
                              functional_form = "constant",
                              time_start = 0,
                              time_stop = n_t,
                              val_start = 0,
                              val_end   = 0)),
  l_phi_inf1_dx = add_period(l_period_def = NULL, 
                             l_period_add = make_period(
                              functional_form = "constant",
                              time_start = 0,
                              time_stop = n_t,
                              val_start = 0,
                              val_end   = 0)),
  l_phi_inf2_dx = add_period(l_period_def = NULL, 
                             l_period_add = make_period(
                              functional_form = "constant",
                              time_start = 0,
                              time_stop = n_t,
                              val_start = 0,
                              val_end   = 0)),
  l_cfr = get_const_multiage_list(n_t, c(0.000013, 0.000087, 0.000374, 0.001460, 0.007725, 0.026200, 0.051450, 0.134000)), # From: Excess_Mortality_Verity.csv
  l_ifr = get_const_multiage_list(n_t, c(0.00000805, 0.00004280, 0.00018925, 0.00084400, 0.00378000, 0.01262500, 0.02517500, 0.07800000)), # From: Excess_Mortality_Verity.csv
  l_interventions = add_intervention(interventions = NULL, 
                                     intervention = make_intervention(
                                                 intervention_type = "StatusQuo",
                                                 time_start = 0,
                                                 time_stop = n_t)),
  l_idx_scale_factor = get_const_multiage_list(n_t, rep(0.2, 8)),
  v_reduced_sus = rep(1, 8),
  v_alpha_dx            = NULL,  # reduction in mortality on detected infectious vs undetcted infectious
  # HOSPITALIZATION PARAMETERS
  # https://www.thelancet.com/action/showPdf?pii=S1473-3099%2820%2930243-7
  # https://www.nejm.org/doi/full/10.1056/NEJMoa2002032:
  l_p_hosp              = get_const_multiage_list(n_t, c(0.00001, 0.0002, 0.005, 0.035, 0.06, 0.095, 0.13, 0.17)),
  v_p_s_hosp            = c(0.1111, 0.1111, 0.1202, 0.1202, (0.1202+0.1746)/2, (0.1746+0.2875)/2, 0.2875, 0.2875),
  v_p_icu_s_hosp        = c(0.191, 0.191, 0.191, 0.191, 0.191, 0.191, 0.191, 0.191),
  v_p_icu_ns_hosp       = c(0.024, 0.024, 0.024, 0.024, 0.024, 0.024, 0.024, 0.024),
  # https://www.ajmc.com/newsroom/intensive-care-unit-usage-for-pneumonia-doubles-length-of-hospital-stay
  # https://www.cdc.gov/nchs/data/nhsr/nhsr116.pdf
  #65+      ICU 7.1 days; Hosp 4.6 days
  #45-64    ICU 7.3 days; Hosp 4.3 days
  #15-44    ICU 7.2 days; Hosp 4.2 days
  #<15      ICU 7.7 days; Hosp 3.1 days
  m_r_exit_hns          = c(1/3.1, 1/3.1, 1/((3.1+4.2)/2), 1/4.2, 1/4.3, 1/4.3, 1/4.6, 1/4.6),
  m_r_exit_icu          = c(1/7.7, 1/7.7, 1/((7.7+7.2)/2), 1/7.2, 1/7.3, 1/7.3, 1/7.1, 1/7.1),
  m_r_exit_hs           = 0.5*m_r_exit_hns + 0.5*m_r_exit_icu,
  m_sigma_hns           = c(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
  m_sigma_hs            = c(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
  m_sigma_icu           = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
  ### Household parameters
  n_hhsize = 3,
  r_tau    = 0,
  r_sigma  = v_sigma[1],
  r_gamma  = v_gamma[1],
  r_omega  = v_omega[1]
){
  
  ### Names of age groups
  v_names_ages <- ordered(c(paste(v_init_age_grps[-length(v_init_age_grps)], 
                                  (v_init_age_grps[-1]-1), sep = "-"), 
                            paste0(v_init_age_grps[length(v_init_age_grps)], "+")),
                          c(paste(v_init_age_grps[-length(v_init_age_grps)], 
                                  (v_init_age_grps[-1]-1), sep = "-"), 
                            paste0(v_init_age_grps[length(v_init_age_grps)], "+")))
  ### Number of age groups
  n_ages <- length(v_names_ages)
  
  ### Number of initial infectious individuals in each age group
  ## Error checking
  if(!is.null(v_inf_init_ages)){ # Check if user added this variable  
    if (length(v_inf_init_ages) != n_ages) {
      stop(paste0("Variable 'v_inf_init_ages' should be of length ", 
                  n_ages, ", same length as 'v_init_age_grps'"))
    }  
  } else{
    v_inf_init_ages <- rep(1, n_ages)
  }
  
  ### Time to start social distancing
  ## Error checking
#  if(!is.null(v_soc_dist_timing)){ # Check if user added this variable 
#    if (length(v_soc_dist_timing) != n_ages) {
#      stop(paste0("Variable 'v_soc_dist_timing' should be of length ", 
#                  n_ages, ", same length as 'v_init_age_grps'"))
#    }  
#  } else{ # Generate default values
#    v_soc_dist_timing <- rep(40, n_ages)
#  }
  
  ### Time to end social distancing
  ## Error checking
#  if(!is.null(v_soc_dist_timing_end)){ # Check if user added this variable 
#    if (length(v_soc_dist_timing_end) != n_ages) { # Check length is equal to the number of age groups
#      stop(paste0("Variable 'v_soc_dist_timing_end' should be of length ", 
#                  n_ages, ", same length as 'v_init_age_grps'"))
#    } else if (!all(v_soc_dist_timing_end > v_soc_dist_timing)){ # Check end date is greater than start date
#      v_ind_err_soc_dist_timing <- v_soc_dist_timing_end <= v_soc_dist_timing
#      stop(paste0("'v_soc_dist_timing_end' should be greater than 'v_soc_dist_timing' in age groups ", 
#                  paste(v_names_ages[v_ind_err_soc_dist_timing], collapse = ", ")))
#    }  
#  } else{ # Generate default values
#    v_soc_dist_timing_end <- rep(154, n_ages)
#  }
  
  ### Social distancing reduction factor
  ## Error checking
#  if(!is.null(v_soc_dist_factor)){
#    if (length(v_soc_dist_factor) != n_ages) {
#      stop(paste0("Variable 'v_soc_dist_factor' should be of length ", 
#                  n_ages, ", same length as 'v_init_age_grps'"))
#    }  else if (!all((v_soc_dist_factor <= 1) & (v_soc_dist_factor >= 0))){
#      v_ind_err_soc_dist_facors <- (v_soc_dist_factor > 1) | (v_soc_dist_factor < 0)
#      stop(paste0("'v_soc_dist_factor' should have values between 0 and 1 in age groups ", 
#                  paste(v_names_ages[v_ind_err_soc_dist_facors], collapse = ", ")))
#    }
#  } else{
#    v_soc_dist_factor <- rep(1, n_ages) 
#  }
  
  ### Reduction factor on COVID deaths from detected infectious (IDX)
  ## Error checking
  if(!is.null(v_alpha_dx)){
    if (length(v_alpha_dx) != n_ages) {
      stop(paste0("Variable 'v_alpha_dx' should be of length ", 
                  n_ages, ", same length as 'v_init_age_grps'"))
    }  else if (!all((v_alpha_dx <= 1) & (v_alpha_dx >= 0))){
      v_ind_err_alpha_dx <- (v_alpha_dx > 1) | (v_alpha_dx < 0)
      stop(paste0("'v_alpha_dx' should have values between 0 and 1 in age groups ", 
                  paste(v_names_ages[v_ind_err_alpha_dx], collapse = ", ")))
    }
  } else{
    v_alpha_dx <- rep(1, n_ages) 
  }
  
  ### Daily birth rate (crude)
  if(r_birth < 0){
    stop("'r_birth' should be a non-negative value")
  } else{
    v_birth <- c(r_birth = r_birth, rep(0, (n_ages-1)))
  }
  
  ### Create list of initial parameters
  l_params_init <- list(
    n_t       = n_t,
    time_step = time_step,
    comp      = comp,
    v_init_age_grps = v_init_age_grps,
    v_names_ages    = v_names_ages,
    n_ages          = n_ages,
    v_inf_init_ages = v_inf_init_ages,
    ctry = ctry,
    ste  = ste, 
    cty  = cty ,
    l_contact_info = l_contact_info,
    v_birth  = v_birth, 
    r_beta   = r_beta,
    l_nu_exp1_dx = l_nu_exp1_dx,
    l_nu_exp2_dx = l_nu_exp2_dx,
    l_nu_inf1_dx = l_nu_inf1_dx,
    l_nu_inf2_dx = l_nu_inf2_dx,
    l_phi_exp1_dx = l_phi_exp1_dx,
    l_phi_exp2_dx = l_phi_exp2_dx,
    l_phi_inf1_dx = l_phi_inf1_dx,
    l_phi_inf2_dx = l_phi_inf2_dx,
    v_sigma  = v_sigma, # infectious
    v_gamma  = v_gamma, # recover
    v_omega  = v_omega, # waning immunity
    l_cfr    = l_cfr,
    l_ifr    = l_ifr,
    ## intervention(s) parameters
    l_interventions        = l_interventions,
    ## Reduced susceptibility 
    v_reduced_sus = v_reduced_sus,
    ## Detection parameters
    l_idx_scale_factor     = l_idx_scale_factor,
    v_alpha_dx             = v_alpha_dx,
    ## hospitalization parameters
    l_p_hosp              = l_p_hosp, 
    v_p_s_hosp            = v_p_s_hosp, 
    v_p_icu_s_hosp        = v_p_icu_s_hosp, 
    v_p_icu_ns_hosp       = v_p_icu_ns_hosp, 
    m_r_exit_hns          = m_r_exit_hns, 
    m_r_exit_icu          = m_r_exit_icu, 
    m_r_exit_hs           = m_r_exit_hs,
    m_sigma_hns           = m_sigma_hns,
    m_sigma_hs            = m_sigma_hs,
    m_sigma_icu           = m_sigma_icu,
    n_hhsize = n_hhsize,
    r_tau    = r_tau,
    r_sigma  = r_sigma,
    r_gamma  = r_gamma,
    r_omega  = r_omega
  )
  return(l_params_init)
}

#' Load all parameters
#'
#' \code{load_all_params} loads all parameters for the SC-COSMO model with 
#' infectious detected states from multiple sources and creates a list.
#'
#' @param l_params_init List with initial set of parameters
#' @param file.init String with the location and name of the file with initial set of parameters
#' @param file.mort String with the location and name of the file with mortality data
#' @return 
#' A list of all parameters used for the decision model.
#' @export
load_all_params <- function(l_params_init = load_params_init(), 
                            file.init = NULL,
                            file.mort = NULL){ # User defined
  #### Load initial set of initial parameters from .csv file ####
  # if(!is.null(file.init)) {
  #   l_params_init <- load(file = file.init)
  if(is.null(l_params_init)) {
    l_params_init <- load_params_init()
  } else{
    l_params_init <- l_params_init
  }
  # list2env(l_params_init, envir = .GlobalEnv)
  #### Error checking ####
  
  #### All-cause age-specific mortality from .csv file ####
  l_params_all <- with(as.list(l_params_init), {
    #### General setup ####
    ### Days at which the model will be evaluated
    v_times <- seq(0, n_t, by = time_step)
    
    ### Names of strategies
    v_names_str  <- c("Do nothing", "Social distancing")  # CEA strategies
    n_str        <- length(v_names_str) # Number of strategies
    
    ### Number of compartments for exposed (E) and infectious (I)
    n_exp_states <- 3 # Might move to load_params_init() function later on
    n_inf_states <- 2 # Might move to load_params_init() function later on
    
    ### Name of states for exposed (E), exposed detected (EDX), infectious (I), and infectious detected (IDX) 
    ### by severity class
    v_names_exp_l1_states    <- paste("E1", letters[seq(1, n_exp_states)], sep = "")
    v_names_exp_l2_states    <- paste("E2", letters[seq(1, n_exp_states)], sep = "")
    v_names_expidx_l1_states <- paste("EDX1", letters[seq(1, n_exp_states)], sep = "")
    v_names_expidx_l2_states <- paste("EDX2", letters[seq(1, n_exp_states)], sep = "")
    v_names_inf_l1_states    <- paste("I1", letters[seq(1, n_inf_states)], sep = "")
    v_names_inf_l2_states    <- paste("I2", letters[seq(1, n_inf_states)], sep = "")
    v_names_infidx_l1_states <- paste("IDX1", letters[seq(1, n_inf_states)], sep = "")
    v_names_infidx_l2_states <- paste("IDX2", letters[seq(1, n_inf_states)], sep = "")
    v_names_exp_states       <- c(v_names_exp_l1_states, v_names_exp_l2_states)
    v_names_expidx_states    <- c(v_names_expidx_l1_states, v_names_expidx_l2_states)
    v_names_exptot_states    <- c(v_names_exp_states, v_names_expidx_states)
    v_names_inf_states       <- c(v_names_inf_l1_states, v_names_inf_l2_states)
    v_names_infidx_states    <- c(v_names_infidx_l1_states, v_names_infidx_l2_states)
    v_names_inftot_states    <- c(v_names_inf_states, v_names_infidx_states)
    v_names_dx_states        <- c(v_names_expidx_states, v_names_infidx_states)
    
    ### Name of age-group-specific states for exposed (E), infectious (I), and infectious detected (IDX)
    v_names_exp_l1_states_ages    <- paste0(rep(v_names_exp_l1_states, each = n_ages), 
                                            "_", v_names_ages)
    v_names_exp_l2_states_ages    <- paste0(rep(v_names_exp_l2_states, each = n_ages), 
                                            "_", v_names_ages)
    v_names_expidx_l1_states_ages <- paste0(rep(v_names_expidx_l1_states, each = n_ages), 
                                            "_", v_names_ages)
    v_names_expidx_l2_states_ages <- paste0(rep(v_names_expidx_l2_states, each = n_ages), 
                                            "_", v_names_ages)
    v_names_inf_l1_states_ages    <- paste0(rep(v_names_inf_l1_states, each = n_ages), 
                                            "_", v_names_ages)
    v_names_inf_l2_states_ages    <- paste0(rep(v_names_inf_l2_states, each = n_ages), 
                                            "_", v_names_ages)
    v_names_infidx_l1_states_ages <- paste0(rep(v_names_infidx_l1_states, each = n_ages), 
                                            "_", v_names_ages)
    v_names_infidx_l2_states_ages <- paste0(rep(v_names_infidx_l2_states, each = n_ages), 
                                            "_", v_names_ages)
    
    v_names_states <- c("S", 
                        v_names_exp_l1_states   ,
                        v_names_exp_l2_states   ,
                        v_names_expidx_l1_states,
                        v_names_expidx_l2_states,
                        v_names_inf_l1_states   ,
                        v_names_inf_l2_states   ,
                        v_names_infidx_l1_states,
                        v_names_infidx_l2_states,
                        "R",
                        "D",
                        "Itot",
                        "IDXtot",
                        "EDXtot",
                        "DXtot",
                        "Dcov",
                        "DcovDX")
    n_states   <- length(v_names_states)
    v_names_states_ages <- paste0(rep(v_names_states, each = n_ages), 
                                  "_", v_names_ages)
    n_states_ages <- length(v_names_states_ages)
    
    ### Alive states
    v_names_states_alive <- v_names_states[1:(1+(n_exp_states*2*2) + (n_inf_states*2*2) + 1)] # exposed times diagnosed (2) and severity (2) + infectious times diagnosed (2) and severity (2)
    n_states_alive       <- length(v_names_states_alive)
    v_names_states_ages_alive <- paste0(rep(v_names_states_alive, each = n_ages), 
                                        "_", v_names_ages)
    
    #### Household  model parameters ####
    ### Vectors names of household model
    v_names_exp_hh    <- paste("E", letters[seq(1, n_exp_states)], sep = "")
    v_names_inf_hh    <- paste("I", letters[seq(1, n_inf_states)], sep = "")
    v_names_rec_hh    <- c("R")
    v_names_states_hh <- c("S", 
                           v_names_exp_hh,
                           v_names_inf_hh,
                           v_names_rec_hh)
    
    n_states_hh      <- length(v_names_states_hh) # Not including the number of IDX
    
    ### Calaculate number of states for household model
    n_hh_mod <- factorial(n_hhsize + (n_states_hh-1))/(factorial(n_hhsize)*factorial(n_states_hh-1)) 
    
    ### Matrix with possible combinations of household members
    df_possibilities <- gen_hh_n(n_hhsize = n_hhsize, 
                                 v_names_states_hh = v_names_states_hh)
    m_possibilities <- as.matrix(df_possibilities)
    
    #### Houshold specific naming vectors ####
    v_hh_names <- as.matrix(tidyr::unite(df_possibilities, col = "names", sep = ""))
    # paste("HH",m_possibilities[, 1:n_states], m_possibilities[,2], sep = "")
    ### Names of household members by class names
    v_names_hh  <- paste("HH", v_hh_names, sep = "")
    ### Derivative names of household members by class
    v_names_dhh <- paste("dHH", v_hh_names, sep = "")
    
    ### All disease states including household states ####
    v_names_states_all <- c(v_names_states_ages, v_names_hh)
    n_states_all       <- length(v_names_states_all)

    
    #### Initial state ####
    ### Array with initial states
    l_init_pop <- obtain_init_pop(v_init_age_grps = v_init_age_grps,
                                  v_inf_init_ages = v_inf_init_ages,
                                  v_names_ages, n_ages, 
                                  v_names_states, n_states,
                                  ctry = ctry, 
                                  ste = ste, 
                                  cty = cty)
    a_states_init <- l_init_pop$a_states
    
    #### WAIFW data ####
    ######### CHANGE THIS #########
    if(ctry == "Italy"){
      WAIFW_data <- read.csv("data-raw/Mossong_2008_ITALY_REDUCED.csv", header=FALSE)
      m_waifw <- as.matrix(WAIFW_data) # CHANGE name of WAIFW!!
      dimnames(m_waifw) <- list(v_names_ages, v_names_ages)
      n_avg_ct <- NULL
      v_r_mort <-  as.matrix(read.csv("data-raw/Italy_WHO_Background_Mortality_2016_REDUCED.csv", 
                                      header=FALSE)) #as.matrix(MORT_data)
      #########
    } else {
      v_r_mort <- l_init_pop$v_r_mort
      ### Overall contact matrix
      m_waifw <- l_contact_info$m_contact
      ### Average number of total contacts
      n_avg_ct <- l_contact_info$n_avg_ct
      ### Setting specific contact matrices
      m_waifw_work            <- l_contact_info$m_contact_work
      m_waifw_school          <- l_contact_info$m_contact_school
      m_waifw_other_locations <- l_contact_info$m_contact_other_locations
      m_waifw_home            <- l_contact_info$m_contact_home
      dimnames(m_waifw) <- dimnames(m_waifw_work) <- dimnames(m_waifw_school) <- dimnames(m_waifw_other_locations) <- dimnames(m_waifw_home)  <- list(v_names_ages, v_names_ages)
      
    }
    ### Calculate total number of within household contacts by age group
    v_hh_contacts <- rowSums(m_waifw_home)
    names(v_hh_contacts) <- v_names_ages
    
    l_waifws_all_internal <- gen_all_waifws(l_interventions, 
                                            l_contact_info, 
                                            v_names_ages, 
                                            n_t)
    l_intervention_effects_all_internal <- gen_all_intervention_effects(l_interventions, 
                                                                        v_names_ages, 
                                                                        n_t)
    m_betas_all_internal <- gen_all_betas(l_interventions, 
                  l_contact_info, 
                  v_names_ages, 
                  n_t,
                  r_beta)
    
    #### Initialize state vector ####
    startingI  <- sum(a_states_init[, "I1a"])/sum(a_states_init)
    ### Household state vector
    v_HH0 <- numeric(length = n_hh_mod)
    names(v_HH0) <- v_names_hh
    v_HH0[1] <- 1 - n_hhsize*startingI
    v_HH0[m_possibilities[, "Ia"] == 1][1] <- n_hhsize*startingI
    v_HH0 <- v_HH0*sum(a_states_init)/n_hhsize
    
    ## Transform to vector
    v_states_init <- c(as.vector(a_states_init), # Community state names and ages
                       v_HH0)                    # Household state names
    names(v_states_init) <- v_names_states_all
    
    ### Indexing vectors for household model
    ## Column index for susceptibles
    v_index_hh_sus <- c(1)
    ## Column index for exposed
    v_index_hh_exp <- which(v_names_states_hh %in% v_names_exp_hh)
    ## Column index for infectious
    v_index_hh_inf <- which(v_names_states_hh %in% v_names_inf_hh)
    ## Column index for recovered
    v_index_hh_rec <- n_states_hh
    
    #### Houshold matrices ####
    l_transition_matrices <- gen_household_matrices_mc_seir(n_hhsize = n_hhsize,
                                                            n_hh_mod = n_hh_mod,
                                                            v_names_states_hh = v_names_states_hh,
                                                            v_names_exp_hh = v_names_exp_hh,
                                                            v_names_inf_hh = v_names_inf_hh,
                                                            v_names_rec_hh = v_names_rec_hh,
                                                            v_index_hh_sus = v_index_hh_sus,
                                                            v_index_hh_exp = v_index_hh_exp,
                                                            v_index_hh_inf = v_index_hh_inf,
                                                            v_index_hh_rec = v_index_hh_rec,
                                                            df_possibilities = df_possibilities,
                                                            r_sigma = r_sigma,
                                                            r_gamma = r_gamma
                                                            )
    m_comm_trans <- l_transition_matrices$m_comm_trans
    m_hh_trans   <- l_transition_matrices$m_hh_trans
    m_hh_prog    <- l_transition_matrices$m_hh_prog
    m_hh_recov   <- l_transition_matrices$m_hh_recov
    m_hh_waning  <- l_transition_matrices$m_hh_waning
    
    #### Detection rates ####
    #### Fill case and screen detection data structures
    v_nu_exp1_dx <- gen_time_varying(l_period_def = l_nu_exp1_dx, 
                                         max_time     = n_t + 1)
    v_nu_exp2_dx <- gen_time_varying(l_period_def = l_nu_exp2_dx, 
                                     max_time     = n_t + 1)
    v_nu_inf1_dx <- gen_time_varying(l_period_def = l_nu_inf1_dx, 
                                     max_time     = n_t + 1)
    v_nu_inf2_dx <- gen_time_varying(l_period_def = l_nu_inf2_dx, 
                                     max_time     = n_t + 1)
    
    v_phi_exp1_dx <- gen_time_varying(l_period_def = l_phi_exp1_dx, 
                                     max_time     = n_t + 1)
    v_phi_exp2_dx <- gen_time_varying(l_period_def = l_phi_exp2_dx, 
                                     max_time     = n_t + 1)
    v_phi_inf1_dx <- gen_time_varying(l_period_def = l_phi_inf1_dx, 
                                     max_time     = n_t + 1)
    v_phi_inf2_dx <- gen_time_varying(l_period_def = l_phi_inf2_dx, 
                                     max_time     = n_t + 1)
    
    #m_p_hosp <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    m_p_hosp <- list()
    for (i in 1:n_ages) {
      #m_p_hosp[i, ] <- gen_time_varying(l_p_hosp[[i]], max_time = n_t + 1)
      m_p_hosp[[i]] <- gen_time_varying(l_p_hosp[[i]], max_time = n_t + 1)
      
    }
    
    #m_cfr <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    m_cfr <- list()
    for (i in 1:n_ages) {
      #m_cfr[i, ] <- gen_time_varying(l_cfr[[i]], max_time = n_t + 1)
      m_cfr[[i]] <- gen_time_varying(l_cfr[[i]], max_time = n_t + 1)
    }

    #m_ifr <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    m_ifr <- list()
    for (i in 1:n_ages) {
      #m_ifr[i, ] <- gen_time_varying(l_ifr[[i]], max_time = n_t + 1)
      m_ifr[[i]] <- gen_time_varying(l_ifr[[i]], max_time = n_t + 1)
    }
    
    #m_idx_scale_factor <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    m_idx_scale_factor <- list()
    for (i in 1:n_ages) {
      #m_idx_scale_factor[i, ] <- gen_time_varying(l_idx_scale_factor[[i]], max_time = n_t + 1)
      m_idx_scale_factor[[i]] <- gen_time_varying(l_idx_scale_factor[[i]], max_time = n_t + 1)
    }
        
    ### Symptom rates
    m_r_exp1_sx <- matrix(0, nrow = n_ages, ncol = n_exp_states)
    m_r_inf1_sx <- matrix(0, nrow = n_ages, ncol = n_inf_states)
    
    sx_gamma_params <- gamma_params(mu = 5, sigma = 2.9, scale = F)
    for (i in 1:(n_exp_states + n_inf_states)) {
      prev_frac <- pgamma(q=i-1, shape=sx_gamma_params$shape, rate=sx_gamma_params$rate)
      end_frac  <- pgamma(q=i, shape=sx_gamma_params$shape, rate=sx_gamma_params$rate)
      frac_sx   <- (end_frac-prev_frac)/(1-prev_frac)
      r_sx <- -1*log(1-frac_sx)
      if (i <= n_exp_states) {
        m_r_exp1_sx[, i] <- r_sx
      } 
      else if ((i>n_exp_states) & (i<(n_exp_states + n_inf_states))) {
        j <- i - n_exp_states
        m_r_inf1_sx[, j] <- r_sx
      } else if (i == (n_exp_states + n_inf_states)) {
        prev_frac <- pgamma(q=i+2, shape=sx_gamma_params$shape, rate=sx_gamma_params$rate)
        end_frac  <- pgamma(q=i+3, shape=sx_gamma_params$shape, rate=sx_gamma_params$rate)
        frac_sx   <- (end_frac-prev_frac)/(1-prev_frac)
        r_sx <- -1*log(1-frac_sx)
        m_r_inf1_sx[, n_inf_states] <- r_sx
      }
    }
    
    #### Create list with all parameters ####
    l_params_all <- list(
      v_times      = v_times,
      v_names_str  = v_names_str,
      n_str        = n_str,
      v_names_ages = v_names_ages, 
      n_ages       = n_ages, 
      n_t          = n_t, 
      n_exp_states = n_exp_states,
      n_inf_states = n_inf_states,
      v_names_exp_l1_states    = v_names_exp_l1_states    ,
      v_names_exp_l2_states    = v_names_exp_l2_states    ,
      v_names_expidx_l1_states = v_names_expidx_l1_states ,
      v_names_expidx_l2_states = v_names_expidx_l2_states ,
      v_names_inf_l1_states    = v_names_inf_l1_states    ,
      v_names_inf_l2_states    = v_names_inf_l2_states    ,
      v_names_infidx_l1_states = v_names_infidx_l1_states ,
      v_names_infidx_l2_states = v_names_infidx_l2_states ,
      v_names_exp_states       = v_names_exp_states       ,      
      v_names_expidx_states    = v_names_expidx_states    ,
      v_names_exptot_states    = v_names_exptot_states    ,
      v_names_inf_states       = v_names_inf_states       ,
      v_names_infidx_states    = v_names_infidx_states    ,
      v_names_inftot_states    = v_names_inftot_states    ,
      v_names_dx_states        = v_names_dx_states        ,
      v_names_states = v_names_states,
      n_states = n_states,
      ### All names of community model
      v_names_states_ages = v_names_states_ages,
      n_states_ages       = n_states_ages,
      ### All states of both community and household model
      v_names_states_all = v_names_states_all,
      n_states_all       = n_states_all,
      ### Alive states in community model
      v_names_states_alive = v_names_states_alive,
      n_states_alive = n_states_alive,
      v_names_states_ages_alive = v_names_states_ages_alive,
      v_states_init = v_states_init,
      v_r_mort = v_r_mort,
      m_waifw                 = m_waifw,
      m_waifw_work            = m_waifw_work            ,
      m_waifw_school          = m_waifw_school          ,
      m_waifw_other_locations = m_waifw_other_locations ,
      m_waifw_home            = m_waifw_home            ,
      l_waifws_all_internal   = l_waifws_all_internal   ,
      v_hh_contacts           = v_hh_contacts, 
      l_intervention_effects_all_internal = l_intervention_effects_all_internal,
      m_betas_all_internal = m_betas_all_internal,
      n_avg_ct = n_avg_ct,
      country = ctry,
      state   = ste,
      county  = cty,
      r_beta  = r_beta,
      v_nu_exp1_dx   = v_nu_exp1_dx,
      v_nu_exp2_dx   = v_nu_exp2_dx,
      v_nu_inf1_dx   = v_nu_inf1_dx,
      v_nu_inf2_dx   = v_nu_inf2_dx,
      v_phi_exp1_dx  = v_phi_exp1_dx,
      v_phi_exp2_dx  = v_phi_exp2_dx,
      v_phi_inf1_dx  = v_phi_inf1_dx,
      v_phi_inf2_dx  = v_phi_inf2_dx,
      m_r_exp1_sx    = m_r_exp1_sx,
      m_r_inf1_sx    = m_r_inf1_sx,
      m_p_hosp       = m_p_hosp,
      m_cfr          = m_cfr,
      m_ifr          = m_ifr,
      m_idx_scale_factor = m_idx_scale_factor,
      ### Household model parameters
      v_names_states_hh = v_names_states_hh,
      v_names_exp_hh    = v_names_exp_hh,
      v_names_inf_hh    = v_names_inf_hh,
      n_states_hh       = n_states_hh, 
      n_hh_mod          = n_hh_mod,
      df_possibilities  = df_possibilities,
      m_possibilities   = m_possibilities,
      v_hh_names        = v_hh_names, 
      v_names_hh        = v_names_hh,  
      v_names_dhh       = v_names_dhh,
      v_index_hh_sus    = v_index_hh_sus,
      v_index_hh_exp    = v_index_hh_exp,
      v_index_hh_inf    = v_index_hh_inf,
      v_index_hh_rec    = v_index_hh_rec,
      m_comm_trans      =  m_comm_trans,
      m_hh_trans        =  m_hh_trans  ,
      m_hh_prog         =  m_hh_prog   ,
      m_hh_recov        =  m_hh_recov  ,
      m_hh_waning       =  m_hh_waning
    )
    return(l_params_all)
  }
  )
  
  l_params_all <- c(l_params_all, 
                    l_params_init) # Add initial set of parameters
  return(l_params_all)
}

#' Update parameters
#'
#' \code{update_param_list} is used to update list of all parameters with new 
#' values for specific parameters.
#'
#' @param l_params_all List with all parameters of decision model
#' @param params_updated Parameters for which values need to be updated
#' @return 
#' A list with all parameters updated.
#' @export
update_param_list <- function(regen_internals_flag = TRUE, l_params_all, params_updated){
  ### Verify if `params_updated` is a list, if not, make it a list
  if (typeof(params_updated) != "list"){
    params_updated <- split(unname(params_updated), names(params_updated)) #converte the named vector to a list
  }
  ### Error check: verify that all elements of `params_updated` are included in `l_params_all`
  ### There are reasons in calibration to pass arguments through so have turned this from
  ### stop to warn and commented out
  ### may be better to add a flag about whether to be strict with this or not
#  if(sum((names(params_updated) %in% names(l_params_all))) != length(params_updated)){
#    warn("Not all variables in `params_updated` are included in `l_params_all`")
#  }
  l_params_all <- modifyList(l_params_all, params_updated) #update the values
  if(regen_internals_flag == TRUE) {
    l_params_all$l_waifws_all_internal <- gen_all_waifws(
                                            l_params_all$l_interventions, 
                                            l_params_all$l_contact_info, 
                                            l_params_all$v_names_ages, 
                                            l_params_all$n_t
                                            )
    l_params_all$l_intervention_effects_all_internal <- gen_all_intervention_effects(
                                                          l_params_all$l_interventions, 
                                                          l_params_all$v_names_ages, 
                                                          l_params_all$n_t
                                                          )
    l_params_all$m_betas_all_internal <- gen_all_betas(
                                          l_params_all$l_interventions, 
                                          l_params_all$l_contact_info, 
                                          l_params_all$v_names_ages, 
                                          l_params_all$n_t,
                                          l_params_all$r_beta
                                        )
    
    
    l_params_all$r_omega  = l_params_all$v_omega[1]
    
    l_params_all$v_nu_exp1_dx <- gen_time_varying(l_period_def = l_params_all$l_nu_exp1_dx, 
                                     max_time     = l_params_all$n_t + 1)
    l_params_all$v_nu_exp2_dx <- gen_time_varying(l_period_def = l_params_all$l_nu_exp2_dx, 
                                     max_time     = l_params_all$n_t + 1)
    l_params_all$v_nu_inf1_dx <- gen_time_varying(l_period_def = l_params_all$l_nu_inf1_dx, 
                                     max_time     = l_params_all$n_t + 1)
    l_params_all$v_nu_inf2_dx <- gen_time_varying(l_period_def = l_params_all$l_nu_inf2_dx, 
                                     max_time     = l_params_all$n_t + 1)
    
    l_params_all$v_phi_exp1_dx <- gen_time_varying(l_period_def = l_params_all$l_phi_exp1_dx, 
                                      max_time     = l_params_all$n_t + 1)
    l_params_all$v_phi_exp2_dx <- gen_time_varying(l_period_def = l_params_all$l_phi_exp2_dx, 
                                      max_time     = l_params_all$n_t + 1)
    l_params_all$v_phi_inf1_dx <- gen_time_varying(l_period_def = l_params_all$l_phi_inf1_dx, 
                                      max_time     = l_params_all$n_t + 1)
    l_params_all$v_phi_inf2_dx <- gen_time_varying(l_period_def =l_params_all$l_phi_inf2_dx, 
                                      max_time     = l_params_all$n_t + 1)  
    
    # l_params_all$m_p_hosp <- matrix(0, nrow = l_params_all$n_ages, ncol = l_params_all$n_t + 1)
    # for (i in 1:l_params_all$n_ages) {
    #   l_params_all$m_p_hosp[i, ] <- gen_time_varying(l_params_all$l_p_hosp[[i]], max_time = l_params_all$n_t + 1)
    # }
    # 
    # l_params_all$m_cfr <- matrix(0, nrow = l_params_all$n_ages, ncol = l_params_all$n_t + 1)
    # for (i in 1:l_params_all$n_ages) {
    #   l_params_all$m_cfr[i, ] <- gen_time_varying(l_params_all$l_cfr[[i]], max_time = l_params_all$n_t + 1)
    # }
    # 
    # l_params_all$m_ifr <- matrix(0, nrow = l_params_all$n_ages, ncol = l_params_all$n_t + 1)
    # for (i in 1:l_params_all$n_ages) {
    #   l_params_all$m_ifr[i, ] <- gen_time_varying(l_params_all$l_ifr[[i]], max_time = l_params_all$n_t + 1)
    # }
    
    #m_p_hosp <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    l_params_all$m_p_hosp <- list()
    for (i in 1:l_params_all$n_ages) {
      #m_p_hosp[i, ] <- gen_time_varying(l_p_hosp[[i]], max_time = n_t + 1)
      l_params_all$m_p_hosp[[i]] <- gen_time_varying(l_params_all$l_p_hosp[[i]], max_time = l_params_all$n_t + 1)
      
    }
    
    #m_cfr <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    l_params_all$m_cfr <- list()
    for (i in 1:l_params_all$n_ages) {
      #m_cfr[i, ] <- gen_time_varying(l_cfr[[i]], max_time = n_t + 1)
      l_params_all$m_cfr[[i]] <- gen_time_varying(l_params_all$l_cfr[[i]], max_time = l_params_all$n_t + 1)
    }
    
    #m_ifr <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    l_params_all$m_ifr <- list()
    for (i in 1:l_params_all$n_ages) {
      #m_ifr[i, ] <- gen_time_varying(l_ifr[[i]], max_time = n_t + 1)
      l_params_all$m_ifr[[i]] <- gen_time_varying(l_params_all$l_ifr[[i]], max_time = l_params_all$n_t + 1)
    }
    
    #m_idx_scale_factor <- matrix(0, nrow = n_ages, ncol = n_t + 1)
    l_params_all$m_idx_scale_factor <- list()
    for (i in 1:l_params_all$n_ages) {
      #m_idx_scale_factor[i, ] <- gen_time_varying(l_idx_scale_factor[[i]], max_time = n_t + 1)
      l_params_all$m_idx_scale_factor[[i]] <- gen_time_varying(l_params_all$l_idx_scale_factor[[i]], max_time = l_params_all$n_t + 1)
    }
  }
  return(l_params_all)
}

#' Obtain initial population
#'
#' \code{obtain_init_pop} generates the number of people at each compartment by
#' age group.
#'
#' @param v_init_age_grps Vector with initial ages of age groups
#' @param v_inf_init_ages Vector with number of individuals in first infectious 
#' class for each age group
#' @param v_names_ages Vector with age groups
#' @param n_ages Number of age groups
#' @param v_names_states Vector with names of compartment states
#' @param n_states NUmber of compartment states
#' @param ctry Country
#' @param ste State
#' @param cty County or municipality
#' @param year Year
#' @return 
#' An array with the initial population at each class by age group.
#' @export
obtain_init_pop <- function(v_init_age_grps = c(0, 5, 15, 25, 45, 55, 65, 70),
                            v_inf_init_ages,
                            v_names_ages, n_ages, 
                            v_names_states, n_states,
                            ctry = "Italy", 
                            ste = NULL, 
                            cty = NULL, 
                            year = 2020){
  ### CHANGE THIS!
  if(ctry == "Italy"){
    POP_data   <- read.csv("data-raw/Italy_2019_REDUCED.csv", header=FALSE)  
    v_pop_init <- as.matrix(POP_data)
  } else{
    df_setting_data <- obtain_setting_data(v_init_age_grps, ctry, ste, cty, year)
    v_pop_init <- df_setting_data$population
    v_r_mort   <- df_setting_data$mort_rate
  }
  
  ### Generate array with initial population
  a_states <- array(data = 0, 
                    dim = list(n_ages,    # Number of age groups
                               n_states), # Number of infectiousness states
                    dimnames = list(v_names_ages, v_names_states))
  a_states[, "S"]   <- v_pop_init - v_inf_init_ages
  a_states[, "I1a"] <- v_inf_init_ages
  return(list(a_states = a_states,
              v_r_mort = v_r_mort))
}

#' Obtain initial population
#'
#' \code{obtain_setting_data} generates the number of people at each compartment by
#' age group.
#'
#' @param v_init_age_grps Vector with initial ages of age groups
#' @param ctry Country
#' @param ste State
#' @param cty County or municipality
#' @param year Year
#' @return 
#' An array with the initial population at each class by age group.
#' @export
obtain_setting_data <- function(v_init_age_grps,
                                ctry = "Mexico", 
                                ste = "Mexico City", 
                                cty = "Mexico City", 
                                year = 2020){
  df_setting_data <- df_pop_state_cty_age %>%
    dplyr::filter(country == ctry & state == ste & county == cty) %>%
    mutate(AgeGrp = cut(age, breaks = c(v_init_age_grps, Inf), 
                        include.lowest = TRUE, right = FALSE)) %>%
    group_by(country, state, county, AgeGrp) %>%
    summarise(population = sum(population), # Population by age group
              tot_pop    = mean(tot_pop),   # Total population for given setting
              deaths     = sum(deaths),     # Total number of deaths
              mort_rate  = (deaths/population)/(365.25) # daily mortality rate
              # land_sqmi  = mean(land_sqmi), 
              # density    = mean(density)
              )
  
  return(df_setting_data)
  
}