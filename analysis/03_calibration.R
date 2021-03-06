#03_calibration.R
#------------------------------------------------------------------------------#
# This script calibrates the Mexican Stanford-CIDE COronavirus Simulation MOdel#
# (SC-COSMO) to time series of symptom detected cases of coronavirus disease   #
# 2019 (Covid-19) as target using a Bayesian approach with by finding the      #
# maximum-a-posteriori (MAP) using optimization algorithms and a sample from   #
# the posterior distribution using the Incremental Mixture Importance          #
# Samping (IMIS) algorithm.                                                    #
#                                                                              # 
# Authors:                                                                     #
#     - Fernando Alarid-Escudero, PhD, <fernando.alarid@cide.edu>              # 
#     - Jeremy D Goldhaber-Fiebert, PhD                                        #
#     - Jason Andrews, MD, MPH                                                 #
#     - Valeria Gracia Olvera, MSc, <valeria.gracia@cide.edu>                  #
#------------------------------------------------------------------------------#

rm(list = ls()) # to clean the workspace

# Load packages, functions and data ---------------------------------------

## Packages ---------------------------------------------------------------
library(sccosmomcma)
library(deSolve)
library(tidyverse)
library(data.table)
library(DEoptim)
library(dampack)
library(reshape2)

## Parallel setup ---------------------------------------------------------
library(doParallel)
library(foreach)

## Functions --------------------------------------------------------------
source("R/03_target_functions.R")
source("R/03_calibration_functions.R")

## Load data --------------------------------------------------------------

load("data/df_age_first_case_mex.rda")    # age group of first case by state
load("data/df_NPIs_2020-12-29.RData")     # mobility mean chagepoints estimates by state
load("data/l_targets_2020-12-13.RData")   # targets generated

# Start calibration -------------------------------------------------------

t0 <- timestamp()

# Lag in the time series of number of infectious
n_lag_inf <- 14
n_lag_conf <- 0

# Save calibration data
save_data <- TRUE

# Run calibration for each type of model ----------------------------------

# Number of age groups
n_ages <- 8

# SUSC
reduced_sus <- 0.25 # based on : https://www.nature.com/articles/s41591-020-0962-9 abstract and figure 1 panel B
l_reduced_sus <- list(SQ = rep(1, n_ages), # SQ
                      SA = c(rep(reduced_sus, 2), rep(1, (n_ages-2)))) # SA


for(n_proj_type in c("SQ", "SA")){ # n_proj_type = "SA"
  
  cat("Running calibration for ", n_proj_type, "\n")

  # State-specific data -----------------------------------------------------
  
  # Set state to calibrate
  state_i <-  "MCMA"
  abbrev_state <- "MCMA"
  
  # Age-specific reduction factor on susceptibility
  v_reduced_sus <- l_reduced_sus[[n_proj_type]]
  
  # Average household size for Mexican states https://www.inegi.org.mx/temas/hogares/
  # Mexico City = 3.3, State of Mexico = 3.7 
  n_hhsize <- 3
  
  # Density
  density <- get_densities_mx() %>%
    filter(state == state_i) %>%
    select(density) %>%
    as.numeric()
  
  # Population 
  df_pop_state_age <- get_population_ages(country = "Mexico", 
                                          state   = state_i,
                                          county  = state_i)
  
  # Deaths
  df_mort_state_cty_age_temp <- get_lifetables(country = "Mexico",
                                               state   = state_i)
  
  df_mort_state_cty_age <- df_mort_state_cty_age_temp %>%
    dplyr::mutate(state  = state_i,
                  county = state_i)
  
  df_pop_state_cty_age_temp <- get_population_ages(country = "Mexico", 
                                                   state   = state_i)
  
  df_pop_state_cty_age <- df_pop_state_cty_age_temp %>%
    dplyr::mutate(state  = state_i,
                  county = state_i)
  
  df_pop_state_cty_age <- dplyr::left_join(df_pop_state_cty_age, 
                                           df_mort_state_cty_age)
  
  df_pop_state_cty_age <- df_pop_state_cty_age %>%
    dplyr::filter(!is.na(dx)) %>%
    dplyr::mutate(deaths = round((dx/lx)*age_pop,0)) 
  
  # Get state-specific contact matrices
  l_contact_matrices <- get_contact_matrix(country = "Mexico",
                                           state   = state_i,
                                           density = density)
  
  
  # Generate and plot targets -----------------------------------------------
  
  # Generate targets
  # l_targets <- gen_targets(v_states_calib = state_i,
  #                          n_time_stamp   = "2020-12-13",
  #                          n_date_last    = "2020-12-07")
  
  # List of start and end dates for outcomes
  l_dates_targets <- list(
    cases      = c(first(l_targets$cases$Date),
                   last(l_targets$cases$Date)),
    cases_inc  = c(first(l_targets$cases_inc$Date),
                   last(l_targets$cases_inc$Date)),
    deaths     = c(first(l_targets$deaths$Date),
                   last(l_targets$deaths$Date)),
    deaths_inc = c(first(l_targets$deaths_inc$Date),
                   last(l_targets$deaths_inc$Date))
  )
  
  # Generate time stamp
  n_time_stamp <- unique(l_targets$cases$time_stamp)
  
  # Generate last day for calibration
  n_date_last <- unique(l_targets$cases$DateLast)
  
  # Plot and save targets
  plot_targets(l_targets, save_plot = F)
  
  
  # Set calibration variables -----------------------------------------------
  
  # First day of calibration
  n_date_ini <- first(l_targets$cases_inc$Date)
  
  # Last day of calibration
  n_date_end_calib <- last(l_targets$cases_inc$Date)
  
  # Number of days to project for selected state
  n_t_calib <- as.numeric(n_date_end_calib - n_date_ini)
  
  ## Calibration parameters -------------------------------------------------
  
  # Vector of parameters
  v_params_calib <- c(r_beta              = 0.20,
                      r_tau               = 0.3,
                      r_soc_dist_factor   = 0.50,
                      r_soc_dist_factor_2 = 0.50,
                      r_soc_dist_factor_3 = 0.50,
                      r_soc_dist_factor_4 = 0.50,
                      r_soc_dist_factor_5 = 0.50,
                      r_nu_exp2_dx_lb     = 0.050,
                      r_nu_exp2_dx_ub     = 0.100,
                      r_nu_exp2_dx_rate   = 0.250,
                      n_nu_exp2_dx_mid    = 35
  )
  
  # Names and number of input parameters to be calibrated
  v_param_names  <- c("r_beta",
                      "r_tau",
                      "r_soc_dist_factor",
                      "r_soc_dist_factor_2",
                      "r_soc_dist_factor_3",
                      "r_soc_dist_factor_4",
                      "r_soc_dist_factor_5",
                      "r_nu_exp2_dx_lb",
                      "r_nu_exp2_dx_ub",
                      "r_nu_exp2_dx_rate",
                      "n_nu_exp2_dx_mid"
  )
  
  # Names and number of input parameters to be calibrated
  n_param <- length(v_param_names)
  
  # Vector with range on input search space
  v_lb <- get_bounds()$v_lb    # lower bound
  v_ub <- get_bounds()$v_ub    # upper bound
  
  # Number of calibration targets
  n_target_names <- "DXIncTot"
  n_target       <- length(n_target_names)
  
  # Interventions -----------------------------------------------------------
  
  ## Define intervention 1 timing
  n_date_NPI <- df_NPIs$Date_INT2[df_NPIs$state == state_i] 
  n_date0_NPI <- as.numeric(n_date_NPI - n_date_ini)
  
  ## Define intervention 2 timing
  n_date_NPI_2 <- df_NPIs$Date_INT3[df_NPIs$state == state_i] 
  n_date0_NPI_2 <- as.numeric(n_date_NPI_2 - n_date_ini)
  n_offset_NPI_2 <- n_date0_NPI_2 - n_date0_NPI
  
  # Define intervention 3 timing
  n_date_NPI_3 <- df_NPIs$Date_INT4[df_NPIs$state == state_i] - 21
  n_date0_NPI_3 <- as.numeric(n_date_NPI_3 - n_date_ini)
  n_offset_NPI_3 <- n_date0_NPI_3 - n_date0_NPI_2
  
  # Define intervention 4 timing
  n_date_NPI_4 <- df_NPIs$Date_INT5[df_NPIs$state == state_i] 
  n_date0_NPI_4 <- as.numeric(n_date_NPI_4 - n_date_ini)
  n_offset_NPI_4 <- n_date0_NPI_4 - n_date0_NPI_3
  
  # Define intervention 5 timing
  n_date_NPI_5 <- df_NPIs$Date_INT6[df_NPIs$state == state_i] 
  n_date0_NPI_5 <- as.numeric(n_date_NPI_5 - n_date_ini)
  n_offset_NPI_5 <- n_date0_NPI_5 - n_date0_NPI_4
  
  # Define status quo interventions
  i1 <- make_intervention(intervention_type = "StatusQuo",
                          time_start = 0,
                          time_stop  = n_date0_NPI + n_lag_inf)
  i2 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = n_date0_NPI + n_lag_inf,
                          time_stop  = n_date0_NPI + n_lag_inf + n_offset_NPI_2,
                          intervention_factor = v_params_calib["r_soc_dist_factor"],
                          intervention_change_rate = 0.5)
  i3 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = n_date0_NPI + n_lag_inf + n_offset_NPI_2,
                          time_stop  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3,
                          intervention_factor = v_params_calib["r_soc_dist_factor_2"],
                          intervention_change_rate = 0.5)
  i4 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3,
                          time_stop   = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + 
                            n_offset_NPI_4,
                          intervention_factor = v_params_calib["r_soc_dist_factor_3"],
                          intervention_change_rate = 0.5)
  i5 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + n_offset_NPI_4,
                          time_stop   = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + n_offset_NPI_4 + n_offset_NPI_5,
                          intervention_factor = v_params_calib["r_soc_dist_factor_4"],
                          intervention_change_rate = 0.5)
  i6 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start  = n_date0_NPI + n_lag_inf + n_offset_NPI_2 + n_offset_NPI_3 + n_offset_NPI_4  + n_offset_NPI_5,
                          time_stop   = n_t_calib + n_lag_inf + 1,
                          intervention_factor = v_params_calib["r_soc_dist_factor_5"],
                          intervention_change_rate = 0.5)
  
  ### Add interventions
  l_interventions <- add_intervention(interventions = NULL, intervention = i1)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i6)
  
  # Load parameters ---------------------------------------------------------
  
  # Age-group vector of first case
  v_inf_init_ages <- as.numeric(df_age_first_case_mex[df_age_first_case_mex$state == "Mexico City", 4:11])
  
  ### Time varying parameters
  l_params_init <- sccosmomcma::load_params_init(n_t = n_t_calib + n_lag_inf, # Number of days
                                             ctry = "Mexico",
                                             ste  = state_i,
                                             cty  = state_i,
                                             v_reduced_sus = v_reduced_sus,
                                             r_beta = v_params_calib["r_beta"],
                                             l_nu_exp2_dx = add_period(l_period_def = NULL,
                                                                       l_period_add = make_period(
                                                                         functional_form = "general logit",
                                                                         time_start = 0,
                                                                         time_stop = n_t_calib + n_lag_inf,
                                                                         val_start = as.numeric(v_params_calib["r_nu_exp2_dx_lb"]),
                                                                         val_end   = as.numeric(v_params_calib["r_nu_exp2_dx_ub"]),
                                                                         v_logit_change_rate = as.numeric(v_params_calib["r_nu_exp2_dx_rate"]),
                                                                         v_logit_change_mid  = as.numeric(v_params_calib["n_nu_exp2_dx_mid"]))),
                                             l_nu_inf2_dx = add_period(l_period_def = NULL,
                                                                       l_period_add = make_period(
                                                                         functional_form = "general logit",
                                                                         time_start = 0,
                                                                         time_stop = n_t_calib + n_lag_inf,
                                                                         val_start = as.numeric(v_params_calib["r_nu_exp2_dx_lb"]),
                                                                         val_end   = as.numeric(v_params_calib["r_nu_exp2_dx_ub"]),
                                                                         v_logit_change_rate = as.numeric(v_params_calib["r_nu_exp2_dx_rate"]),
                                                                         v_logit_change_mid  = as.numeric(v_params_calib["n_nu_exp2_dx_mid"]))),
                                             v_inf_init_ages  = v_inf_init_ages,
                                             l_contact_info  = l_contact_matrices,
                                             l_interventions = l_interventions,
                                             n_hhsize = n_hhsize,
                                             r_tau   = v_params_calib["r_tau"],
                                             r_omega = 0
  )
  
  # Load all parameter values
  l_params_all <- sccosmomcma::load_all_params(l_params_init = l_params_init)
  

  # Save image --------------------------------------------------------------
  
  if(save_data){
    save.image(file = paste0("data/MCMA_calibration_data_",n_proj_type,"_",n_time_stamp,".RData"))
  }
  
  # Find a maximum-a-posteriori: IMIS ---------------------------------------
  
  # Specify seed (for reproducible sequence of random numbers)
  set.seed(072218)
  
  # Number of random samples to obtain from the posterior distribution 
  n_resamp <- 1000 
  
  # log_lik_par(sample.prior(n_samp = 4),
  #             log_lik_offset = 0,#-1481.144,
  #             l_params_all = l_params_all,
  #             n_lag_inf  = n_lag_inf,
  #             n_lag_conf = n_lag_conf,
  #             l_dates_targets = l_dates_targets)
  # likelihood(sample.prior(n_samp = 4))
  # m_prior <- sample.prior(n_samp = 4)
  # system.time(
  #   v_llik <- likelihood(m_prior)
  # )
  # v_weights <- v_llik/sum(v_llik) 
  
  # Run IMIS
  system.time(
    l_fit_imis <- IMIS::IMIS(B        =  1000,       # incremental sample size at each iteration of IMIS
                             B.re     =  n_resamp,   # desired posterior sample size
                             number_k =  55,         # maximum number of iterations in IMIS
                             D        =  0)
  )
  
  ### Exploring posterior distribution --------------------------------------
  
  # Obtain posterior
  m_calib_post <- l_fit_imis$resample
  
  # Prior vs posterior distribution
  m_calib_post_eff <- m_calib_post
  m_calib_post_eff[, 3:7] <- 1 - m_calib_post[, 3:7]
  
  m_samp_prior <- sample.prior(n_resamp)
  m_samp_prior_eff <- m_samp_prior
  m_samp_prior_eff[, 3:7] <- 1-m_samp_prior_eff[, 3:7]
  
  # Ordered prior vs posterior distribution in terms of NPI effectiveness
  df_samp_prior_eff <- melt(cbind(PDF = "Prior", 
                                          as.data.frame(m_samp_prior_eff)), 
                                    variable.name = "Parameter")
  df_samp_post_imis_eff  <- melt(cbind(PDF = "Posterior IMIS", 
                                               as.data.frame(m_calib_post_eff)), 
                                         variable.name = "Parameter")
  df_samp_prior_post_eff <- bind_rows(df_samp_prior_eff, 
                                              df_samp_post_imis_eff)
  df_samp_prior_post_eff$PDF <- ordered(df_samp_prior_post_eff$PDF, 
                                                levels = c("Prior", "Posterior IMIS")) # "Posterior SIR", 
  gg_post_imis_eff <- ggplot(df_samp_prior_post_eff, 
                                     aes(x = value, y = ..density.., fill = PDF)) +
    facet_wrap(~Parameter, scales = "free", ncol = 3) +
    scale_x_continuous(breaks = number_ticks(6)) +
    geom_density(alpha=0.5) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom")
  
  #Save plot  
  ggsave(plot = gg_post_imis_eff,
         paste0("figs/03_posterior_vs_prior_marginals_",abbrev_state,"_", n_time_stamp,"_eff.jpg"), 
         width = 10, height = 8, dpi = 300)
  
  # Prior vs posterior distribution
  df_samp_prior <- melt(cbind(PDF = "Prior", 
                              as.data.frame(m_samp_prior)), 
                        variable.name = "Parameter")
  df_samp_post_imis  <- melt(cbind(PDF = "Posterior IMIS", 
                                   as.data.frame(m_calib_post)), 
                             variable.name = "Parameter")
  df_samp_prior_post <- bind_rows(df_samp_prior, 
                                  df_samp_post_imis)
  df_samp_prior_post$PDF <- ordered(df_samp_prior_post$PDF, 
                                    levels = c("Prior", "Posterior IMIS")) # "Posterior SIR", 
  gg_post_imis <- ggplot(df_samp_prior_post, 
                         aes(x = value, y = ..density.., fill = PDF)) +
    facet_wrap(~Parameter, scales = "free", ncol = 3) +
    scale_x_continuous(breaks = number_ticks(6)) +
    geom_density(alpha=0.5) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom")
  
  # Save plot
  ggsave(plot = gg_post_imis,
         paste0("figs/03_posterior_vs_prior_marginals_",abbrev_state,"_",n_proj_type,
                "_", n_time_stamp,".jpg"), 
         width = 10, height = 8, dpi = 300)
  
  
  ### Summary statistics of posterior distribution --------------------------
  
  # Compute posterior mean
  v_calib_post_mean <- colMeans(m_calib_post)
  
  # Compute posterior SE
  v_calib_post_se <- matrixStats::colSds(m_calib_post)
  
  # Compute posterior median and 95% credible interval
  m_calib_post_95cr <- matrixStats::colQuantiles(m_calib_post, 
                                                 probs = c(0.025, 0.5, 0.975))
  
  # Compute maximum-a-posteriori (MAP) as the mode of the sampled values
  v_calib_post_mode <- apply(l_fit_imis$resample, 2,
                             function(x) as.numeric(modeest::mlv(x, method = "shorth")[1]))
  
  v_map <- v_calib_post_mean
  
  # Summary statistics
  df_posterior_summ <- data.frame(
    Parameter = v_param_names,
    Mean      = v_calib_post_mean,
    m_calib_post_95cr,
    MAP       = v_calib_post_mode,
    check.names = FALSE)
  
  n_log_post_IMIS <- log_post(v_calib_post_mode,
                              l_params_all = l_params_all, 
                              n_lag_inf = n_lag_inf, 
                              n_lag_conf = n_lag_conf,
                              l_dates_targets = l_dates_targets)
  
  df_calib_post_map_IMIS <- data.frame(`Time stamp` = n_time_stamp,
                                       LastDay      = n_date_last,
                                       state = state_i,
                                       Method = "IMIS",
                                       value = n_log_post_IMIS,
                                       map = v_calib_post_mean,
                                       lb = m_calib_post_95cr[, 1],
                                       ub = m_calib_post_95cr[, 3],
                                       se = v_calib_post_se, 
                                       check.names = FALSE)
  
  rownames(df_calib_post_map_IMIS) <- v_param_names
  
  #### Save summary statistics ----------------------------------------------
  
  # As .RData
  save(df_posterior_summ, df_calib_post_map_IMIS, 
       file = paste0("output/03_summary_posterior_IMIS_",abbrev_state,"_",n_proj_type,"_",n_time_stamp,".RData"))
  
  # As .csv
  write.csv(df_posterior_summ, 
            file = paste0("tables/03_summary_posterior_IMIS_",abbrev_state,"_",n_proj_type,"_",n_time_stamp,".csv"),
            row.names = FALSE)
  
  # Save posterior and MAP
  save(l_fit_imis,
       df_calib_post_map_IMIS,
       m_calib_post,
       v_map,
       v_calib_post_mode,
       df_samp_prior_post,
       gg_post_imis,
       file = paste0("output/03_map_output_IMIS_",abbrev_state,"_",n_proj_type,"_",n_time_stamp,".RData"))
  
  t1 <- timestamp()
  
}
