rm(list = ls()) # to clean the workspace

#### 02.1 Load packages and functions ####
#### 02.1.1 Load packages and functions ####
library(deSolve)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(grid) # For watermark
# library(sccosmo)
devtools::load_all(".")
library(sccosmoData)
library(stringi)
library(usmap)
library(countrycode)
library(dampack)

######## Initial setup ######
choose_country <- "Mexico"
#choose_state   <- "Ciudad de Mexico"
choose_state   <- "Hidalgo"
df_density <- get_densities(country = choose_country, 
              state = choose_state)

l_contact_matrices <- get_contact_matrix(country = choose_country, 
                                         state = choose_state, 
                                         density = df_density$density)
df_mort_state_cty_age_temp <- get_lifetables(country = choose_country,
                                             state = choose_state)
df_mort_state_cty_age <- df_mort_state_cty_age_temp %>%
  dplyr::mutate(state = choose_state,
                county = choose_state)
df_pop_state_cty_age_temp <- get_population_ages(country = choose_country, 
                                                 state   = choose_state)
df_pop_state_cty_age <- df_pop_state_cty_age_temp %>%
  dplyr::mutate(state = choose_state,
                county = choose_state)

df_pop_state_cty_age <- dplyr::left_join(df_pop_state_cty_age, 
                                         df_mort_state_cty_age)

df_pop_state_cty_age <- df_pop_state_cty_age %>%
  dplyr::filter(!is.na(dx)) %>%
  dplyr::mutate(deaths = (dx/lx)*age_pop) #%>%
#  dplyr::rename(population = age_pop)

v_params_calib <- c(r_beta                = 0.20, 
                    r_soc_dist_factor     = 0.75, 
                    r_nu_exp2_dx_lb       = 0.100, # 0.900,  
                    r_nu_exp2_dx_ub       = 0.005, # 0.900,   
                    r_nu_exp2_dx_rate     = 0.250, 
                    n_nu_exp2_dx_mid      = 16
)

test_sq = FALSE
if (test_sq == FALSE ) {
  i1 <- make_intervention(intervention_type = "StatusQuo",
                          time_start = 0,
                          time_stop = 30)
  i2 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = 30,
                          time_stop = 40,
                          intervention_factor = v_params_calib["r_soc_dist_factor"],
                          intervention_change_rate = 0.5)
  i3 <- make_intervention(intervention_type = "SocialDistancingLinear",
                          time_start = 40,
                          time_stop  = 55,
                          intervention_factor = v_params_calib["r_soc_dist_factor"],
                          intervention_factor_end = 0.80)
  i4 <- make_intervention(intervention_type = "SocialDistancing",
                          time_start = 55,
                          time_stop = 240,
                          intervention_factor = v_params_calib["r_soc_dist_factor"],
                          intervention_change_rate = 0.5,
                          resume_school = TRUE,
                          school_intervention_factor = 1)
  #                        resume_school = FALSE)
  i5 <- make_intervention(intervention_type = "StatusQuo",
                          time_start = 240,
                          time_stop = 1000)
  
  l_interventions <- add_intervention(interventions = NULL, intervention = i1)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i2)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i3)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i4)
  l_interventions <- add_intervention(interventions = l_interventions, intervention = i5)
} else {
  i1 <- make_intervention(intervention_type = "StatusQuo",
                          time_start = 0,
                          time_stop = 1000)
  
  l_interventions <- add_intervention(interventions = NULL, intervention = i1)
  
}
n_t <- 200
l_params_init <- load_params_init(n_t = n_t, # Number of days
                                  ctry = choose_country,
                                  ste  = choose_state,
                                  cty  = choose_state,
                                  l_contact_info = l_contact_matrices,
                                  # v_omega = rep(1/10, 8),
                                  # r_birth = (18/1000)/365.25,
                                  r_beta = v_params_calib["r_beta"],
                                  l_nu_exp2_dx = add_period(l_period_def = NULL, 
                                                            l_period_add = make_period(
                                                              functional_form = "general logit",
                                                              time_start = 0,
                                                              time_stop = n_t,
                                                              val_start = as.numeric(v_params_calib["r_nu_exp2_dx_lb"]),
                                                              val_end   = as.numeric(v_params_calib["r_nu_exp2_dx_ub"]),
                                                              v_logit_change_rate = as.numeric(v_params_calib["r_nu_exp2_dx_rate"]),
                                                              v_logit_change_mid  = as.numeric(v_params_calib["n_nu_exp2_dx_mid"]))),
                                  l_nu_inf2_dx = add_period(l_period_def = NULL, 
                                                            l_period_add = make_period(
                                                              functional_form = "general logit",
                                                              time_start = 0,
                                                              time_stop = n_t,
                                                              val_start = as.numeric(v_params_calib["r_nu_exp2_dx_lb"]),
                                                              val_end   = as.numeric(v_params_calib["r_nu_exp2_dx_ub"]),
                                                              v_logit_change_rate = as.numeric(v_params_calib["r_nu_exp2_dx_rate"]),
                                                              v_logit_change_mid  = as.numeric(v_params_calib["n_nu_exp2_dx_mid"]))),
                                  v_inf_init_ages   = c(0, 0, 0, 1, 0, 0, 0, 0),
                                  l_idx_scale_factor = get_const_multiage_list(n_t, rep(0, 8)),
                                  l_interventions = l_interventions,
                                  r_tau = v_params_calib["r_beta"]*1.5, 
                                  # v_cfr = rep(0, 8),
                                  # v_ifr = rep(0, 8),
                                  n_hhsize = 3,
                                  m_r_exit_hns          = rep(1/20, 8),
                                  m_r_exit_icu          = rep(1/30, 8),
                                  m_r_exit_hs           = 0.5*rep(1/20, 8) + 0.5*rep(1/30, 8),
                                  m_sigma_hns           = rep(10, 8),
                                  m_sigma_hs            = rep(10, 8),
                                  m_sigma_icu           = rep(7,  8)
                                  #,
#                                  v_omega = rep(1/180, 8)
)
# l_params_init$v_birth
# l_params_init$v_omega
# v_init_age_grps <- l_params_init$v_init_age_grps


l_params_all <- load_all_params(l_params_init = l_params_init)
# l_params_all$v_r_mort <- rep(0, 8)
# l_params_all$m_r_exp1_sx <- matrix(0, nrow = 8, ncol = 3)
# l_params_all$m_r_inf1_sx <- matrix(0, nrow = 8, ncol = 2)
# list2env(l_params_init, envir = .GlobalEnv)
# list2env(l_params_all, envir = .GlobalEnv)
sum(l_params_all$v_states_init[1:15])

system.time(
  l_out_cosmo <- cosmo(l_params_all = l_params_all)
)
# diagnostics(deSolve::ode(y = l_params_all$v_states_init, 
#                          times = l_params_all$v_times, 
#                          func = ifelse(l_params_all$comp, 
#                                        yes = cosmo_dXdt_comp, 
#                                        no = cosmo_dxdt), 
#                          parms = l_params_all))
View(l_out_cosmo$df_out_cosmo)
# v_pop <- as.numeric(l_out_cosmo$df_out_cosmo[17, -1])
rowSums(l_out_cosmo$df_out_cosmo[1:n_t, 2:(l_out_cosmo$l_params_all$n_states_ages+1)] < 0)
rowSums(l_out_cosmo$df_out_cosmo[1:n_t, (l_out_cosmo$l_params_all$n_states_ages+2):ncol(l_out_cosmo$df_out_cosmo)]<0.0)
# View(l_out_cosmo$df_out_cosmo[23:26, 2:(l_out_cosmo$l_params_all$n_states_ages+1)])
# View(l_out_cosmo$df_out_cosmo[, 2:(l_out_cosmo$l_params_all$n_states_ages+1)])
# View(l_out_cosmo$df_out_cosmo[1:40, (l_out_cosmo$l_params_all$n_states_ages+2):ncol(l_out_cosmo$df_out_cosmo)])


colnames(l_out_cosmo$df_out_cosmo)[l_out_cosmo$df_out_cosmo[23, ]<0]

apply(l_out_cosmo$df_out_cosmo, 1, 
      function(x) colnames(l_out_cosmo$df_out_cosmo)[x<0])

### Plot population size over time
plot_popsize_totals(l_out_cosmo,
                    only_all = TRUE, print_plot = FALSE)
plot_popsize_totals(l_out_cosmo, 
                    only_all = FALSE, print_plot = FALSE)
# ggsave(filename = "figs/02_cumulative_infections_total_mex_nothing.png", width = 8, height = 6)
### Plot total cumulative number of diagnosed infections over time
# df_InfIDXcumtot <- data.frame(Outcome = "Cumulative diagnosed infections",
#                               Intervention = "Do Nothing", 
#                               time = calc_infidxcum_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                               value = calc_infidxcum_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                               check.names = F)
plot_dxcum_totals(l_out_cosmo,
                  only_all = TRUE, print_plot = FALSE)
plot_dxcum_totals(l_out_cosmo,
                  only_all = FALSE, print_plot = FALSE)
plot_dxcum_totals(l_out_cosmo,
                  proportion = TRUE,
                  only_all = TRUE, print_plot = FALSE)
plot_dxcum_totals(l_out_cosmo,
                  proportion = TRUE, age_proportion = TRUE,
                  only_all = FALSE, print_plot = FALSE)
# plot_infidxcum_totals(l_out_cosmo,
#                    proportion = TRUE, age_proportion = TRUE,
#                    only_all = TRUE, print_plot = FALSE)

### Plot total cumulative number of diagnosed infections over time
# df_InfIDXInctot <- data.frame(Outcome = "Incidental diagnosed infections",
#                               Intervention = "Do Nothing", 
#                               time = calc_infidxinc_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                               value = calc_infidxinc_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                               check.names = F)
plot_dxinc_totals(l_out_cosmo,
                  only_all = TRUE, print_plot = FALSE)
plot_dxinc_totals(l_out_cosmo,
                  only_all = FALSE, print_plot = FALSE)
plot_dxinc_totals(l_out_cosmo,
                  proportion = TRUE,
                  only_all = TRUE, print_plot = FALSE)
plot_dxinc_totals(l_out_cosmo,
                  proportion = TRUE, age_proportion = TRUE,
                  only_all = FALSE, print_plot = FALSE)
# plot_infidxcum_totals(l_out_cosmo,
#                    proportion = TRUE, age_proportion = TRUE,
#                    only_all = TRUE, print_plot = FALSE)
### Plot total cumulative number of infections over time
# df_Infcumtot <- data.frame(Outcome = "Cumulative infections",
#                            Intervention = "Do Nothing", 
#                            time = calc_infcum_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                            value = calc_infcum_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                            check.names = F)
plot_infcum_totals(l_out_cosmo,
                   only_all = TRUE, print_plot = FALSE)
plot_infcum_totals(l_out_cosmo,
                   only_all = FALSE, print_plot = FALSE)
plot_infcum_totals(l_out_cosmo,
                   proportion = TRUE,
                   only_all = TRUE, print_plot = FALSE)
plot_infcum_totals(l_out_cosmo,
                   proportion = TRUE, age_proportion = TRUE,
                   only_all = FALSE, print_plot = FALSE)
# plot_infcum_totals(l_out_cosmo,
#                    proportion = TRUE, age_proportion = TRUE,
#                    only_all = TRUE, print_plot = FALSE)
# ggsave(filename = "figs/02_cumulative_infections_total_mex_nothing.png", width = 8, height = 6)

### Plot total incidental number of infections over time
# df_InfInctot <- data.frame(Outcome = "Incident infections",
#                            Intervention = "Do Nothing", 
#                            time = calc_infinc_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                            value = calc_infinc_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                            check.names = F)
plot_infinc_totals(l_out_cosmo,
                   only_all = TRUE, print_plot = FALSE)
plot_infinc_totals(l_out_cosmo,
                   only_all = FALSE, print_plot = FALSE)
plot_infinc_totals(l_out_cosmo,
                   proportion = TRUE,
                   only_all = TRUE, print_plot = FALSE)
plot_infinc_totals(l_out_cosmo,
                   proportion = TRUE, age_proportion = TRUE,
                   only_all = FALSE, print_plot = FALSE)
# plot_infcum_totals(l_out_cosmo,
#                    proportion = TRUE, age_proportion = TRUE,
#                    only_all = TRUE, print_plot = FALSE)
# ggsave(filename = "figs/02_cumulative_infections_total_mex_nothing.png", width = 8, height = 6)


### Plot total number of infections over time
# df_Inftot <- data.frame(Outcome = "Infections",
#                         Intervention = "Do Nothing", 
#                         time = calc_inf_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                         value = calc_inf_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                         check.names = F)
plot_inf_totals(l_out_cosmo,
                only_all = TRUE, print_plot = FALSE)
plot_inf_totals(l_out_cosmo,
                only_all = FALSE, print_plot = FALSE)
plot_inf_totals(l_out_cosmo,
                proportion = TRUE, 
                only_all = FALSE, print_plot = FALSE)
plot_inf_totals(l_out_cosmo,
                proportion = TRUE, age_proportion = TRUE,
                only_all = FALSE, print_plot = FALSE)
# plot_inf_totals(l_out_cosmo,
#                 proportion = TRUE, age_proportion = TRUE,
#                 only_all = TRUE, print_plot = FALSE)
# ggsave(filename = "figs/02_infections_total_mex_nothing.png", width = 8, height = 6)

### Plot total number of diagnosed infections over time
# df_InfIDXtot <- data.frame(Outcome = "Prevalent diagnosed Infections",
#                         Intervention = "Do Nothing", 
#                         time = calc_infidx_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                         value = calc_infidx_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                         check.names = F)
plot_dx_totals(l_out_cosmo,
               only_all = TRUE, print_plot = FALSE)
plot_dx_totals(l_out_cosmo,
               only_all = FALSE, print_plot = FALSE)
plot_dx_totals(l_out_cosmo,
               proportion = TRUE, 
               only_all = FALSE, print_plot = FALSE)
plot_dx_totals(l_out_cosmo,
               proportion = TRUE, age_proportion = TRUE,
               only_all = FALSE, print_plot = FALSE)
# plot_inf_totals(l_out_cosmo,
#                 proportion = TRUE, age_proportion = TRUE,
#                 only_all = TRUE, print_plot = FALSE)
# ggsave(filename = "figs/02_infections_total_mex_nothing.png", width = 8, height = 6)

### Plot total number of COVID19 deaths over time
# df_Dcov <- data.frame(Outcome = "COVID deaths",
#                       Intervention = "Do Nothing", 
#                       time = calc_deaths_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                       value = calc_deaths_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                       check.names = F)
plot_deaths_totals(l_out_cosmo = l_out_cosmo,
                   only_all = TRUE, print_plot = FALSE)
plot_deaths_totals(l_out_cosmo = l_out_cosmo,
                   only_all = FALSE, print_plot = FALSE)
plot_deaths_totals(l_out_cosmo = l_out_cosmo,
                   proportion = TRUE, 
                   only_all = FALSE, print_plot = FALSE)
plot_deaths_totals(l_out_cosmo = l_out_cosmo,
                   proportion = TRUE, age_proportion = TRUE, 
                   only_all = FALSE, print_plot = FALSE)

### Plot total number of COVID19 deaths from diagnosed infections over time
# df_DcovDX <- data.frame(Outcome = "COVID deaths from diagnosed infections",
#                       Intervention = "Do Nothing", 
#                       time = calc_deathsdx_totals(l_out_cosmo = l_out_cosmo)[, 1],
#                       value = calc_deathsdx_totals(l_out_cosmo = l_out_cosmo)[, 10], 
#                       check.names = F)
plot_deathsdx_totals(l_out_cosmo = l_out_cosmo,
                     only_all = TRUE, print_plot = FALSE)
plot_deathsdx_totals(l_out_cosmo = l_out_cosmo,
                     only_all = FALSE, print_plot = FALSE)
plot_deathsdx_totals(l_out_cosmo = l_out_cosmo,
                     proportion = TRUE, 
                     only_all = FALSE, print_plot = FALSE)
plot_deathsdx_totals(l_out_cosmo = l_out_cosmo,
                     proportion = TRUE, age_proportion = TRUE, 
                     only_all = FALSE, print_plot = FALSE)


### Plot symptomatics
plot_sx_totals(l_out_cosmo = l_out_cosmo,
                     only_all = TRUE, print_plot = FALSE)
plot_sx_totals(l_out_cosmo = l_out_cosmo,
                     only_all = FALSE, print_plot = FALSE)
plot_sx_totals(l_out_cosmo = l_out_cosmo,
                     proportion = TRUE, 
                     only_all = FALSE, print_plot = FALSE)
plot_sx_totals(l_out_cosmo = l_out_cosmo,
                     proportion = TRUE, age_proportion = TRUE, 
                     only_all = FALSE, print_plot = FALSE)

### Plot E and Is
plot_expinf_totals(l_out_cosmo = l_out_cosmo,
                   only_all = TRUE, print_plot = FALSE)
plot_expinf_totals(l_out_cosmo = l_out_cosmo,
               only_all = FALSE, print_plot = FALSE)
plot_expinf_totals(l_out_cosmo = l_out_cosmo,
               proportion = TRUE, 
               only_all = FALSE, print_plot = FALSE)
plot_expinf_totals(l_out_cosmo = l_out_cosmo,
               proportion = TRUE, age_proportion = TRUE, 
               only_all = FALSE, print_plot = FALSE)

#### Reproduction number ####
calc_basic_reproduction_number_kr(l_out_cosmo = l_out_cosmo, 
                                  v_time = 2:100)
calc_reproduction_number_wt(l_out_cosmo = l_out_cosmo, 
                            v_time = 2:100,
                            nsim_chosen = 200)
#calc_reproduction_number_wt(l_out_cosmo = l_out_cosmo, 
#                            v_time = 2:(l_out_cosmo$l_params_all$n_t-3))

#### Hospitalizations ####
l_out_hosp <- prep_dx_hospitalizations(l_out_cosmo = l_out_cosmo)
df_hosp_prev <- calc_dx_hosp(l_hosp = l_out_hosp, l_out_cosmo = l_out_cosmo)
plot(All ~ time, df_hosp_prev)
