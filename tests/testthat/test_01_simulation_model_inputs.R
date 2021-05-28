context("testing 02_simulation_model_functions.R")

# library(dplyr)    # For data manipulation
library(sccosmo)
library(dplyr)

#### Unit tests start ####
test_that("invalid initial inputs", {
  ### Initial infectious individuals in each age group
  expect_error(load_params_init(v_inf_init_ages = rep(1, 10)), 
               regexp = "Variable 'v_inf_init_ages' should be of length 8, same length as 'v_init_age_grps'"
               )
  expect_error(load_params_init(v_inf_init_ages = rep(1, 7)), 
               regexp = "Variable 'v_inf_init_ages' should be of length 8, same length as 'v_init_age_grps'"
  )
  expect_error(load_params_init(v_init_age_grps = c(0, 50),
                                v_inf_init_ages = rep(1, 8)), 
               regexp = "Variable 'v_inf_init_ages' should be of length 2, same length as 'v_init_age_grps'"
  )
  expect_error(load_params_init(v_init_age_grps = c(0, 50),
                                v_inf_init_ages = rep(1, 1)), 
               regexp = "Variable 'v_inf_init_ages' should be of length 2, same length as 'v_init_age_grps'"
  )
  ### Time to start social distancing
  # expect_error(load_params_init(v_soc_dist_timing = rep(30, 10)), 
  #              regexp = "Variable 'v_soc_dist_timing' should be of length 8, same length as 'v_init_age_grps'"
  # )
  # expect_error(load_params_init(v_soc_dist_timing = rep(30, 7)), 
  #              regexp = "Variable 'v_soc_dist_timing' should be of length 8, same length as 'v_init_age_grps'"
  # )
  
  ### Time to end social distancing
  # expect_error(load_params_init(v_soc_dist_timing_end = rep(30, 10)), 
  #              regexp = "Variable 'v_soc_dist_timing_end' should be of length 8, same length as 'v_init_age_grps'"
  # )
  # expect_error(load_params_init(v_soc_dist_timing_end = rep(30, 7)), 
  #              regexp = "Variable 'v_soc_dist_timing_end' should be of length 8, same length as 'v_init_age_grps'"
  # )
  # expect_error(load_params_init(v_soc_dist_timing     = rep(30, 8),
  #                               v_soc_dist_timing_end = c(rep(35, 4), rep(25, 4))), 
  #              regexp = "'v_soc_dist_timing_end' should be greater than 'v_soc_dist_timing' in age groups 45-54, 55-64, 65-69, 70+"
  # )
  ### Social distancing reduction factor
  # expect_error(load_params_init(v_soc_dist_factor = rep(1, 10)), 
  #              regexp = "Variable 'v_soc_dist_factor' should be of length 8, same length as 'v_init_age_grps'"
  # )
  # expect_error(load_params_init(v_soc_dist_factor = rep(1, 10)), 
  #              regexp = "Variable 'v_soc_dist_factor' should be of length 8, same length as 'v_init_age_grps'"
  # )
  # expect_error(load_params_init(v_soc_dist_factor = rep(1.3, 8)), 
  #              regexp = "'v_soc_dist_factor' should have values between 0 and 1 in age groups 0-4, 5-14, 15-24, 25-44, 45-54, 55-64, 65-69, 70+"
  # )
  # expect_error(load_params_init(v_soc_dist_factor = rep(-1, 8)), 
  #              regexp = "'v_soc_dist_factor' should have values between 0 and 1 in age groups 0-4, 5-14, 15-24, 25-44, 45-54, 55-64, 65-69, 70+"
  # )
  # expect_error(load_params_init(v_soc_dist_factor = c(rep(-1, 4), rep(0.05, 4))), 
  #              regexp = "'v_soc_dist_factor' should have values between 0 and 1 in age groups 0-4, 5-14, 15-24, 25-44"
  # )
})

# load inputs
# l_params_all <- load_all_params()

#### Unit tests start ####
# test_that("invalid inputs all", {
#   ### Test Mexico parameters 
#   n_pop_mex <- df_pop_state_cty_age_mx %>%
#     filter(country == "Mexico", state == "National", county == "National") %>%
#     summarise(tot_pop = tot_pop[1]) %>%
#     select(tot_pop) %>%
#     as.numeric()
#   l_params_all_mex <- load_all_params(l_params_init = load_params_init(ctry = "Mexico",
#                                                                        ste = "National",
#                                                                        cty = "National"))
#   expect_that(sum(l_params_all_mex$v_states_init[1:l_params_all_mex$n_ages]), 
#               is_equivalent_to(n_pop_mex))
#   
  # ### Test Mexico City parameters
  # n_pop_cdmx <- df_pop_state_cty_age_mx %>%
  #   filter(country == "Mexico", state == "Mexico City", county == "Mexico City") %>%
  #   summarise(tot_pop = tot_pop[1]) %>%
  #   select(tot_pop) %>%
  #   as.numeric()
  # l_params_all_cdmx <- load_all_params(l_params_init = load_params_init(ctry = "Mexico",
  #                                                                      ste = "Mexico City",
  #                                                                      cty = "Mexico City"))
  # expect_that(sum(l_params_all_cdmx$v_states_init[1:l_params_all_cdmx$n_ages]),
  #             is_equivalent_to(n_pop_cdmx))

  # ### Test Brazil parameters on São Paulo
  # n_pop_bra_sp <- df_pop_state_cty_age_bra %>%
  #   filter(country == "Brazil", state == "São Paulo", county == "São Paulo") %>%
  #   summarise(tot_pop = tot_pop[1]) %>%
  #   select(tot_pop) %>%
  #   as.numeric()
  # l_params_all_bra_sp <- load_all_params(l_params_init = load_params_init(ctry = "Brazil",
  #                                                                      ste = "São Paulo",
  #                                                                      cty = "São Paulo"))
  # expect_that(sum(l_params_all_bra_sp$v_states_init[1:l_params_all_bra_sp$n_ages]), 
  #             is_equivalent_to(n_pop_bra_sp))
# })
