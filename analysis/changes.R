reduced_sus <- 0.25 # based on : https://www.nature.com/articles/s41591-020-0962-9 abstract and figure 1 panel B
v_reduced_sus <- rep(1, n_ages) # SQ
v_reduced_sus <- c(rep(reduced_sus, 2), rep(1, (n_ages-2))) # SA

reduced_hosp <- 0.1 # based on : https://www.nature.com/articles/s41591-020-0962-9/figures/1 panel C
v_reduced_hosp <- rep(1, n_ages) # SQ
v_reduced_hosp <- c(rep(reduced_hosp, 2), rep(1, (n_ages-2))) # SA

# l_p_hosp <- get_const_multiage_list(n_t, c(0.00001, 0.0002, 0.005, 0.035, 0.06, 0.095, 0.13, 0.17)) * v_reduced_hosp # whatever we are doing, doing it there

reduced_mort <- 0.1
v_reduced_mort <- rep(1, n_ages) # SQ
v_reduced_mort <- c(rep(reduced_mort, 2), rep(1, (n_ages-2))) # SA