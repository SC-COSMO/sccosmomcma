#' The Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{cosmo} implements the Stanford-CIDE COronavirus Simulation MOdel 
#' (SC-COSMO), which is an age-structure Susceptible-Exposed-Infectious-Recovered
#' (SEIR) model of the coronavirus disease (COVID-19).
#' @param l_params_all List with all parameters of decision model
#' @return 
#' A data.frame with population for each class by age group over time.
#' @export
cosmo <- function(l_params_all){
  df_out_cosmo <- as.data.frame(deSolve::ode(y = l_params_all$v_states_init, 
                                             times = l_params_all$v_times, 
                                             func = ifelse(l_params_all$comp, 
                                                           yes = cosmo_dXdt_comp, 
                                                           no = cosmo_dXdt), 
                                             parms = l_params_all))
  l_out_cosmo <- list(df_out_cosmo = df_out_cosmo,
                      l_params_all = l_params_all)
  return(l_out_cosmo)
}

#' Derivatives of the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{cosmo_dXdt} computes the derivatives of the Stanford-CIDE COronavirus 
#' Simulation MOdel (SC-COSMO).
#'
#' @param time Time at which the derivatives will be computed
#' @param v_pop Vector with the population by class
#' @param l_params_all List with all parameters of SC-COSMO model
#' @return 
#' A list with a vector of derivatives for each class by age group.
#' @export
cosmo_dXdt <- function(time, v_pop, l_params_all) {
  # with(as.list(l_params_all), {
  ## Getting local variables off of l_params_all to improve performance
  n_ages                    <- l_params_all$n_ages
  n_states                  <- l_params_all$n_states
  v_names_ages              <- l_params_all$v_names_ages
  v_names_states            <- l_params_all$v_names_states
  n_states_ages             <- l_params_all$n_states_ages
  n_states_all              <- l_params_all$n_states_all
  v_states_init             <- l_params_all$v_states_init
  v_names_states_ages_alive <- l_params_all$v_names_states_ages_alive
  n_hhsize                  <- l_params_all$n_hhsize
  n_exp_states              <- l_params_all$n_exp_states
  n_inf_states              <- l_params_all$n_inf_states
  v_names_exp_l1_states     <- l_params_all$v_names_exp_l1_states
  v_names_exp_l2_states     <- l_params_all$v_names_exp_l2_states
  v_names_expidx_l1_states  <- l_params_all$v_names_expidx_l1_states
  v_names_expidx_l2_states  <- l_params_all$v_names_expidx_l2_states
  v_names_inf_l1_states     <- l_params_all$v_names_inf_l1_states
  v_names_inf_l2_states     <- l_params_all$v_names_inf_l2_states
  v_names_infidx_l1_states  <- l_params_all$v_names_infidx_l1_states
  v_names_infidx_l2_states  <- l_params_all$v_names_infidx_l2_states
  v_reduced_sus             <- l_params_all$v_reduced_sus
  v_names_hh                <- l_params_all$v_names_hh
  r_tau                     <- l_params_all$r_tau
  v_hh_contacts             <- l_params_all$v_hh_contacts
  v_index_hh_sus            <- l_params_all$v_index_hh_sus
  v_index_hh_exp            <- l_params_all$v_index_hh_exp
  v_index_hh_inf            <- l_params_all$v_index_hh_inf
  v_index_hh_rec            <- l_params_all$v_index_hh_rec
  m_possibilities           <- l_params_all$m_possibilities
  v_birth                   <- l_params_all$v_birth
  v_omega                   <- l_params_all$v_omega
  v_r_mort                  <- l_params_all$v_r_mort
  v_sigma                   <- l_params_all$v_sigma
  v_gamma                   <- l_params_all$v_gamma
  v_alpha_dx                <- l_params_all$v_alpha_dx
  m_comm_trans              <- l_params_all$m_comm_trans
  m_hh_trans                <- l_params_all$m_hh_trans
  m_hh_prog                 <- l_params_all$m_hh_prog
  m_hh_recov                <- l_params_all$m_hh_recov
  m_hh_waning_ones          <- l_params_all$m_hh_waning
  m_r_exp1_sx               <- l_params_all$m_r_exp1_sx
  m_r_inf1_sx               <- l_params_all$m_r_inf1_sx
  r_omega                   <- l_params_all$r_omega
  
    ## Convert the vector of incoming state values into a matrix, where 
    ## each column in the matrix contains those values for a common state 
    ## Each column now represents a model state for each cohort
    # To test: v_pop = v_states_init
    a_states <- array(data = v_pop, 
                      dim = list(n_ages,    # Number of age groups
                                 n_states), # Number of states
                      dimnames = list(v_names_ages, v_names_states))
    v_HH <- (v_pop[(n_states_ages+1):n_states_all]/sum(v_states_init[v_names_states_ages_alive]))*n_hhsize
    # if(sum(a_states<0)){
    #   warning("Entries in community model cannot be negative")
    #   print(time)
    #   stop("Entries in community model cannot be negative")
    # }
    # if(sum(v_HH < 0) == 0){
    #   print(time)
    #   # print(v_HH)
    #   # stop("Entries in community model cannot be negative")
    # }
    ### Detection rates
    #r_nu_exp1_dx <- general_logit(r_nu_exp1_dx_lb, r_nu_exp1_dx_ub, r_nu_exp1_dx_rate, n_nu_exp1_dx_mid, time)
    #r_nu_exp2_dx <- general_logit(r_nu_exp2_dx_lb, r_nu_exp2_dx_ub, r_nu_exp2_dx_rate, n_nu_exp2_dx_mid, time)
    r_nu_exp1_dx <- get_parameters(param = "v_nu_exp1_dx", time = time, l_params_all = l_params_all)
    r_nu_exp2_dx <- get_parameters(param = "v_nu_exp2_dx", time = time, l_params_all = l_params_all)
    
    # we assume that severity level (not whether someone is exposed vs. infectious) determines detectability w/o active screening
    #r_nu_inf1_dx <- r_nu_exp1_dx
    #r_nu_inf2_dx <- r_nu_exp2_dx 
    r_nu_inf1_dx <- get_parameters(param = "v_nu_inf1_dx", time = time, l_params_all = l_params_all)
    r_nu_inf2_dx <- get_parameters(param = "v_nu_inf2_dx", time = time, l_params_all = l_params_all) 
    #r_nu_inf1_dx <- general_logit(r_nu_inf1_dx_lb, r_nu_inf1_dx_ub, r_nu_inf1_dx_rate, n_nu_inf1_dx_mid, time)
    #r_nu_inf2_dx <- general_logit(r_nu_inf2_dx_lb, r_nu_inf2_dx_ub, r_nu_inf2_dx_rate, n_nu_inf2_dx_mid, time)
    
    r_phi_exp1_dx <- get_parameters(param = "v_phi_exp1_dx", time = time, l_params_all = l_params_all)
    r_phi_exp2_dx <- get_parameters(param = "v_phi_exp2_dx", time = time, l_params_all = l_params_all)
    r_phi_inf1_dx <- get_parameters(param = "v_phi_inf1_dx", time = time, l_params_all = l_params_all)
    r_phi_inf2_dx <- get_parameters(param = "v_phi_inf2_dx", time = time, l_params_all = l_params_all)
    
    
    m_nu_exp1_dx  <-  diag(rep(r_nu_exp1_dx, n_exp_states))
    m_nu_exp2_dx  <-  diag(rep(r_nu_exp2_dx, n_exp_states))
    m_nu_inf1_dx  <-  diag(rep(r_nu_inf1_dx, n_inf_states))
    m_nu_inf2_dx  <-  diag(rep(r_nu_inf2_dx, n_inf_states))
    
    m_phi_exp1_dx  <-  diag(rep(r_phi_exp1_dx, n_exp_states))
    m_phi_exp2_dx  <-  diag(rep(r_phi_exp2_dx, n_exp_states))
    m_phi_inf1_dx  <-  diag(rep(r_phi_inf1_dx, n_inf_states))
    m_phi_inf2_dx  <-  diag(rep(r_phi_inf2_dx, n_inf_states))
    
    ### Compute the force of infection including any intervention effects
    #m_Beta <- get_betas(time, l_params_all)
    m_Beta <- get_parameters(param = "m_Beta", time = time, l_params_all = l_params_all)
    
    ### Pull the time-varying case fatality rate
    v_cfr <- get_parameters(param = "v_cfr", time = time, l_params_all = l_params_all)
    
    ### Pull the time-varying reduction factor on detected infectious
    v_idx_scale_factor <- get_parameters(param = "v_idx_scale_factor", 
                                         time = time, 
                                         l_params_all = l_params_all)
    # print(v_idx_scale_factor)
    # if(sum(m_Beta<0)){
    #   stop("Negative m_Beta")
    # }
    
    ### Extract three one-column matrices from each of the matrix columns 
    ### to obtain all of the state values for each model states. Conveniently 
    ### organized by stock type
    
    ### Susceptibles
    S <- a_states[, "S"]
    
    #### Exposed
    ### Exposed l=1
    m_E1 <- a_states[, v_names_exp_l1_states]
    ### Exposed l=2
    m_E2 <- a_states[, v_names_exp_l2_states]
    
    #### Exposed detected (EDX)
    ### Exposed detected (EDX) l=1
    m_EDX1 <- a_states[, v_names_expidx_l1_states]
    ### Exposed detected (EDX) l=2
    m_EDX2 <- a_states[, v_names_expidx_l2_states]
    
    #### Infectious
    ### Infectious l=1
    m_I1 <- a_states[, v_names_inf_l1_states]
    ### Infectious l=2
    m_I2 <- a_states[, v_names_inf_l2_states]
    
    #### Infectious detected (IDX)
    ### Infectious detected (IDX) l=1
    m_IDX1 <- a_states[, v_names_infidx_l1_states]
    ### Infectious detected (IDX) l=2
    m_IDX2 <- a_states[, v_names_infidx_l2_states]
    
    ### Recovered
    R <- a_states[, "R"]
    
    ### Dead
    D <- a_states[, "D"]
    
    ### Aggregated exposed
    E <- rowSums(m_E1 + m_E2) #m_E[drop = FALSE]
    
    ### Aggregated infectious detected
    EDX <- rowSums(m_EDX1 + m_EDX2) #m_E[drop = FALSE]
    
    ### Aggregated infectious 
    I <- rowSums(m_I1 + m_I2) #m_I[drop = FALSE]
    
    ### Aggregated infectious detected
    IDX <- rowSums(m_IDX1 + m_IDX2) #m_I[drop = FALSE]
    
    ### Population alive by age
    N <- S + E + EDX + I + IDX + R
    
    ### Total population alive
    Ntot <- sum(N)
    
    ### Proportion of population in each age group
    v_prop_ages_curr <- N/Ntot
    
    ### Proportion of susceptibles in each age group
    v_prop_ages_curr_sus <- rep(0, n_ages)
    if(sum(S) > 0){
      v_prop_ages_curr_sus[S > 0] <- S[S > 0]/sum(S[S > 0])
    }
    
    ### Proportion of infections that are DX
    if(sum((I + IDX) == 0) == n_ages){
      p_dx <- 0
      v_p_dx <- rep(0, n_ages)
    } else{
      v_p_dx <- rep(0, n_ages)
      v_p_dx <- ifelse(IDX > 0,
                       IDX/(I + IDX),
                       rep(0, n_ages))
      # v_dx[(I + IDX) > 0] <- ((IDX[(I + IDX) > 0]/(I[(I + IDX) > 0] + IDX[(I + IDX) > 0]))) 
      # v_prop_ages_curr_inf <- rep(0, n_ages)
      # 
      # v_prop_ages_curr_inf[(I + IDX) > 0] <- (v_prop_ages_curr[(I + IDX) > 0]/sum(v_prop_ages_curr[(I + IDX) > 0]))
      # p_dx <- as.numeric(v_dx[(I + IDX) > 0]  %*% 
      #                      v_prop_ages_curr_inf[(I + IDX) > 0])  
      p_dx <- sum(IDX)/sum(IDX+I)
    }
    
    
    ### Adjust transmission matrices to account for detection
    # v_prop_ages_curr_idx          <- rep(0, n_ages)
    # v_prop_ages_curr_idx[IDX > 0] <- (v_prop_ages_curr[IDX > 0]/sum(v_prop_ages_curr[IDX > 0]))
    # p_idx_scale_factor <- as.numeric(v_idx_scale_factor[IDX > 0] %*% 
    #                                    v_prop_ages_curr_idx[IDX > 0])
    p_idx_scale_factor <- ifelse(sum(IDX) == 0,
                                 0,
                                 as.numeric(v_idx_scale_factor[IDX > 0] %*% IDX[IDX > 0]/sum(IDX)))
    r_tau_avg     <- r_tau*(1 - p_dx) + r_tau*p_idx_scale_factor*p_dx
    # print(p_idx_scale_factor)
    
    ### Household infection
    n_contacts_hh <- as.numeric(v_hh_contacts %*% v_prop_ages_curr)
    household_infection_rate <- gen_household_transmission_mc_seir(r_tau = r_tau_avg,
                                                                   n_hhsize = n_hhsize,
                                                                   n_contacts_hh = n_contacts_hh,
                                                                   v_HH = v_HH,
                                                                   v_index_hh_sus = v_index_hh_sus,
                                                                   v_index_hh_exp = v_index_hh_exp,
                                                                   v_index_hh_inf = v_index_hh_inf,
                                                                   v_index_hh_rec = v_index_hh_rec,
                                                                   m_possibilities = m_possibilities)
    n_household_infection_rate <- as.numeric(household_infection_rate)
    ## Force of infection from household to community
    v_lambda_hh_community <- n_household_infection_rate * v_prop_ages_curr_sus
    
    ### Compute one-column matrix of forces of infection lambda
    v_lambda <- (m_Beta %*% (I/Ntot) + (m_Beta %*% (IDX/Ntot))*v_idx_scale_factor) * v_reduced_sus    
    # cbind(S, round(v_lambda*S, 0), round(v_lambda_hh_community*N, 0))
    
    ### ODE model
    ## S
    dS_dt  <- v_birth*Ntot +
      v_omega * R - 
      (v_lambda + v_r_mort) * S - # Leaving S from community transmission
      v_lambda_hh_community * S   # Leaving S from household transmission
    
    ### E
    ## E1
    m_dE1_dt <- cbind(S = (v_lambda* S) + v_lambda_hh_community * S, 
                      n_exp_states*v_sigma*m_E1[, -n_exp_states]) - # Incoming E1
      matrix((v_r_mort + n_exp_states*v_sigma), nrow = n_ages, ncol = n_exp_states) * m_E1 - 
      m_E1 %*% (m_nu_exp1_dx + m_phi_exp1_dx) - 
      m_r_exp1_sx * m_E1
    ## E2
    m_dE2_dt <- cbind(E2 = 0, n_exp_states*v_sigma*m_E2[, -n_exp_states]) +          # Incoming E2
      m_r_exp1_sx * m_E1 -                                                           # Incoming E2
      matrix((v_r_mort + n_exp_states*v_sigma), nrow = n_ages, ncol = n_exp_states) * m_E2 - 
      m_E2 %*% (m_nu_exp2_dx + m_phi_exp2_dx)
    
    ### EDX 
    ## EDX1
    m_dEDX1_dt <- m_E1 %*% (m_nu_exp1_dx + m_phi_exp1_dx)  +                                  # Incoming EDX1
      cbind(EDX1 = 0, n_exp_states*v_sigma*m_EDX1[, -n_exp_states]) -                         # Incoming EDX1
      matrix((v_r_mort + n_exp_states*v_sigma), nrow = n_ages, ncol = n_exp_states) * m_EDX1
    ## EDX2
    m_dEDX2_dt <- m_E2 %*% (m_nu_exp2_dx + m_phi_exp2_dx)  +                                  # Incoming EDX2
      cbind(EDX2 = 0, n_exp_states*v_sigma*m_EDX2[, -n_exp_states]) -                         # Incoming EDX2
      matrix((v_r_mort + n_exp_states*v_sigma), nrow = n_ages, ncol = n_exp_states) * m_EDX2
    
    ### I
    ## I1
    m_dI1_dt <- cbind(E1 = n_exp_states*v_sigma*m_E1[, n_exp_states], n_inf_states*v_gamma*m_I1[, -n_inf_states]) - # Incoming I1
      matrix((v_r_mort + n_inf_states*v_gamma), nrow = n_ages, ncol = n_inf_states) * m_I1 -
      m_I1 %*% (m_nu_inf1_dx + m_phi_inf1_dx) - # Detected infectious 
      m_r_inf1_sx * m_I1 
    ## I2
    m_dI2_dt <- cbind(E2 = n_exp_states*v_sigma*m_E2[, n_exp_states], n_inf_states*v_gamma*m_I2[, -n_inf_states]) + # Incoming I2
      m_r_inf1_sx * m_I1 -                                                                                        # Incoming I2
      matrix((v_r_mort + n_inf_states*v_gamma), nrow = n_ages, ncol = n_inf_states) * m_I2 -
      m_I2 %*% (m_nu_inf2_dx + m_phi_inf2_dx) # Detected infectious 
    
    ### IDX
    ## IDX1
    m_dIDX1_dt <- m_I1 %*% (m_nu_inf1_dx + m_phi_inf1_dx)  +                                                    # Incoming IDX1
      cbind(EDX1 = n_exp_states*v_sigma*m_EDX1[, n_exp_states], n_inf_states*v_gamma*m_IDX1[, -n_inf_states]) - # Incoming IDX1
      matrix((v_r_mort + n_inf_states*v_gamma), nrow = n_ages, ncol = n_inf_states) * m_IDX1
    ## IDX2
    m_dIDX2_dt <- m_I2 %*% (m_nu_inf2_dx + m_phi_inf2_dx)  +                                                    # Incoming IDX2
      cbind(EDX2 = n_exp_states*v_sigma*m_EDX2[, n_exp_states], n_inf_states*v_gamma*m_IDX2[, -n_inf_states]) - # Incoming IDX2
      matrix((v_r_mort + n_inf_states*v_gamma), nrow = n_ages, ncol = n_inf_states) * m_IDX2
    
    ### R
    dR_dt  <- (1-v_cfr)*n_inf_states*v_gamma*(m_I1[, n_inf_states] + m_I2[, n_inf_states]) +          # Incoming R
      (1-(v_alpha_dx*v_cfr))*n_inf_states*v_gamma*(m_IDX1[, n_inf_states] + m_IDX2[, n_inf_states]) - # Incoming R
      v_omega * R - # Leaving R through waning immunity 
      v_r_mort * R # Leaving R through background mortality
    
    ### D
    dD_dt  <- v_r_mort * N + 
      v_cfr*n_inf_states*v_gamma*(m_I1[, n_inf_states] + m_I2[, n_inf_states]) + 
      (v_alpha_dx*v_cfr)*n_inf_states*v_gamma*(m_IDX1[, n_inf_states] + m_IDX2[, n_inf_states]) 
    
    ### Total infections
    dItot_dt <- n_exp_states*v_sigma*(m_E1[, n_exp_states] + m_E2[, n_exp_states])
    
    ### Total infectious detected
    dIDXtot_dt <- rowSums((m_I1 %*% (m_nu_inf1_dx + m_phi_inf1_dx))) + rowSums((m_I2 %*% (m_nu_inf2_dx + m_phi_inf2_dx)))
    
    ### Total exposed detected
    dEDXtot_dt <- rowSums((m_E1 %*% (m_nu_exp1_dx + m_phi_exp1_dx))) + rowSums((m_E2 %*% (m_nu_exp2_dx + m_phi_exp2_dx)))
    
    ### Total exposed and infectious detected (i.e., detected cases)
    dDXtot_dt <- dIDXtot_dt + dEDXtot_dt
    
    ### Deaths from coronavirus
    dDcov_dt <- v_cfr*n_inf_states*v_gamma*(m_I1[, n_inf_states] + m_I2[, n_inf_states]) + 
      (v_alpha_dx*v_cfr)*n_inf_states*v_gamma*(m_IDX1[, n_inf_states] + m_IDX2[, n_inf_states])
    ### Deaths from coronavirus detected
    dDcovDX_dt <- (v_alpha_dx*v_cfr)*n_inf_states*v_gamma*(m_IDX1[, n_inf_states] + m_IDX2[, n_inf_states])
    
    #### Household model derivatives ####
    ### Household FOI from community
    # v_lambda_hh <- ((t(m_Beta %*% (I/Ntot)) %*% 
    #                     v_prop_ages_curr_inf)*
    #                    (1 - p_dx)) + 
    #   (((t(m_Beta %*% (IDX/Ntot))*v_idx_scale_factor) %*% 
    #       v_prop_ages_curr_idx)*
    #      p_dx)
    
    v_lambda_I_hh   <- t(m_Beta %*% (I/Ntot))
    v_lambda_IDX_hh <- t(m_Beta %*% (IDX/Ntot))*v_idx_scale_factor
    
    v_lambda_hh     <- as.numeric((v_lambda_I_hh*(1-v_p_dx) + v_lambda_IDX_hh*v_p_dx) %*% v_prop_ages_curr_sus)
    
    # v_lambda_hh     <- as.numeric((v_lambda_I_hh + v_lambda_IDX_hh) %*% v_prop_ages_curr_sus)
    
    m_comm_trans_dx <- m_comm_trans*v_lambda_hh
    ### Household within transmission rate
    m_hh_trans_dx   <- m_hh_trans*r_tau_avg
    ### Waning rate composition
    m_hh_waning     <- r_omega * m_hh_waning_ones
    
    ### Total death rate
    r_death_tot <- sum(dD_dt)/Ntot
    ### Births
    ## Total birth rate
    r_birth_tot <- r_death_tot
    ## Bring offsprings to susceptible states
    v_index_keep_sus    <- (m_possibilities[, v_index_hh_sus] > 0)
    v_tot_births        <- rep(0, length(v_HH))
    names(v_tot_births) <- v_names_hh
    if(sum(v_HH[v_index_keep_sus]) == 0){
      v_tot_births[v_names_hh[1]] <- r_birth_tot
    } else{
      v_tot_births[v_index_keep_sus] <- (v_HH[v_index_keep_sus]/sum(v_HH[v_index_keep_sus]))*r_birth_tot
    }
    
    v_index_keep_alive  <- v_HH > 0 
    v_tot_deaths        <- rep(0, length(v_HH))
    names(v_tot_deaths) <- v_names_hh
    if(sum(v_HH[v_index_keep_alive]) == 0){
      v_tot_deaths[v_names_hh[1]] <- r_death_tot
    } else{
      v_tot_deaths[v_index_keep_alive] <- (v_HH[v_index_keep_alive]/sum(v_HH[v_index_keep_alive]))*r_death_tot
    }
    
    #### Rates of change in within household epidemics ####
    v_dHH  <- v_tot_births        +
      t(m_comm_trans_dx) %*% v_HH + 
      t(m_hh_trans_dx)   %*% v_HH + 
      t(m_hh_prog)       %*% v_HH + 
      t(m_hh_recov)      %*% v_HH +
      t(m_hh_waning)     %*% v_HH +
      - v_tot_deaths

    v_dHH <- v_dHH*sum(v_states_init[v_names_states_ages_alive])/n_hhsize
    
    # print(sum(v_dHH))
    
    ### Return list of derivatives
    return(list(c(dS_dt,  
                  as.vector(m_dE1_dt), 
                  as.vector(m_dE2_dt), 
                  as.vector(m_dEDX1_dt),
                  as.vector(m_dEDX2_dt),
                  as.vector(m_dI1_dt),
                  as.vector(m_dI2_dt),
                  as.vector(m_dIDX1_dt),
                  as.vector(m_dIDX2_dt),
                  dR_dt, 
                  dD_dt,
                  dItot_dt,
                  dIDXtot_dt,
                  dEDXtot_dt,
                  dDXtot_dt, 
                  dDcov_dt,
                  dDcovDX_dt,
                  v_dHH)))
  # }
  # )
}

#' Compiled Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO) with detected
#' infectious states
#'
#' \code{cosmo_dXdt_comp} Compiled version of the Stanford-CIDE COronavirus 
#' Simulation MOdel (SC-COSMO).
#' @export
cosmo_dXdt_comp <- compiler::cmpfun(cosmo_dXdt)
