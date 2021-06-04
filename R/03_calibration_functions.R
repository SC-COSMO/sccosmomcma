#' Generate model outputs for calibration from a parameter set
#'
#' \code{calibration_out} computes model outputs for the SC-COSMO
#' to be used for calibration routines.
#'
#' @param v_params_calib Vector of parameters that need to be calibrated.
#' @param l_params_all List with all parameters of the decision model.
#' @param n_lag_inf Lag in time series of infectious individuals.
#' @param n_lag_conf Lag in time series to account for the difference from first symptomatic.
#' @return 
#' A list with model-predicted outcomes for each target.
#' @export
calibration_out <- function(v_params_calib, 
                            l_params_all,
                            n_lag_inf = NULL,
                            n_lag_conf = 12,
                            l_dates_targets){ # User defined
  l_params_all$l_interventions[[2]]$intervention_factor <- v_params_calib["r_soc_dist_factor"]
  l_params_all$l_interventions[[3]]$intervention_factor <- v_params_calib["r_soc_dist_factor_2"]
  l_params_all$l_interventions[[4]]$intervention_factor <- v_params_calib["r_soc_dist_factor_3"]
  l_params_all$l_interventions[[5]]$intervention_factor <- v_params_calib["r_soc_dist_factor_4"]
  l_params_all$l_interventions[[6]]$intervention_factor <- v_params_calib["r_soc_dist_factor_5"]
  l_params_all$l_nu_exp2_dx[[1]]$val_start <- as.numeric(v_params_calib["r_nu_exp2_dx_lb"])     # 1 refers to period
  l_params_all$l_nu_exp2_dx[[1]]$val_end <- as.numeric(v_params_calib["r_nu_exp2_dx_ub"])
  l_params_all$l_nu_exp2_dx[[1]]$v_logit_change_rate <- as.numeric(v_params_calib["r_nu_exp2_dx_rate"])
  l_params_all$l_nu_exp2_dx[[1]]$v_logit_change_mid <- as.numeric(v_params_calib["n_nu_exp2_dx_mid"])
  l_params_all$l_nu_inf2_dx[[1]]$val_start <- as.numeric(v_params_calib["r_nu_exp2_dx_lb"])
  l_params_all$l_nu_inf2_dx[[1]]$val_end <- as.numeric(v_params_calib["r_nu_exp2_dx_ub"])
  l_params_all$l_nu_inf2_dx[[1]]$v_logit_change_rate <- as.numeric(v_params_calib["r_nu_exp2_dx_rate"])
  l_params_all$l_nu_inf2_dx[[1]]$v_logit_change_mid <- as.numeric(v_params_calib["n_nu_exp2_dx_mid"])
  
  #print(v_params_calib)
  ## Upadate parameter values to be used for calibration
  l_params_all <- sccosmomcma::update_param_list(l_params_all = l_params_all, 
                                             params_updated = v_params_calib)
  # print(l_params_all$v_soc_dist_factor)
  # Run model with updated calibrated parameters
  l_out_cosmo <- sccosmomcma::cosmo(l_params_all = l_params_all)
  
  ####### Epidemiological Output ###########################################
  v_dates  <- n_date_ini + 0:n_t_calib
  v_dates0 <- 0:sum(n_t_calib)
  # v_dates  <- df_n_t_mex_states_project_i$Date_init + 0:sum(df_n_t_mex_states_project_i$n_t)
  # v_dates0 <- 0:sum(df_n_t_mex_states_project_i$n_t)
  #### Cumulative DX Cases ####
  df_dxcum_ages <- calc_dxcum_totals(l_out_cosmo)
  
  if(is.null(n_lag_inf)){
    df_DX_cum_tot <- data.frame(Outcome = "Cumulative diagnosed infections",
                                time = df_dxcum_ages$time[-c(1:(n_lag_conf))],
                                Date = v_dates,
                                Date0 = v_dates0,
                                value = df_dxcum_ages[-c(1:(n_lag_conf)), (l_params_all$n_ages + 2)])
  } else {
    df_DX_cum_tot <- data.frame(Outcome = "Cumulative diagnosed infections", 
                                time = df_dxcum_ages$time[-c(1:(n_lag_inf + n_lag_conf))],
                                Date = v_dates,
                                Date0 = v_dates0,
                                value = df_dxcum_ages[-c(1:(n_lag_inf + n_lag_conf)), (l_params_all$n_ages + 2)])
  }
  ### Age-specific
  if(is.null(n_lag_inf)){
    df_DX_cum_ages <- data.frame(Outcome = "Cumulative diagnosed infections", 
                                 time = df_dxcum_ages$time[-c(1:(n_lag_conf))],
                                 Date = v_dates,
                                 Date0 = v_dates0,
                                 df_dxcum_ages[-c(1:(n_lag_conf)), as.character(l_params_all$v_names_ages)],
                                 check.names = FALSE)
    
  } else {
    df_DX_cum_ages <- data.frame(Outcome = "Cumulative diagnosed infections", 
                                 time = df_dxcum_ages$time[-c(1:(n_lag_inf + n_lag_conf))],
                                 Date = v_dates,
                                 Date0 = v_dates0,
                                 df_dxcum_ages[-c(1:(n_lag_inf + n_lag_conf)), as.character(l_params_all$v_names_ages)],
                                 check.names = FALSE)
  }
  
  #### Incident DX Cases ####
  df_dxinc_ages <- calc_dxinc_totals(l_out_cosmo)
  
  if(is.null(n_lag_inf)){
    df_DX_inc_tot <- data.frame(Outcome = "Incident confirmed infections",
                                time = df_dxinc_ages$time[-c(1:(n_lag_conf))],
                                Date = v_dates,
                                Date0 = v_dates0,
                                value = df_dxinc_ages[-c(1:(n_lag_conf)), (l_params_all$n_ages + 2)])
  } else {
    df_DX_inc_tot <- data.frame(Outcome = "Incident confirmed infections", 
                                time = df_dxinc_ages$time[-c(1:(n_lag_inf + n_lag_conf))],
                                Date = v_dates,
                                Date0 = v_dates0,
                                value = df_dxinc_ages[-c(1:(n_lag_inf + n_lag_conf)), (l_params_all$n_ages + 2)])
  }
  
  #### Cumulative DX COVID Deaths  ####
  df_DCov_ages <- calc_deathsdx_totals(l_out_cosmo) # calc_deaths_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_Dcov_total <- data.frame(#county =  v_states_calib, 
      Outcome = "Cumulative COVID19 deaths infections", 
      time = df_DCov_ages$time[-c(1:(n_lag_conf))],
      Date = v_dates,
      Date0 = v_dates0,
      value = df_DCov_ages[-c(1:(n_lag_conf)), (l_params_all$n_ages + 2)])
    
  } else {
    df_Dcov_total <- data.frame(#county =  v_states_calib, 
      Outcome = "Cumulative COVID19 deaths infections", 
      time = df_DCov_ages$time[-c(1:(n_lag_inf + n_lag_conf))],
      Date = v_dates,
      Date0 = v_dates0,
      value = df_DCov_ages[-c(1:(n_lag_inf + n_lag_conf)), (l_params_all$n_ages + 2)]  )
    
  }
  ### Age-specific
  if(is.null(n_lag_inf)){
    df_Dcov_cum_ages <- data.frame(Outcome = "Cumulative COVID19 deaths infections",
                                   time = df_DCov_ages$time[-c(1:(n_lag_conf))],
                                   Date = v_dates,
                                   Date0 = v_dates0,
                                   df_DCov_ages[-c(1:(n_lag_conf)), as.character(l_params_all$v_names_ages)],
                                   check.names = FALSE)
    
  } else {
    df_Dcov_cum_ages <- data.frame(Outcome = "Cumulative COVID19 deaths infections",
                                   time = df_DCov_ages$time[-c(1:(n_lag_inf + n_lag_conf))],
                                   Date = v_dates,
                                   Date0 = v_dates0,
                                   df_DCov_ages[-c(1:(n_lag_inf + n_lag_conf)), as.character(l_params_all$v_names_ages)],
                                   check.names = FALSE)
  }
  
  #### Incident COVID Deaths  ####
  df_DCov_inc_ages <- calc_incdeathsdx_totals(l_out_cosmo)
  if(is.null(n_lag_inf)){
    df_Dcov_inc_total <- data.frame(#county =  v_states_calib, 
      Outcome = "Incident COVID19 deaths infections", 
      time = df_DCov_inc_ages$time[-c(1:(n_lag_conf))],
      Date = v_dates,
      Date0 = v_dates0,
      value = df_DCov_inc_ages[-c(1:(n_lag_conf)), (l_params_all$n_ages + 2)])
    
  } else {
    df_Dcov_inc_total <- data.frame(#county =  v_states_calib, 
      Outcome = "Incident COVID19 deaths infections", 
      time = df_DCov_inc_ages$time[-c(1:(n_lag_inf + n_lag_conf))],
      Date = v_dates,
      Date0 = v_dates0,
      value = df_DCov_inc_ages[-c(1:(n_lag_inf + n_lag_conf)), (l_params_all$n_ages + 2)]  )
    
  }
  
  ### Filter by dates
  df_DX_cum_tot <- df_DX_cum_tot %>%
    filter(Date >= l_dates_targets$cases[1] & Date <= l_dates_targets$cases[2])
  
  df_DX_inc_tot <- df_DX_inc_tot %>%
    filter(Date >= l_dates_targets$cases_inc[1] & Date <= l_dates_targets$cases_inc[2])
  
  df_Dcov_total <- df_Dcov_total %>%
    filter(Date >= l_dates_targets$deaths[1] & Date <= l_dates_targets$deaths[2])
  
  df_Dcov_inc_total <- df_Dcov_inc_total %>%
    filter(Date >= l_dates_targets$deaths_inc[1] & Date <= l_dates_targets$deaths_inc[2])
  
  
  ####### Return Output ###########################################
  l_out <- list(DXCumTot   = df_DX_cum_tot,
                DXIncTot   = df_DX_inc_tot,
                DcovTot    = df_Dcov_total,
                DcovIncTot = df_Dcov_inc_total,
                DXCumAges  = df_DX_cum_ages,
                DcovAges   = df_Dcov_cum_ages)
  return(l_out)
}

#' Function to plot decision model output
#' \code{plot_model_out_vs_targets} plots the decision model outputs vs. 
#' calibration targets.
#'
#' @param l_model_out List with decision model outputs.
#' @param l_targets List with calibration targets.
#' @param print_plot Logical. Prints plot if TRUE.
#' @param return_plot Logical. Returns plot if TRUE.
#' @return 
#' A list with a ggplot object per decision model output.
#' @export
plot_model_out_vs_targets <- function(l_model_out,
                                      print_plot = TRUE, 
                                      return_plot = TRUE){
  with(as.list(l_model_out, l_targets), {
    #### Cumulative deaths
    gg.out.adeno <- ggplot(cases, aes(x = Date, y = value,
                                      ymin = lb, ymax = ub)) +
      geom_point(shape = 1, size = 2) +
      geom_errorbar() +
      geom_line(data = df.out.adeno,
                aes(x = Age, y = value), col = "blue") +
      geom_ribbon(data = df.out.adeno,
                  aes(ymin = lb, ymax = ub), alpha = 0.4, fill = "blue") +
      facet_wrap(~Target) +
      scale_x_continuous(limits = c(NA, 90)) +
      scale_y_continuous(breaks = number_ticks(6)) +
      ylab("Proportion per 100") +
      theme_bw(base_size = 16)
    # #### CRC Incidence
    # gg.out.crc <- ggplot(df.targets_crc, aes(x = Age, y = value, 
    #                                          ymin = lb, ymax = ub)) +
    #   geom_point(shape = 1, size = 2) +
    #   geom_errorbar() +
    #   geom_line(data = df.out.crc, 
    #             aes(x = Age, y = value), col = "red") +
    #   geom_ribbon(data = df.out.crc,
    #               aes(ymin = lb, ymax = ub), alpha = 0.4, fill = "red") +
    #   facet_wrap(~Target) +
    #   scale_x_continuous(limits = c(NA, 85)) +
    #   scale_y_continuous(breaks = number_ticks(6)) +
    #   ylab("Incidence per 100,000") +
    #   theme_bw(base_size = 16)
    # #### CRC Mortality
    # gg.out.crc.mort <- ggplot(df.target_crc_mort, aes(x = Age, y = value, 
    #                                                   ymin = lb, ymax = ub)) +
    #   geom_point(shape = 1, size = 2) +
    #   geom_errorbar() +
    #   geom_line(data = df.out.crc.mort, 
    #             aes(x = Age, y = value), col = "red") +
    #   geom_ribbon(data = df.out.crc.mort,
    #               aes(ymin = lb, ymax = ub), alpha = 0.4, fill = "red") +
    #   facet_wrap(~Target) +
    #   scale_x_continuous(limits = c(NA, 85)) +
    #   scale_y_continuous(breaks = number_ticks(6)) +
    #   ylab("Mortality rate per 100,000") +
    #   theme_bw(base_size = 16)
    # 
    # if(print.plot == T){
    #   print(gg.out.adeno)
    #   print(gg.out.crc)
    #   print(gg.out.crc.mort)
    # }
    # if(return.plot){
    #   return(list(gg.out.adeno    = gg.out.adeno,
    #               gg.out.crc      = gg.out.crc,
    #               gg.out.crc.mort = gg.out.crc.mort))  
    # }
  }
  )
  
}

#' Log-likelihood function for a parameter set
#'
#' \code{log_lik} computes a log-likelihood value for one (or multiple) 
#' parameter set(s) to be used in optimization procedures
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @param l_params_all List with all parameters of the decision model. 
#' @param comp Logil. Should the compiled versio of SC-COSMO be used?
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @export
log_lik_opt <- function(v_params,
                        ...){ # User defined
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  v_target_names <- c("DXIncTot"#, "DcovIncTot" #, "DXCumAges", "DcovCumAges"
  )
  n_target       <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=2
    jj <- tryCatch( { 
      ###   Run model for parametr set "v_params" ###
      l_model_res <- calibration_out(v_params_calib = v_params[j, ],...)
      l_model_res$DXIncTot <- l_model_res$DXIncTot %>% 
        mutate(value = replace(value, value == 0, 1e-20))
      # l_model_res$DcovIncTot <- l_model_res$DcovIncTot %>% 
      #   mutate(value = replace(value, value == 0, 1e-20))
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      ## TARGET 1: Incident DX Cases ("DXCumTot")
      ## Normal log-likelihood  
      # v_llik[j, "DXIncTot"] <- sum(dnorm(x = l_targets$cases_inc$value,# 
      #                                    mean = l_model_res$DXIncTot$value, 
      #                                    # sd = l_targets$cases$se,
      #                                    log = T), na.rm = T)
      ## Poisson log-likelihood
      # v_llik[j, "DXIncTot"] <- sum(dpois(x = l_targets$cases_inc$value,# 
      #                                    lambda = l_model_res$DXIncTot$value, 
      #                                    log = T), na.rm = T)
      # print("Poisson")
      ## Negative-Binomial log-likelihood
      v_llik[j, "DXIncTot"] <- sum(dnbinom(x = l_targets$cases_inc$value, 
                                           size = 1,
                                           mu = l_model_res$DXIncTot$value,
                                           # sd = l_targets$cases$se,
                                           log = T), na.rm = T)
      # print("Negative-Binomial")
      # ## Multivariate-Normal log-likelihood  
      # v_llik[j, "DXIncTot"] <- mvnfast::dmvn(X = l_targets$cases_inc$value,
      #                                        mu = l_model_res$DXIncTot$value,
      #                                        sigma = l_targets$cov_cases_inc,
      #                                        log = TRUE)
      # print("MVN")
      
      ## TARGET 2: DX COVID Deaths ("DcovTot")
      ## Poisson log-likelihood
      # v_llik[j, "DcovIncTot"] <- sum(dpois(x = l_targets$deaths_inc$value,
      #                                   lambda = l_model_res$DcovIncTot$value, 
      #                                   log = T), na.rm = T)
      # Negative-Binomial log-likelihood
      # v_llik[j, "DcovIncTot"] <- sum(dnbinom(x = l_targets$deaths_inc$value,
      #                                        size = 1,
      #                                        mu= l_model_res$DcovIncTot$value,
      #                                        log = T), na.rm = T)
      # ## Multivariate-Normal log-likelihood  
      # v_llik[j, "DcovIncTot"] <- mvnfast::dmvn(X = l_targets$deaths_inc$value,
      #                                        mu = l_model_res$DcovIncTot$value,
      #                                        sigma = l_targets$cov_deaths_inc,
      #                                        log = TRUE)
      # print("MVN")
      
      # ## TARGET 3: Age-specific Cumulative DX Cases ("DXCumAges")
      # n_calib_t <- length(l_model_res$DXCumTot$value)
      # v_dates_calib_ages <- c(n_calib_t) # c(floor(n_calib_t/2), n_calib_t)
      # # Multinomial log-likelihood
      # v_llik[j, "DXCumAges"] <- sum(mc2d::dmultinomial(x = as.matrix(l_targets$cases_age[v_dates_calib_ages, -c(1:4)]), 
      #                                                  prob =  as.matrix(l_model_res$DXCumAges[v_dates_calib_ages, -c(1:4)])/l_model_res$DXCumTot$value[v_dates_calib_ages],
      #                                                  log = T), na.rm = T)
      # 
      # ## TARGET 4: Age-specific DX COVID Deaths ("DcovAges")
      # # Multinomial log-likelihood
      # v_llik[j, "DcovCumAges"] <- sum(mc2d::dmultinomial(x = as.matrix(l_targets$deaths_age[v_dates_calib_ages, -c(1:4)]), 
      #                                                  prob =  as.matrix(l_model_res$DcovAges[v_dates_calib_ages, -c(1:4)])/l_model_res$DcovTot$value[v_dates_calib_ages],
      #                                                  log = T), na.rm = T)
      ## OVERALL
      ## can give different targets different weights (user must change this)
      v_weights <- rep(1, n_target)
      ## weighted sum
      v_llik_overall[j] <- v_llik[j, ] %*% v_weights
    }, error = function(e) NA) 
    if(is.na(jj)) { v_llik_overall <- -Inf }
  } ## End loop over sampled parameter sets
  
  ## return GOF
  return(v_llik_overall)
}

#' Log-likelihood function for a parameter set
#'
#' \code{log_lik} computes a log-likelihood value for one (or multiple) 
#' parameter set(s).
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @param l_params_all List with all parameters of the decision model. 
#' @param comp Logil. Should the compiled versio of SC-COSMO be used?
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @export
log_lik <- function(v_params,
                    ...){ # User defined
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  v_target_names <- c("DXIncTot" #, "DcovIncTot" #, "DXCumAges", "DcovCumAges"
  )
  n_target       <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=2
    jj <- tryCatch( { 
      if(sum(v_params[j, ] < get_bounds()$v_lb | v_params[j, ] > get_bounds()$v_ub) > 0){
        v_llik_overall[j] <- NA
      } else{
        ###   Run model for parametr set "v_params" ###
        l_model_res <- calibration_out(v_params_calib = v_params[j, ],
                                       #l_params_all = l_params_all,  n_lag_inf  = n_lag_inf, n_lag_conf = n_lag_conf, l_dates_targets = l_dates_targets)
                                       ...)
        
        l_model_res$DXIncTot <- l_model_res$DXIncTot %>% 
          mutate(value = replace(value, value == 0, 1e-20))
        # l_model_res$DcovIncTot <- l_model_res$DcovIncTot %>% 
        #   mutate(value = replace(value, value == 0, 1e-20))
        
        ###  Calculate log-likelihood of model outputs to targets  ###
        ## TARGET 1: Cumulative DX Cases ("DXCumTot")
        ## Normal log-likelihood  
        # v_llik[j, "DXIncTot"] <- sum(dnorm(x = l_targets$cases_inc$value,# 
        #                                    mean = l_model_res$DXIncTot$value, 
        #                                    # sd = l_targets$cases$se,
        #                                    log = T), na.rm = T)
        ## Poisson log-likelihood
        # v_llik[j, "DXIncTot"] <- sum(dpois(x = l_targets$cases_inc$value,# 
        #                                    lambda = l_model_res$DXIncTot$value, 
        #                                    # sd = l_targets$cases$se,
        #                                    log = T), na.rm = T)
        # print("POisson")
        ## Negative-Binomial log-likelihood
        v_llik[j, "DXIncTot"] <- sum(dnbinom(x = l_targets$cases_inc$value, 
                                             size = 1,
                                             mu = l_model_res$DXIncTot$value,
                                             log = T), na.rm = T)
        print(v_llik[j, "DXIncTot"])
        # print("Negative-Binomial")
        ## Multivariate-Normal log-likelihood  
        # v_llik[j, "DXIncTot"] <- mvnfast::dmvn(X = l_targets$cases_inc$value,
        #                                        mu = l_model_res$DXIncTot$value,
        #                                        sigma = l_targets$cov_cases_inc,
        #                                        log = TRUE)
        # print("MVN")
        
        ## TARGET 2: DX COVID Deaths ("DcovTot")
        ## Poisson log-likelihood
        # v_llik[j, "DcovIncTot"] <- sum(dpois(x = l_targets$deaths_inc$value,
        #                                      lambda = l_model_res$DcovIncTot$value, 
        #                                      # sd = l_targets$deaths$se,
        #                                      log = T), na.rm = T)
        # Negative-Binomial log-likelihood
        # v_llik[j, "DcovIncTot"] <- sum(dnbinom(x = l_targets$deaths_inc$value,
        #                                        size = 1,
        #                                        mu= l_model_res$DcovIncTot$value,
        #                                        log = T), na.rm = T)
        ## Multivariate-Normal log-likelihood  
        # v_llik[j, "DcovIncTot"] <- mvnfast::dmvn(X = l_targets$deaths_inc$value,
        #                                          mu = l_model_res$DcovIncTot$value,
        #                                          sigma = l_targets$cov_deaths_inc,
        #                                          log = TRUE)
        # print("MVN")
        
        # ## TARGET 3: Age-specific Cumulative DX Cases ("DXCumAges")
        # n_calib_t <- length(l_model_res$DXCumTot$value)
        # v_dates_calib_ages <- c(n_calib_t) # c(floor(n_calib_t/2), n_calib_t)
        # # Multinomial log-likelihood
        # v_llik[j, "DXCumAges"] <- sum(mc2d::dmultinomial(x = as.matrix(l_targets$cases_age[v_dates_calib_ages, -c(1:4)]), 
        #                                                  prob =  as.matrix(l_model_res$DXCumAges[v_dates_calib_ages, -c(1:4)])/l_model_res$DXCumTot$value[v_dates_calib_ages],
        #                                                  log = T), na.rm = T)
        # 
        # ## TARGET 4: Age-specific DX COVID Deaths ("DcovAges")
        # # Multinomial log-likelihood
        # v_llik[j, "DcovCumAges"] <- sum(mc2d::dmultinomial(x = as.matrix(l_targets$deaths_age[v_dates_calib_ages, -c(1:4)]), 
        #                                                    prob =  as.matrix(l_model_res$DcovAges[v_dates_calib_ages, -c(1:4)])/l_model_res$DcovTot$value[v_dates_calib_ages],
        #                                                    log = T), na.rm = T)
        
        ## OVERALL
        ## can give different targets different weights (user must change this)
        v_weights <- rep(1, n_target)
        ## weighted sum
        v_llik_overall[j] <- v_llik[j, ] %*% v_weights
      }
    }, error = function(e) NA) 
    if(is.na(jj)) { v_llik_overall[j] <- -Inf }
  } ## End loop over sampled parameter sets
  
  ## return GOF
  return(v_llik_overall)
}

#' Parallel evaluation of log-likelihood function for a sets of parameters
#'
#' \code{log_lik_par} computes a log-likelihood value for one (or multiple) 
#' parameter set(s) using parallel computation.
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @param l_params_all List with all parameters of the decision model. 
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @export
log_lik_par <- function(v_params, 
                        log_lik_offset,
                        ...) { 
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  
  n_samp <- nrow(v_params)
  
  ### Get OS
  os <- get_os()
  print(paste0("Parallelized Likelihood calculations on ", os))
  
  no_cores <- detectCores() - 1
  
  n_time_init_likpar <- Sys.time()
  
  if(os == "macosx"){
    # Initialize cluster object
    cl <- parallel::makeForkCluster(no_cores) 
    doParallel::registerDoParallel(cl)
    llk <- foreach(i = 1:n_samp, .combine = c) %dopar% {
      log_lik(v_params[i, ], ...)
    }
    n_time_end_likpar <- Sys.time()
  }
  if(os == "windows"){
    # Initialize cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
    opts <- list(attachExportEnv = TRUE)
    llk <- foreach(i = 1:n_samp, .combine = c,
                   .export = ls(globalenv()),
                   .packages=c("sccosmomcma",
                               "ggplot2",
                               # "tictoc",
                               "tidyverse",
                               "dplyr",
                               # "lubridate",
                               "epitools"
                               # "dampack", 
                               # "sccosmoData"
                   ),
                   .options.snow = opts) %dopar% {
                     log_lik(v_params[i, ], ...)
                   }
    n_time_end_likpar <- Sys.time()
  }
  if(os == "linux"){
    # Initialize cluster object
    cl <- makeCluster(no_cores)
    registerDoMC(cl)
    llk <- foreach(i = 1:n_samp, .combine = c) %dopar% {
      log_lik(v_params[i, ], ...)
    }
    n_time_end_likpar <- Sys.time()
  }
  
  stopCluster(cl)
  n_time_total_likpar <- difftime(n_time_end_likpar, n_time_init_likpar, 
                                  units = "hours")
  print(paste0("Runtime: ", round(n_time_total_likpar, 2), " hrs."))
  #-# Try this: # PO
  rm(cl)        # PO
  gc()          # PO
  #-#           # PO
  return(llk - (log_lik_offset))
}

#' Likelihood
#'
#' \code{likelihood} computes a likelihood value for one (or multiple) 
#' parameter set(s).
#'
#' @param v_params Vector (or matrix) of model parameters. 
#' @return 
#' A scalar (or vector) with likelihood values.
#' @export
likelihood <- function(v_params){ 
  v_like <- exp(log_lik_par(v_params, 
                            log_lik_offset = -2162.375,
                            l_params_all,
                            n_lag_inf, 
                            n_lag_conf, 
                            l_dates_targets)) # 504 + 
  return(v_like)
}

#' Evaluate log-posterior of calibrated parameters
#'
#' \code{log_post} Computes a log-posterior value for one (or multiple) 
#' parameter set(s) based on the simulation model, likelihood functions and 
#' prior distributions.
#' @param v_params Vector (or matrix) of model parameters 
#' @return 
#' A scalar (or vector) with log-posterior values.
#' @export
log_post <- function(v_params, ...) { 
  v_lpost <- log_prior(v_params) + log_lik(v_params, ...)
  return(v_lpost) 
}

#' Evaluate log-posterior of calibrated parameters for OPTIMIZATION purposes
#'
#' \code{log_post_opt} Computes a log-posterior value for one (or multiple) 
#' parameter set(s) based on the simulation model, likelihood functions and 
#' prior distributions and it's used for OPTIMIZATION purposes, NOT Bayesian
#' estimation.
#' @param v_params Vector (or matrix) of model parameters 
#' @return 
#' A scalar (or vector) with log-posterior values.
#' @export
log_post_opt <- function(v_params, ...) { 
  
  v_param_names = names(get_bounds()$v_lb) #make sure this is actually same name as in model params
  v_lb = get_bounds()$v_lb # get rid of defaults and then user can pass in
  v_ub = get_bounds()$v_ub
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_param <- length(v_param_names)
  n_samp <- nrow(v_params)
  #  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  oob <- FALSE
  for (i in 1:n_param){
    if (v_params[,i]<v_lb[i] | v_params[, i] > v_ub[i]) {
      lprior <- lprior + 
        -1000000000*ifelse(v_params[,i]<v_lb[i], abs(v_lb[i]-v_params[,i]), abs(v_ub[i]-v_params[,i])) +
        -10000000
      oob <- TRUE
    }
  }
  if (oob==TRUE) {
    return(lprior)
  }
  v_lpost <- log_prior_opt(v_params = v_params,
                           v_param_names = v_param_names, #make sure this is actually same name as in model params
                           v_lb = v_lb, # get rid of defaults and then user can pass in
                           v_ub = v_ub
  ) +
    log_lik_opt(v_params, ...)
  return(v_lpost) 
}

#' Evaluate posterior of calibrated parameters
#'
#' \code{posterior} computes a posterior value for one (or multiple) parameter 
#' set(s).
#' @param v_params Vector (or matrix) of model parameters 
#' @return 
#' A scalar (or vector) with posterior values.
#' @export
posterior <- function(v_params, ...) { 
  v_posterior <- exp(log_post(v_params)) 
  return(v_posterior)
}

#' Get bounds for parameters to calibrate (and their standard errors)
#'
#' \code{get_bounds} defines the bounds for each of the calibrated parameters.
#'
#' @return 
#' A list with lower and upper bounds and their standard errors
#' @export
get_bounds <- function() {
  
  # lower bounds
  v_lb <- c(r_beta              = 0.100,
            r_tau               = 0.150,
            r_soc_dist_factor   = 0.250,
            r_soc_dist_factor_2 = 0.250,
            r_soc_dist_factor_3 = 0.250,
            r_soc_dist_factor_4 = 0.250,
            r_soc_dist_factor_5 = 0.250,
            r_nu_exp2_dx_lb     = 0.005,
            r_nu_exp2_dx_ub     = 0.005,
            r_nu_exp2_dx_rate   = 0.010,
            n_nu_exp2_dx_mid    = 30
            
  )
  
  ## upper bounds
  v_ub <- c(r_beta              = 0.300,
            r_tau               = 0.400,
            r_soc_dist_factor   = 0.750,
            r_soc_dist_factor_2 = 0.750,
            r_soc_dist_factor_3 = 0.750,
            r_soc_dist_factor_4 = 0.750,
            r_soc_dist_factor_5 = 0.750,
            r_nu_exp2_dx_lb     = 0.120,
            r_nu_exp2_dx_ub     = 0.250,
            r_nu_exp2_dx_rate   = 1.000,
            n_nu_exp2_dx_mid    = 100 # 60
            
  )
  
  ## standard errors based on bounds
  v_se <- (v_ub - v_lb)/(2*2) # even more conservative an 1.96
  
  return(list(v_lb = v_lb, 
              v_ub = v_ub,
              v_se = v_se))
}

#' Sample from prior distributions of calibrated parameters
#'
#' \code{sample.prior} generates a sample of parameter sets from their prior 
#' distribution.
#' @param n_samp Number of samples.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A matrix with 3 rows and \code{n_samp} rows. Each row corresponds to a 
#' parameter set sampled from their prior distributions
#' @export
sample.prior <- function(n_samp,
                         v_param_names = names(get_bounds()$v_lb), #make sure this is actually same name as in model params
                         v_lb = get_bounds()$v_lb, # get rid of defaults and then user can pass in 
                         v_ub = get_bounds()$v_ub){
  n_param <- length(v_param_names)
  m_lhs_unit   <- lhs::randomLHS(n = n_samp, k = n_param)
  m_param_samp <- matrix(nrow = n_samp, ncol = n_param)
  colnames(m_param_samp) <- v_param_names
  for (i in 1:n_param){
    m_param_samp[, i] <- qunif(m_lhs_unit[,i],
                               min = v_lb[i],
                               max = v_ub[i])
    # ALTERNATIVE prior using beta (or other) distributions
    # m_param_samp[, i] <- qbeta(m_lhs_unit[,i],
    #                            min = 1,
    #                            max = 1)
  }
  return(m_param_samp)
}

#' Evaluate log-prior of calibrated parameters
#'
#' \code{log_prior} computes a log-prior value for one (or multiple) parameter 
#' set(s) based on their prior distributions
#' @param v_params Vector (or matrix) of model parameters.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A scalar (or vector) with log-prior values.
#' @export
log_prior <- function(v_params, 
                      v_param_names = names(get_bounds()$v_lb), #make sure this is actually same name as in model params
                      v_lb = get_bounds()$v_lb, # get rid of defaults and then user can pass in 
                      v_ub = get_bounds()$v_ub){
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_param <- length(v_param_names)
  n_samp <- nrow(v_params)
  #  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  for (i in 1:n_param){
    lprior <- lprior + dunif(v_params[, i],
                             min = v_lb[i],
                             max = v_ub[i], 
                             log = T)
    # ALTERNATIVE prior using beta distributions
    # lprior <- lprior + dbeta(v_params[, i],
    #                          min = 1,
    #                          max = 1, 
    #                          log = T)
  }
  return(lprior)
}

#' Evaluate log-prior of calibrated parameters for OPTIMIZATION purposes
#'
#' \code{log_prior} computes a log-prior value for one (or multiple) parameter 
#' set(s) based on their prior distributions and it's used for OPTIMIZATION 
#' purposes, NOT Bayesian estimation.
#' @param v_params Vector (or matrix) of model parameters.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A scalar (or vector) with log-prior values.
#' @export
log_prior_opt <- function(v_params, 
                          v_param_names = names(get_bounds()$v_lb), #make sure this is actually same name as in model params
                          v_lb = get_bounds()$v_lb, # get rid of defaults and then user can pass in 
                          v_ub = get_bounds()$v_ub){
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params) 
  }
  n_param <- length(v_param_names)
  n_samp <- nrow(v_params)
  #  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  for (i in 1:n_param){
    if ((v_params[, i] < v_lb[i]) | (v_params[, i] > v_ub[i])) {
      lprior <- lprior + -1000000000*ifelse(v_params[, i] < v_lb[i], 
                                            abs(v_lb[i] - v_params[,i]), 
                                            abs(v_ub[i] - v_params[,i]))
    } else{
      lprior <- lprior + dunif(v_params[, i],
                               min = v_lb[i],
                               max = v_ub[i], 
                               log = T)
      # ALTERNATIVE prior using beta distributions
      # lprior <- lprior + dbeta(v_params[, i],
      #                          min = 1,
      #                          max = 1, 
      #                          log = T)  
    }
    
  }
  return(lprior)
}

#' Evaluate prior of calibrated parameters
#'
#' \code{prior} computes a prior value for one (or multiple) parameter set(s).
#' @param v_params Vector (or matrix) of model parameters 
#' @return 
#' A scalar (or vector) with prior values.
#' @export
prior <- function(v_params) { 
  v_prior <- exp(log_prior(v_params)) 
  return(v_prior)
}

#' Get operating system
#' 
#' @return 
#' A string with the operating system.
#' @export
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "MacOSX"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

#' Helper convenience function that defines age-specific lists of periods for values that 
#' change over the duration of the simulation.
#' 
#' \code{get_non_const_multiage_list} Adds lists of time periods as new entries 
#' in list of age group lists
#' @param v_time_stop vector with time period at which values change (in days).
#' Last entry should be n_t.
#' @param m_ageval matrix of values for each age group (n_ages) at each time period
#' @return 
#' List (l_l_period) list of lists of periods (1 per age group) with 
#' non-constant values
#' @export
get_non_const_multiage_list <- function(v_time_stop, m_ageval) {
  n_ages     <- nrow(m_ageval)
  n_periods  <- length(v_time_stop)
  l_l_period <- NULL
  
  for(curr_age in 1:n_ages) { # curr_age <- 1
    curr_l_period <- NULL
    time_start    <- 0
    for(period in 1:n_periods){ # period <- 1
      time_stop <- v_time_stop[period] 
      curr_period <- make_period(functional_form = "constant", 
                                 time_start      = time_start, 
                                 time_stop       = time_stop,
                                 val_start       = m_ageval[curr_age, period],
                                 val_end         = m_ageval[curr_age, period]
      )
      if(period==1){
        curr_l_period <- add_period(NULL, curr_period) 
      } else{
        curr_l_period <- add_period(curr_l_period, curr_period)
      }
      time_start <- time_stop
    }
    l_l_period <- add_list_periods(l_l_period, curr_l_period)
  }
  return(l_l_period)
} 

