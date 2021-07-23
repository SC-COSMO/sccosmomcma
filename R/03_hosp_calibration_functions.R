#' Generate model outputs for calibration from a hospitalization parameter set
#'
#' \code{hosp_calibration_out} computes model hospitalization outputs
#' to be used for calibration routines.
#'
#' @param v_params_calib Vector of hospitalization parameters that need to be 
#' calibrated.
#' @param l_params_all List with all parameters of the decision model.
#' @param n_lag_inf Lag in time series of infectious individuals.
#' @param n_lag_conf Lag in time series to account for the difference from 
#' first symptomatic.
#' @param l_dates_hosp_targets List of initial and final date of target series.
#' @return 
#' A list with model-predicted outcomes for each target.
#' @export
hosp_calibration_out <- function(v_params_calib, 
                                 l_params_all,
                                 n_lag_inf = NULL,
                                 n_lag_conf = 12,
                                 l_dates_hosp_targets){ # User defined

  # Parameters to be used for calibration
  l_params_all$m_r_exit_tot    <- v_params_calib["m_r_exit_tot"]
  l_params_all$m_r_exit_nonicu <- v_params_calib["m_r_exit_nonicu"]
  l_params_all$m_r_exit_icu    <- v_params_calib["m_r_exit_icu"]
  l_params_all$m_sigma_tot     <- v_params_calib["m_sigma_tot"]
  l_params_all$m_sigma_nonicu  <- v_params_calib["m_sigma_nonicu"]
  l_params_all$m_sigma_icu     <- v_params_calib["m_sigma_icu"]
  
  ## Update parameter values to be used for calibration
  l_params_all <- sccosmomcma::update_param_list(l_params_all   = l_params_all, 
                                                 params_updated = v_params_calib)
  
  
  # Epidemiological Output --------------------------------------------------
  
  v_dates  <- n_date_ini + 0:n_t_calib
  v_dates0 <- 0:sum(n_t_calib)
  v_time   <- n_lag_inf:(n_t_calib + n_lag_inf)
  
  ## Hospitalized incidence and prevalence ----------------------------------

  l_hosp_out <- prep_dx_hosp_calib(l_params_all)
  
  df_hosp_out_inc  <- l_hosp_out$m_hosp_inc
  df_hosp_out_prev <- l_hosp_out$m_hosp_prev


  ### HOSPITALIZED PREVALENCE -----------------------------------------------
  
  df_hosp_prev <- data.frame(Outcome = "Hospitalized prevalence",
                             time    = v_time,
                             Date    = v_dates,
                             Date0   = v_dates0,
                             value   = df_hosp_out_prev[1 , c(l_params_all$n_ages+1), ]
  )
  
  # Age specific
  df_hosp_prev_ages <- data.frame(Outcome = "Hospitalized prevalence",
                                  time    = v_time,
                                  Date    = v_dates,
                                  Date0   = v_dates0,
                                  t(df_hosp_out_prev[1 , -c(l_params_all$n_ages+1), ]),
                                  check.names = FALSE
  )


  ## HOSPITALIZED NON ICU PREVALENCE ----------------------------------------

  df_hosp_novent_prev <- data.frame(Outcome = "Hospitalized non ICU prevalence",
                                    time    = v_time,
                                    Date    = v_dates,
                                    Date0   = v_dates0,
                                    value   = df_hosp_out_prev[2, c(l_params_all$n_ages+1), ])
  
  # Age-specific
  df_hosp_novent_prev_ages <- data.frame(Outcome = "Hospitalized non ICU prevalence",
                                         time    = v_time,
                                         Date    = v_dates,
                                         Date0   = v_dates0,
                                         t(df_hosp_out_prev[2, -c(l_params_all$n_ages+1), ]),
                                         check.names = FALSE
    )
  

  ## HOSPITALIZED ICU PREVALENCE --------------------------------------------
  
  df_hosp_vent_prev <- data.frame(Outcome = "Hospitalized ICU prevalence",
                                  time    = v_time,
                                  Date    = v_dates,
                                  Date0   = v_dates0,
                                  value   = df_hosp_out_prev[3, c(l_params_all$n_ages+1), ])
  
  # Age-specific
  df_hosp_vent_prev_ages <- data.frame(Outcome = "Hospitalized ICU prevalence",
                                       time    = v_time,
                                       Date    = v_dates,
                                       Date0   = v_dates0,
                                       t(df_hosp_out_prev[3, -c(l_params_all$n_ages+1), ]),
                                       check.names = FALSE)
  

  ## HOSPITALIZED INCIDENCE -------------------------------------------------

  df_hosp_inc <- data.frame(Outcome = "Hospitalized incidence",
                            time    = v_time,
                            Date    = v_dates,
                            Date0   = v_dates0,
                            value   = df_hosp_out_inc[1, c(l_params_all$n_ages+1), ])
  
  # Age-specific
  df_hosp_inc_ages <- data.frame(Outcome = "Hospitalized incidence",
                                 time    = v_time,
                                 Date    = v_dates,
                                 Date0   = v_dates0,
                                 t(df_hosp_out_inc[1, -c(l_params_all$n_ages+1), ]),
                                 check.names = FALSE)


  ## HOSPITALIZED NON ICU INCIDENCE -----------------------------------------

  df_hosp_novent_inc <- data.frame(Outcome = "Hospitalized non ICU incidence",
                                     time    = v_time,
                                     Date    = v_dates,
                                     Date0   = v_dates0,
                                     value   = df_hosp_out_inc[2, c(l_params_all$n_ages+1), ])

  # Age-specific
  df_hosp_novent_inc_ages <- data.frame(Outcome = "Hospitalized non ICU incidence",
                                        time    = v_time,
                                        Date    = v_dates,
                                        Date0   = v_dates0,
                                        t(df_hosp_out_inc[2, -c(l_params_all$n_ages+1), ]),
                                        check.names = FALSE)
  

  ## HOSPITALIZED ICU INCIDENCE ---------------------------------------------

  df_hosp_vent_inc <- data.frame(Outcome = "Hospitalized ICU incidence",
                                 time    = v_time,
                                 Date    = v_dates,
                                 Date0   = v_dates0,
                                 value   = df_hosp_out_inc[3, c(l_params_all$n_ages+1), ])
  

  # Age-specific
  df_hosp_vent_inc_ages <- data.frame(Outcome = "Hospitalized ICU incidence",
                                      time    = v_time,
                                      Date    = v_dates,
                                      Date0   = v_dates0,
                                      t(df_hosp_out_inc[3, -c(l_params_all$n_ages+1), ]),
                                      check.names = FALSE)


  ### Filter by dates -------------------------------------------------------

  df_H_prev <- df_hosp_prev %>%
    filter(Date >= l_dates_hosp_targets$hosp[1] & Date <= l_dates_hosp_targets$hosp[2])
  
  df_H_inc <- df_hosp_inc %>%
    filter(Date >= l_dates_hosp_targets$hosp[1] & Date <= l_dates_hosp_targets$hosp[2])

  df_HNovent_prev <- df_hosp_novent_prev %>%
    filter(Date >= l_dates_hosp_targets$novent[1] & Date <= l_dates_hosp_targets$novent[2])
  
  df_HNovent_inc <- df_hosp_novent_inc%>%
    filter(Date >= l_dates_hosp_targets$novent[1] & Date <= l_dates_hosp_targets$novent[2])
  
  df_Hvent_prev <- df_hosp_vent_prev %>%
    filter(Date >= l_dates_hosp_targets$vent[1] & Date <= l_dates_hosp_targets$vent[2])
  
  df_Hvent_inc <- df_hosp_vent_inc%>%
    filter(Date >= l_dates_hosp_targets$vent[1] & Date <= l_dates_hosp_targets$vent[2])

  # Return Output -----------------------------------------------------------
  
  l_out <- list(H_prev        = df_H_prev,
                H_vent_prev   = df_Hvent_prev,
                H_novent_prev = df_HNovent_prev,
                H_inc         = df_H_inc,
                H_vent_inc    = df_Hvent_inc,
                H_novent_inc  = df_HNovent_inc)
  
  return(l_out)
}

#' Log-likelihood function for a hospitalization parameter set 
#'
#' \code{hosp_log_lik_opt} computes a log-likelihood value for one (or multiple) 
#' hospitalization parameter set(s) to be used in optimization procedures.
#'
#' @param v_params Vector (or matrix) of hospitalization parameters.
#' @param ... Further arguments to be passed to.
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @export
hosp_log_lik_opt <- function(v_params,
                             ...){ # User defined
  
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  
  n_samp <- nrow(v_params)
  v_target_names <- c("H_prev", "H_vent_prev", "H_novent_prev"#, "H_inc", "H_vent_inc"
  )
  n_target       <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j = 1
    jj <- tryCatch( { 
      ###   Run model for parameter set j ###
      # m_params_samp <- c()
      # for(i in 1:length(l_params)){
      #   m_params_samp <- rbind(m_params_samp, as.vector(l_params[[i]][j,]))
      # }
      # row.names(m_params_samp) <- names(l_params)
      
      l_model_res <- hosp_calibration_out(v_params_calib = v_params[j,], 
                                          # l_params_all,
                                          # n_lag_inf = 14,
                                          # n_lag_conf = 0,
                                          # l_dates_hosp_targets)
                                          ...)
      l_model_res$H_prev <- l_model_res$H_prev %>% 
        mutate(value = replace(value, value == 0, 1e-20))
      
      l_model_res$H_vent_prev <- l_model_res$H_vent_prev %>% 
        mutate(value = replace(value, value == 0, 1e-20))
      
      l_model_res$H_novent_prev <- l_model_res$H_novent_prev %>% 
        mutate(value = replace(value, value == 0, 1e-20))

      ###  Calculate log-likelihood of model outputs to targets  ###
      #### TARGET 1: Hospitalized prevalence ("H_prev") ####
      ## Normal log-likelihood  
      v_llik[j, "H_prev"] <- sum(dnorm(x   = l_hosp_targets$hosp$value,#
                                      mean = l_model_res$H_prev$value,
                                      sd   = l_hosp_targets$hosp$se,
                                      log  = T),
                                na.rm = T)
      # print("Log-normal")
      
      ## Poisson log-likelihood
      # v_llik[j, "H_prev"] <- sum(dpois(x     = l_hosp_targets$hosp$value,#
      #                                 lambda = l_model_res$H_prev$value,
      #                                 log    = T), 
      #                           na.rm = T)
      # print("Poisson")
      
      ## Negative-Binomial log-likelihood
      # v_llik[j, "H_prev"] <- sum(dnbinom(x   = l_hosp_targets$hosp$value, 
      #                                   size = 1,
      #                                   mu   = l_model_res$H_prev$value,
      #                                   # sd = l_hosp_targets$hosp$se,
      #                                   log  = T), 
      #                           na.rm = T)
      
      # print("Negative-Binomial")
      
      # ## Multivariate-Normal log-likelihood  
      # v_llik[j, "H_prev"] <- mvnfast::dmvn(X    = l_hosp_targets$hosp$value,
      #                                     mu    = l_model_res$H_prev$value,
      #                                     sigma = l_hosp_targets$cov_hosp,
      #                                     log   = TRUE)
      # print("MVN")
      
      #### TARGET 2: Hospitalized prevalence with ventilator ("H_vent_prev") ####
      ## Normal log-likelihood  
      v_llik[j, "H_vent_prev"] <- sum(dnorm(x    = l_hosp_targets$vent$value,#
                                            mean = l_model_res$H_vent_prev$value,
                                            sd   = l_hosp_targets$vent$se,
                                            log  = T),
                                      na.rm = T)
      # print("Log-normal")
      
      # Poisson log-likelihood
      # v_llik[j, "H_vent_prev"] <- sum(dpois(x     = l_hosp_targets$vent$value,
      #                                      lambda = l_model_res$H_vent_prev$value,
      #                                      log    = T), 
      #                                na.rm = T)
      
      # Negative-Binomial log-likelihood
      # v_llik[j, "H_vent_prev"] <- sum(dnbinom(x   = l_hosp_targets$vent$value,
      #                                        size = 1,
      #                                        mu   = l_model_res$H_vent_prev$value,
      #                                        log  = T),
      #                                na.rm = T)
      
      # Multivariate-Normal log-likelihood  
      # v_llik[j, "H_vent_prev"] <- mvnfast::dmvn(X    = l_hosp_targets$vent$value,
      #                                          mu    = l_model_res$H_vent_prev$value,
      #                                          sigma = l_hosp_targets$cov_ven,
      #                                          log   = TRUE)
      # print("MVN")
      
      #### TARGET 3: Hospitalized prevalence without ventilator ("H_vent_prev") ####
      ## Normal log-likelihood  
      v_llik[j, "H_novent_prev"] <- sum(dnorm(x    = l_hosp_targets$novent$value,#
                                              mean = l_model_res$H_novent_prev$value,
                                              sd   = l_hosp_targets$novent$se,
                                              log  = T),
                                        na.rm = T)
      # print("Log-normal")
      
      # Poisson log-likelihood
      # v_llik[j, "H_vent_prev"] <- sum(dpois(x     = l_hosp_targets$novent$value,
      #                                      lambda = l_model_res$H_vent_prev$value,
      #                                      log    = T), 
      #                                na.rm = T)
      
      # Negative-Binomial log-likelihood
      # v_llik[j, "H_vent_prev"] <- sum(dnbinom(x   = l_hosp_targets$icu$value,
      #                                        size = 1,
      #                                        mu   = l_model_res$H_icu$value,
      #                                        log  = T),
      #                                na.rm = T)
      
      # Multivariate-Normal log-likelihood  
      # v_llik[j, "H_vent_prev"] <- mvnfast::dmvn(X    = l_hosp_targets$icu$value,
      #                                          mu    = l_model_res$H_vent_prev$value,
      #                                          sigma = l_hosp_targets$cov_icu,
      #                                          log   = TRUE)
      # print("MVN")

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

#' Log-likelihood function for a hospitalization parameter set
#'
#' \code{hosp_log_lik} computes a log-likelihood value for one (or multiple) 
#' hospitalization parameter set(s).
#'
#' @param v_params Vector (or matrix) of hospitalization parameters.
#' @param ... Further arguments to be passed to.
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @export
hosp_log_lik <- function(v_params,
                         ...){ # User defined

  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  
  n_samp <- nrow(v_params)
  v_target_names <- c("H_prev", "H_vent_prev", "H_novent_prev"
  )
  n_target       <- length(v_target_names)
  v_llik <- matrix(0, nrow = n_samp, ncol = n_target) 
  colnames(v_llik) <- v_target_names
  v_llik_overall <- numeric(n_samp)
  for(j in 1:n_samp) { # j=2
    jj <- tryCatch( { 

      # m_params_samp <- c()
      # for(i in 1:length(l_params)){
      #   m_params_samp <- rbind(m_params_samp, as.vector(l_params[[i]][j,]))
      # }
      # row.names(m_params_samp) <- names(l_params)

      if(sum(v_params[j,] < get_hosp_bounds()$v_lb | v_params[j,]  > get_hosp_bounds()$v_ub ) > 0){
        v_llik_overall[j] <- NA
      }else{
      
        ###   Run model for parameter set j ###
        l_model_res <- hosp_calibration_out(v_params_calib = v_params[j,] , 
                                            # l_params_all,
                                            # n_lag_inf = 14,
                                            # n_lag_conf = 0,
                                            # l_dates_hosp_targets)
                                            ...)
        l_model_res$H_prev <- l_model_res$H_prev %>% 
          mutate(value = replace(value, value == 0, 1e-20))
        
        l_model_res$H_vent_prev <- l_model_res$H_vent_prev %>% 
          mutate(value = replace(value, value == 0, 1e-20))
        
        l_model_res$H_novent_prev <- l_model_res$H_novent_prev %>% 
          mutate(value = replace(value, value == 0, 1e-20))
        
        ###  Calculate log-likelihood of model outputs to targets  ###
        #### TARGET 1: Hospitalized prevalence ("H_prev") ####
        ## Normal log-likelihood  
        v_llik[j, "H_prev"] <- sum(dnorm(x   = l_hosp_targets$hosp$value,#
                                        mean = l_model_res$H_prev$value,
                                        sd   = l_hosp_targets$hosp$se,
                                        log  = T),
                                  na.rm = T)
        # print("Log-normal")
        
        ## Poisson log-likelihood
        # v_llik[j, "H_prev"] <- sum(dpois(x     = l_hosp_targets$hosp$value,#
        #                                 lambda = l_model_res$H_prev$value,
        #                                 log    = T), 
        #                           na.rm = T)
        # print("Poisson")
        
        ## Negative-Binomial log-likelihood
        # v_llik[j, "H_prev"] <- sum(dnbinom(x   = l_hosp_targets$hosp$value, 
        #                                   size = 1,
        #                                   mu   = l_model_res$H_prev$value,
        #                                   # sd = l_hosp_targets$hosp$se,
        #                                   log  = T), 
        #                           na.rm = T)
        
        # print("Negative-Binomial")
        
        # ## Multivariate-Normal log-likelihood  
        # v_llik[j, "H_prev"] <- mvnfast::dmvn(X    = l_hosp_targets$hosp$value,
        #                                     mu    = l_model_res$H_prev$value,
        #                                     sigma = l_hosp_targets$cov_hosp,
        #                                     log   = TRUE)
        # print("MVN")
        
        #### TARGET 2: Hospitalized prevalence with ventilator ("H_vent_prev") ####
        ## Normal log-likelihood  
        v_llik[j, "H_vent_prev"] <- sum(dnorm(x    = l_hosp_targets$vent$value,#
                                              mean = l_model_res$H_vent_prev$value,
                                              sd   = l_hosp_targets$vent$se,
                                              log  = T),
                                        na.rm = T)
        # print("Log-normal")
        
        # Poisson log-likelihood
        # v_llik[j, "H_vent_prev"] <- sum(dpois(x     = l_hosp_targets$vent$value,
        #                                      lambda = l_model_res$H_vent_prev$value,
        #                                      log    = T), 
        #                                na.rm = T)
        
        # Negative-Binomial log-likelihood
        # v_llik[j, "H_vent_prev"] <- sum(dnbinom(x   = l_hosp_targets$vent$value,
        #                                        size = 1,
        #                                        mu   = l_model_res$H_vent_prev$value,
        #                                        log  = T),
        #                                na.rm = T)
        
        # Multivariate-Normal log-likelihood  
        # v_llik[j, "H_vent_prev"] <- mvnfast::dmvn(X    = l_hosp_targets$vent$value,
        #                                          mu    = l_model_res$H_vent_prev$value,
        #                                          sigma = l_hosp_targets$cov_vent,
        #                                          log   = TRUE)
        #print("MVN")
        
        #### TARGET 3: Hospitalized prevalence without ventilator ("H_novent_prev") ####
        ## Normal log-likelihood  
        v_llik[j, "H_novent_prev"] <- sum(dnorm(x    = l_hosp_targets$novent$value,#
                                              mean = l_model_res$H_novent_prev$value,
                                              sd   = l_hosp_targets$novent$se,
                                              log  = T),
                                        na.rm = T)
        # print("Log-normal")
        
        # Poisson log-likelihood
        # v_llik[j, "H_novent_prev"] <- sum(dpois(x     = l_hosp_targets$novent$value,
        #                                      lambda = l_model_res$H_novent_prev$value,
        #                                      log    = T), 
        #                                na.rm = T)
        
        # Negative-Binomial log-likelihood
        # v_llik[j, "H_novent_prev"] <- sum(dnbinom(x   = l_hosp_targets$novent$value,
        #                                        size = 1,
        #                                        mu   = l_model_res$H_novent_prev$value,
        #                                        log  = T),
        #                                na.rm = T)
        
        # Multivariate-Normal log-likelihood  
        # v_llik[j, "H_novent_prev"] <- mvnfast::dmvn(X    = l_hosp_targets$novent$value,
        #                                          mu    = l_model_res$H_novent_prev$value,
        #                                          sigma = l_hosp_targets$cov_novent,
        #                                          log   = TRUE)
        #print("MVN")
        
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

#' Parallel evaluation of log-likelihood function for a sets of hospitalization
#' parameters
#'
#' \code{hosp_log_lik_par} computes a log-likelihood value for one (or multiple) 
#' hospitalization parameter set(s) using parallel computation.
#'
#' @param v_params Vector (or matrix) of model parameters.
#' @param log_lik_offset Offset applied to log-likelihood values. 
#' @param ... Further arguments to be passed to.
#' @return 
#' A scalar (or vector) with log-likelihood values.
#' @importFrom stats dnorm dunif quantile qunif rbeta rgamma sd
#' @export
hosp_log_lik_par <- function(v_params, 
                             log_lik_offset,
                             ...) { 
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  
  n_samp <- nrow(v_params)
  
  ### Get OS
  os <- get_os()
  print(paste0("Parallelized Likelihood calculations on ", os))

  no_cores <- detectCores() - 10
  
  n_time_init_likpar <- Sys.time()
  
  if(os == "macosx"){
    # Initialize cluster object
    cl <- makeForkCluster(no_cores) 
    registerDoParallel(cl)
    llk <- foreach(i = 1:n_samp, .combine = c) %dopar% {
      
      # l_temp <- list()
      # for(j in 1:length(l_params)){ # j = 1
      #   l_temp[[j]] <- t(as.matrix(l_params[[j]][i,]))
      # }
      # names(l_temp) <- names(l_params)
      
      hosp_log_lik(v_params = v_params[i,] ,
              # l_params_all = l_params_all,
              # n_lag_inf = 14,
              # n_lag_conf = 0,
              # l_dates_hosp_targets)
              ...)
    }
    n_time_end_likpar <- Sys.time()
  }
  if(os == "windows"){
    # Initialize cluster object
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    opts <- list(attachExportEnv = TRUE)
    llk <- foreach(i = 1:n_samp, .combine = c,
                   .export = ls(globalenv()),
                   .packages=c("sccosmomcma",
                               "ggplot2",
                               "tidyverse",
                               "dplyr",
                               "lubridate",
                               "epitools"
                   ),
                   .options.snow = opts) %dopar% { # i = 1
                     
                     # l_temp <- list()
                     # for(j in 1:length(l_params)){ # j = 1
                     #   l_temp[[j]] <- t(as.matrix(l_params[[j]][i,]))
                     # }
                     # names(l_temp) <- names(l_params)
                     
                     hosp_log_lik(v_params = v_params[i,] ,
                             # l_params_all = l_params_all,
                             # n_lag_inf = 14,
                             # n_lag_conf = 0,
                             # l_dates_hosp_targets)
                             ...)
                   }
    n_time_end_likpar <- Sys.time()
  }
  if(os == "linux"){
    # Initialize cluster object
    cl <- makeCluster(no_cores)
    registerDoMC(cl)
    llk <- foreach(i = 1:n_samp, .combine = c) %dopar% {
      
      # l_temp <- list()
      # for(j in 1:length(l_params)){ # j = 1
      #   l_temp[[j]] <- t(as.matrix(l_params[[j]][i,]))
      # }
      # names(l_temp) <- names(l_params)
      
      hosp_log_lik(v_params = v_params[i,],
                   # l_params_all = l_params_all,
                   # n_lag_inf = 14,
                   # n_lag_conf = 0,
                   # l_dates_hosp_targets)
                   ...)
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
#' \code{hosp_likelihood} computes a likelihood value for one (or multiple) 
#' hospitalization parameter set(s).
#'
#' @param v_params Vector (or matrix) of hospitalization parameters. 
#' @return 
#' A scalar (or vector) with likelihood values.
#' @export
hosp_likelihood <- function(v_params){ 
  v_like <- exp(hosp_log_lik_par(v_params, 
                            log_lik_offset = 0,
                            l_params_all,
                            n_lag_inf, 
                            n_lag_conf, 
                            l_dates_hosp_targets))
  return(v_like)
}

#' Evaluate log-posterior of calibrated hospitalization parameters
#'
#' \code{hosp_log_post} computes a log-posterior value for one (or multiple) 
#' hospitalization parameter set(s) based on the simulation model, 
#' likelihood functions and prior distributions.
#' 
#' @param v_params Vector (or matrix) of model parameters.
#' @param ... Further arguments to be passed to.
#' @return 
#' A scalar (or vector) with log-posterior values.
#' @export
hosp_log_post <- function(v_params, ...) { 
  v_lpost <- hosp_log_prior(v_params) + hosp_log_lik(v_params, ...)
  return(v_lpost) 
}

#' Evaluate log-posterior of calibrated hospitalization parameters for OPTIMIZATION purposes
#'
#' \code{hosp_log_post_opt} computes a log-posterior value for one (or multiple) 
#' hospitalization parameter set(s) based on the simulation model, likelihood 
#' functions and prior distributions. Used for OPTIMIZATION purposes, 
#' NOT Bayesian estimation.
#' 
#' @param v_params Vector (or matrix) of hospitalization parameters.
#' @param ... Further arguments to be passed to.
#' @return 
#' A scalar (or vector) with log-posterior values.
#' @export
hosp_log_post_opt <- function(v_params, ...) { 
  
  v_param_names = names(get_hosp_bounds()$v_lb) #make sure this is actually same name as in model params
  v_lb = get_hosp_bounds()$v_lb # get rid of defaults and then user can pass in
  v_ub = get_hosp_bounds()$v_ub
  if(is.null(dim(v_params))) { # If vector, change to matrix
    v_params <- t(v_params)
  }
  n_param <- length(v_param_names)
  n_samp <- nrow(v_params)
  #  colnames(v_params) <- v_param_names
  lprior <- rep(0, n_samp)
  oob <- FALSE
  for (i in 1:n_param){ # i = 1
    if (v_params[,i]<v_lb[i] | v_params[, i] > v_ub[i]) {
      lprior <- lprior + 
        -1000000000*ifelse(v_params[,i]<v_lb[i], 
                           abs(v_lb[i]-v_params[,i]), 
                           abs(v_ub[i]-v_params[,i])) + -10000000
      oob <- TRUE
    }
  }
  if (oob==TRUE) {
    return(lprior)
  }
  v_lpost <- hosp_log_prior_opt(v_params = v_params,
                                v_param_names = v_param_names, #make sure this is actually same name as in model params
                                v_lb = v_lb, # get rid of defaults and then user can pass in
                                v_ub = v_ub) + 
    hosp_log_lik_opt(v_params             = v_params,
                     l_params_all         = l_params_all,
                     n_lag_inf            = 14,
                     n_lag_conf           = 0,
                     l_dates_hosp_targets = l_dates_hosp_targets)
                # ...)
  return(v_lpost) 
}

#' Evaluate posterior of calibrated hospitalization parameters
#'
#' \code{hosp_posterior} computes a posterior value for one (or multiple) hospitalization
#' parameter set(s).
#' 
#' @param v_params Vector (or matrix) of model parameters.
#' @param ... Further arguments to be passed to. 
#' @return 
#' A scalar (or vector) with posterior values.
#' @export
hosp_posterior <- function(v_params, ...) { 
  v_posterior <- exp(hosp_log_post(v_params)) 
  return(v_posterior)
}

#' Get bounds for hospitalization parameters to calibrate
#'
#' \code{get_hosp_bounds} defines the bounds for each of the calibrated 
#' hospitalization parameters.
#' 
#' @return 
#' A list with lower and upper bounds and their standard errors.
#' @export
get_hosp_bounds <- function() {
  
  # lower bounds
  v_lb <- c(m_r_exit_tot    = 0,
            m_r_exit_nonicu = 0,
            m_r_exit_icu    = 0,
            m_sigma_tot     = 0,
            m_sigma_nonicu  = 0,
            m_sigma_icu     = 0
            )
  
  ## upper bounds
  v_ub <- c(m_r_exit_tot    = 1,
            m_r_exit_nonicu = 1,
            m_r_exit_icu    = 1,
            m_sigma_tot     = 10,
            m_sigma_nonicu  = 10,
            m_sigma_icu     = 30
  )
  
  ## standard errors based on bounds
  v_se <- (v_ub - v_lb)/(2*2) # even more conservative an 1.96
  
  return(list(v_lb = v_lb, 
              v_ub = v_ub,
              v_se = v_se))
}

#' Sample from prior distributions of calibrated parameters
#'
#' \code{hosp_sample_prior} generates a sample of hospitalization parameter sets 
#' from their prior distribution.
#' 
#' @param n_samp Number of samples.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A matrix with 6 columns and \code{n_samp} rows. Each row corresponds to a 
#' parameter set sampled from their prior distributions.
#' @export
hosp_sample_prior <- function(n_samp,
                         v_param_names = names(get_hosp_bounds()$v_lb), #make sure this is actually same name as in model params
                         v_lb = get_hosp_bounds()$v_lb, # get rid of defaults and then user can pass in 
                         v_ub = get_hosp_bounds()$v_ub){
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

#' Evaluate log-prior of calibrated hospitalization parameters
#'
#' \code{hosp_log_prior} computes a log-prior value for one (or multiple) 
#' hospitalization parameter set(s) based on their prior distributions.
#' 
#' @param v_params Vector (or matrix) of model parameters.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A scalar (or vector) with log-prior values.
#' @export
hosp_log_prior <- function(v_params, 
                           v_param_names = names(get_hosp_bounds()$v_lb), #make sure this is actually same name as in model params
                           v_lb = get_hosp_bounds()$v_lb, # get rid of defaults and then user can pass in 
                           v_ub = get_hosp_bounds()$v_ub){
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

#' Evaluate log-prior of calibrated hospitalization parameters for OPTIMIZATION purposes
#'
#' \code{hosp_log_prior_opt} computes a log-prior value for one (or multiple)
#' hospitalization parameter set(s) based on their prior distributions. Used for 
#' OPTIMIZATION purposes, NOT Bayesian estimation.
#' 
#' @param v_params Vector (or matrix) of model parameters.
#' @param v_param_names Vector with parameter names.
#' @param v_ub Vector with lower bounds for each parameter.
#' @param v_lb Vector with upper bounds for each parameter.
#' @return 
#' A scalar (or vector) with log-prior values.
#' @export
hosp_log_prior_opt <- function(v_params, 
                               v_param_names = names(get_hosp_bounds()$v_lb), #make sure this is actually same name as in model params
                               v_lb = get_hosp_bounds()$v_lb, # get rid of defaults and then user can pass in 
                               v_ub = get_hosp_bounds()$v_ub){
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

#' Evaluate prior of calibrated hospitalization parameters
#'
#' \code{hosp_prior} computes a prior value for one (or multiple) hospitalization
#' parameter set(s).
#' @param v_params Vector (or matrix) of model parameters. 
#' @return 
#' A scalar (or vector) with prior values.
#' @export
hosp_prior <- function(v_params) { 
  v_prior <- exp(hosp_log_prior(v_params)) 
  return(v_prior)
}

#' Produce hospitalization occupation over time
#' 
#' \code{prep_dx_hosp_calib} computes both incident hospitalizations and prevalence 
#' in hospitalizations (Tot, Non ICU and ICU) by age groups to be used for 
#' calibration routines.
#' 
#' @param l_out_cosmo List with output from SC-COSMO and all parameters.
#' @param use_prevalence Flag (default is FALSE) of whether to use prevalence as 
#' input to hospitalization model. If FALSE, incidence is used instead of 
#' prevalence.
#' @return 
#' A list of incident hospitalizations by type (Tot , Non ICU, or ICU), 
#' age group (including All), and time
#' and of prevalent hospitalizations by type (Tot , Non ICU, or ICU), 
#' age group (including All), and time.
#' @export
prep_dx_hosp_calib <- function(l_params_all, use_prevalence = FALSE) {
  with(as.list(l_params_all), {  
    
    df_temp <- l_targets$cases_inc_age %>%
      ungroup()

    m_inc_dx <- as.matrix(df_temp[14:21])
    m_inc_dx <- t(m_inc_dx)
    n_time <- length(df_temp$Date)
    
    ## Obtain time-varying gamma parameters of exit hospitalization
    v_gamma_params_tot <- get_gamma_params(m_r_exit     = m_r_exit_tot, 
                                           m_sigma_exit = m_sigma_tot, 
                                           gamma_times  = n_time)
    
    v_gamma_params_nonicu  <- get_gamma_params(m_r_exit    = m_r_exit_nonicu, 
                                               m_sigma_exit = m_sigma_nonicu, 
                                               gamma_times  = n_time)
    
    v_gamma_params_icu <- get_gamma_params(m_r_exit     = m_r_exit_icu, 
                                           m_sigma_exit = m_sigma_icu, 
                                           gamma_times  = n_time)

    # Age groups plus "All"
    n_age_groups <- n_ages + 1
    
    # Empty matrices to fill
    m_hosp_inc  <- array(0, dim = c(3, n_age_groups, n_time),
                         dimnames = list(c("tot_hosp", "novent_hosp", "vent_hosp"), 
                                         c(rownames(m_inc_dx), "All"), 
                                         1:n_time))
    m_hosp_prev <- array(0, dim = c(3, n_age_groups, n_time), 
                         dimnames = list(c("tot_hosp", "novent_hosp", "vent_hosp"), 
                                         c(rownames(m_inc_dx), "All"), 
                                         1:n_time)) 
    
    max_forward_flow <- 3*ceiling(1/min(m_r_exit_tot, m_r_exit_nonicu, m_r_exit_icu))
    
    m_hosp_inc[1, 1:n_ages, ] <- m_inc_dx * m_p_tot_hosp
    m_hosp_inc[2, 1:n_ages, ] <- m_inc_dx * m_p_tot_hosp * m_p_nonicu_hosp
    m_hosp_inc[3, 1:n_ages, ] <- m_inc_dx * m_p_tot_hosp * m_p_icu_hosp
    
    m_hosp_inc[1, n_age_groups, ] <-  colSums(m_hosp_inc[1, 1:n_ages, ])
    m_hosp_inc[2, n_age_groups, ] <-  colSums(m_hosp_inc[2, 1:n_ages, ])
    m_hosp_inc[3, n_age_groups, ] <-  colSums(m_hosp_inc[3, 1:n_ages, ])  
    
    for(t in 1:(n_time)) { # t = 1
      
      ## COMPUTE HOSPITALIZED PREVALENCE BASED ON DISTRIBUTION OF LENGTH OF STAY/GAMMA
      for(k in 0:max_forward_flow) { # k = 0
        ## IF FUTURE TIME IS NOT AFTER END OF SIMULATION TOTAL TIME
        if (t+k<=n_time) {
          fraction_remaining_tot <- 1 - pgamma(q=k, shape=v_gamma_params_tot[[t]]$shape, rate=v_gamma_params_tot[[t]]$rate)
          fraction_remaining_nonicu  <- 1 - pgamma(q=k, shape=v_gamma_params_nonicu[[t]]$shape,  rate=v_gamma_params_nonicu[[t]]$rate)
          fraction_remaining_icu <- 1 - pgamma(q=k, shape=v_gamma_params_icu[[t]]$shape, rate=v_gamma_params_icu[[t]]$rate)

          m_hosp_prev[1, 1:n_ages, t+k] <- m_hosp_prev[1, 1:n_ages, t+k] + (fraction_remaining_tot)*m_hosp_inc[1, 1:n_ages, t]
          m_hosp_prev[2, 1:n_ages, t+k] <- m_hosp_prev[2, 1:n_ages, t+k] + (fraction_remaining_nonicu)*m_hosp_inc[2, 1:n_ages, t]
          m_hosp_prev[3, 1:n_ages, t+k] <- m_hosp_prev[3, 1:n_ages, t+k] + (fraction_remaining_icu)*m_hosp_inc[3, 1:n_ages, t]        
        }
      }
    }
    
    ## Compute total prevalence across age groups
    m_hosp_prev[1, n_age_groups, ] <-  colSums(m_hosp_prev[1, 1:n_ages, ])
    m_hosp_prev[2, n_age_groups, ] <-  colSums(m_hosp_prev[2, 1:n_ages, ])
    m_hosp_prev[3, n_age_groups, ] <-  colSums(m_hosp_prev[3, 1:n_ages, ])
    
    ret_list <- list(m_hosp_inc, m_hosp_prev)
    names(ret_list) <- c("m_hosp_inc", "m_hosp_prev")
    return(ret_list)
  })
}
