#' Get parameters for gamma distribution 
#'
#' \code{get_gamma_params} Get parameters for gamma distribution based on mean 
#' rates of exiting and standard deviation of time to exiting over time
#' @param m_r_exit Matrix or vector with mean rates of exiting by age groups 
#' over time (matrix) or for a given day (vector)
#' @param m_sigma_exit Matrix or vector with standard deviations of time to 
#' exiting by age groups over time (matrix) or for a given day (vector)
#' @param gamma_times Number of times for which time-varying gamma parameters
#' should be produced.
#' @return 
#' List with time-specific shape and rate parameters for the gamma distribution
#' @export
get_gamma_params <- function(m_r_exit, m_sigma_exit, gamma_times){
  if(is.null(dim(m_r_exit)) & is.null(dim(m_sigma_exit))){
    m_r_exit     <- matrix(rep(m_r_exit, gamma_times), ncol = gamma_times)
    m_sigma_exit <- matrix(rep(m_sigma_exit, gamma_times), ncol = gamma_times)
  }
  if(!is.null(dim(m_r_exit)) & !is.null(dim(m_sigma_exit))){
    if(sum(dim(m_r_exit) == dim(m_sigma_exit)) != length(dim(m_r_exit)))
      stop("Dimensions of matrix of rates of exiting and matrix of variances are not the same")
  }
  if(!is.null(dim(m_r_exit)) & is.null(dim(m_sigma_exit))){
    stop("Variances is a vector and rates of exit is a matrix")
  }
  if(is.null(dim(m_r_exit)) & !is.null(dim(m_sigma_exit))){
    stop("Rates of exit is a vector and variances is a matrix")
  }
  
  n_times_hosp <- ncol(m_r_exit)
  
  l_gamma_params_exit <- vector(mode = "list", length = gamma_times)
  
  for(gtime in 1:gamma_times){
    gtime_index <- min(gtime, n_times_hosp)
    l_gamma_params_exit[[gtime]] <- dampack::gamma_params(mu = (1/m_r_exit[, gtime_index]), 
                                                          sigma = m_sigma_exit[, gtime_index], 
                                                          scale = F)  
  }
  return(l_gamma_params_exit)
}

#' Get proportion of cases that go into hospitalization over time
#'
#' \code{get_prop_hosp} pulls age-specific proportion of incident detected cases 
#' hospitalized for a given time
#' @param time the time (numeric, in days) at which p_hosp is evaluated
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' Proportion hospitalized (age-specific) at a given point in time.
#' @export
get_prop_hosp <- function(time, l_params_all) {
#  with(as.list(l_params_all), {
    
  v_p_hosp_time <- as.vector(unlist(lapply(l_params_all$m_p_hosp, do.call, list(time))))
  
    #v_p_hosp_time <- as.vector(l_params_all$m_p_hosp[,floor(time)+1])
    
    return(v_p_hosp_time)
    
#  }
#  )
}

#' Produce hospitalization occupation over time
#' 
#' Computes both incident hospitalizations and prevalence in hospital
#' also by Non ICU and ICU, and age groups.
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @param use_prevalence Indicator to use prevalence as input to hospitalization
#' model. If FALSE, incidence is used instead of prevalence (Default == FALSE).
#' @return 
#' A list of incident hospitalizations by type (Tot , Non ICU, or ICU), 
#' age group (including All), and time
#' and of prevalent hospitalizations by type (Tot , Non ICU, or ICU), 
#' age group (including All), and time
#' @export
prep_dx_hospitalizations <- function(l_out_cosmo, use_prevalence = FALSE) {
  with(as.list(l_out_cosmo$l_params_all), {  
    
    ## GET DETECTED INCIDENCE (PREVALENCE??? THE LANCET ARTICLE IS AMBIGUOUS)
    if(use_prevalence){
      df_dxcases <- calc_dx_totals(l_out_cosmo)
    } else {
      df_dxcases <- calc_dxinc_totals(l_out_cosmo)
    }
    
    dx <- array(df_dxcases)
    ## chopping off last column (All)
    ## chopping off first column (time)
    dx <- dx[, -10]
    dx <- dx[, -1]
    
    n_time <- length(df_dxcases$time)
    
    ## Obtain time-varying gamma parameters of exit hospitalization
    v_gamma_params_tot <- get_gamma_params(m_r_exit = m_r_exit_tot, 
                                           m_sigma_exit = m_sigma_tot, 
                                           gamma_times = n_time)
    
    v_gamma_params_nonicu  <- get_gamma_params(m_r_exit = m_r_exit_nonicu, 
                                               m_sigma_exit = m_sigma_nonicu, 
                                               gamma_times = n_time)
    
    v_gamma_params_icu <- get_gamma_params(m_r_exit = m_r_exit_icu, 
                                           m_sigma_exit = m_sigma_icu, 
                                           gamma_times = n_time)
    
    ## Create matrices to store prevalent and incident hospitalization by 
    ## severity and for each age groups
    # extra row for All
    n_age_groups <- n_ages + 1
    m_hosp_inc  <- array(0, dim = c(3, n_age_groups, n_time))
    m_hosp_prev <- array(0, dim = c(3, n_age_groups, n_time))  
    
    inc_dx <- t(dx)
    
    max_forward_flow <- 3*ceiling(1/min(m_r_exit_tot, m_r_exit_nonicu, m_r_exit_icu))
    #max_forward_flow <- n_time-t
    
    ## COMPUTE HOSPITALIZED INCIDENCE AND PREVALENCE
    m_hosp_inc[1, 1:n_ages, ] <- inc_dx * m_p_tot_hosp_proj
    m_hosp_inc[2, 1:n_ages, ] <- inc_dx * m_p_tot_hosp_proj * m_p_nonicu_hosp_proj
    m_hosp_inc[3, 1:n_ages, ] <- inc_dx * m_p_tot_hosp_proj * m_p_icu_hosp_proj
    
    m_hosp_inc[1, n_age_groups, ] <-  colSums(m_hosp_inc[1, 1:n_ages, ])
    m_hosp_inc[2, n_age_groups, ] <-  colSums(m_hosp_inc[2, 1:n_ages, ])
    m_hosp_inc[3, n_age_groups, ] <-  colSums(m_hosp_inc[3, 1:n_ages, ])  
    
    for(t in 1:n_time) { # t = 4
      
      ## COMPUTE HOSPITALIZED PREVALENCE BASED ON DISTRIBUTION OF LENGTH OF STAY/GAMMA
      for(k in 0:max_forward_flow) {
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


#' Total Incident DX hospitalizations
#' 
#' Calculate total number of incident diagnosed hospitalizations.
#' 
#' @param l_hosp List with output from the hospitalization simulation 
#' (prep_hospitalizations)
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of total number of incident diagnosed 
#' hospitalizations by age group and overall as columns over time.
#' @export
calc_dx_inchosp<-function(l_hosp, l_out_cosmo) {
  m_hosp_inc   <- l_hosp$m_hosp_inc
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  df_HospInctot <- data.frame(time = l_out_cosmo$df_out_cosmo$time, 
                              t(m_hosp_inc[1,,]),
                              check.names = FALSE)
  colnames(df_HospInctot)[-1] <- c(levels(v_names_ages), "All")
  return(df_HospInctot)  
}

#' Total Incident DX ICU hospitalizations
#' 
#' Calculate total number of incident diagnosed ICU hospitalizations.
#' 
#' @param l_hosp List with output from the hospitalization simulation 
#' (prep_hospitalizations)
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of total number of incident diagnosed 
#' ICU hospitalizations by age group and overall as columns over time.
#' @export
calc_dx_inchosp_icu<-function(l_hosp, l_out_cosmo) {
  m_hosp_inc   <- l_hosp$m_hosp_inc
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  df_ICUInctot <- data.frame(time = l_out_cosmo$df_out_cosmo$time, 
                             t(m_hosp_inc[3,,]),
                             check.names = FALSE)
  colnames(df_ICUInctot)[-1] <- c(levels(v_names_ages), "All")
  return(df_ICUInctot)
}

#' Total Incident DX Non-ICU hospitalizations
#' 
#' Calculate total number of incident diagnosed Non-ICU hospitalizations.
#' 
#' @param l_hosp List with output from the hospitalization simulation 
#' (prep_hospitalizations)
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of total number of incident diagnosed 
#' Non-ICU hospitalizations by age group and overall as columns over time.
#' @export
calc_dx_inchosp_nonicu<-function(l_hosp, l_out_cosmo) {
  m_hosp_inc   <- l_hosp$m_hosp_inc
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  df_NonICUInctot <- data.frame(time = l_out_cosmo$df_out_cosmo$time, 
                                t(m_hosp_inc[2,,]),
                                check.names = FALSE)
  colnames(df_NonICUInctot)[-1] <- c(levels(v_names_ages), "All")
  return(df_NonICUInctot)    
}

#' Total Prevalent DX hospitalizations
#' 
#' Calculate total number of prevalent diagnosed hospitalizations.
#' 
#' @param l_hosp List with output from the hospitalization simulation 
#' (prep_hospitalizations)
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of total number of prevalent diagnosed 
#' hospitalizations by age group and overall as columns over time.
#' @export
calc_dx_hosp<-function(l_hosp, l_out_cosmo) {
  m_hosp_prev   <- l_hosp$m_hosp_prev
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  df_Hosptot <- data.frame(time = l_out_cosmo$df_out_cosmo$time, 
                           t(m_hosp_prev[1,,]),
                           check.names = FALSE)
  colnames(df_Hosptot)[-1] <- c(levels(v_names_ages), "All")
  return(df_Hosptot)   
}

#' Total Prevalent DX ICU hospitalizations
#' 
#' Calculate total number of prevalent diagnosed ICU hospitalizations.
#' 
#' @param l_hosp List with output from the hospitalization simulation 
#' (prep_hospitalizations)
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of total number of prevalent diagnosed 
#' ICU hospitalizations by age group and overall as columns over time.
#' @export
calc_dx_hosp_icu<-function(l_hosp, l_out_cosmo) {
  m_hosp_prev   <- l_hosp$m_hosp_prev
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  df_ICUtot <- data.frame(time = l_out_cosmo$df_out_cosmo$time, 
                          t(m_hosp_prev[3,,]),
                          check.names = FALSE)
  colnames(df_ICUtot)[-1] <- c(levels(v_names_ages), "All")
  return(df_ICUtot)   
}

#' Total Prevalent DX Non-ICU hospitalizations
#' 
#' Calculate total number of prevalent diagnosed Non-ICU hospitalizations.
#' 
#' @param l_hosp List with output from the hospitalization simulation 
#' (prep_hospitalizations)
#' @param l_out_cosmo List with output from SC-COSMO and all parameters
#' @return 
#' A data.frame with the total number of total number of prevalent diagnosed 
#' Non-ICU hospitalizations by age group and overall as columns over time.
#' @export
calc_dx_hosp_nonicu<-function(l_hosp, l_out_cosmo) {
  m_hosp_prev   <- l_hosp$m_hosp_prev
  l_params_all <- l_out_cosmo$l_params_all
  
  v_names_ages <- l_params_all$v_names_ages
  
  df_NonICUtot <- data.frame(time = l_out_cosmo$df_out_cosmo$time, 
                             t(m_hosp_prev[2,,]),
                             check.names = FALSE)
  colnames(df_NonICUtot)[-1] <- c(levels(v_names_ages), "All")
  return(df_NonICUtot)   
}
