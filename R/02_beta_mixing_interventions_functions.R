#' Computing force of infection at each time step for the Stanford-CIDE 
#' COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{get_betas} provides a wrapper to compute the force of infection
#' at each time step for the SC-COSMO model. It handles
#' applying intervention effects and timing to various WAIFW matrices.
#' 
#' @param time Current time at which the model is evaluated.
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' A matrix of the beta components of the forces of infection for each subgroup 
#' applied to the overall (composite) mixing matrix in the model.
#' 
#' @export
get_betas <-function(time, l_params_all) {
  with(as.list(l_params_all), {
    m_betas_intervention <- m_betas_all_internal[, , floor(time)+1]
    return(m_betas_intervention)
  }
  )
}

#' Gets correct mixing matrices for the Stanford-CIDE COronavirus Simulation 
#' MOdel (SC-COSMO)
#'
#' \code{get_waifw} provides a wrapper to determine which WAIFW matrices
#' are the correct ones given intervention timings. 
#' 
#' @param time Time at which the model is evaluated.
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' An overall WAIFW matrix appropriate for the model given the 
#' current time and interventions being used.
#' 
#' @export
get_waifw <- function(time, l_params_all) {
  with(as.list(l_params_all), {
    m_curr_waifw <- m_waifws_all_internal[, , floor(time)+1]
    m_curr_waifw <- as.matrix(m_curr_waifw)
    return(m_curr_waifw)
  }
  )
}

#' Gets correct intervention effect multipliers for the Stanford-CIDE COronavirus
#' Simulation MOdel (SC-COSMO)
#'
#' \code{get_intervention_effects} provides a wrapper to determine
#' the intervention effects given intervention timings. 
#' 
#' @param time Current time at which the model is evaluated.
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' A vector of intervention effects appropriate for the model given the 
#' current time and interventions being used.
#' 
#' @export
get_intervention_effects <- function(time, l_params_all) {
  with(as.list(l_params_all), {
    v_intervention_current <- m_intervention_effects_all_internal[, floor(time)+1]
    return(v_intervention_current)
  }
  )  
}

#' Efficiency by precomputing WAIFW matrices at each time step for
#' the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{gen_all_waifws} provides a way to precompute all WAIFW matrices
#' given the pattern of interventions at each time step for the 
#' SC-COSMO model. 
#' 
#' @param l_interventions List of interventions and timing.
#' @param l_contact_info List of WAIFW matrix inputs.
#' @param v_names_ages Vector of names for rows and columns of WAIFW matrices.
#' @param max_time Number of time steps to precompute WAIFW matrices for.
#' 
#' @return 
#' A 3D matrix of appropriate WAIFW matrices for each time step in the model.
#' 
#' @export
gen_all_waifws <- function(l_interventions, l_contact_info, v_names_ages, max_time) {

  n_ages <- length(v_names_ages)
  m_all_waifws <- array(0, 
                        dim = c(n_ages, n_ages, max_time+1), 
                        dimnames = list(v_names_ages, v_names_ages, seq(0, max_time)))
  
  m_zero_contacts <- array(0, dim = c(n_ages, n_ages), dimnames(list(v_names_ages, v_names_ages)))
  m_school_waifws <- array(0, 
                        dim = c(n_ages, n_ages, max_time+1), 
                        dimnames = list(v_names_ages, v_names_ages, seq(0, max_time)))
  
  idx = 1
  for(time_idx in 1:(max_time+1)) {
    if (time_idx > l_interventions[[idx]]$stop) {
      idx <- idx+1
      if (idx > length(l_interventions)) {
        break
      }
    }
       
    if (l_interventions[[idx]]$type == "SocialDistancing" | 
        l_interventions[[idx]]$type == "SocialDistancingLinear") {
        m_all_waifws[, , time_idx] <- l_contact_info$m_contact_work + 
                                      l_contact_info$m_contact_other_locations
        m_school_waifws[, , time_idx] <- m_zero_contacts
        if ("resume_school" %in% names(l_interventions[[idx]])) {
          if (l_interventions[[idx]]$resume_school == TRUE) {
            m_school_waifws[, , time_idx] <- l_contact_info$m_contact_school
          } 
        }
    } else {
      m_all_waifws[, , time_idx] <- l_contact_info$m_contact_work + 
                                    l_contact_info$m_contact_other_locations 
      
      m_school_waifws[, , time_idx] <- l_contact_info$m_contact_school
    }
  }    
  
  l_all_waifws <- list(m_all_waifws    = m_all_waifws,
                       m_school_waifws = m_school_waifws)
  
  return(l_all_waifws)
}

#' Efficiency by precomputing intervention effects at each time step for
#' the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{gen_all_intervention_effects} provides a way to precompute all 
#' intervention effects at each time step for the SC-COSMO model. 
#' 
#' @param l_interventions List of interventions and timing.
#' @param v_names_ages Vector of names for rows and columns of WAIFW matrices.
#' @param max_time Number of time steps to precompute intervention effects for.
#' 
#' @return 
#' A 2D matrix of appropriate intervention vectors for each time step 
#' in the model.
#' 
#' @export
gen_all_intervention_effects <- function(l_interventions, v_names_ages, max_time) {
  
  n_ages <- length(v_names_ages)
  m_all_intervention_effects <- array(0, dim=c(n_ages,max_time+1))
  m_school_intervention_effects <- array(0, dim=c(n_ages,max_time+1))
  
  #    print(l_interventions)
  ### find index of the intervention that is active
  ### if it is SocialDistancing, 
  ### also find the index of the intervention that immediately preceded it

  idx = 1
  #  for(idx in 1:n_intervention_periods) {
  for(time_idx in 1:(max_time+1)) {
    if (time_idx > l_interventions[[idx]]$stop) {
      idx <- idx+1
      if (idx > length(l_interventions)) {
        break
      }
    }

    active_intervention_idx         <- idx
    active_intervention_type        <- l_interventions[[idx]]$type
    active_intervention_start       <- l_interventions[[idx]]$start
    active_intervention_end         <- l_interventions[[idx]]$stop
    active_intervention_effect      <- 1
    active_intervention_change_rate <- 0.5
    
    if (active_intervention_type == "SocialDistancing") {
      active_intervention_effect <- l_interventions[[idx]]$intervention_factor
      active_intervention_change_rate <- l_interventions[[idx]]$intervention_change_rate
    } else if (active_intervention_type == "SocialDistancingLinear") {
      active_intervention_effect_start <- l_interventions[[idx]]$intervention_factor
      active_intervention_effect_end   <- l_interventions[[idx]]$intervention_factor_end
      active_intervention_day_slope    <- (active_intervention_effect_end - active_intervention_effect_start) / 
                                          ((active_intervention_end - active_intervention_start) + 1)
    }
    
    prev_intervention_effect <- 1
    if(active_intervention_type == "SocialDistancing") {
      if (active_intervention_idx > 1) {
        if(l_interventions[[active_intervention_idx-1]]$type == "SocialDistancing") {
          prev_intervention_effect <- l_interventions[[active_intervention_idx-1]]$intervention_factor
        }
        else if (l_interventions[[active_intervention_idx-1]]$type == "SocialDistancingLinear") {
          prev_intervention_effect <- l_interventions[[active_intervention_idx-1]]$intervention_factor_end
        }
      }
    }
    
    if (active_intervention_type == "SocialDistancingLinear") {
      v_soc_dist_factor_time_dep <- active_intervention_effect_start + 
                                    active_intervention_day_slope * (time_idx - active_intervention_start)
    } else {
      v_soc_dist_factor_time_dep <- m_general_logit(v_logit_lb = rep(prev_intervention_effect, n_ages), 
                                                  v_logit_ub = rep(active_intervention_effect, n_ages), 
                                                  v_logit_change_rate = rep(active_intervention_change_rate, n_ages), 
                                                  v_logit_change_mid  = rep(active_intervention_start, n_ages),
                                                  logit_t = time_idx)
     
    }
    m_all_intervention_effects[,time_idx] <- v_soc_dist_factor_time_dep

    v_school_factor_time_dep <- v_soc_dist_factor_time_dep    
    if (active_intervention_type == "SocialDistancing" | 
        active_intervention_type == "SocialDistancingLinear") {
      if ("resume_school" %in% names(l_interventions[[idx]])) {
        if (l_interventions[[idx]]$resume_school == TRUE) {
          if ("school_intervention_factor" %in% names(l_interventions[[idx]])) {
            v_school_factor_time_dep <- rep(l_interventions[[idx]]$school_intervention_factor, n_ages)
          } else{
            v_school_factor_time_dep <- rep(1, n_ages)
          }
        } 
      }
    }
    m_school_intervention_effects[, time_idx] <- v_school_factor_time_dep
  }
  
  l_all_intervention_effects <- list(m_all_intervention_effects    = m_all_intervention_effects,
                                     m_school_intervention_effects = m_school_intervention_effects)
  
  return(l_all_intervention_effects)
}

#' Efficiency by precomputing beta component of force of infection 
#' at each time step for the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{gen_all_betas} provides a way to precompute all beta matrices
#' given the pattern of interventions and WAIFW matrices
#' at each time step for the SC-COSMO model. 
#' 
#' @param l_interventions List of interventions and timing.
#' @param l_contact_info List of WAIFW matrix inputs.
#' @param v_names_ages Vector of names for rows and columns of WAIFW matrices.
#' @param max_time Number of time steps to precompute WAIFW matrices for.
#' @param r_beta Beta component of force of infection w/o intervention.
#' 
#' @return 
#' A 3D matrix of appropriate beta matrices for each time step in the model.
#' 
#' @export
gen_all_betas <- function(l_interventions, l_contact_info, 
                          v_names_ages, max_time, r_beta) {

  n_ages <- length(v_names_ages)
  n_t    <- max_time
  m_all_betas <- array(0, 
                       dim = c(n_ages, n_ages, n_t+1), 
                       dimnames = list(v_names_ages, v_names_ages, seq(0, n_t)))
  
  l_intervention_effects_all_internal <- gen_all_intervention_effects(l_interventions = l_interventions,
                                                                      v_names_ages = v_names_ages,
                                                                      max_time = n_t)
  l_waifws_all_internal <- gen_all_waifws(l_interventions = l_interventions,
                                          l_contact_info = l_contact_info,
                                          v_names_ages = v_names_ages,
                                          max_time = n_t)

  m_intervention_effects_all_internal    <- l_intervention_effects_all_internal$m_all_intervention_effects
  m_intervention_effects_school_internal <- l_intervention_effects_all_internal$m_school_intervention_effects
  m_waifws_all_internal                  <- l_waifws_all_internal$m_all_waifws   
  m_waifws_school_internal               <- l_waifws_all_internal$m_school_waifws
  
  ### Transmission rate by age
  v_beta <- rep(r_beta, n_ages)   
  for(time_idx in 1:(n_t+1)) {
    m_all_betas[, , time_idx] <- (
                                  (
                                    as.numeric(v_beta * m_intervention_effects_all_internal[, time_idx]) *
                                    (m_waifws_all_internal[, , time_idx])
                                  ) +
                                  (
                                    as.numeric(v_beta * m_intervention_effects_school_internal[, time_idx]) *
                                    (m_waifws_school_internal[, , time_idx])
                                  )
                               )
  }
  
  return(m_all_betas)
}