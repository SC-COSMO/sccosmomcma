#' The Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{get_cfr} pulls age-specific fraction of detected cases that die at a given time
#' @param time the time (numeric, in days) at which v_cfr is evaluated
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' Proportion (age-specific) of detected cases that die at a given point in time.
#' @export
get_cfr <- function(time, l_params_all) {
  #with(as.list(l_params_all), {
    
    v_cfr_time <- as.vector(unlist(lapply(l_params_all$m_cfr, do.call, list(time))))
    #v_cfr_time <- as.vector(l_params_all$m_cfr[,floor(time)+1])
    
    return(v_cfr_time)
    
#  }
#  )
}

#' The Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{get_ifr} pulls age-specific fraction of undetected cases that die at a given time
#' @param time the time (numeric, in days) at which v_ifr is evaluated
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' Proportion (age-specific) of undetected cases that die at a given point in time.
#' @export
get_ifr <- function(time, l_params_all) {
#  with(as.list(l_params_all), {
    
    v_ifr_time <- as.vector(unlist(lapply(l_params_all$m_ifr, do.call, list(time))))
  
    #v_ifr_time <- as.vector(l_params_all$m_ifr[,floor(time)+1])
    
    return(v_ifr_time)
    
#  }
#  )
}