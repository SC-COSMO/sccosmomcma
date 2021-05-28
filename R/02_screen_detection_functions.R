#' The Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{get_screen_detection} pulls screen detection parameters for a given time
#' @param time the time (numeric, in days) at which screen detection rate is evaluated
#' @param l_params_all List with all parameters of decision model.
#' Screen detection rates should be stored as a list within l_params_all$l_screen_detection
#' @param severity defines the severity level and state (exp1, exp2, inf1, inf2)
#' @return 
#' Screen detection rate for a given severity (1, 2) or state (exp, inf) at a
#' given point in time.
#' @export

get_screen_detection <- function(time, l_params_all, severity) {
#  with(as.list(l_params_all), {
    
    ## TO DO USE v_screen_detection
    
    # Logic flow to confirm that the screen detection parameter exists in l_params_all  
    if (severity %in%c("v_phi_exp1_dx", "v_phi_exp2_dx", "v_phi_inf1_dx", "v_phi_inf2_dx") &
        !(is.null(l_params_all[[severity]]))) {
      
      ## Pull case detection for a given severity level, indexed by time
      screen_detection_param <- l_params_all[[severity]](time)
      
      #screen_detection_param <- l_params_all[[severity]][[floor(time)+1]]
    }
    else {
      stop("Parameter ", severity, " does not exist in l_params_all")
    }
    
    return(screen_detection_param)
#  }
#  )
}