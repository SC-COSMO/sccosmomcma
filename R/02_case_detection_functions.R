#' Case detection parameters
#'
#' \code{get_case_detection} pulls case detection parameters for a given time.
#' 
#' @param time Time (numeric, in days) at which case detection rate is evaluated.
#' @param l_params_all List with all parameters of decision model.
#' Case detection rates should be stored as a list within 
#' \code{l_params_all$l_case_detection}.
#' @param severity Severity level and state (exp1, exp2, inf1, inf2).
#' @return 
#' Case detection rate for a given severity (1, 2) or state (exp, inf) at a
#' given point in time.
#' @export
get_case_detection <- function(time, l_params_all, severity) {
#with(as.list(l_params_all), {
    
    # Logic flow to confirm that the case detection parameter exists in l_params_all  
    if (severity %in% c("v_nu_exp1_dx", "v_nu_exp2_dx", "v_nu_inf1_dx", "v_nu_inf2_dx") &
        !(is.null(l_params_all[[severity]]))) {
        
      case_detection_param <- l_params_all[[severity]](time)
      
        ## Pull case detection for a given severity level, indexed by time
      # if (severity %in% c( "v_nu_exp2_dx", "v_nu_inf2_dx")) {
      #   case_detection_param <- l_params_all$jgf_magic(time)
      # } else {
      #   case_detection_param <- l_params_all[[severity]][[floor(time)+1]] 
      # }
      
          # general_logit(logit_lb = l_params_all$l_nu_exp2_dx[[1]]$val_start,
          #             logit_ub = l_params_all$l_nu_exp2_dx[[1]]$val_end,
          #             logit_change_rate = l_params_all$l_nu_exp2_dx[[1]]$v_logit_change_rate,
          #             logit_change_mid = l_params_all$l_nu_exp2_dx[[1]]$v_logit_change_mid,
          #             logit_t = time)
        
    } 
    else {
      stop("Parameter ", severity, " does not exist in l_params_all")
    }
      
    return(case_detection_param)
#  }
#  )
}