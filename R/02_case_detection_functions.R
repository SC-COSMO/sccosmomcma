#' The Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{get_case_detection} pulls case detection parameters for a given time
#' @param time the time (numeric, in days) at which case detection rate is evaluated
#' @param l_params_all List with all parameters of decision model.
#' Case detection rates should be stored as a list within l_params_all$l_case_detection
#' @param severity defines the severity level and state (exp1, exp2, inf1, inf2)
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



## Testing graveyard

# v_case_detection_r_nu_exp1_dx <- general_logit(0.1, 0.05, .2, 40, 0:100)
# v_case_detection_r_nu_exp2_dx <- general_logit(0.2, 0.25, .2, 40, 101:200)
# 
# l_params_all$l_case_detection[["r_nu_exp1_dx"]][["values"]] <- v_case_detection_r_nu_exp1_dx
# l_params_all$l_case_detection[["r_nu_exp2_dx"]][["values"]] <- v_case_detection_r_nu_exp1_dx
# 
# 
# cdr1 <- make_period(functional_form = "general_logit", time_start = 0, time_stop = 100, val_start = .1, val_end = .05,
#                             v_logit_change_rate = .2, v_logit_change_mid = 40)
# cdr2 <- make_period(functional_form = "general_logit", time_start = 100, time_stop = 200, val_start = .05, val_end = .2,
#                             v_logit_change_rate = .2, v_logit_change_mid = 150)
# cdr3 <- make_period(functional_form = "linear", time_start = 200, time_stop = 300, val_start = .2, val_end = .1)
# cdr4 <- make_period(functional_form = "constant", time_start = 300, time_stop = 400, val_start = .1, val_end = .1)
# 
# l_case_detection <- add_period(l_period_def = NULL, l_period_add = cdr1)
# l_case_detection <- add_period(l_period_def = l_case_detection, l_period_add = cdr2)
# l_case_detection <- add_period(l_period_def = l_case_detection, l_period_add = cdr3)
# l_case_detection <- add_period(l_period_def = l_case_detection, l_period_add = cdr4)

# test <- approxfun(x = 1:length(gen_time_varying(l_period_def = l_case_detection, max_time = 450)), y = gen_time_varying(l_period_def = l_case_detection, max_time = 450))
# 
# plot(gen_time_varying(l_period_def = l_case_detection, max_time = 450))
