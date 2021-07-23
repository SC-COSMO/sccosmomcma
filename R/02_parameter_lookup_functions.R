#' Get time-varying parameters for the Stanford-CIDE COronavirus Simulation 
#' MOdel (SC-COSMO)
#'
#' \code{get_parameters} wraps functions that pull time-varying parameters.
#' 
#' @param param Single character string defining which parameter to pull.
#' Currently available parameters: r_nu_exp1_dx, r_nu_exp2_dx, r_nu_inf1_dx, 
#' r_nu_inf2_dx, r_phi_exp1_dx, r_phi_exp2_dx, r_phi_inf1_dx, r_phi_inf2_dx,
#' m_Beta, v_cfr, v_p_hosp.
#' @param time Time (numeric, in days) since model start.
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' Parameter values for the parameter specified in "param", either a number, 
#' vector, or matrix.
#' @export
get_parameters <- function(param, time, l_params_all) {

  ## Nu
  if (param %in% c("v_nu_exp1_dx", "v_nu_exp2_dx", "v_nu_inf1_dx", "v_nu_inf2_dx")) {
    param_val <- get_case_detection(time, l_params_all, severity = param)
  }
  
  ## Phi
  if (param %in% c("v_phi_exp1_dx", "v_phi_exp2_dx", "v_phi_inf1_dx", "v_phi_inf2_dx")) {
    param_val <- get_screen_detection(time, l_params_all, severity = param)
  }
  
  ## Beta
  if (param == "m_Beta") {
    param_val <- get_betas(time, l_params_all)
  }
  
  ## Case Fatality
  if (param == "v_cfr") {
    param_val <- get_cfr(time, l_params_all)
  }
  
  ## Reduction factor on transmission from detected infectious (IDX)
  if (param == "v_idx_scale_factor") {
    param_val <- get_idx_scale_factor(time, l_params_all)
  }
  
  ## Proportion Hospitalized
  if (param == "v_p_hosp") {
    param_val <- get_prop_hosp(time, l_params_all)
  }
  
  return(param_val)
}


