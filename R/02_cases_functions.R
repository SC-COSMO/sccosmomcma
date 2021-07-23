#' Get reduction factor on transmission from detected infectious (IDX)
#'
#' \code{get_idx_scale_factor} pulls age-specific reduction factor on 
#' on transmission from detected infectious.
#' @param time Time (numeric, in days) at which \code{v_idx_scale_factor} is 
#' evaluated.
#' @param l_params_all List with all parameters of decision model.
#' @return 
#' Reduction (age-specific) on infectiousness of detected cases at a given point 
#' in time.
#' @export
get_idx_scale_factor <- function(time, l_params_all) {
  
  v_idx_scale_factor_time <- as.vector(unlist(lapply(l_params_all$m_idx_scale_factor, do.call, list(time))))
  
  return(v_idx_scale_factor_time)
}