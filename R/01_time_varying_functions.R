#' Generalized logistic function
#' 
#' \code{general_logit} computes a generalized logistic function whose starting (lb) 
#' and ending values (ub) are determined by parameters. Rate of change from lb to ub
#' and timing of change are also controlled by parameters.
#' @param logit_lb Lower asymptote.
#' @param logit_ub Upper asymptote.
#' @param logit_change_rate Rate of change from lower to upper asymptote.
#' @param logit_change_mid Point in time (t) where function has gone halfway 
#' from lb to ub.
#' @param logit_t Time point(s) at which function should be evaluated.
#' @return 
#' Value(s) of generalized logit function | inputs.
#' @export
general_logit <- function(logit_lb, 
                          logit_ub, 
                          logit_change_rate, 
                          logit_change_mid, 
                          logit_t) {
  
  if (logit_lb<0 | 
      logit_ub<0 
  ) {
    stop("Lower and upperbound parameters should be non-negative")
  }
  simple_logit <- 1/(1+exp(-1*logit_change_rate*(logit_t - logit_change_mid)))
  gen_logit <- logit_lb+(logit_ub-logit_lb)*simple_logit
  
  return(gen_logit)
}

#' Generalized logistic function (matrix version)
#' 
#' \code{m_general_logit} computes a generalized logistic function whose 
#' starting (lb) and ending values (ub) are determined by parameters. 
#' Rate of change from lb to ub and timing of change are also controlled by 
#' parameters.
#' @param v_logit_lb Age-specific lower asymptote.
#' @param v_logit_ub Age-specific upper asymptote.
#' @param v_logit_change_rate Age-specific rate of change from lower to upper 
#' asymptote.
#' @param v_logit_change_mid Age-specific point in time (t) where function has 
#' gone halfway from lb to ub.
#' @param logit_t time point(s) at which function should be evaluated.
#' @return 
#' Matrix with value(s) of generalized logit function | inputs.
#' @export
m_general_logit <- function(v_logit_lb, 
                            v_logit_ub, 
                            v_logit_change_rate, 
                            v_logit_change_mid, 
                            logit_t) {
  
  if (sum(v_logit_lb < 0) > 0 | 
      sum(v_logit_ub < 0) > 0 
  ) {
    stop("Lower and upperbound parameters should be non-negative")
  }
  if(length(v_logit_lb) != length(v_logit_ub) |
     length(v_logit_ub) != length(v_logit_change_rate) | 
     length(v_logit_change_rate) != length(v_logit_change_mid)){
    stop("All vectors should be of the same length")
  }
  n_ages <- length(v_logit_change_mid)
  m_logit_t <- matrix(rep(logit_t, n_ages), 
                      nrow = n_ages, byrow = T)
  m_simple_logit <- 1/(1+exp(-1*v_logit_change_rate*(m_logit_t - v_logit_change_mid)))
  m_gen_logit <- v_logit_lb + (v_logit_ub - v_logit_lb)*m_simple_logit
  
  return(m_gen_logit)
}



#' Helper function to define a time chunk for a time varying parameter
#' 
#' \code{make_period} defines a time chunk for a generic time varying parameter.
#' Includes some required arguments, but also allows additional arguments.
#' @param functional_form Character string of "linear", "constant", or 
#' "general_logit".
#' @param time_start Numeric time period start (in days).
#' @param time_stop Numeric time period end (in days).
#' @param val_start Starting value for time varying parameter at \code{time_start}.
#' @param val_end Ending value for time varying parameter at \code{time_end}.
#' @param ... Further arguments to be passed to.
#' @return 
#' A list with parameters initialized.
#' @export
make_period <- function(functional_form, 
                        time_start, 
                        time_stop, 
                        val_start, 
                        val_end, 
                        ...) {
  
  l_input_list <- list(...)
  
  l_period_def <- list(type  = functional_form, 
                       start = time_start,
                       stop  = time_stop,
                       val_start = val_start,
                       val_end = val_end
  )
  
  # Ensure that general logit midpoint is within the start and end times
  # if (functional_form %in% c("General Logit", "general_logit", "General_Logit", "general logit",
  #                 "Generalized Logit", "Generalized_Logit", "generalized_logit", "generalized logit",
  #                 "General Logistic", "general_logisitc", "General_Logistic", "general_logistic",
  #                 "Generalized Logistic", "Generalized_Logistic", "generalized_logistic", "generalizedlogistic")) {
  #   if (v_logit_change_mid > time_stop | v_logit_change_mid < time_start) {
  #     stop("Midpoint provided is outside of the time range specified")
  #   }
  # }
  
  l_period_def <- append(l_period_def,
                         l_input_list)
  
  return(l_period_def)
}


#' Helper function add periods together into a list containing all periods
#' 
#' \code{add_period} adds time periods as new entries in list of time periods
#' @param l_period_def Existing list of period definitions.
#' @param l_period_add Period to add to the full list, should be defined
#' using \code{make_period}.
#' @return 
#' List (\code{l_period_def}) appended with \code{l_period_add} as a new entry.
#' @export
add_period <- function(l_period_def = NULL, l_period_add) {
  if (is.null(l_period_def)) {
    l_period_def <- list()
  }
  
  current_length = length(l_period_def)
  
  ## Check for overlap
  ## Really we want to have start == end of previous
  if (current_length > 0) {
    if (l_period_def[[current_length]][["stop"]] != l_period_add[["start"]]) {
      stop("Stop date of previous period should match start date of current period.")
    }
    # if (l_period_def[[current_length]][["val_end"]] != l_period_add[["val_start"]]) {
    #   stop("End val of previous period should match the start val of current period.")
    # }
  }
  
  l_period_def[[current_length + 1]] <- l_period_add
  return(l_period_def)
}


#' Helper function add lists of periods into single containing list which for 
#' each age group has values for all periods
#' 
#' \code{add_list_periods} adds lists of time periods as new entries in list of 
#' age group lists.
#' @param l_l_period_def Existing list of lists.
#' @param l_period_add List of periods to add to the containing list, should 
#' be defined using add_period calls.
#' @return 
#' List (\code{l_l_period_def}) appended with \code{l_period_def} as a new entry.
#' @export
add_list_periods <- function(l_l_period_def = NULL, l_period_def) {
  if (is.null(l_l_period_def)) {
    l_l_period_def <- list()
  }
  
  current_length = length(l_l_period_def)
  l_l_period_def[[current_length + 1]] <- l_period_def
  return(l_l_period_def)
}

#' Helper convenience function that defines age-specific lists of periods for values that 
#' will remain constant over the duration of the simulation.
#' 
#' \code{get_const_multiage_list} adds lists of time periods as new entries in 
#' list of age group lists.
#' @param time_stop Time period end (numeric, in days) (should be \code{n_t}).
#' @param v_ageval Vector of constant values for each each group.
#' @return 
#' List (\code{l_l_period}) list of lists of periods (1 per age group) 
#' with constant values.
#' @export
get_const_multiage_list <- function(time_stop, v_ageval) {
  n_ages <- length(v_ageval)
  
  l_l_period <- NULL
  for(curr_age in 1:n_ages) {
    curr_period <- make_period(functional_form = "constant", 
                               time_start      = 0, 
                               time_stop       = time_stop,
                               val_start       = v_ageval[curr_age],
                               val_end         = v_ageval[curr_age]
                              )
    curr_l_period <- add_period(NULL, curr_period)
    
    l_l_period <- add_list_periods(l_l_period, curr_l_period)
  }
  
  return(l_l_period)
}


#' Helper convenience function that defines age-specific lists of periods for values that 
#' change over the duration of the simulation.
#' 
#' \code{get_non_const_multiage_list} adds lists of time periods as new entries 
#' in list of age group lists.
#' @param v_time_stop Vector with time period at which values change (in days).
#' Last entry should be \code{n_t}.
#' @param m_ageval Matrix of values for each age group (\code{n_ages}) at each
#' time period.
#' @return 
#' List (\code{l_l_period}) list of lists of periods (1 per age group) with 
#' non-constant values.
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

#' Helper function for creating time-varying vector of values
#' 
#' \code{gen_time_varying} creates a time-varying vector of values based on 
#' entries in \code{l_period_def}.
#' @param l_period_def Existing list of period definitions, created using 
#' \code{add_period}.
#' @param max_time Number of days that should be estimated, will fill if not 
#' fully covered by final stop in \code{l_period_def}.
#' @return 
#' Vector of time-varying parameter values at minimum of length max_time 
#' (could be longer).
#' @export
gen_time_varying <- function(l_period_def, max_time) {
  
  ## Loop through the periods that have been defined and generate the time-varying parameter vector
  tvp_out <- NULL
  for (i in 1:length(l_period_def)) {
    
    ## Linear / Constant
    if (l_period_def[[i]]$type %in% c("Linear", "linear", "Constant", "constant")) {
      temp_out <- approx(x = c(l_period_def[[i]]$start, l_period_def[[i]]$stop), y = c(l_period_def[[i]]$val_start, l_period_def[[i]]$val_end),
                         xout = seq(l_period_def[[i]]$start, l_period_def[[i]]$stop, 1), method = "linear")$y
      if (is.null(tvp_out)) {
        tvp_out <- temp_out
      } else {
        tvp_out <- c(tvp_out, temp_out[-1]) # Since start and end overlap, drop the first entry when binding
      }
    }
    
    ## General Logit
    if (l_period_def[[i]]$type %in% c("General Logit", "general_logit", "General_Logit", "general logit",
                                      "Generalized Logit", "Generalized_Logit", "generalized_logit", "generalized logit",
                                      "General Logistic", "general_logisitc", "General_Logistic", "general logistic",
                                      "Generalized Logistic", "Generalized_Logistic", "generalized_logistic", "generalizedlogistic")) {
      temp_out <- general_logit(logit_lb = l_period_def[[i]]$val_start,
                                logit_ub = l_period_def[[i]]$val_end,
                                logit_change_rate = l_period_def[[i]]$v_logit_change_rate,
                                logit_change_mid = l_period_def[[i]]$v_logit_change_mid,
                                logit_t = l_period_def[[i]]$start:l_period_def[[i]]$stop)
      if (is.null(tvp_out)) {
        tvp_out <- temp_out
      } else {
        tvp_out <- c(tvp_out, temp_out[-1]) # Since start and end overlap, drop the first entry when binding
      }
    }
  }
  # Make sure that the output is at least as long as the time requested by max_out (copy forward last value (linear))
  if (length(tvp_out) < max_time) {
    warning("Periods provided are shorter than max_time requested, copying forward last value to fill")
    tvp_out <- c(tvp_out, rep.int(x = last(tvp_out), times = max_time - length(tvp_out)))
  }
  
  tvp_out <- approxfun(x = 0:(max_time-1), y = tvp_out, rule = 2)
  
  return(tvp_out)
}