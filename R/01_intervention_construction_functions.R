#' Constructing intervention descriptions for each intervention in the
#' the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{make_intervention} provides a generic constructor to build
#' intervention "objects".
#' 
#' @param intervention_type One of the currently implemented intervention
#' types in the SC-COSMO model: 1) "StatusQuo"; 2) "SocialDistancing";
#' 3) "SocialDistancingLinear".
#' @param time_start Days since model start when intervention begins.
#' @param time_stop Days since model start when intervention ends.
#' @param ... Named values specific to a given intervention type.
#' 
#' @return 
#' A list object containing appropriately named values for intervention 
#' objects.
#' 
#' @export
make_intervention <- function(intervention_type, 
                              time_start, 
                              time_stop,
                              ...) {
  l_input_list <- list(...)

  l_intervention <- list(type  = intervention_type, 
                        start = time_start,
                        stop  = time_stop
                       )

  l_intervention <- append(l_intervention,
                           l_input_list)

  return(l_intervention)

}

#' Building list of interventions to be simulated in
#' the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{add_intervention} constructs a list to hold the
#' intervention "objects". By convention, it is best to add them in order 
#' of their start times.
#' 
#' @param interventions List of the current interventions to append to.
#' If NULL, we construct a new intervention list and then append to it. For
#' first add to list should be called omitting (or setting to NULL) this
#' param.
#' @param intervention Intervention "object" to add to intervention list.
#' 
#' @return 
#' A list object containing interventions including the one just added.
#' 
#' @export
add_intervention <- function(interventions = NULL, intervention) {
  if (is.null(interventions)) {
    interventions <- list()
  }

  current_length = length(interventions)
  interventions[[current_length + 1]] <- intervention
  return(interventions)

}