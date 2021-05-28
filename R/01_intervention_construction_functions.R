#' Constructing intervention descriptions for each intervention in the
#' the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)
#'
#' \code{make_intervention} provides a generic constructor to build
#' intervention "objects".
#' 
#' @param intervention_type one of the currently implemented intervention
#' types in the SC-COSMO model: 1) "StatusQuo"; 2) "SocialDistancing";
#' 3) "SocialDistancingLinear"
#' @param time_start days since model start when intervention begins
#' @param time_stop days since model start when intervention ends
#' @param ... named values specific to a given intervention type
#' 
#' @return 
#' a list object containing appropriately named values for intervention 
#' objects
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
#' @param interventions a list of the current interventions to append to.
#' If NULL, we construct a new intervention list and then append to it. For
#' first add to list should be called omitting (or setting to NULL) this
#' param
#' @param intervention an intervention "object" to add to intervention list
#' 
#' @return 
#' a list object containing interventions including the one just added
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

### SOME CODE TO TEST THESE FUNCTIONS
# i1 <- make_intervention(intervention_type = "StatusQuo",
#                   time_start = 0,
#                   time_stop = 10)
# i2 <- make_intervention(intervention_type = "SocialDistancing",
#                   time_start = 10,
#                   time_stop = 30,
#                   intervention_factor = 0.9,
#                   intervention_rate = 1.2)
# i3 <- make_intervention(intervention_type = "SocialDistancing",
#                   time_start = 30,
#                   time_stop = 90,
#                   intervention_factor = 0.6,
#                   intervention_rate = 0.9)
# i4 <- make_intervention(intervention_type = "StatusQuo",
#                   time_start = 90,
#                   time_stop = 1000)
# 
# #interventions <- add_intervention(interventions = NULL, intervention = i1)
# interventions <- add_intervention(, intervention = i1)
# interventions <- add_intervention(interventions = interventions, intervention = i2)
# interventions <- add_intervention(interventions = interventions, intervention = i3)
# interventions <- add_intervention(interventions = interventions, intervention = i4)
# 
# print(i1)
# print(i2)
# print(i3)
# print(i4)
# print(interventions)
# length(interventions)