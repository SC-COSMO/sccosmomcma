#' Generate household transmission
#' 
#' \code{gen_household_transmission} generates household transmission (i.e., 
#' force of infection).
#' 
#' @param r_tau Transmission rate in the household model.
#' @param n_hhsize  Household size (integer).
#' @param n_contacts_hh Average number of within household contacts.
#' @param v_HH State vector of within household epidemics.
#' @param v_index_hh_sus Indexing column vector of susceptible.
#' @param v_index_hh_inf Indexing column vector of infectious.
#' @param v_index_hh_rec Indexing column vector of recovered.
#' @param m_possibilities Matrix with possible combinations of household 
#' members within each epidemic compartment.
#' @return 
#' A scalar with the force of infection.
#' @export
gen_household_transmission_mc_seir <- function(r_tau,
                                               n_hhsize,
                                               n_contacts_hh,
                                               v_HH,
                                               v_index_hh_sus,
                                               v_index_hh_exp,
                                               v_index_hh_inf,
                                               v_index_hh_rec,
                                               m_possibilities){
  #### Indexing vectors ####
  ## Indices with active susceptibles and infectious
  v_index_keep_tau <- (m_possibilities[, v_index_hh_sus] > 0 & 
                         rowSums(m_possibilities[, v_index_hh_inf, 
                                                 drop = FALSE]) > 0)
  #### Estimate within HH transmission ####
  ### Vector with combinations that can infect and get infected
  v_hh_size_inf <- m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
    rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop=FALSE])
  
  n_hh_inf_rate <- (r_tau * n_contacts_hh) * 
    (v_hh_size_inf %*% v_HH[v_index_keep_tau])/n_hhsize ### CHECK IN FUTURE
  return(n_hh_inf_rate)
}