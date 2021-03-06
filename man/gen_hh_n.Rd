% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_household_model_inputs_functions.R
\name{gen_hh_n}
\alias{gen_hh_n}
\title{Generate household possibilities}
\usage{
gen_hh_n(n_hhsize, v_names_states_hh)
}
\arguments{
\item{n_hhsize}{Household size.}

\item{v_names_states_hh}{Vector with state names of household model.}
}
\value{
A data.frame with possible combinations of household members within each
epidemic compartment.
}
\description{
\code{gen_hh_n} generates household transmission (i.e., force of infection).
}
