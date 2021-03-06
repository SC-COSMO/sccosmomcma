% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_intervention_construction_functions.R
\name{add_intervention}
\alias{add_intervention}
\title{Building list of interventions to be simulated in
the Stanford-CIDE COronavirus Simulation MOdel (SC-COSMO)}
\usage{
add_intervention(interventions = NULL, intervention)
}
\arguments{
\item{interventions}{List of the current interventions to append to.
If NULL, we construct a new intervention list and then append to it. For
first add to list should be called omitting (or setting to NULL) this
param.}

\item{intervention}{Intervention "object" to add to intervention list.}
}
\value{
A list object containing interventions including the one just added.
}
\description{
\code{add_intervention} constructs a list to hold the
intervention "objects". By convention, it is best to add them in order
of their start times.
}
