% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_calibration_functions.R
\name{log_post}
\alias{log_post}
\title{Evaluate log-posterior of calibrated parameters}
\usage{
log_post(v_params, ...)
}
\arguments{
\item{v_params}{Vector (or matrix) of model parameters.}

\item{...}{Further arguments to be passed to.}
}
\value{
A scalar (or vector) with log-posterior values.
}
\description{
\code{log_post} computes a log-posterior value for one (or multiple)
parameter set(s) based on the simulation model, likelihood functions and
prior distributions.
}
