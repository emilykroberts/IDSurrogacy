% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_traceplots.R
\name{plot_traceplots}
\alias{plot_traceplots}
\title{traceplots}
\usage{
plot_traceplots(params_matrix, variable)
}
\arguments{
\item{params_matrix}{matrix of parameters from mcmc}

\item{variable}{variable of which to plot}
}
\value{
traceplots
}
\description{
creates traceplots
}
\examples{
example(plot_traceplots(params_matrix = params_matrix, variable = "int"))
}
