% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/like12_shape.R
\name{like12_shape}
\alias{like12_shape}
\title{likelihood function}
\usage{
like12_shape(shape12_0, scale12_0, dat0, c12_0, omega12_z0)
}
\arguments{
\item{shape12_0}{shape of Weibull distribution}

\item{scale12_0}{scale of Weibull distribution}

\item{dat0}{data}

\item{c12_0}{frailty coefficient}

\item{omega12_z0}{frailty term}
}
\value{
likelihood
}
\description{
evaluate likelihood ...
}
\examples{
example(like12_shape(shape12_0 = 1, scale12_0 = 1, dat0 = 1, 
c12_0 = 1, omega12_z0 = 1))
}
