% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hbic.R
\name{hbic}
\alias{hbic}
\title{High Dimensional BIC.}
\usage{
hbic(coefs, tau, x, y, largep)
}
\arguments{
\item{coefs}{matrix of coefficients.}

\item{tau}{scalar between 0 and 1 corresponding to the conditional quantile.}

\item{x}{design matrix used in fitting the model.}

\item{y}{vector containing the values of the outcome used in fitting the model.}
}
\value{
a vector of HBIC values
}
\description{
High Dimensional BIC.
}
