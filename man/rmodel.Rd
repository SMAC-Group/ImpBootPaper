% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap.R
\name{rmodel}
\alias{rmodel}
\title{Simulation from a parametric model}
\usage{
rmodel(n = NULL, x = NULL, theta, model, seed)
}
\arguments{
\item{n}{number of observations}

\item{x}{matrix of design}

\item{theta}{a list of model parameters (see 'Details')}

\item{model}{either "\code{lomax}" for Lomax distribution,
"\code{treg}" for t regression, "\code{betareg}" for beta regression}

\item{seed}{seed for random number generation}
}
\value{
a \code{matrix} of bootstrap estimates
}
\description{
Simulation from a parametric model
}
\details{
The list of \code{theta} depends on the model. The names of the element of the list
must follow this convention:
\itemize{
   \item \code{lomax}: positive parameters \code{b}
   and \code{q};
   \item \code{treg}: the regression coefficients \code{beta},
   the variance \code{sig2}\eqn{>0} and the
   degrees of freedom \code{nu}\eqn{>0};
   \item \code{betareg}: the regression coefficients \code{beta}
   and the precision parameter (positive) \code{phi}.
}
}
