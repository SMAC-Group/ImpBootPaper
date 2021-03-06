% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confidence_interval.R
\name{parboot_interval}
\alias{parboot_interval}
\title{Percentile parametric bootstrap interval}
\usage{
parboot_interval(
  y,
  x = NULL,
  initial,
  boot = NULL,
  model,
  alpha = c(0.025, 0.975),
  control = boot.control()
)
}
\arguments{
\item{y}{observations}

\item{x}{matrix of design}

\item{initial}{MLE on observations}

\item{boot}{a \code{matrix} of parametric bootstrap estimates,
if \code{NULL} it is estimated internally}

\item{model}{either "\code{lomax}" for Lomax distribution,
"\code{treg}" for t regression, "\code{betareg}" for beta regression}

\item{alpha}{levels}

\item{control}{see \code{\link{boot.control}}}
}
\value{
a \code{matrix} levels in rows and MLE in columns
}
\description{
Percentile parametric bootstrap interval
}
