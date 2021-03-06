% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPfunctions.R
\name{traceProbMP}
\alias{traceProbMP}
\title{Traceline function for an item from a fitted MP model.}
\usage{
traceProbMP(x, j, theta = seq(-5, 5, 0.1))
}
\arguments{
\item{x}{Fitted mxModel (e.g., from \code{\link{fitMP}}.}

\item{j}{Item index.}

\item{theta}{Grid for theta, the latent trait.}
}
\description{
Traceline function for an item from a fitted MP model.
}
\details{
Extracts the given item from the fitted model, and then computes values of the traceline along the grid for theta.
  Each row of the output corresponds to a grid point, and each column to a response option category (from smallest to largest category).
}
