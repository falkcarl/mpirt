% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPfunctions.R
\name{EAPscoresMP}
\alias{EAPscoresMP}
\title{Compute EAP scores based on a fitted MP model.}
\usage{
EAPscoresMP(x, items, dat)
}
\arguments{
\item{x}{Fitted mxModel (e.g., from \code{\link{fitMP}}.}

\item{items}{Vector of indices corresponding to items to include in EAP computations.}

\item{dat}{Raw data to use in computing EAP scores.}
}
\description{
Compute EAP scores based on a fitted MP model.
}
\details{
Computes EAP score for the respondents in \code{dat} based on the fitted item response model. Currently assumes standard normal prior. An argument could be added to change this as the function calls
  \code{\link[rpf]{EAPscores}}. It is possible to only compute EAP scores based on a subset of items.
}
