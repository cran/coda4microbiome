% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coda4microbiome_data.R
\docType{data}
\name{Event_time}
\alias{Event_time}
\title{data_survival}
\format{
The dataset contains three objects:
\describe{
\item{x}{microbiome abundance matrix for 150 individuals (rows) and 48 genera (columns)}
\item{Event}{a \code{numeric}, event occurrence. Vector (type: numeric or logical) specifying 0 or FALSE for no event occurrence, and 1 or TRUE for event occurrence.}
\item{Event_time}{a \code{numeric}, time to event or follow up time for right censored data. Must be a vector (type:numeric) specifying time to event for each sample for right censored data.}
}
}
\description{
Survival Data simulated from the Crohn's disease original study: 48 taxa and 150 individuals
}
\references{
\doi{10.1016/j.chom.2014.02.005}
}
\keyword{data}
