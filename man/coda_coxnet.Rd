% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coda_coxnet.R
\name{coda_coxnet}
\alias{coda_coxnet}
\title{coda_coxnet}
\usage{
coda_coxnet(
  x,
  time,
  status,
  covar = NULL,
  lambda = "lambda.1se",
  nvar = NULL,
  alpha = 0.9,
  nfolds = 10,
  showPlots = TRUE,
  coef_threshold = 0
)
}
\arguments{
\item{x}{abundance matrix or data frame (rows are samples, columns are variables (taxa))}

\item{time}{time to event or follow up time for right censored data. Must be a numericvector.}

\item{status}{event occurrence. Vector (type: numeric or logical) specifying 0, or FALSE, for no event occurrence, and 1, or TRUE, for event occurrence.}

\item{covar}{data frame with covariates (default = NULL)}

\item{lambda}{penalization parameter (default = "lambda.1se")}

\item{nvar}{number of variables to use in the glmnet.fit function (default = NULL)}

\item{alpha}{elastic net parameter (default = 0.9)}

\item{nfolds}{number of folds}

\item{showPlots}{if TRUE, shows the plots (default = TRUE)}

\item{coef_threshold}{coefficient threshold, minimum absolute value of the coefficient for a variable to be included in the model (default =0)}
}
\value{
list with "taxa.num","taxa.name","log-contrast coefficients","risk.score","apparent Cindex","mean cv-Cindex","sd cv-Cindex","risk score plot","signature plot".
}
\description{
Microbial signatures in survival studies
The algorithm performs variable selection through an elastic-net penalized Cox regression conveniently adapted to CoDA.
The result is expressed as the (weighted) balance between two groups of taxa.
It allows the use of non-compositional covariates.
}
\examples{

data(data_survival, package = "coda4microbiome")
time <- Event_time
status <- Event
set.seed(12345)
coda_coxnet(x = x,
           time = Event_time,
           status = Event,
           covar = NULL,
           lambda = "lambda.1se", nvar = NULL,
           alpha = 0.9, nfolds = 10, showPlots = TRUE, coef_threshold = 0)



}
\author{
M. Calle, M. Pujolassos, T. Susin
}
