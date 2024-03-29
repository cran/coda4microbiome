% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coda4microbiome_functions.R
\name{coda_glmnet_null}
\alias{coda_glmnet_null}
\title{coda_glmnet_null}
\usage{
coda_glmnet_null(
  x,
  y,
  niter = 100,
  covar = NULL,
  lambda = "lambda.1se",
  alpha = 0.9,
  sig = 0.05
)
}
\arguments{
\item{x}{abundance matrix or data frame (rows are samples, columns are variables (taxa))}

\item{y}{outcome (binary or continuous); data type: numeric, character or factor vector}

\item{niter}{number of iterations (default = 100)}

\item{covar}{data frame with covariates (default = NULL)}

\item{lambda}{penalization parameter (default = "lambda.1se")}

\item{alpha}{elastic net parameter (default = 0.9)}

\item{sig}{significance level for the confidence interval (default = 0.05)}
}
\value{
a list with "accuracy" and "confidence interval"
}
\description{
Performs a permutational test for the coda_glmnet() algorithm:
It provides the distribution of results under the null hypothesis by
implementing the coda_glmnet() on different rearrangements of the response variable.
}
\examples{

data(Crohn, package = "coda4microbiome")

coda_glmnet_null(x=x_Crohn[,(1:10)], y=y_Crohn, niter=2,covar=NULL,lambda="lambda.1se",
                                                alpha=0.9,sig=0.05)


}
\author{
M. Calle - T. Susin
}
