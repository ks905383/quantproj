% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_quantile_surfaces.R
\name{make.quantile.surfaces}
\alias{make.quantile.surfaces}
\title{Construct quantile estimates}
\usage{
make.quantile.surfaces(params, bulk.x, tail.x, norm.x)
}
\arguments{
\item{params}{the list output from \code{\link{get.quantiles}}
giving the fit quantiles and their corresponding
coefficients}

\item{bulk.x, tail.x, norm.x}{the basis functions used to
estimate the coefficients; i.e. loaded from file
or recalculated using \code{\link{get.predictors}}}
}
\value{
a \code{[365 x num_years x num_quantiles]} numeric array
       giving the estimated quantiles (\code{[day of year,year,quantile]})
}
\description{
\code{make.quantile.surfaces} takes the quantile regression
parameters estimated by \code{\link{get.quantiles}} and
constructs the resultant quantile estimates using the basis
functions used as predictors in the quantile regression.
}
