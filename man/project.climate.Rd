% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_climate.R
\name{project.climate}
\alias{project.climate}
\title{Project climate distribution changes with resampling}
\usage{
project.climate(defaults, params, base.data, output.years, norm.x = numeric(),
  bulk.x = numeric(), tail.x = numeric(), norm.x.base = numeric(),
  output = "full", index.type = "resampling",
  resampling.timescale = "day", rand.seed.set = 42)
}
\arguments{
\item{defaults}{the output from \code{\link{various_defaults}};
for \code{base.norm.x.df}, \code{aux.dir}, and \code{base.year.range}}

\item{params}{the calculated quantile fit parameters of the
projecting model, taken directly from the output of
\code{\link{get.quantiles}}}

\item{base.data}{the time series to be projected, MUST BE
AN "XTS" OBJECT (this is to reduce the number of
parameters that have to input into this function -
the xts object carries the time information as well).}

\item{output.years}{the desired time frame (in years) for the
output, projected climate time series. THIS MUST
INCLUDE EVERY DESIRED FUTURE YEAR - not just the range -
\code{output.years=c(2011,2099)} will result in an
output time series, with two years, \code{2011} and
\code{2099}.}

\item{norm.x, bulk.x, tail.x, norm.x.base}{For each basis function,
if empty, the code first attempts to load it from file. If
the file is missing, the basis function is calculated through
\code{\link{get.predictors}}. Degrees of freedom (and whether
or not to include the volcanic CO2 forcing) for the basis
functions are taken from the \code{norm.x.df}, \code{bulk.x.df},
\code{tail.x.df}, and \code{get.volc} parameters in the
\code{params} input (to ensure the quantile surfaces are
constructed using the same inputs as they were estimated);
\code{base.norm.x.df} is taken from \code{defaults}.}

\item{output}{set the output options as explained the "Output
options" section below.}

\item{index.type, resampling.timescale}{set the temporal resampling
options as explained in the "Temporal resampling" section
below.}

\item{rand.seed.set}{set random seed for resampling (default=\code{42})}
}
\value{
see "Output options" section above.
}
\description{
\code{project.climate} projects a climate timeseries by the
distributional changes estimated by \code{\link{get.quantiles}},
with several temporal projection options, including
resampling an existing time period by day or year.
}
\section{Projection}{

 The \code{base.data} time series is projected by the
 distributional changes from a model, estimated through
 \code{\link{get.quantiles}}, using the following p
 \enumerate{
     \item The estimated model quantiles are unpacked by
         evaluating the quantile regression model using
         the quantile fit parameters (\code{params}) and
         the bases (either loaded or calculated)
     \item The data to be projected is normalized by
         fitting the normalizing quantiles given by
         \code{defaults$base.norm.x.df} (subtracting
         the second element, and dividing by the
         difference between the third and first element;
         by default, \code{defaults$base.norm.x.df=c(0.1,0.5,0.9)})
     \item The projecting indices are found - for each
         day in the projected time series, which day is
         to be used as a base for comparison (see section
         on "Temporal resampling", below)
     \item Match base period data points in the time series
         to be projected to quantiles in the projecting
         model data
     \item For each projected day, apply the change in the
         matched quantile in the (normalized) projecting
         data time series between the projected day and the
         set base day to the time series to be projected
     \item Un-normalize the now projected time series
 }
}

\section{Temporal resampling}{

 To generate a new, projected, climate time series over
 the time period given by \code{output.years}, the base
 inputted climate time series must be projected in some
 temporal mapping - say, if the base climate time series
 spans 1979 - 2010, and the desired \code{output.years}
 are 2011 - 2099, for each of the future years, a
 corresponding 'base' year must be chosen to compare against.
 \code{project.climate} offers three options for this
 projection:
 \describe{
     \item{resample_rep}{resampling with replacement (the
         default). For \code{resampling.timescale="year"}
         (the default), for each projected year in
         \code{output.years}, a random year is chosen from
         the base input time frame, and projected by the
         change in the model distribution between the
         corresponding output year and that chosen base year.
         For \code{resampling.timescale="day"}, for each
         day of year in \code{output.years}, that same
         day of year is chosen from a random year in the
         base input time frame.}
     \item{raw}{If \code{output.years} spans the same length
         as the base input time frame, it is just wholesale
         projected forward (i.e. 1979-2010 is projected to
         2011-2042 year-by-year). If \code{output.years} is
         shorter, just the first years are chosen up to the
         span of \code{output.years}. If \code{output.years}
         is longer, the projection wraps around - i.e.
         2011-2050 is projected from [1979-2010 1979-1987].}
     \item{raw_mean}{The behvavior is identical to that of the
         \code{raw} option above, but the model quantiles are 
         averaged across the time periods before the change is 
         applied - in other words, each matched base quantile is
         projected by the change in the average quantile between
         the base and projecting time periods.}
     \item{resample}{NOT YET IMPLEMENTED}
 }
}

\section{Output options}{

 \describe{
 \item{full}{(default) a list object giving the projected time
     series as an \code{xts} object (\code{proj.data}), the
     coefficients of the quantiles used for normalizing the
     base climate time series (\code{base.norm.coef}),
     the \code{lat} and \code{lon} of the pixel, the
     \code{output.years}, and the linear indices (in the base
     period) used a the reference year/days in the projected
     period (\code{base.idxs})}
 \item{xts}{an xts object with just the projected time series}
 \item{vec}{a numeric vector with just the projected time series
     (generated from \code{as.numeric([xts])})}
 }
}

