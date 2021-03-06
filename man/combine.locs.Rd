% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_locs.R
\name{combine.locs}
\alias{combine.locs}
\title{Re-merge quantile map projections}
\usage{
combine.locs(defaults, comments = character(), output.type = "linear",
  lat_vec = numeric(), lon_vec = numeric())
}
\arguments{
\item{defaults}{the output from \code{\link{various_defaults}}}

\item{comments}{if desired, a 'comments' global attribute will be
added to the NetCDF file, with the string input to this
option as the body}

\item{output.type}{how the time series is exported into the netcdf file; either 
    "linear" (default), as an [idx x time] array, or "grid," as a 
[lon x lat x time] array}

\item{lat_vec, lon_vec}{if \code{output.type=="grid"}, this allows the specification
    of the output grid (useful if the projected data doesn't have values for 
the whole global/etc. grid, or if you just want a subset)}
}
\description{
Concatenate the individual projection files calculated through separate calls
to \code{build.projections} into a single file, removing duplicates along the
way.
}
\section{Output Files}{

The file is saved as both an .RData file and a NetCDF file. The final output
in the NetCDF file is a \code{[loc x time]} array (variable name
\code{[defaults$filevar]}), with variables \code{lat} and \code{lon} being
\code{[loc x 1]} vectors giving the location. The file also includes a variable
\code{base_idx} giving, for each projected day, the linear index in the
original base data that was used as the base day for the projection.
}

\section{Beware}{

Files are searched for in the directory code{[base.data.dir]/output}, with
\code{[base.data.dir]} set by \code{defaults}. As a result, this code depends
on the desired temporary files being the only files in that directory with
the filename pattern searched for (from the \code{file.list} variable in the
code).
}

