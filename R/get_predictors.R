#' Build spline basis functions
#'
#' \code{get.predictors} builds basis functions as predictors for
#' the quantile regression in the quantile mapping climate
#' projection method of \emph{Haugen et al. 2018}. The basis
#' functions are constructed such that the resultant quantile
#' estimate is smoothly varying and depends on:
#' \enumerate{
#'    \item the seasonal cycle (dependent only on month-of-year),
#'          estimated by a periodic b-spline
#'    \item the long-term trend (dependent only on the year index),
#'          estimated by a natural cubic spline
#'    \item an interaction between the seasonal cycle and long-term
#'          trend (to capture changes in the seasonal cycle)
#' }
#' The degrees of freedom for each of these basis functions are set through the
#' parameters \code{df.x}, \code{df.t}, and \code{df.xt}, respectively, or through a \code{3 x 1} vector \code{dfs}.
#'
#' @section Multiple runs/files:
#' \code{get.predictors} is specifically designed to allow for basis
#' functions for quantile regression across multiple data sets that
#' span the same time period. In the context of the quantile mapping
#' climate projection method, these multiple data sets are multiple
#' runs of the same model over the same time period. In this case, the
#' basis vectors are repeated end-to-end \code{n_files} times.
#'
#' @section Volcanic forcing:
#' The \code{get.volc} input additionally allows for the addition
#' of a historical volcanic CO2 forcing as a predictor for the
#' quantile regression. This volcanic forcing is taken from the
#' XXXXXX hindcast from CCSM4XXXX. If \code{get.volc=TRUE}, then a
#' latitude must be specified through the \code{lat} parameter.
#' Since the volcanic forcing is only available 1850-2008, this
#' is discouraged for values beyond 2008...
#'
#' @param dfs degrees of freedom - a \code{[3 x 1]} vector giving the
#'    seasonal, long term change, and interaction degrees of
#'    freedom, respectively
#' @param df.x,df.t,df.xt alternatively, input the dfs separately,
#'     with seasonal \code{df.x}, long-term \code{df.t}, and
#'     interaction \code{df.xt} (kept for backwards compatibility
#'     with some older functions)
#' @param n_files number of runs (if \code{> 1},
#'            this ensures that Jan 1, 1971, Run
#'            1 has the same basis function value
#'            as Jan 1, 1971, Run 2, etc.). Basis
#'            functions are repeated end-to-end
#'            \code{n_files} times.
#' @param year.range \code{c(start_year,end_year)}, or, if
#'          \code{get.volc=FALSE}, can also just be a number giving
#'          the length in years of the used time series.
#' @param get.volc whether to add a volcanic degree of freedom
#' @param lat only necessary if \code{get.volc=TRUE}; This is the
#'         latitude of the pixel for which basis functions
#'         should be calculated.
#' @param save.predictors if true, the basis functions are saved
#'         to an .RData file (def: \code{FALSE})
#' @param save.fn if \code{save.predictors=TRUE}, this sets the name of the
#'         file to save the basis functions to. By default, the
#'         filename is: "\code{[aux.dir]/bases/spline_basis_functions_[nyears]years_[n_files]runs_[df.x]-[df.t]-[df.xt]df(_[lat.volc]volc).RData}"
#'
#' @return a numerical matrix
#' 
#' @importFrom pbs pbs
#' @importFrom splines ns

get.predictors <- function(n_files=1,
  dfs=numeric(),df.x=numeric(),df.t=numeric(),df.xt=numeric(),
  year.range=c(1850,2099), get.volc=FALSE, lat=numeric(),
  save.predictors=FALSE,save.fn=character()) {

  # Extract from dfs
  if (length(df.x)==0&&length(dfs)>0) {
    df.x <- dfs[1]; df.t <- dfs[2]; df.xt <- dfs[3]; rm(dfs)
  }

  # Allow for "length" inputs of year.range
  if (length(year.range)==1) {
    year.range <- c(1,year.range)
  }

  day_of_year <- 1:365
  nyears <- diff(year.range) + 1

  # GET BASIS VECTORS FOR SEASONAL TREND --------------------------------------------
  if (df.x != 0) {
    # Create a basis matrix for *periodic* B-splines for the seasonal (sub-yearly) (this
    # uses [pbs], which creates a b-spline basis for polynomial splines, with df.x-1
    # knots).
    # [x] in this case is constructed by scaling a vector giving the day-of-year index
    # for every day in the time series (so 1:365 == 1:365, 366:720 == 1:365, etc.) by
    # the standard deviation of that vector
    x <- scale(rep(day_of_year, n_files*nyears), center=FALSE)
    x.basis <- as.matrix(pbs(x, df=df.x))
    rm(x)
  }

  # GET BASIS VECTORS FOR YEARLY TREND ----------------------------------------------
  if (df.t != 0) {
    # Create a basis matrix for (not periodic) natural cubic splines for the long term
    # (yearly) temporal trend with df.t degrees of freedom (this uses [ns], which
    # creates a basis for cubic splines, with df.t-1 knots).
    # [t] in this case is constructed by scaling a vector giving the year index for
    # every day in the time series (so 1:365 == 1, 366:720 == 2, etc.) by the standard
    # deviation of that vector
    t <- scale(rep(rep((1:nyears),each=365),n_files),center=FALSE)
    t.basis <- ns(t, df=df.t)
    rm(t)
  }

  # GET BASIS VECTORS FOR INTERACTION ----------------------------------------------
  if (df.xt != 0) {
    # Create a basis matrix for the interaction between the periodic seasonal and non-
    # periodic long-term trend (allowing for changes in the seasonal cycle). (this uses
    # [pbs], which creates a b-spline basis for polynomial splines, with df.xt-1 knots).
    # [xt] in this case is constructed by
    xt <- scale(rep(day_of_year, n_files*nyears), center=FALSE)
    xt.basis <- scale(as.matrix(pbs(xt, df=df.xt)), center=FALSE)
    rm(xt)
  }

  # MAKE DESIGN / MODEL MATRIX (the [X] in Y <- Xa+b) ------------------------------
  if (df.xt != 0) {
    # add interaction, if desired
    X <- model.matrix(~x.basis + t.basis + xt.basis:t.basis)
    rm(list=c("x.basis","t.basis","xt.basis"))
  } else {
    X <- model.matrix(~x.basis + t.basis)
    rm(list=c("x.basis","t.basis"))
  }

  # ADD VOLCANIC VARIABLE IF DESIRED -----------------------------------------------
  if (get.volc) {
    if (length(lat)==0) {
      stop("pixel latitude must be sepcified (through 'lat=#') to determine volcanic forcing.")
    }
    if (min(year.range)<1850 || max(year.range)>2008) {
      warning(paste0("year.range (",paste0(year.range,collapse=", "),") lies outside the bounds of ",
        "the volcanic forcing dataset (1850-2008). Values outside of this range will be interpolated ",
        "using [na.approx]."))
    }

    volc.time <- c(t(matrix(rep(seq(1850, 2099), 365), 250, 365)))
    volc.data <- get.volcanic(lat)
    lat.volc <- volc.data$lat.volc
    volc <- scale(rep(volc.data$volc.daily[volc.time %in% seq(year.range[1], year.range[2])], n_files), center=FALSE)
    # If the forcing is 0 for all years in the timeframe, then
    # [scale] produces only NaNs (dividing 0 by sd(0)), so just
    # don't add the variable in that case (also that means that
    # there's no variation in the volcanic forcing, and therefore
    # no point in including it as a predictor)
    if (!all(is.na(volc))) {
      X <- cbind(X, volc)
    }
    rm(list=c("volc.time","volc","volc.data"))
  }

  # REMOVE UNNEEDED BLOATING THINGS -----------------------------------------------
  attr(X,"dimnames") <- NULL
  attr(X,"assign") <- NULL

  # SAVE IF DESIRED ---------------------------------------------------------------
  # Don't need the "NORM"/"BULK"/etc. Suffix bc all you need is the df in the filename and whether the volcanic thing is added... should change the rest to
  if (save.predictors) {
    if (length(save.fn)==0) {
      save.fn <- paste0(defaults$aux.dir,"bases/","spline_basis_functions_",
                                                          diff(year.range)+1,"years_",n_files,"runs_",
                                                          paste0(c(df.x,df.t,df.xt),collapse="-"),"df")
      if (get.volc) {
        save.fn <- paste0(save.fn,"_volc",gsub("\\.","-",round(lat.volc[which.min(abs(as.vector(lat.volc)-as.vector(lat)))],1)))
      }
      save.fn <- paste0(save.fn,".RData")
    }
    save(list=c("X"),file=save.fn)
    print(paste0(save.fn," saved!"))
  }

  # RETURN ------------------------------------------------------------------------
  invisible(X)
}


# Function to access/interpolate the (monthly) volcanic data
# (THe volcanic dataset only spans 1850-2008, so there's no
# volcanic activity and therefore no predictor past 2008,
# so it probably shouldn't be used)
get.volcanic <- function(mlat) {
  # volc.data is internal to the package

  # Find the latitude in the volcanic dataset closest to the
  # latitude of the inputted pixel (mlat)
  volc.lat.idx <- which.min(abs(as.vector(volc.data$lat) - as.vector(mlat)))

  # Get mean across pressure level dimension (so mean column MMR)
  # The first element is ignored, because the .nc file starts at
  # December 1849.
  mean.volc <- apply(volc.data$volc[volc.lat.idx, ,-1], 2, mean)
  rm(list=c("volc"))

  # Get date sequence of the desired daily output time series
  out.days <- seq.Date(as.Date("1850-01-01"),as.Date("2099-12-31"),"days")
  out.days <- out.days[strftime(out.days,"%j")!=366]
  # Create daily NA time series over the 250 years
  volc.daily <- rep(NA, 365*250)
  # Assign the monthly values from the NetCDF file to the middle of each month
  volc.daily[as.numeric(strftime(out.days,format="%Y%m%d"))%in%volc.data$volc.date] <- mean.volc

  # Interpolate the NAs from the once-a-month values, and set further NAs to 0
  volc.daily <- na.approx(volc.daily, na.rm=FALSE)
  volc.daily[is.na(volc.daily)] <- 0

  # Normalize
  volc.daily <- (volc.daily) / sd(volc.daily)

  # Return the interpolated, daily volcanic time series
  return(list(volc.daily=volc.daily,lat.volc=volc.data$lat[volc.lat.idx]))
}


