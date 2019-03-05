#' Estimate quantiles of a normalized time series
#'
#' \code{get.quantiles} performs quantile regressions on the time series
#' \code{model.y} to estimate the quantiles given by \code{q_bulk} and
#' \code{q_tail}, after normalizing using \code{q_norm}. Specifically,
#' \enumerate{
#'    \item the quantiles \code{q_norm} are estimated on \code{model.y}
#'    \item \code{model.y} is normalized, by subtracting the estimated quantile
#'      \code{q_norm[2]} and dividing by the IQR of estimated quantiles
#'      \code{q_norm[3]-q_norm[1]}
#'    \item the quantiles \code{q_bulk} are estimated on the normalized \code{model.y}
#'    \item the high tail and low tail are calculated by subtracting the max and min
#'      estimated quantile (from 3.) from the normalized \code{model.y}
#'    \item the quantiles \code{q_tail} are estimated from the high tail and low tail
#'      exceedences (so from \code{exceedence[exceedence>0]} for the high and low)
#' }
#'
#' @section Notes on inputted quantiles to be estimated:
#' \describe{
#'   \item{q_norm}{quantiles used for normalization. \code{model.y} is normalized by
#'           subtracting the estimated \code{q_norm[2]} and dividing by the
#'           estimated IQR \code{q_norm[3]-q_norm[1]}. We therefore suggest to keep
#'           \code{q_norm[2] == 0.5} (the median).}
#'   \item{q_bulk}{primary, non-tail quantiles to estimate (post-normalization)}
#'   \item{q_tail}{tail quantiles to estimate; based on exceedences of the post-
#'           normalized \code{model.y} beyond the max and min estimated \code{q_bulk}.
#'           Note that \code{q_tail} is based on the exceedence, so if you set
#'           \code{q_tail[1]} to \code{0.75}, that 'real' value of that quantile is actually
#'           \code{max(q_bulk)+(1-max(q_bulk))*0.75}, and so forth.}
#'  }
#'
#' @section Notes on estimation:
#'   Quantiles are estimated using quantile regression, with cubic spline basis
#'   functions, with degrees of freedom set by \code{norm.x.df},
#'   \code{bulk.x.df}, and \code{tail.x.df}, for the normalization, bulk, and
#'   tail exceedence quantile calculations, respectively. The basis functions
#'   are either loaded (if they exist in the \code{bases.dir}, this requires
#'   bases to be calculated using \code{\link{get.predictors}}), calculated from
#'   scratch using \code{\link{get.predictors}}, or directly inputted using the
#'   function parameters \code{norm.x}, \code{bulk.x}, and \code{tail.x}.
#'   Loading pre- calculated bases is generally the fastest method.
#'
#'   Each of \code{norm.x.df}, \code{bulk.x.df}, and \code{tail.x.df} is a
#'   \code{3 x 1} vector giving the degrees of freedom for the seasonal cycle
#'   (based only on the month-of-year), the long-term change (based only on the
#'   year), and the interaction (changing seasonal cycle, based on both
#'   month-of-year and year), respectively.
#'
#' @param model.y the temperature time series, as either a numeric vector or an
#'   \code{xts} object. If multiple runs are analyzed, they should be end-to-end
#'   (so [run 1 1979-2099, run 2 1979-2099, etc...]). The number of runs is
#'   inferred from the length of the time series and the \code{year.range}
#'   parameter below.
#' @param norm.x.df,bulk.x.df,tail.x.df degrees of freedom for the basis
#'   functions, as \code{3 x 1} vectors giving the dfs for the seasonal cycle,
#'   long-term trend, and interaction, respectively
#' @param q_norm a \code{3 x 1} vector giving the quantiles used to normalize
#'   \code{model.y}; the estimated \code{q_norm[2]} is subtracted from
#'   \code{model.y}, and the resultant time series is divided by the IQR
#'   \code{q_norm[3]-q_norm[1]} (i.e. \code{c(0.1,0.5,0.9)} - \code{q_norm}[2]
#'   should \emph{probably} always be \code{0.5}, and \code{q_norm}[c(1,3)]
#'   should \emph{probably} always be symmetric)
#' @param q_bulk a vector giving the non-tail quantiles to be estimated
#' @param q_tail a vector giving the tail quantiles to be estimated - these are
#'   calculated with respect to the min and max of \code{q_bulk}, so if you set
#'   \code{q_tail[1]} to \code{0.75}, that 'real' value of that quantile is
#'   actually \code{max(q_bulk)+(1-max(q_bulk))*0.75}, and so forth.
#' @param year.range manual specification of the year range (as a \code{2 x 1}
#'   vector giving the first and last year of the inputted time series. This is
#'   to allow for a general treatment of inputted time series with multiple runs
#'   - xts objects are ordered, which is problematic when you have 40 data
#'   points for the same date in what should be separate time series).
#' @param lat,lon the latitude/longitude of the pixel. If \code{get.volc=TRUE},
#'   \code{lat} is used to determine which base function to load. \code{lon} (and
#'   \code{lat}, if \code{get.volc=FALSE}) is merely outputted directly into the
#'   output list. Can be left empty.
#' @param bases.dir directory where the basis functions are stored. Within the
#'   normal file structure, this would be "\code{[defaults$aux.dir],bases/}".
#' @param norm.x,bulk.x,tail.x directly input bases (calculated with
#'   \code{\link{get.predictors}}) if desired. Otherwise, the code will first
#'   attempt to load them in \code{bases.dir}, then recalculate using
#'   \code{\link{get.predictors}} if not found. Can be left empty.
#' @param get.volc if loading or calculating basis functions, sets whether or
#'   not to include the volcanic CO2 fit in the normalization basis function (by
#'   default \code{FALSE}).
#'
#' @return A list is returned, with list members giving the fit coefficients for
#'   the normalization, bulk/primary, and tail fits (\code{coef_norm},
#'   \code{coef_bulk}, and \code{coef_tail}, respectively), in addition to the
#'   inputted \code{q_bulk}, \code{q_tail}, \code{q_norm}, \code{q_all} (a
#'   listing of all the estimated quantiles), \code{norm.x.df},
#'   \code{bulk.x.df}, and \code{tail.x.df}, and \code{lat}, \code{lon}, and
#'   \code{year.range}.
#'
#' @importFrom tictoc tic toc
#' @importFrom quantreg rq.fit.pfn
#' @importFrom xts coredata 

get.quantiles <- function(model.y,norm.x.df, bulk.x.df,tail.x.df,
                          q_norm, q_bulk, q_tail, year.range=c(1920,2099),
                          lat=numeric(), lon=numeric(),
                          bases.dir,
                          norm.x=numeric(),bulk.x=numeric(),tail.x=numeric(),
                          get.volc=TRUE) {

  #----- SETUP ----------------------------------------------------------------
  tic("All quantile fits")
  # Gather quantiles, set threshold for tail quantiles (?)
  q_low = q_tail
  q_high = rev(1-q_tail)
  q_all = c(q_low*q_bulk[1], q_bulk, q_bulk[length(q_bulk)] + q_high*(1-q_bulk[length(q_bulk)])) # All quantiles

  # Get number of different ensemble members
  n_files <- length(model.y)/365/(diff(year.range)+1)

  #----- NORMALIZING FIT ------------------------------------------------------
  tic("total normalization")
  # Get basis vectors for spline fit
  tic("loading bases")
  # Get spline basis functions for normalization fit
  if (length(norm.x)==0) {
    if (get.volc) {
      if (length(lat)==0) {
        stop(paste0("The option get.volc=TRUE needs an inputted latitude to provide the correct forcing."))
      }
      lat.volc <- volc.data$lat
      # Load basis file
      basis.fn <- paste0(bases.dir,"spline_basis_functions_",diff(year.range)+1,"years_",
                         n_files,"runs_",paste0(norm.x.df,collapse="-"),
                         "df_volc",gsub("\\.","-",round(lat.volc[which.min(abs(as.vector(lat.volc)-lat))],1)),".RData")
      if (file.exists(basis.fn)) {load(basis.fn);norm.x<-X;rm(x)} else {norm.x <- get.predictors(n_files=n_files,dfs=norm.x.df,year.range=mod.year.range,get.volc=TRUE,lat=lat)}
      rm(list=c("basis.fn","lat.volc"))
    } else {
      basis.fn <- paste0(bases.dir,"spline_basis_functions_",diff(year.range)+1,"years_",n_files,"runs_",paste0(norm.x.df,collapse="-"),"df.RData")
      if (file.exists(basis.fn)) {load(basis.fn);norm.x<-X;rm(x)} else {norm.x <- get.predictors(n_files=n_files,dfs=norm.x.df,year.range=year.range)}
    }
  }
  toc()

  # Get quantile fits
  tic("fitting quantiles")
  coef_norm = matrix(0, dim(norm.x)[2], length(q_norm))
  try({for (i in 1:length(q_norm)) {
    coef_norm[,i] = rq.fit.pfn(norm.x, y=model.y, tau=q_norm[i], max.bad.fixup=20)$coefficients
  }})
  yqhat = norm.x %*% coef_norm

  # Normalize by subtracting the median and scaling by the IQR set by q_norm[c(1,3)]
  y_norm = (coredata(model.y) - yqhat[,q_norm==0.5]) / (yqhat[,3]  - yqhat[,1])
  toc()

  toc()
  print("Normalization complete")
  rm(list=c("norm.x","yqhat","model.y"))

  #----- BULK FIT (central quantiles after normalizing) -----------------------
  tic("total bulk")
  # Get basis vectors for spline fit
  tic("loading bases")
  if (length(bulk.x)==0) {
    # Load bulk basis function
    basis.fn <- paste0(bases.dir,"spline_basis_functions_",diff(year.range)+1,"years_",
      n_files,"runs_",paste0(bulk.x.df,collapse="-"),"df.RData")
    if (file.exists(basis.fn)) {
      load(basis.fn); bulk.x <- X; rm(list=c("X","basis.fn"))
    } else {
      bulk.x <- get.predictors(n_files=n_files,dfs=bulk.x.df,year.range=year.range)
    }
  }
  toc()

  # Get quantile fits
  tic("quantile fit")
  coef_bulk = matrix(0, dim(bulk.x)[2], length(q_bulk))
  try({for (i in 1:length(q_bulk)) {
    coef_bulk[,i] = rq.fit.pfn(bulk.x, y=y_norm, tau=q_bulk[i], max.bad.fixup=20)$coefficients
  }})
  # Calculate threshold for exceedences (for tail fit below)
  yhat_low = bulk.x %*% coef_bulk[,1]
  yhat_high = bulk.x %*% coef_bulk[,length(q)]
  toc()

  toc()
  print("Bulk spline fit complete")
  rm(list=c("bulk.x"))


  #----- TAIL FIT (exceedences past max/min bulk fit) -------------------------
  tic("total tail")
  # Get basis vectors for spline fit
  tic("loading bases")
  if (length(tail.x)==0) {
    if (load.bases) {
      load(paste0(bases.dir,"spline_basis_functions_",diff(year.range)+1,"years_",
        n_files,"runs_",paste0(tail.x.df,collapse="-"),"df.RData"))
    } else {
      tail.x = get.predictors(n_files=n_files, dfs=tail.x.df, year.range=year.range)
      attr(tail.x,"dimnames") <- NULL
    }
  }
  toc()

  # For the tails use a model based on exceedences
  tic("quantile fit")
  # Calculate the exceedence past the estimated top and bottom quantiles of [q_bulk]
  lowtail = y_norm - yhat_low
  hightail = y_norm - yhat_high

  coef_low = matrix(0, dim(tail.x)[2], length(q_low))
  coef_high = matrix(0, dim(tail.x)[2], length(q_high))
  for (i in 1:length(q_tail)) {
    coef_low[,i] = rq.fit.pfn(tail.x[lowtail<0,], y=lowtail[lowtail<0], tau=q_low[i])$coef
    coef_high[,i] = rq.fit.pfn(tail.x[hightail>0,], y=hightail[hightail>0], tau=q_high[i])$coef
  }
  toc()

  toc()
  print("Tail fit complete")
  rm(list=c("tail.x","y_norm","yhat_low","yhat_high"))

  print("All fits complete; exporting results")

  #----- OUTPUT ---------------------------------------------------------------
  toc()
  list(coef_norm=coef_norm, coef_bulk=coef_bulk, coef_tail=list(coef_low, coef_high),
    q_bulk=q_bulk, q_tail=q_tail, q_norm=q_norm, q_all=q_all, norm.x.df=norm.x.df,
    bulk.x.df=bulk.x.df, tail.x.df=tail.x.df, lat=lat, lon=lon, year.range=year.range)
}
