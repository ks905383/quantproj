#' Project climate distribution changes with resampling
#'
#' \code{project.climate} projects a climate timeseries by the
#' distributional changes estimated by \code{\link{get.quantiles}},
#' with several temporal projection options, including
#' resampling an existing time period by day or year.
#'
#' @section Projection:
#'  The \code{base.data} time series is projected by the
#'  distributional changes from a model, estimated through
#'  \code{\link{get.quantiles}}, using the following p
#'  \enumerate{
#'      \item The estimated model quantiles are unpacked by
#'          evaluating the quantile regression model using
#'          the quantile fit parameters (\code{params}) and
#'          the bases (either loaded or calculated)
#'      \item The data to be projected is normalized by
#'          fitting the normalizing quantiles given by
#'          \code{defaults$base.norm.x.df} (subtracting
#'          the second element, and dividing by the
#'          difference between the third and first element;
#'          by default, \code{defaults$base.norm.x.df=c(0.1,0.5,0.9)})
#'      \item The projecting indices are found - for each
#'          day in the projected time series, which day is
#'          to be used as a base for comparison (see section
#'          on "Temporal resampling", below)
#'      \item Match base period data points in the time series
#'          to be projected to quantiles in the projecting
#'          model data
#'      \item For each projected day, apply the change in the
#'          matched quantile in the (normalized) projecting
#'          data time series between the projected day and the
#'          set base day to the time series to be projected
#'      \item Un-normalize the now projected time series
#'  }
#'
#' @section Temporal resampling:
#'  To generate a new, projected, climate time series over
#'  the time period given by \code{output.years}, the base
#'  inputted climate time series must be projected in some
#'  temporal mapping - say, if the base climate time series
#'  spans 1979 - 2010, and the desired \code{output.years}
#'  are 2011 - 2099, for each of the future years, a
#'  corresponding 'base' year must be chosen to compare against.
#'  \code{project.climate} offers three options for this
#'  projection:
#'  \describe{
#'      \item{resample_rep}{resampling with replacement (the
#'          default). For \code{resampling.timescale="year"}
#'          (the default), for each projected year in
#'          \code{output.years}, a random year is chosen from
#'          the base input time frame, and projected by the
#'          change in the model distribution between the
#'          corresponding output year and that chosen base year.
#'          For \code{resampling.timescale="day"}, for each
#'          day of year in \code{output.years}, that same
#'          day of year is chosen from a random year in the
#'          base input time frame.}
#'      \item{raw}{If \code{output.years} spans the same length
#'          as the base input time frame, it is just wholesale
#'          projected forward (i.e. 1979-2010 is projected to
#'          2011-2042 year-by-year). If \code{output.years} is
#'          shorter, just the first years are chosen up to the
#'          span of \code{output.years}. If \code{output.years}
#'          is longer, the projection wraps around - i.e.
#'          2011-2050 is projected from [1979-2010 1979-1987].}
#'      \item{resample}{NOT YET IMPLEMENTED}
#'  }
#'
#' @section Output options:
#'  \describe{
#'  \item{full}{(default) a list object giving the projected time
#'      series as an \code{xts} object (\code{proj.data}), the
#'      coefficients of the quantiles used for normalizing the
#'      base climate time series (\code{base.norm.coef}),
#'      the \code{lat} and \code{lon} of the pixel, the
#'      \code{output.years}, and the linear indices (in the base
#'      period) used a the reference year/days in the projected
#'      period (\code{base.idxs})}
#'  \item{xts}{an xts object with just the projected time series}
#'  \item{vec}{a numeric vector with just the projected time series
#'      (generated from \code{as.numeric([xts])})}
#'  }
#'
#' @param defaults the output from \code{\link{various_defaults}};
#'          for \code{base.norm.x.df}, \code{aux.dir}, and \code{base.year.range}
#' @param params the calculated quantile fit parameters of the
#'          projecting model, taken directly from the output of
#'          \code{\link{get.quantiles}}
#' @param base.data the time series to be projected, MUST BE
#'          AN "XTS" OBJECT (this is to reduce the number of
#'          parameters that have to input into this function -
#'          the xts object carries the time information as well).
#' @param output.years the desired time frame (in years) for the
#'          output, projected climate time series. THIS MUST
#'          INCLUDE EVERY DESIRED FUTURE YEAR - not just the range -
#'          \code{output.years=c(2011,2099)} will result in an
#'          output time series, with two years, \code{2011} and
#'          \code{2099}.
#' @param norm.x,bulk.x,tail.x,norm.x.base For each basis function,
#'          if empty, the code first attempts to load it from file. If
#'          the file is missing, the basis function is calculated through
#'          \code{\link{get.predictors}}. Degrees of freedom (and whether
#'          or not to include the volcanic CO2 forcing) for the basis
#'          functions are taken from the \code{norm.x.df}, \code{bulk.x.df},
#'          \code{tail.x.df}, and \code{get.volc} parameters in the
#'          \code{params} input (to ensure the quantile surfaces are
#'          constructed using the same inputs as they were estimated);
#'          \code{base.norm.x.df} is taken from \code{defaults}.
#' @param output set the output options as explained the "Output
#'          options" section below.
#' @param index.type,resampling.timescale set the temporal resampling
#'          options as explained in the "Temporal resampling" section
#'          below.
#' @param rand.seed.set set random seed for resampling (default=\code{42})
#'
#' @return see "Output options" section above.
#' 
#' @importFrom tictoc tic toc
#' @importFrom xts xts timeBasedSeq
#' @importFrom quantreg rq
#' @importFrom abind abind

# unlike with the wrappers, here the dfs are taken from [params], not [defaults] (except for base.norm.x.df)

project.climate <- function(defaults,params,base.data,output.years,
  norm.x=numeric(),bulk.x=numeric(),tail.x=numeric(),norm.x.base=numeric(),
  output="full",
  index.type="resampling",resampling.timescale="day",rand.seed.set=42) {

if (all(is.nan(base.data))) {
    # If the input is all NaNs, don't bother with all
    # this computation, just send some NaNs back
    t.out <- do.call("c",lapply(output.years,function(yr) {timeBasedSeq(paste0(yr,'/',yr,'/d'))}))
    t.out <- t.out[strftime(t.out,format="%j")!='366'] #(removing leap years)
    proj.data <- xts(rep(NaN,length(output.years)*365),order.by=t.out)

    if (output=="full") {
		list(proj.data=proj.data,base.coef=NaN,lat=NaN,lon=NaN,
     	 n_tails=c(NaN,NaN),base.norm.x.df=defaults$base.norm.x.df, output.years=output.years,base.idxs=NA*numeric(length(output.years)*365))
    } else if (output=="xts") {
    	proj.data
    } else if (output=="vec") {
    	rep(NaN,length(output.years)*365)
    }

} else {

	# ----- SETUP ----------------------------------------------------------------
    tic("Total processing")

    # Get some temporal features of the reanalysis data
    t = time(base.data)
    base.year.range = as.numeric(strftime(t, format = "%Y"))

    # Make sure the projecting model data spans at least
    # the output.years desired (with a warning if this
    # is not the case)
    if (any(output.years>params$year.range[2])) {
    	warning(paste0('The desired [output.years] (max: ',max(output.years),') span farther than the projecting model year range (',paste0(params$year.range,collapse="-"),'); only years ',output.years[1],'-',params$year.range[2],' will be outputted.'))
    	output.years <- output.years[output.years<=params$year.range[2]]
    }

    # ----- GET (CALCULATED) SPLINE FITS --------------------------------------
    # Load bases and calculate spline fits one-by-one
    # to minimize memory burden (one of these for 180
    # years/40 runs can run close to 1GB in memory)
    tic("Evaluating quantile fits")

    # Get spline basis functions for normalization fit
    if (length(norm.x)==0) {
        if (defaults$get.volc) {
            lat.volc <- volc.data$lat
            # Load basis file
            basis.fn <- paste0(defaults$aux.dir,"bases/spline_basis_functions_",diff(defaults$base.year.range)+1,"years_1runs_",
                paste0(defaults$base.norm.x.df,collapse="-"),
                "df_volc",gsub("\\.","-",round(lat.volc[which.min(abs(as.vector(lat.volc)-as.vector(process.inputs.tmp$lat)))],1)),".RData")
            if (file.exists(basis.fn)) {load(basis.fn)} else { norm.x.base <- get.predictors(n_files=1,dfs=params$norm.x.df,year.range=params$year.range,get.volc=TRUE,lat=process.inputs.tmp$lat,save.predictors=T)}
            rm(list=c("basis.fn","lat.volc"))
        } else {
            basis.fn <- paste0(defaults$aux.dir,"bases/spline_basis_functions_",diff(defaults$base.year.range)+1,"years_1runs_",paste0(params$norm.x.df,collapse="-"),"df.RData")
            if (file.exists(basis.fn)) {load(basis.fn)} else {norm.x.base <- get.predictors(n_files=1,dfs=params$norm.x.df,year.range=-params$year.range,save.predictors=T)}
        }
    }
    # Calculate spline fit for normalization
    yqhat_norm = norm.x %*% params$coef_norm
    rm(norm.x)

    # Get spline basis for bulk post-normalization fit
    base.fn <- paste0(defaults$aux.dir,"bases/spline_basis_functions_",diff(params$year.range)+1,"years_",
        "1","runs_",paste0(params$bulk.x.df,collapse="-"),"df.RData")
    if (file.exists(base.fn)) {
        load(base.fn); bulk.x <- X; rm(X)
    } else {
        bulk.x <- get.predictors(n_files=1,dfs=params$bulk.x.df,year.range=params$year.range,save.predictors=T)
    }
    rm(base.fn)
    # Calculate spline post-normalization bulk fit
    yqhat_bulk = bulk.x %*% params$coef_bulk
    rm(bulk.x)

    # Get spline basis functions for tail fit
    base.fn <- paste0(defaults$aux.dir,"bases/spline_basis_functions_",diff(params$year.range)+1,"years_",
        "1","runs_",paste0(params$tail.x.df,collapse="-"),"df.RData")
    if (file.exists(base.fn)) {
        load(base.fn); tail.x <- X; rm(X)
    } else {
        tail.x <- get.predictors(n_files=1,dfs=params$tail.x.df,year.range=params$year.range,save.predictors=T)
    }
    rm(base.fn)
    # Calculate spline tail exceedence fits
    yqhat_low = tail.x %*% params$coef_tail[[1]] + c(yqhat_bulk[,1])
    yqhat_high = tail.x %*% params$coef_tail[[2]] + c(yqhat_bulk[,length(params$q_bulk)])
    rm(tail.x)

    # Concatenate all post-normalization quantiles together
    q_full <- abind(yqhat_low,yqhat_bulk,yqhat_high,along=2)

    # Clean house
    rm(list=c("yqhat_bulk","yqhat_low","yqhat_high"))
    toc()

    # ----- NORMALIZE REANALYSIS DATA -----------------------------------------
    # Get basis functions for reanalysis data (with reanalysis dfs set)
    # (Loading instead of recalculating saves around 0.25 secs/run)
    tic("Normalizing base data")
    if (length(norm.x.base)==0) {
        if (defaults$get.volc) {
            lat.volc <- volc.data$lat
            # Load basis file
            basis.fn <- paste0(defaults$aux.dir,"bases/spline_basis_functions_",diff(defaults$base.year.range)+1,"years_1runs_",
                paste0(defaults$base.norm.x.df,collapse="-"),
                "df_volc",gsub("\\.","-",round(lat.volc[which.min(abs(as.vector(lat.volc)-as.vector(process.inputs.tmp$lat)))],1)),".RData")
            if (file.exists(basis.fn)) {
                # Make sure not to overwrite
                tmp.norm.x <- norm.x
                load(basis.fn)
                # Rename to not conflict with LENS basis loading below
                norm.x.base <- norm.x; rm(norm.x);norm.x <- tmp.norm.x; rm(tmp.norm.x)
            } else {
                norm.x.base <- get.predictors(n_files=1,dfs=defaults$base.norm.x.df,year.range=defaults$base.year.range,get.volc=TRUE,lat=process.inputs.tmp$lat)
            }
            rm(list=c("basis.fn","lat.volc"))
        } else {
            basis.fn <- paste0(defaults$aux.dir,"bases/spline_basis_functions_",diff(defaults$base.year.range)+1,"years_1runs_",paste0(defaults$base.norm.x.df,collapse="-"),"df.RData")
            if (file.exists(basis.fn)) {load(basis.fn)} else {norm.x.base <- get.predictors(n_files=1,dfs=defaults$base.norm.x.df,year.range=defaults$base.year.range)}
        }
    }
    # Fit normalizing quantiles (spline fit) to reanalysis data
    base.coef <- rq(as.numeric(base.data)~norm.x.base-1,tau=params$q_norm)$coef
    # Calculate out the fit from the calculated coefficients
    yqhat_base <- norm.x.base %*% base.coef

    # Normalize reanalysis data
    base.data <- (base.data-yqhat_base[,2]) / (yqhat_base[,3]-yqhat_base[,1])
    toc()

    # Clean house
    rm(norm.x.base)

    # ----- GET PROJECTING TEMPORAL INDICES -----------------------------------
    tic("Get projecting temporal indices")
    # Preallocate index vector
    base.idxs <- numeric(length=length(output.years)*365)

    if (index.type=="resampling_rep") {
	    	if (resampling.timescale=="day") {
		      # For each day in the projection (in [output.years]), take a random day from
		      # the original [base.data] in the same position in the year (so for Jan 1., one
		      # of the Jan 1sts in [base.data]), with replacement, and project it by the quantile
		      # change for that quantile-day between that random day's original year and the
		      # projected day's future year.

		    	# So, create a vector of indices of the original data giving the days to project
		    	# for each future day, randomly chosen with resampling, of the length of the
		    	# desired output data (length(output.years)*365), making sure to keep Jan 1sts
		    	# Jan 1st, and so forth.
		    	for (doy.idx in seq(1,365)) {
		    		set.seed(rand.seed.set)
		    		base.idxs[(seq(1,length(output.years))-1)*365+doy.idx] <- (sample.int((diff(range(base.year.range))+1),size=length(output.years),replace=TRUE)-1)*365+doy.idx
		    	}
		    	rm(doy.idx)
		    } else if (resampling.timescale=="year") {
		      # For each year in the projection (in [output.years]), take a random year from
		      # the original [base.data] with replacement, and project it day-by-day by the
		      # quantile change for that quantile-day between the random year and the projected
		      # future year

		    	# Get indices year-by-year
		    	set.seed(rand.seed.set)
		    	base.idxs <- as.vector(t(kronecker((sample.int((diff(range(base.year.range))+1),size=length(output.years),replace=TRUE)-1)*365,matrix(1,1,365))) +
		    						kronecker(seq(1:365),matrix(1,1,length(output.years))))

			} else {
				stop(paste0('The chosen resampling.timescale, ',resampling.timescale,', is not supported. Please choose "year", or "day"'))
			}
	} else if (index.type=="resampling") {
			# For resampling without replacement, count up to the indices to the lenght of the
			# output.year time series - if the original timeseries isn't long enough, restart
			# the sampling again
			if (resampling.timescale=="day") {
				if (length(output.years)<=(diff(range(base.year.range))+1)) {
					base.idxs <- seq(1,length(output.years)*365)
				} else if (length(output.years)>(diff(range(base.year.range))+1)) {
					base.idxs <- c(rep(seq(1,(diff(range(base.year.range))+1)*365),floor(length(output.years)/(diff(range(base.year.range))+1))),
						seq(1,(length(output.years)*365)%%((diff(range(base.year.range))+1)*365)))
				}
			} else if (resmapling.timescale=="year") {

			}

	} else if (index.type=="raw") {
		# For non-resampling indexing (just taking data points in the original order),
		# just count up the indices to the length of the output.year time series - if
		# the original timeseries isn't long enough, wrap around.
		if (length(output.years)<=(diff(range(base.year.range))+1)) {
			base.idxs <- seq(1,length(output.years)*365)
		} else if (length(output.years)>(diff(range(base.year.range))+1)) {
			base.idxs <- c(rep(seq(1,(diff(range(base.year.range))+1)*365),floor(length(output.years)/(diff(range(base.year.range))+1))),
				seq(1,(length(output.years)*365)%%((diff(range(base.year.range))+1)*365)))
		}
	}

    # Now, match those indices to those same days in the projecting model data
    # (by merely accounting for the differents start dates of the two)
    base.idxs.model <- base.idxs+365*(min(base.year.range)-params$year.range[1])

    # Now, get the indices of the future days to be used as the comparison when
    # scaling by the quantile change (indices of the projecting model)
    proj.idxs <- rep((output.years-params$year.range[1])*365,each=365) + rep(seq(1,365),length(output.years))

    # Split up present and future projecting model quantiles, for legibility
    q_full_i <- q_full[base.idxs.model,]
    q_full_f <- q_full[proj.idxs,]
    # Remove the full vector, to save space
    rm(q_full)
    toc()

    # ----- MAP QUANTILES -----------------------------------------------------
    tic("Project distribution changes")
    # Get match of (post-normalization) reanalysis value to their corresponding
    # (post-normalization) model quantiles by getting the quantile index of
    # the largest model quantile value smaller than the reanalysis value
    q_idxs <- apply(kronecker(matrix(1,1,length(params$q_all)),
                              as.numeric(base.data)[base.idxs]) > matrix(q_full_i,nrow=length(base.idxs),),
                    1,
                    function(q_excs) max(which(q_excs),1))

    # Get pre-allocated matrix for reanalysis output
    base.dq <- vector("numeric",length(base.idxs))

    # For the days in the data to be projected that are beyond the lowest
    # or highest extreme quantiles, scale the data by the change in that
    # lowest or highest quantile from the first year, to the final year
    edge.idxs <- is.element(q_idxs,c(1,length(params$q_all)))
    base.dq[edge.idxs] <- as.numeric(base.data)[base.idxs[edge.idxs]] + q_full_f[cbind(which(edge.idxs),q_idxs[edge.idxs])] -
                                                                q_full_i[cbind(which(edge.idxs),q_idxs[edge.idxs])]

    # For the non-edge/extreme quantiles, scale the reanalysis data point by
    # the (linearly) interpolated quantile change between the two closest
    # quantiles:
    # T_f = Q_f(idx) + [(T_i - Q_i(idx)) * (Q_f(idx+1) - Q_f(idx)) / (Q_i(idx+1) - Q_i(idx))]
    # for output (projected) T T_f, input (original) T T_i, year_i quantiles
    # Q_i, and year_f quantiles Q_f at the indices idx and idx+1 from [q_idxs]
    frac_q_loc <- (as.numeric(base.data)[base.idxs[!edge.idxs]] -
                    q_full_i[cbind(which(!edge.idxs),q_idxs[!edge.idxs])]) /
                  (q_full_i[cbind(which(!edge.idxs),q_idxs[!edge.idxs]+1)] -
                    q_full_i[cbind(which(!edge.idxs),q_idxs[!edge.idxs])])
    base.dq[!edge.idxs] <- q_full_f[cbind(which(!edge.idxs),q_idxs[!edge.idxs])] +
                frac_q_loc * (q_full_f[cbind(which(!edge.idxs),q_idxs[!edge.idxs]+1)] -
                q_full_f[cbind(which(!edge.idxs),q_idxs[!edge.idxs])])
    rm(list=c("frac_q_loc","edge.idxs"))

    # Get final transformed reanalysis data by un-normalizing the scaled data:
    # T_out = IQR_rea * (IQR_model_f / IQR_model_i) * T_f + Med_model_f - Med_model_i + Med_rea
    # for inter-quartile range (or whatever quantiles used to normalize data,
    # assumed here to be the 1st and 3rd columns of yqhat_norm) IQR, medians Med,
    # and T_f from the scaling above
    proj.data <-  (yqhat_base[base.idxs,3] - yqhat_base[base.idxs,1]) *
                  (yqhat_norm[proj.idxs, 3] - yqhat_norm[proj.idxs, 1]) / (yqhat_norm[base.idxs.model, 3] - yqhat_norm[base.idxs.model, 1]) *
                  base.dq +
                 yqhat_norm[proj.idxs,params$q_norm==0.5] - yqhat_norm[base.idxs.model,params$q_norm==0.5] + yqhat_base[base.idxs,2]

    # ----- OUTPUT ------------------------------------------------------------
    # Change back into xts object
    t.out <- do.call("c",lapply(output.years,function(yr) {timeBasedSeq(paste0(yr,'/',yr,'/d'))}))
    t.out <- t.out[strftime(t.out,format="%j")!='366'] #(removing leap years)
    proj.data <- xts(as.vector(proj.data),order.by=(t.out))


    toc()
    toc()

    # Return based on output options
    if (output=="full") {
    	list(proj.data=proj.data, base.coef=base.coef,lat=params$lat,lon=params$lon,
      		base.norm.x.df=defaults$base.norm.x.df,
      		output.years=output.years,base.idxs=base.idxs)
    } else if (output=="xts") {
    	proj.data
    } else if (output=="vec") {
    	as.numeric(proj.data)
    }


}
}

