#' Re-merge quantile map projections
#'
#' Concatenate the individual projection files calculated through separate calls
#' to \code{build.projections} into a single file, removing duplicates along the
#' way.
#'
#' @section Output Files:
#' The file is saved as both an .RData file and a NetCDF file. The final output
#' in the NetCDF file is a \code{[loc x time]} array (variable name
#' \code{[defaults$filevar]}), with variables \code{lat} and \code{lon} being
#' \code{[loc x 1]} vectors giving the location. The file also includes a variable
#' \code{base_idx} giving, for each projected day, the linear index in the
#' original base data that was used as the base day for the projection.
#'
#'
#' @section Beware:
#' Files are searched for in the directory code{[base.data.dir]/output}, with
#' \code{[base.data.dir]} set by \code{defaults}. As a result, this code depends
#' on the desired temporary files being the only files in that directory with
#' the filename pattern searched for (from the \code{file.list} variable in the
#' code).
#'
#' @param defaults the output from \code{\link{various_defaults}}
#' @param comments if desired, a 'comments' global attribute will be
#' 			added to the NetCDF file, with the string input to this
#' 			option as the body
#' @param output.type how the time series is exported into the netcdf file; either 
#' 			"linear" (default), as an [idx x time] array, or "grid," as a 
#'			[lon x lat x time] array
#' @param lat_vec,lon_vec if \code{output.type=="grid"}, this allows the specification
#' 			of the output grid (useful if the projected data doesn't have values for 
#'			the whole global/etc. grid, or if you just want a subset)
#' 
#' @importFrom rlist list.select
#' @importFrom ncdf4 ncdim_def ncvar_def nc_create ncvar_put ncatt_put nc_close


combine.locs <- function(defaults,comments=character(),output.type="linear",
	lat_vec=numeric(),lon_vec=numeric()) {

	#------ Set Output Filename -------------------------------------------------
	# Get better string filename summary for the projection base year choosing
	if (defaults$index.type=="resampling_rep") {
		proj.method <- paste0("resample",defaults$resampling.timescale)
		proj.desc <- paste0("resampling ",defaults$resampling.timescale," with replacement")
	} else if (defaults$index.type=="resampling") {
		proj.method <- paste0("NoRepResample",defaults$resampling.timescale)
		proj.desc <- paste0("resampling ",defaults$resampling.timescale," without replacement")
	} else if (defaults$index.type=="raw") {
		proj.method <- "projectraw"
		proj.desc <- paste0("projecting ",defaults$resampling.timescale,"-by-",defaults$resampling.timescale," in order")
	} else if (defaults$index.type=="raw_mean") {
	  proj.method <- "projectbyavg"
	  proj.desc <- paste0("projecting by the change in the quantiles averaged over the time periods")
	}
	# Set output filename
	fn <- paste0(defaults$base.data.dir,"output/",
		defaults$filevar,"_day_",defaults$base.name,
		"_",defaults$mod.name,"proj_",proj.method,"_",
		paste0(paste0(defaults$proj.year.range,collapse="0101-"),"1231"))
	if (defaults$bootstrapping) {
		fn <- paste0(fn,defaults$block.size,"block",defaults$nboots,"runs")
	}


	#------ Find Output Files ----------------------------------------------------
	if (defaults$bootstrapping) {
		file.list <- dir(path=paste0(defaults$base.data.dir,"output/"),
					 paste0(defaults$filevar,"_day_",defaults$base.name,".*",
					 		paste0(defaults$proj.year.range,collapse="-"),".*_",
					 		defaults$block.size,"block",defaults$nboots,"runs.RData")
					 )
	} else {
		file.list <- dir(path=paste0(defaults$base.data.dir,"output/"),
					 paste0(defaults$filevar,"_day_",defaults$base.name,".*",
					 		paste0(defaults$proj.year.range,collapse="-"),".*RData")
					 )
	}

	#------ Load and Concatenate Files -------------------------------------------
	# Initialize list
	output.all <- list()

	# Load and concatenate files
	for (filen in file.list) {

		# Load output
		load(paste0(defaults$base.data.dir,"output/",filen))

		# Separate out resampling indices
		if (filen==file.list[1]) {
			base.idxs <- output.map[[1]]$base.idxs
		}

		# Attach
		output.all <- c(output.all,list.select(output.map,lat,lon,global_loc,reg,proj.data))

	}

	# Remove duplicated latlon combinations
	dup.idxs <- duplicated(cbind(list.select(output.all,lat),list.select(output.all,lon)))
	output.all <- output.all[!dup.idxs]
	# Rename back to output.map
	output.map <- output.all; rm(list=c("output.all"))

	# Save as .RData file
	save(file=paste0(fn,".RData"),output.map)


	#------ Save as NetCDF file! -------------------------------------------------
	# Get variables (all the messiness below is to make sure NaNs stay NaNs. I don't know why they're there, but that's a question for soon after)
	lon <- unlist(lapply(list.select(output.map,lon),function(x) as.numeric(as.character(x))))
	lat <- unlist(lapply(list.select(output.map,lat),function(x) as.numeric(as.character(x))))
	global_loc <- unlist(lapply(list.select(output.map,global_loc),function(x) as.numeric(as.character(x))))
	Raw <- t(matrix(unlist(list.select(output.map,proj.data)),ncol=length(global_loc),byrow=FALSE))
	rm(output.map)

	if (output.type=="linear") {
		# Dimensions
		latdim <- ncdim_def("lat","deg",lat)
		londim <- ncdim_def("lon","deg",lon)
		locdim <- ncdim_def("global_loc","idx",global_loc,longname="global location index")
		timedim <- ncdim_def("time","days",1:ncol(Raw),calendar="noleap")

		# Variables
		tas <- ncvar_def("tas","K",list(locdim,timedim),NaN,longname="Near-Surface Air Temperature",prec="double")
		latv <- ncvar_def("lat","deg",locdim,NaN,longname="latitude",prec="double")
		lonv <- ncvar_def("lon","deg",locdim,NaN,longname="longitude",prec="double")
		idxv <- ncvar_def("base_idx","days",timedim,NaN,longname="day index of base climate",prec="integer")
		locv <- ncvar_def("global_loc","idx",locdim,NaN,longname="global location index",prec="integer")

		# Create netcdf file
		if (defaults$bootstrapping) {
			# Uhhhhhh, i don't think this should be [locdim]...
			#runid <- ncvar_def("run_id","id",locdim,NaN,longname="bootstrap run index",prec="integer")
			ncout <- nc_create(paste0(fn,".nc"),list(lonv,latv,idxv,runid,tas))
		} else {
			#ncout <- nc_create(paste0(fn,".nc"),list(latv,lonv,idxv,locv,tas))
			ncout <- nc_create(paste0(fn,".nc"),list(lonv,latv,idxv,tas))
		}

		# Add variables
		ncvar_put(ncout,latv,lat)
		ncvar_put(ncout,lonv,lon)
		ncvar_put(ncout,idxv,base.idxs)
		ncvar_put(ncout,tas,Raw)
		if (defaults$bootstrapping) {
			ncvar_put(ncout,runid,rep(1:defaults$nboots,length(unique(lat))))
		}


	} else if (output.type=="grid") {
		# Get the grid vectors (ignores gaps - so if there are pixels with lat = 1 2 5 6, 
		# the vector will be c(1 2 5 6))
		if (length(lon_vec)==0) {lon_vec <- sort(unique(lon))}
		if (length(lat_vec)==0) {lat_vec <- sort(unique(lat))}
		

		if (max(diff(lat_vec))/median(unique(diff(lat_vec)))>2 || max(diff(lon_vec))/median(unique(diff(lon_vec)))>2) {
			warning("grid may be discontinuous (not every lat or lon value is sequential); suggest inputting target grid using the lon_vec and lat_vec options")
		}

		# Assign the values from the original [Raw] vector to the grid, by longitude band
		Raw_out <- array(data=NA,c(length(lon_vec),length(lat_vec),dim(Raw)[2]))
		global_loc_out <- array(data=NA,c(length(lon_vec),length(lat_vec)))
		for (lon_idx in 1:length(lon_vec)) {
			match.idx <- match(lat_vec,lat[which(lon==lon_vec[lon_idx])])
			# This includes support for the ordering in lat_vec being different than in lat
			Raw_out[lon_idx,which(!is.na(match.idx)),] <- Raw[which(lon==lon_vec[lon_idx])[match.idx[!is.na(match.idx)]],]
			global_loc_out[lon_idx,which(!is.na(match.idx))] <- global_loc[which(lon==lon_vec[lon_idx])[match.idx[!is.na(match.idx)]]]
		}
		Raw <- Raw_out; rm(Raw_out)

		# Dimensions
		latdim <- ncdim_def("lat","deg",lat_vec)
		londim <- ncdim_def("lon","deg",lon_vec)
		locdim <- ncdim_def("global_loc","idx",global_loc,longname="global location index")
		timedim <- ncdim_def("time","days",1:dim(Raw)[3],calendar="noleap")

		# Variables
		tas <- ncvar_def("tas","K",list(londim,latdim,timedim),NaN,longname="Near-Surface Air Temperature",prec="double")
		latv <- ncvar_def("lat","deg",latdim,NaN,longname="latitude",prec="double")
		lonv <- ncvar_def("lon","deg",londim,NaN,longname="longitude",prec="double")
		idxv <- ncvar_def("base_idx","days",timedim,NaN,longname="day index of base climate",prec="integer")
		locv <- ncvar_def("global_loc","idx",list(londim,latdim),NaN,longname="global location index",prec="integer")

		# Create netcdf file
		if (defaults$bootstrapping) {
			# Uhhhhhh, i don't think this should be [locdim]...
			#runid <- ncvar_def("run_id","id",locdim,NaN,longname="bootstrap run index",prec="integer")
			ncout <- nc_create(paste0(fn,".nc"),list(idxv,runid,locv,tas))
		} else {
			#ncout <- nc_create(paste0(fn,".nc"),list(latv,lonv,idxv,locv,tas))
			ncout <- nc_create(paste0(fn,".nc"),list(idxv,locv,tas))
		}

		# Add variables to netcdf file
		ncvar_put(ncout,idxv,base.idxs)
		ncvar_put(ncout,locv,global_loc_out)
		ncvar_put(ncout,tas,Raw)
		if (defaults$bootstrapping) {
			ncvar_put(ncout,runid,rep(1:defaults$nboots,length(unique(lat))))
		}
	} else {
		error(paste0("output type '",output.type,"' unrecognized"))
	}

	# add global attributes
	ncatt_put(ncout,0,"variable_short",defaults$filevar)
	ncatt_put(ncout,0,"variable_long","Near-Surface Air Temperature")
	ncatt_put(ncout,0,"frequency","day")
	ncatt_put(ncout,0,"experiment",paste0(defaults$mod.name,"-projected"))
	ncatt_put(ncout,0,"range",paste0(defaults$proj.year.range,collapse="-"))
	ncatt_put(ncout,0,"model_id",defaults$base.name)
	ncatt_put(ncout,0,"projection",paste0("Projected from ",paste0(defaults$base.year.range,collapse="-"),
		" using day-by-day quantile maps from the ",defaults$mod.name," ensemble. "))
	if (length(comments)>0) {
	  ncatt_put(ncout,0,"comments",comments)
	}
	ncatt_put(ncout,0,"projection_method",proj.desc)
	ncatt_put(ncout,0,"attribution","This file was created using the 'quantproj' code package by Kevin Schwarzwald based on methodology and code by Matz Haugen and Haugen et al. (2018).")

	# close the file, writing data to disk
	nc_close(ncout)

	# Fin
	print(paste0(fn," saved!"))
	invisible()
}
