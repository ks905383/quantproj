#' Get File Information to Allow Processing of Subsets
#'
#' Get information about the climate data files in a folder,
#' and split them up by region-by-latitude chunk to allow for
#' processing subsections of the (very large) files at a time.
#'
#' @section Expected File Structure:
#'  In general, most common forms of climate file structures are supported, 
#' 	especially the CMIP5 structure. Variables can either be on a \code{lon x lat} 
#'	grid or stored by linear location. Files can either contain all runs of a model
#' 	or can be saved by run. Files can either contain the whole timeframe of a model
#'	run or be split up in consecutive temporal chunks. Furthermore:
#'  \describe{
#'    \item{filename}{the code searches for \emph{NetCDF} files using the search
#'      string "\code{[defaults$filevar]_day_.*nc}" (by default; this can be 
#'		changed by setting \code{defaults$search.str}). Make sure no other NetCDF
#'      files with that pattern are present in the search directory (by default
#'      \code{defaults$mod.data.dir}).}
#'    \item{variable setup}{Currently, the code expects the primary variable to 
#' 		have either a location dimension (giving the linear index of a location),
#'		or a lon x lat grid. These are all identified by name - the search terms
#' 		used can be set in \code{defaults$varnames} - out-of-the-box, the package 
#'		for example supports "lat", "latitude", "Latitude", and "latitude_1" as 
#'		possible names for the "lat" dimension.}
#'    \item{locations}{The code expects there to be two location variables,
#'      \code{lat} and \code{lon} (CMIP5 syntax), giving the lat/lon location of
#'       every pixel in the file. The names of those variables can be any of the
#' 		 alternatives given by \code{defaults$varnames}.}
#'    \item{multiple runs}{If there are multiple runs in the file, there should
#'      be a \code{run} variable/dimension in the file giving the run id as an
#'      integer}
#'  }
#'
#' @param defaults the output from \code{\link{set.defaults}}. The defaults used
#'   are \code{filevar} and, if \code{search.dir=numeric()} (by default),
#'   \code{mod.data.dir}.
#' @param save.output whether to save the file information as a file in the
#'   search directory (as "\code{process_inputs.RData}"), by default
#'   \code{FALSE}
#' @param search.dir by default, \code{get.process.chunks} searches in the
#'   \code{mod.data.dir} in \code{defaults}. If a different directory should be
#'   examined, use this to set the path.
#' @param show.messages whether to show some useful descriptions of the search
#'	 procedure (by default \code{TRUE})
#'
#' @return a list giving, for each region-by-latitude chunk the subset suffix
#'   (\code{reg}) if available, the filename (\code{fn}), the latitude and
#'   longitude coordinates of the pixels in the subset (\code{lat}, one element;
#'   and \code{lon}), the global pixel id (the variable \code{global_loc} if it
#'   exists in the \emph{NetCDF} file, otherwise a new one is created, counting
#'   pixels up by files alphabetically), the local id (\code{local_idxs}) linear
#'   index within the file of each pixel), the number of experiment runs in
#'   the file (the variable \code{run} in the \emph{NetCDF} file; 1 otherwise), 
#' 	 and the within-file indices along each location dimension (either just 
#'	 location or lon x lat) \code{dim_idxs}. 
#'
#' @importFrom ncdf4 nc_open ncvar_get


# SHOULD A FUTURE VERSION PUT IN WILDCARDS IN THE MODEL / EXPERIMENT FIELDS BY DEFAULT? 
# LIKE WHEN CHECKING REGIONS/RUNS, TO LOOK FOR "tas_day_.*_.*_[run]_[year]_([suffix]).nc"?
# ALSO SHOULD IT BE ABLE TO DETECT MULTIPLE CONVENTIONS IN THE SAME FOLDER? 


get.process.chunks <- function(defaults,save.output=FALSE,search.dir=character(),show.messages=TRUE) {
	if (show.messages) {cat("Beginning file characterization...\n",fill=T)}

	# Start list
	process.inputs <- list()

	# The default search directory is defaults$mod.data.dir
	if (length(search.dir)==0) {
		search.dir <- defaults$mod.data.dir
	}
	if (show.messages) {cat(paste0("Scanning ",search.dir," for files following the search string: '",defaults$search.str,"'"),fill=T)}
	# Get all the regionally-separated LENS files
	fns <- dir(path=search.dir, pattern=defaults$search.str)
	if (show.messages) {cat(paste0(length(fns)," files found.\n(Sample filename: ",fns[1],")\n"),fill=T)}

	# ----- FIGURE OUT HOW THE FILES ARE STRUCTURED ---------------------------
	nc <- nc_open(paste0(search.dir,fns[1]))
	# Figure out whether the files are by region or not (by whether there's a 
	# region suffix in the filename - assuming CMIP5 standard, there are 
	# 7 parts to the filename without the region suffix, and 8 parts with)
	if (length(strsplit(fns[1],"\\_|\\.")[[1]])>7) {
		files.subset <- TRUE
		subset.fns <- dir(path=search.dir,pattern=paste0(paste0(strsplit(fns[1],"\\_|\\.")[[1]][1:6],collapse='_'),".*","\\.nc"))
		if (show.messages) {
			cat(paste0("Files seem to be saved by regions/geographic subsets. ",length(subset.fns)," regions found."),
			"Some sample regions:",paste0(sapply(subset.fns,function(fn) strsplit(fn,"\\_|\\.")[[1]][7])[sample(1:length(subset.fns),min(3,length(subset.fns)))]),
			fill=T)
		}
		rm(subset.fns)
	} else {
		files.subset <- FALSE
		if (show.messages) {cat(paste0("Each file seems to span all analyzed locations (no specific regional subset files found)."),fill=T)}
	}

	# Get number of dimensions of the primary variable
	var.ndims <- nc$var[[defaults$filevar]]$ndims
	# Identify the dimensions of the primary variable by setting them to 
	# standardized variable names (lat, lon, time, run, loc)
	dim.list.var <- character(length=var.ndims)
	dim.list.var.f <- as.character(unlist(list.select(nc$var[[defaults$filevar]]$dim,name)))
	for (dim.idx in 1:length(dim.list.var.f)) {
		dim.list.var[dim.idx] <- names(defaults$varnames)[unlist(lapply(defaults$varnames,function(x) dim.list.var.f[dim.idx]%in%x))]
	}
	# If a dim name has failed to assign, this is a problem. Currently an error
	if (any(sapply(dim.list.var,nchar)==0)) {
		# Maybe see if you can safely switch this to a warning instead?
		stop(paste0("Could not identify what the dimension named '",dim.list.var.f[which(sapply(dim.list.var,nchar)==0)],"'' is."))
	}
	cat(paste0("The primary variable ",defaults$filevar," seems to have ",var.ndims," dimensions: ",paste0(dim.list.var.f,collapse=", ")),fill=T)
	rm(dim.idx)
	# Do the same thing, but with all listed dimensions (ncdf4 in R has a weird
	# way of dealing with variables/dimensions, so sometimes stuff gets lost).
	# This is to find the lat/lon variables in the file if they're just listed
	# as dimensions, but the primary variable is loc-based so they're not stored
	# there. Don't need the 'stop' call as above. 
	dim.list.all <- character(length(nc$dim))
	dim.list.all.f <- unique(c(names(nc$dim),names(nc$var)))
	for (dim.idx in 1:length(dim.list.all.f)) {
		try(dim.list.all[dim.idx] <- names(defaults$varnames)[unlist(lapply(defaults$varnames,function(x) dim.list.all.f[dim.idx]%in%x))],silent=T)
	}
	#dim.list.all <- dim.list.all[!is.na(dim.list.all)]
	rm(dim.idx)


	# Figure out whether the files are by run, or each contains all runs (by 
	# whether the run index changes for files of the same region)
	# Basically, remove the run field from the filename, and see how many files
	# there are with a wildcard instead of that run field
	run.fns <- dir(path=search.dir,
		pattern=paste0(paste0(strsplit(fns[1],"\\_|\\.")[[1]][1:4],collapse='_'),".*",
		paste0(strsplit(fns[1],"\\_|\\.")[[1]][6:(6+files.subset)],collapse='_'),"\\.nc"))
	if (length(run.fns)==1) {
		files.by.run <- FALSE
	} else {
		files.by.run <- TRUE
	}

	# Determine number of runs
	if (!files.by.run) {
		if ("run"%in%dim.list.var) {
			# The keeping of dim.list.var.f is for the case that 
			nruns <- nc$dim[[dim.list.var.f[dim.list.var=="run"]]]$len
			if (show.messages) {cat(paste0("Files seem to be saved with multiple experiment runs per file. ",nruns," experiment runs found."),fill=T)}
		} else {
			nruns <- 1
			if (show.messages) {cat(paste0("One experiment run found."),fill=T)}
		}
	} else {
		nruns <- length(run.fns)
		if (nruns>1) {
			if (show.messages) {cat(paste0("Files seem to be saved separately by experiment run. ",nruns," experiment runs found."),fill=T)}
		} else {
			if (show.messages) {cat(paste0("One experiment run found."),fill=T)}
		}
	}
	rm(run.fns)

	# Figure out whether the files are split up by time, or each contain the 
	# whole timeframe. Basically, remove the year field from the filename, and
	# see how many files there are with a wildcard instead of that year field, 
	# *with* the added check to see whether the files are continuous.
	if (files.subset) {
		#year.fns <- dir(path=search.dir,
		#	pattern=paste0(paste0(strsplit(fns[1],"\\_|\\.")[[1]][1:5],collapse='_'),"\\_[^\\_]*\\_",
		#				   strsplit(fns[1],"\\_|\\.")[[1]][7],"\\.nc"))
		year.fns <- dir(path=search.dir,
						pattern=paste0(paste0(strsplit(fns[1],"\\_|\\.")[[1]][1:5],collapse='_'),
									   "_[[:digit:]\\-]*_",
					                   paste0(strsplit(fns[1],"\\_|\\.")[[1]][7:(length(strsplit(fns[1],"\\_|\\.")[[1]])-1)],collapse="_"),
					                   ".nc")
						)

	} else {
		#year.fns <- dir(path=search.dir,
        #        pattern=paste0(paste0(strsplit(fns[1],"\\_|\\.")[[1]][1:5],collapse='_'),"\\_[^\\_]*\\_",
        #                        "\\.nc"))
		year.fns <- dir(path=search.dir,
						pattern=paste0(paste0(strsplit(fns[1],"\\_|\\.")[[1]][1:5],collapse='_'),
									   "_[[:digit:]\\-]*.nc")
						)
	}	 
	if (length(year.fns)==1) {
		files.by.timesub <- FALSE
		if (show.messages) {cat(paste0("Files seem to be saved in one temporal chunk, spanning ",strsplit(fns[1],"\\_|\\.")[[1]][6],"."),fill=T)}
	} else {
		if (nchar(strsplit(strsplit(year.fns[1],"\\_|\\.")[[1]][6],'-')[[1]][1])==8) {dtform <- "%Y%m%d"
		} else if (nchar(strsplit(strsplit(year.fns[1],"\\_|\\.")[[1]][6],'-')[[1]][1])==6) {dtform <- "%Y%m"
		} else if (nchar(strsplit(strsplit(year.fns[1],"\\_|\\.")[[1]][6],'-')[[1]][1])==4) {dtform <- "%Y"
		} else {stop(paste0("Can't identify what format the date, [",strsplit(strsplit(year.fns[1],"\\_|\\.")[[1]][6],'-')[[1]][1],"] is... Please make sure the 6th filename position is in the format %Y(%m(%d))."))}
		# Now see if the files are temporally continuous
		time.frame <- lapply(year.fns,function(fn) as.Date(strsplit(strsplit(fn,"\\_|\\.")[[1]][6],'-')[[1]],format=dtform))		
		time.cont <-sapply(1:(length(time.frame)-1),function(idx) time.frame[[idx+1]][1]-time.frame[[idx]][2]==1)

		if (all(time.cont)) {
			files.by.timesub <- TRUE
			if (show.messages) {cat(paste0("Files seem to be saved in several temporal chunks, spanning ",
				paste0(sapply(year.fns,function(fn) paste0(sapply(strsplit(strsplit(fn,"\\_|\\.")[[1]][6],"-")[[1]],function(dt) substr(dt,1,4)),collapse="-")),collapse=", "),"."),fill=T)}
		} else {
			files.by.timesub <- FALSE
			if (show.messages) {cat(paste0("Multiple temporal chunks found, but they don't seem to be continguous... (spanning ",
				paste0(sapply(year.fns,function(fn) paste0(sapply(strsplit(strsplit(fn,"\\_|\\.")[[1]][6],"-")[[1]],function(dt) substr(dt,1,4)),collapse="-")),
					collapse=", "),") The files will be treated as separate load chunks completely."),fill=T)}

		}
		rm(list=c("time.frame","time.cont","dtform","year.fns"))
	}

	# If files.by.run or files.by.timesub, reduce the filenames in fns; each process chunk will 
	# contain nrun*ntimesubs files. They're different search patterns because the time listing can
	# either be right before the .nc or separated by the region field; this ensures that the search
	# string ends in .nc either way. 
	if (files.by.run) {fns <- unique(as.character(sapply(fns,function(fn) gsub(strsplit(fn,"\\_|\\.")[[1]][5],".*",fn))))}
	if (files.by.timesub) {fns <- unique(as.character(sapply(fns,function(fn) gsub(strsplit(gsub("\\.nc","",fn),"\\_")[[1]][6],".*",fn))))}

	nc_close(nc)
	rm(nc)

	if (show.messages) {
		if (files.by.run || files.by.timesub) {
			cat(paste0("\n",length(fns)," file group(s) found, each containing a unique combination of models, experiments, locations, and time segments. ",
				"Data in separate files but containing continuous time chunks and/or different runs of the same experiment will be concatenated in [get.ncdf]."),fill=T)
			cat(paste0("\nNow analyzing each of ",length(fns)," file group(s) for load chunks by latitude band..."),fill=T)
		} else {
			cat(paste0("\nNow analyzing each of ",length(fns)," file(s) for load chunks by latitude band...\n"),fill=T)
		}
	}

	# ----- GENERATE PROCESS LIST ---------------------------------------------
	# Set global id start in case no global ids are heading in the netcdf files
	global_loc_idxs_tmp <- 0

	# Get list of inputs necessary to run qmapping procedure
	for (fn in fns) {
		# If files.by.run, this automatically gets all the different run 
		# filenames; if !files.by.run, this should just get the original 
		# filename. 
		fns.tmp <- dir(path=search.dir,pattern=fn)

		# Open netcdf file
		nc <- nc_open(paste0(search.dir,fns.tmp[1]))

		# Get geo grid vars (allowing for multiple naming options of lat/lon)
		# (which() call to avoid issues when dim.list.all has NAs/unidentified 
		# dimensions/variables)
		lons <- ncvar_get(nc,dim.list.all.f[which(dim.list.all=="lon")])
		lats <- ncvar_get(nc,dim.list.all.f[which(dim.list.all=="lat")])
		
		# Switch lons/lats to linear indices from grid if lon/lat grid
		if ("lon" %in% dim.list.var && "lat" %in% dim.list.var) {
			lat_vec <- lats; lon_vec <- lons;
			lats <- matrix(rep(lats,each=length(lons)),ncol=length(lats))
			lons <- matrix(rep(lons,ncol(lats)),ncol=ncol(lats))
		}


		# Get indices
		global_loc_idxs <- tryCatch(ncvar_get(nc,dim.list.all.f[which(dim.list.all=="loc")]),
			error=function(e) {
				# If no indices provided, set manual ones (important for 
				# temporary filenames) by just using the linear indices
				global_loc_idxs_tmp + 1:(length(lons)) #(length works like numel() in matlab, not length() in matlab)
			})
		global_loc_idxs_tmp <- max(global_loc_idxs)

		# Get unique latitudes
		lats_unique <- as.vector(unique(lats))

		# Make list elements per lat band
		process.inputs.tmp <- list()
		for (lat.idx in 1:length(lats_unique)) {
			# Get which indices (counted from inside the file) are in that latitude band
			idxs.tmp <- which(lats==lats_unique[lat.idx])
			lons.out <- lons[idxs.tmp]

			# Get indices of the location dimension(s) <- 2 if lat/lon, 1 if loc
			if ("lon" %in% dim.list.var && "lat" %in% dim.list.var) {
				idx.dims <- list(lon=which(lon_vec==lons.out),lat=which(lat_vec==lats_unique[lat.idx]))
			} else {
				idx.dims <- list(loc=idxs.tmp)
			}

			# Throw error if indexing is not linear (BUT WE CAN FIX THIS BY LOADING CHUNKS, RIGHT?)
			if (length(idxs.tmp)>1) {
				if (length(unique(diff(idxs.tmp)))>1 || unique(diff(idxs.tmp))!=1) {
					if (files.subset) {
						warning(paste0('The linear indexing is not monotonically increasing by one in the lat band: ',lats_unique[lat.idx]),". This slows down NetCDF file loading.")
					} else {
						warning(paste0('The linear indexing is not monotonically increasing by one in region: ',strsplit(fns.tmp[1],"\\_|\\.")[[1]][7],', lat band: ',lats_unique[lat.idx]),". This slows down NetCDF file loading.")
					}
				}
			}

			# Build list
			process.inputs.tmp[[lat.idx]] <- list(
				reg=if (strsplit(fns.tmp[1],"\\_|\\.")[[1]][7]=="nc") {"N/A"} else {strsplit(fns.tmp[1],"\\_|\\.")[[1]][7]},
				fn=fns.tmp,
				fn_path=search.dir,
				lat=lats_unique[lat.idx],
				lon=lons.out,
				global_loc=global_loc_idxs[idxs.tmp],
				local_idxs=idxs.tmp,
				dim_idxs=idx.dims,
				nruns=nruns,
				dim_list=dim.list.var
				)
		}

		# Close netcdf file
		nc_close(nc)

		# Concatenate with existing list
		process.inputs <- c(process.inputs,process.inputs.tmp)
		# ...and clean house
		rm(list=c("nc","global_loc_idxs","lons","lats","lats_unique","lat.idx",
			"process.inputs.tmp","idxs.tmp"))
	}
	if (show.messages) {
		if (files.by.timesub || files.by.run) {
			cat(paste0("Returning information for ",length(process.inputs)," latitude band(s) in ",length(fns)," file group(s)."),fill=T)
		} else {	
			cat(paste0("Returning information for ",length(process.inputs)," latitude band(s) in ",length(fns)," file(s)."),fill=T)
		}
	}
	rm(list=c("fn","fns"))

	# Save
	if (save.output) {save(file=paste0(search.dir,"process_inputs.RData"),process.inputs);print(paste0(search.dir,"process_inputs.RData saved!"))}

	# Return
	invisible(process.inputs)
}
