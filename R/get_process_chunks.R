#' Get File Information to Allow Processing of Subsets
#'
#' Get information about the climate data files in a folder,
#' and split them up by region-by-latitude chunk to allow for
#' processing subsections of the (very large) files at a time.
#'
#' @section Expected file structure:
#'  \code{get.process.chunks} and the \code{quantproj} packages as a whole
#'  expect the climate data to follow a few conventions:
#'  \describe{
#'    \item{filename}{the code searches for \emph{NetCDF} files using the search
#'      string "\code{[defaults$filevar]_day_.*nc}". Make sure no other NetCDF
#'      files with that pattern are present in the search directory (by default
#'      \code{defaults$mod.data.dir}).}
#'    \item{variable setup}{Currently, the code expects the primary variable to
#'      be a \code{[loc x time (x run)]} \emph{2-dimensional} variable named
#'      \code{defaults$filevar}. Future versions will support \code{[lon x lat x
#'      time (x run)]} grids as well. }
#'    \item{locations}{The code expects there to be two location variables,
#'      \code{lat} and \code{lon} (CMIP5 syntax), giving the lat/lon location of
#'       every pixel in the file}
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
#'
#' @return a list giving, for each region-by-latitude chunk the subset suffix
#'   (\code{reg}) if available, the filename (\code{fn}), the latitude and
#'   longitude coordinates of the pixels in the subset (\code{lat}, one element;
#'   and \code{lon}), the global pixel id (the variable \code{global_loc} if it
#'   exists in the \emph{NetCDF} file, otherwise a new one is created, counting
#'   pixels up by files alphabetically), the local id (\code{local_idxs}) linear
#'   index within the file of each pixel), and the number of experiment runs in
#'   the file (the variable \code{run} in the \emph{NetCDF} file; 1 otherwise).

#should explain the "not monotonically increasing" error...

get.process.chunks <- function(defaults,save.output=FALSE,search.dir=character()) {
	if (length(search.dir)==0) {
		search.dir <- defaults$mod.data.dir
	}
	# Get all the regionally-separated LENS files
	fns <- dir(path=search.dir, pattern=paste0(defaults$filevar,"_day_.*nc"))

	# Start list
	process.inputs <- list()

	# Set global id start in case no global ids are heading in the netcdf files
	global_loc_idxs_tmp <- 0

	# Get list of inputs necessary to run qmapping procedure
	for (fn in fns) {
		# Open netcdf file
		nc <- nc_open(paste0(search.dir,fn))
		# Get indices
		tryCatch(global_loc_idxs <- ncvar_get(nc,"global_loc"),
			error=function(e) {
			# If no indices provided, set manual ones (important for temporary filenames)
			global_loc_idxs <- global_loc_idxs_tmp + seq(1,nc$dim$global_loc$len)
			global_loc_idxs_tmp <- global_loc_idxs_tmp + nc$dim$global_loc$len
			})

		# Get geo grid vars
		lons <- ncvar_get(nc,'lon')
		lats <- ncvar_get(nc,'lat')

		# Get unique latitudes
		lats_unique <- unique(lats)

		# Make list elements per lat band
		process.inputs.tmp <- list()
		for (lat.idx in 1:length(lats_unique)) {
			# Get which indices (counted from inside the file) are in that latitude band
			idxs.tmp <- which(lats==lats_unique[lat.idx])

			# Throw error if indexing is not linear
			if (length(idxs.tmp)>1) {
				if (1!=unique(diff(idxs.tmp))) {
					stop(paste0('The linear indexing is not monotonically increasing by one in the lat band: ',idxs.tmp))
				}
			}

			# Build list
			process.inputs.tmp[[lat.idx]] <- list(
				reg=if (strsplit(fn,"\\_|\\.")[[1]][7]=="nc") {"NA"} else {strsplit(fn,"\\_|\\.")[[1]][7]},
				fn=fn,
				lat=lats_unique[lat.idx],
				lon=lons[idxs.tmp],
				global_loc=global_loc_idxs[idxs.tmp],
				local_idxs=idxs.tmp,
				nruns=if ("run" %in% names(nc$dim)) {length(ncvar_get(nc,'run'))} else {1})
		}

		# Close netcdf file
		nc_close(nc)

		# Concatenate with existing list
		process.inputs <- c(process.inputs,process.inputs.tmp)
		# ...and clean house
		rm(list=c("nc","global_loc_idxs","lons","lats","lats_unique","lat.idx",
			"process.inputs.tmp","idxs.tmp"))
	}
	rm(list=c("fn","fns"))

	# Save
	if (save.output) {save(file=paste0(search.dir,"process_inputs.RData"),process.inputs);print(paste0(search.dir,"process_inputs.RData saved!"))}

	# Return
	invisible(process.inputs)
}
