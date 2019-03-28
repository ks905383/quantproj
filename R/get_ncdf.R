#' Load Subsets of NetCDF Files By Location, Run, Time
#'
#'  \code{get.ncdf} loads the data stored in NetCDF files, as determined by 
#'  the information in \code{process.inputs.tmp}, output from 
#'  \code{\link{get.process.chunks}}. Data is loaded by region, latitude 
#'  band, run, and time range, and is concatenated across files if necessary.
#'  Data is returned as a list by pixel, with every list member containing 
#'  the \code{lat} and \code{lon} of the pixel, and the raw data 
#'  ("\code{Raw}").
#
#'
#' @section Expected File Structure:
#'  In general, most common forms of climate file structures are supported, 
#'  especially the CMIP5 structure. Variables can either be on a \code{lon x lat} 
#'  grid or stored by linear location. Files can either contain all runs of a model
#'  or can be saved by run. Files can either contain the whole timeframe of a model
#'  run or be split up in consecutive temporal chunks. Furthermore:
#'  \describe{
#'    \item{filename}{the code searches for \emph{NetCDF} files using the search
#'      string "\code{[defaults$filevar]_day_.*nc}" (by default; this can be 
#'    changed by setting \code{defaults$search.str}). Make sure no other NetCDF
#'      files with that pattern are present in the search directory (by default
#'      \code{defaults$mod.data.dir}).}
#'    \item{variable setup}{Currently, the code expects the primary variable to 
#'    have either a location dimension (giving the linear index of a location),
#'    or a lon x lat grid. These are all identified by name - the search terms
#'    used can be set in \code{defaults$varnames} - out-of-the-box, the package 
#'    for example supports "lat", "latitude", "Latitude", and "latitude_1" as 
#'    possible names for the "lat" dimension.}
#'    \item{locations}{The code expects there to be two location variables,
#'      \code{lat} and \code{lon} (CMIP5 syntax), giving the lat/lon location of
#'       every pixel in the file. The names of those variables can be any of the
#'     alternatives given by \code{defaults$varnames}.}
#'    \item{multiple runs}{If there are multiple runs in the file, there should
#'      be a \code{run} variable/dimension in the file giving the run id as an
#'      integer}
#'  }
#'
#' @section Warnings and Notes on File Structure: 
#'  \itemize{
#'    \item Files have to begin on 01/01 and end on 12/31 of a year. 
#'    \item Loading speeds up a lot if as much of the data as possible is stored
#'      in a single file - an example is storing files by region, but with all
#'      runs and years in the same file for the pixels in that region.
#'  }
#' 
#'
#' @param defaults the output from \code{\link{set.defaults}}. The defaults 
#'  used are \code{filevar} and \code{varnames} to search for grid/run  
#'  variablesin the files. 
#' @param process.inputs.tmp the output from \code{\link{get.process.chunks}},
#'  giving the information needed to load a chunk of data from NetCDF file(s), 
#'  including filename(s), locations within files, etc.
#' @param year.range which year range should be output, in a vector 
#'  \code{c(year_i,year_f)}. If nothing set, the full possible year range of 
#'  the file is loaded. 
#'
#' @return A list; each element reperesenting a pixel, containing named
#'  subelements \code{lat}, \code{lon}, and the raw data \code{Raw}.  
#'
#' @importFrom ncdf4 nc_open ncvar_get
#' @importFrom abind abind


get.ncdf <- function(defaults,process.inputs.tmp,year.range=numeric()) {

  #----- SETUP --------------------------------------------------------------

  # Define the array that the climate data will be appended onto. Basically Raw should be [time(x run) x loc].
  Raw <- array(numeric(),c(0,length(process.inputs.tmp$global_loc)))

  # Figure out if there are discontinuities (from the point of view of linear
  # indexing) in the indices (remember that the dim.idxs can be a list too...)
  split.idxs <- lapply(process.inputs.tmp$dim_idxs,function(idx) {
    splits <- c(0,which(diff(idx)>1),length(idx))
    out.tmp <- list()
    for (split.idx in 1:(length(splits)-1)) {
      out.tmp[[split.idx]] <- idx[seq(splits[split.idx]+1,splits[split.idx+1])]
    }
    return(out.tmp)
  })

  # Need to transform split.idxs into a list containing all combinations of the component lists,
  # each with a corresponding idxs for each dimension
  # Get all combinations of subsets between dimension
  subset.idxs <- expand.grid(lapply(unlist(lapply(split.idxs,length)),function(x) seq(1:x)))
  subset.loads <- list()
  for (subset.id in 1:nrow(subset.idxs)) {
    subset.loads[[subset.id]] <- list()
    for (dim.idx in 1:ncol(subset.idxs)) {
      subset.loads[[subset.id]][colnames(subset.idxs)[dim.idx]] <- split.idxs[[dim.idx]][subset.idxs[subset.id,dim.idx]]
    }
  }
  rm(list=c("split.idxs","dim.idx","subset.id"))


  # See if the files are saved in split up time chunks
  time.cleared.fns <- unique(as.character(sapply(process.inputs.tmp$fn,function(fn) gsub(strsplit(fn,'\\_')[[1]][6],'',fn))))
  if (length(time.cleared.fns)<length(process.inputs.tmp$fn)) {
    fn.groups <- list()
    for (fn.idx in 1:length(time.cleared.fns)) {
      fn.groups[[fn.idx]] <- process.inputs.tmp$fn[sapply(process.inputs.tmp$fn,
                                                          function(fn) gsub(strsplit(fn,'\\_')[[1]][6],'',fn))==time.cleared.fns[fn.idx]]
    }
  } else {
    fn.groups <- as.list(process.inputs.tmp$fn)
  }

  #----- LOAD DATA BY CONTINUOUS LINEAR CHUNK -------------------------------
  for (fns in fn.groups) {

    # First determine whether the desired year.range is actually feasible
    start <- array(1,length(process.inputs.tmp$dim_list))
    count <- array(-1,length(process.inputs.tmp$dim_list))
    names(start) <- names(count) <- process.inputs.tmp$dim_list

    starts <- lapply(seq(1:length(fns)),function(x) start)
    counts <- lapply(seq(1:length(fns)),function(x) count)
    #starts <- rep(as.data.frame(start),length(fns)); counts <- rep(as.data.frame(count),length(fns))

    fn.year.ranges <- sapply(fns,function(fn) strtoi(substr(strsplit(strsplit(fn,'\\_|\\.')[[1]][6],'\\-')[[1]],1,4)))

    # If no desired output year range set, just make it the full range of the files
    if (length(year.range)==0) {
      year.range.tmp <- c(min(fn.year.ranges),max(fn.year.ranges))
    } else {
      year.range.tmp <- year.range
    }

    # Determine which time subsections of the file(s) to load (location
    # subdivisions are determined later, by file)
    if (min(year.range.tmp)>=min(fn.year.ranges) && max(year.range.tmp)<=max(fn.year.ranges)) {
      for (fn.idx in 1:length(fns)) {
        if (min(year.range.tmp)>max(fn.year.ranges[,fn.idx]) || max(year.range.tmp)<min(fn.year.ranges[,fn.idx])) {
          # If the year range is completely outside of the range of that file, skip it
          starts[[fn.idx]]["time"] <- 0 # This will be taken as a sign to skip this file
        } else if (min(year.range.tmp)<=min(fn.year.ranges[,fn.idx]) && max(year.range.tmp)>=max(fn.year.ranges[,fn.idx])) {
          # If the year range spans beyond both the edges of that file, load it wholly
          # Don't have to do anything, because starts is already at 1, counts already at -1
        } else if (min(year.range.tmp)>=min(fn.year.ranges[,fn.idx]) && max(year.range.tmp)<=max(fn.year.ranges[,fn.idx])) {
          # If the year range *starts* *and* *ends* in this file...
          starts[[fn.idx]]["time"] <- (year.range.tmp[1]-fn.year.ranges[1,fn.idx])*365+1
          counts[[fn.idx]]["time"] <- (year.range.tmp[2]-year.range.tmp[1]+1)*365
        } else if (min(year.range.tmp)>=min(fn.year.ranges[,fn.idx])) {
          # If the year range *starts* in this file
          starts[[fn.idx]]["time"] <- (year.range.tmp[1]-fn.year.ranges[1,fn.idx])*365+1
          counts[[fn.idx]]["time"] <- (fn.year.ranges[2,fn.idx]-year.range.tmp[1]+1)*365
        } else if (max(year.range.tmp)<=max(fn.year.ranges[,fn.idx])) {
          # If the year range *ends* in this file, starts at 1, so don't need to change starts
          counts[[fn.idx]]["time"] <- (year.range.tmp[2]-fn.year.ranges[1,fn.idx]+1)*365
        }
      }
    } else {
      stop(paste0("The desired output year range, ",paste0(year.range.tmp,collapse='-'),
                  ", is not feasible in the files, which span ",min(fn.year.ranges),"-",max(fn.year.ranges)))
    }

    for (fn.idx in 1:length(fns)) {
      fn <- fns[fn.idx]; start <- starts[[fn.idx]]; count <- counts[[fn.idx]]
      # If the desired year range is in this file, get the data
      if (starts[[fn.idx]]["time"]!=0) {

        # Load data
        ncdata <- nc_open(paste0(process.inputs.tmp$fn_path,fn))

        # Populate an empty array, with dimension length 0 in the
        # location/lon direction, to be used in binding possible
        # multiple loads together below
        out.dims <- count
        out.dims[out.dims==-1] <- ncdata$var$tas$size[out.dims==-1]
        out.dims[process.inputs.tmp$dim_list%in%c("loc","lon","lat")] <- 0
        out.dims <- out.dims[!names(out.dims)%in%"lat"]
        Raw_tmp <- array(numeric(),out.dims); dimnames(Raw) <- names(process.inputs.tmp$dim_list); rm(out.dims)

        if (fn==process.inputs.tmp$fn[1]) {lat <- numeric(); lon <- numeric()}

        # Load by individual chunk of indices
        for (subset.id in 1:nrow(subset.idxs)) {
          # Set location info into start/count
          for (dimx in names(process.inputs.tmp$dim_idxs)) {
            start[[dimx]] <- subset.loads[[subset.id]][[dimx]][1]
            count[[dimx]] <- length(subset.loads[[subset.id]][[dimx]])
          }

          # Load actual data
          Raw_file <- ncvar_get(ncdata,varid=defaults$filevar,
                                start=start,
                                count=count)
          # Concatenate data along location dimension
          Raw_tmp <- abind(Raw_tmp,Raw_file,along=which(process.inputs.tmp$dim_list%in%c("loc","lon")))

          rm(list=c("Raw_file","dimx"))

          # Get lat lon (only need to do it once)
          if (fn==process.inputs.tmp$fn[1] && subset.id==1) {
            # Allow for multiple conventions for each of the dimensions (set by the defaults$varnames list)
            dim.list.all.f <- unique(c(names(ncdata$dim),names(ncdata$var)))

            lat <- c(lat,ncvar_get(ncdata,
                                   dim.list.all.f[dim.list.all.f%in%defaults$varnames$lat],
                                   start=start[c("loc","lat")[c("loc","lat")%in%process.inputs.tmp$dim_list]],
                                   count=count[c("loc","lat")[c("loc","lat")%in%process.inputs.tmp$dim_list]]))
            lon <- c(lon,ncvar_get(ncdata,
                                   dim.list.all.f[dim.list.all.f%in%defaults$varnames$lon],
                                   start=start[c("loc","lon")[c("loc","lat")%in%process.inputs.tmp$dim_list]],
                                   count=count[c("loc","lon")[c("loc","lat")%in%process.inputs.tmp$dim_list]]))

            rm(dim.list.all.f)
          }

        }
        rm(list=c("start","count"))

        # If by run, concatenate end-to-end by run
        # First, needs to be in the form time x loc x run
        dim.list.tmp <- process.inputs.tmp$dim_list; dim.list.tmp <- dim.list.tmp[!dim.list.tmp=="lat"]
        Raw_tmp <- aperm(Raw_tmp,
                         c(which(dim.list.tmp=="time"),
                           which(!dim.list.tmp%in%c("time","run")),
                           which(dim.list.tmp=="run")))
        if (length(dim.list.tmp)>2) {
          Raw_tmp <- apply(Raw_tmp,which(dim.list.tmp%in%c("loc","lon")),cbind)
        }

        nc_close(ncdata); rm(ncdata)

        # Add to the main output Raw
        Raw <- abind(Raw,Raw_tmp,along=1)
        rm(Raw_tmp)
      }
    }
  }

  #----- SETUP INTO LISTS GIVING DATA, lAT, LON; RETURN ---------------------

  # If length(lat) == 1, make it to fit the longitude
  if (length(lat)==1) {
    lat <- rep(lat,each=length(lon))
  }

  # Make into list (the [if] is because it otherwise breaks for 1-D
  # vectors since they're stored explicitly as 1-D in R instead of 2-D
  # with a singleton)
  raw.list <- list()
  if (length(lat)>1) {
    for (x in seq(1,length(lat))) {
      raw.list[[x]] <- list(Raw=Raw[,x],lat=lat[x],lon=lon[x])
    }
    rm(x)
  } else {
    raw.list[[1]] <- list(Raw=Raw,lat=lat,lon=lon)
  }
  Raw <- raw.list
  rm(raw.list)


  # Return
  return(Raw)
}

