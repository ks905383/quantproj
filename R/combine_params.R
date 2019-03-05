#' Re-merge Quantile Regression Coefficients
#'
#' Concatenate the individual parameter files calculated through separate calls
#' to \code{\link{estimate.quantiles}} into a single file, removing duplicates along the
#' way. This is only really necessary for \code{\link{get.quantile.changes}}.
#'
#' @param defaults the output from \code{\link{various_defaults}}
#'
#' @return a merged file containing all calculated regression coefficients in a
#'   large list, \code{params}, and saved in \code{[defaults$mod.data.dir]/params/}
#'   as \code{[filevar]_day_[mod.name]_quantfit_params_[mod.year.range[1]-mod.year.range[2]]_alllocs.RData}.
#'
#' @importFrom rlist list.select

combine.params <- function(defaults) {
  # Set output filename
  output.fn <- paste0(defaults$mod.data.dir,"params/",
                      defaults$filevar,"_day_",defaults$mod.name,
                      "_quantfit_params_",
                      paste0(defaults$mod.year.range,collapse="-"),
                      "_alllocs.RData")

  # Find coefficient files
  fn.search.pattern <- paste0(defaults$filevar,
                      "_day_",defaults$mod.name,"_quantfit_params_",
                      paste0(defaults$mod.year.range,collapse="-"),"_locs.*")
  if (defaults$bootstrapping) {
    fn.search.pattern <- paste0(fn.search.pattern,"_",defaults$block.size,"block",defaults$nboots,"runs")
  }
  fn.search.pattern <- paste0(fn.search.pattern,".RData")

  file.list <- dir(path=paste0(defaults$mod.data.dir,"params/"),
                   pattern=fn.search.pattern)

  # Load and concatenate files
  
  for (filen in file.list) {
    # Load output
    load(paste0(defaults$mod.data.dir,"params/",filen))

    # Attach
    params.all <- c(params.all,params)
    rm(params)
  }

  # Remove duplicated latlon combinations
  dup.idxs <- duplicated(cbind(unlist(list.select(params,lat)),unlist(list.select(params,lon))))
  params.all <- params.all[!dup.idxs]
  # Rename back to output.map
  params <- params.all; rm(list=c("params.all"))

  # Save as .RData file
  save(file=paste0(output.fn,".RData"),params)
}



