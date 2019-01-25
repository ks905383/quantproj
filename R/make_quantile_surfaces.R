  #' Construct quantile estimates
  #'
  #' \code{make.quantile.surfaces} takes the quantile regression
  #' parameters estimated by \code{\link{get.quantiles}} and
  #' constructs the resultant quantile estimates using the basis
  #' functions used as predictors in the quantile regression.
  #'
  #' @param params the list output from \code{\link{get.quantiles}}
  #'            giving the fit quantiles and their corresponding
  #'            coefficients
  #' @param bulk.x,tail.x,norm.x the basis functions used to
  #'            estimate the coefficients; i.e. loaded from file
  #'            or recalculated using \code{\link{get.predictors}}
  #'
  #' @return a \code{[365 x num_years x num_quantiles]} numeric array
  #'        giving the estimated quantiles (\code{[day of year,year,quantile]})

  make.quantile.surfaces <- function(params,bulk.x,tail.x,norm.x) {
    # Build normalizing estimates
    yqhat = norm.x %*% params$coef_norm

    # Build post-normalized bulk fit estimates
    yqhat_bulk = bulk.x %*% params$coef_bulk
    # Un-normalize bulk fit estimates
    yqhat_bulk_real = yqhat_bulk * (yqhat[,3] - yqhat[,1]) + yqhat[,2]

    # Build exceedence-based tail fit estimates and add them back to the bulk fit
    yqhat_low = tail.x %*% params$coef_tail[[1]]
    yqhat_low = (yqhat_low + yqhat_bulk[,1]) * (yqhat[,3] - yqhat[,1]) + yqhat[,2]

    yqhat_high = tail.x %*% params$coef_tail[[2]]
    yqhat_high = (yqhat_high + yqhat_bulk[,length(params$q_bulk)]) * (yqhat[,3] - yqhat[,1]) + yqhat[,2]

    # Concatenate all quantile estimates
    yqhat_real = cbind(yqhat_low, yqhat_bulk_real, yqhat_high)
    # Reshape into [day of year x year x quantile]
    yqhat_real = array(yqhat_real, dim=c(365,diff(params$year.range)+1,length(params$q_all)))
    # Return
    return(yqhat_real)
  }
