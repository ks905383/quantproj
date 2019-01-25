#' Set default inputs for quantile mapping
#'
#' \code{set.defaults} creates a list object that stores the default
#' input parameters for the package \code{quantproj}, to be
#' used as an input for most other functions. It also creates the
#' directories set as defaults (and their subdirectories used in
#' various bits of code in this package) if they do not yet exist.
#'
#' @section Directory structure:
#' Generally, this package uses three main directories (though
#' \code{aux.dir} may be the same as either of the data directories,
#' the user is discouraged from setting the two data directories equal
#' or using pre-existing directories, for file dependency reasons):
#' \describe{
#'		\item{\code{aux_dir}}{a directory containing two sub-directories,
#'			\code{bases} for saved spline basis functions (the
#'			predictors in the quantile regression), and \code{run_logs}
#'			for logging progress on the bigger scripts}
#'		\item{\code{base.data.dir}}{the directory in which to place the raw data to
#'		be projected. Contains a sub-directory, \code{output}, for the resultant
#'		projections.}
#'		\item{\code{mod.data.dir}}{the directory in which to place the model data
#'		whose distributional changes are used to project the data in
#'		\code{base.data.dir}. Contains a sub-directory, \code{params}, in which the
#'		quantile fit parameters are stored.}
#' }
#'
#' @section Note on file structure:
#'  This package assumes all files follow CMIP5 standard file naming
#'  conventions; in other words, it assumes that variable, data frequency, model
#'  name, and file year range are all encoded in the filenames themselves (see
#'  the
#'  \href{https://cmip.llnl.gov/cmip5/docs/cmip5_data_reference_syntax.pdf}{CMIP5
#'  syntax guide} for more information). The filename format is:
#'  \code{[variable]_[frequency]_[model]_[experiment]_[run/ensemble
#'  member]_[timeframe](_[suffix])}
#'
#' @param base.data.dir the directory with raw data to be projected
#' @param mod.data.dir the directory with raw model ensemble data used to project
#' @param aux.dir directory for auxiliary files - saved basis functions, etc.
#'
#' @param filevar variable shorthand (CMIP5 syntax, def: \code{"tas"} for near-surface air temperature)
#' @param freq data frequency shorthand (CMIP5 syntax, def: \code{"day"} for daily)
#' @param mod.name name of model ensemble (def: \code{"LENS"})
#' @param base.name name of data product to be projected (def: \code{"ERA-INTERIM"})
#'
#' @param base.year.range desired time frame of the "base" time period to project (def: \code{c(1979,2010)})
#' @param mod.year.range desired time frame of the time frame to be processed (quantiles calculated) for the model ensmeble; make sure this range includes \code{base.year.range} (def: \code{c(1979,2099)})
#' @param proj.year.range desired time frame to project to (def: \code{c(2011,2099)})
#'
#' @param q_norm normalizing quantiles (must be a \code{3 x 1} vector, i.e. \code{c(0.1,0.5,0.9)})
#' @param q_bulk bulk quantiles (must be a \code{[n_bulk_quantiles x 1]} vector, i.e. \code{c(0.25,0.3,0.5,0.6,0.75)})
#' @param q_tail tail quantiles (must be a \code{[n_bulk_quantiles x 1]} vector, i.e. \code{c(0.5,0.9,0.95)}). These are fractions of the difference between \code{max(q_bulk)} and \code{0}/\code{1}!
#' @param norm.x.df,bulk.x.df,tail.x.df degrees of freedom for normalizing, bulk, and tail quantile regressions, in a \code{3 x 1} vector for \code{c([seasonal],[long-term],[interaction])}
#' @param get.volc include a volcanic predictor in the normalization (def: \code{FALSE})
#'
#' @param nboots number of bootstrap runs (def: \code{50})
#' @param block.size bootstrap block size (def: \code{40}) (ACROSS RUNS - SO 20 means each bootstrap block is across 20 randomly chosen runs)
#'
#' @return a list object containing all of the parameters listed above, to be used
#'  as an input in thd other functions in this package.


set.defaults <- function(
#----- FILE PATHS -----
# The directory in which code is stored
aux.dir="aux/",
# The directory with weather data to be projected
base.data.dir="base_data/",
# The directory with raw model ensemble data
mod.data.dir="model_data/",

#----- VARIABLE -----
# Variable shorthand (use CMIP5)
filevar="tas",

#----- MODEL NAMES AND IDENTIFIERS -----
# Name of the model ensmeble
mod.name="LENS",
# Name of the weather product
base.name="ERA-INTERIM",

#------ QUANTILE MAPPING PARAMETERS -----
# Base data process year range (this is the desired
# length of the processing - if the reanalysis files
# spread farther than this range, only a subset will
# be used)
base.year.range=c(1979,2010),
# Model process year range (this is the desired
# length of the processing, as above. Make sure this
# range includes the base.year.range above!)
mod.year.range=c(1979,2099),

# Years to project to
proj.year.range=c(2011,2099),

# Normalizing quantiles
q_norm=c(0.1,0.5,0.9),
# Normalizing degrees of freedom (model)
norm.x.df=c(14,6,3),
# Normalizing degrees of freedom (weather)
base.norm.x.df=c(10,1,0),
# Include a volcanic predictor in the normalization
get.volc=F,

# Bulk quantiles
q_bulk=c(0.10,0.18,0.25,0.35,0.42,0.50,0.58,0.65,0.75,0.82,0.90),
# Bulk degrees of freedom
bulk.x.df=c(14,6,3),

# Tail quantiles
q_tail=c(0.01,0.10,0.25,0.50,0.75),
# Tail degrees of freedom
tail.x.df=c(3,1,0),

#------ PROJECTION OPTIONS -----
# How to determine 'base' days for projection
index.type="resampling_rep",
resampling.timescale="year",

#------ BOOTSTRAPPING -----
bootstrapping=F,
nboots=50,
block.size=40
) {

	# Warning if the years don't overlap.
	if (min(base.year.range)<min(mod.year.range) || max(base.year.range)>max(mod.year.range)) {
		warning(paste0('The year range chosen for the model processing (',
						paste0(mod.year.range,collapse='-'),') does not ',
						'include the years chosen for the weather data processing (',
						paste0(base.year.range,collapse='-'),'), please reconsider.'))
	}

	# Create directories if they don't yet exist ---------
	# Auxiliary directory and subdirectories
	if (!dir.exists(aux.dir)) {dir.create(aux.dir); print(paste0(aux.dir," created for auxiliary files!"))}
	if (!dir.exists(paste0(aux.dir,"bases/"))) {dir.create(paste0(aux.dir,"bases/")); print(paste0(aux.dir,"bases/","created for predictor basis files!"))}
	if (!dir.exists(paste0(aux.dir,"run_logs/"))) {dir.create(paste0(aux.dir,"run_logs/")); print(paste0(aux.dir,"run_logs/ created for run logs!"))}

	# Bsse data directory
	if (!dir.exists(base.data.dir)) {dir.create(base.data.dir); print(paste0(base.data.dir," created for data to be projected!"))}
	if (!dir.exists(paste0(base.data.dir,"output/"))) {dir.create(paste0(base.data.dir,"output/")); print(paste0(base.data.dir,"output/ created for projected data!"))}

	# Model/projecting data directory
	if (!dir.exists(mod.data.dir)) {dir.create(mod.data.dir); print(paste0(mod.data.dir," created for projecting data!"))}
	if (!dir.exists(paste0(mod.data.dir,"params/"))) {dir.create(paste0(mod.data.dir,"params/")); print(paste0(mod.data.dir,"params/ created for projecting data quantile fits!"))}

	# Save -----------------------------------------------
	defaults <- list()
	for (var.idx in ls()[ls()!="defaults"]) {
		defaults[[var.idx]] <- eval(parse(text=var.idx))
	}
	rm(var.idx)

	return(defaults)
}


