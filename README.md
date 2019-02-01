# quantproj

quantproj is a package allowing the projection of climate time series using the estimated day-of-year distributional changes of different climate time series, based on the methodologies developed by _Haugen et al. (2018)_. In its original intent, the package can be used to project weather time series (reanalysis, station data, etc.) using the estimated changes in quantiles of climate model time series at the same locations. Specifically, the climate model data used in projection can come from models with multiple runs of data over the same timeframe (so-called 'large ensembles', such as CESM-LE), allowing much finer estimations of quantiles than possible with one-run climate models. Normalizing the model and weather data before projecting means no assumptions about the distribution of the weather data have to be made. 

The final outcome is a set of 'future' time series based on weather data, but reflecting the relative distributional change in a model - in other words, a "relatively cold January 1st" in the historical record is projected by the difference in what a "relatively cold January 1st" looks like in the model between the present in the future, which may be different from what the change in "somewhat warm January 1st" or a "relatively cold March 12th," etc. could be. 

This package is optimized for speed and allows for file management, batch processing, and the calculation of some diagnostic metrics. 

This package is based on code originally developed by Matz Haugen (@matzhaugen - [original code here](https://github.com/matzhaugen/future-climate-emulations-analysis)) of which segments are used here with permission.

## Base Code Run
Do you just want to quickly get started building a distributionally-scaled projection of weather data? Here's the basic code run that will get you there:

1. **Set Inputs and Build Data Structure:** Run `defaults <- set.defaults(...)`, which will build the directories needed (if no base directories are specified in the input, `set.defaults()` will build them in your current working directory), and allows you to set which quantiles to estimate, how many degrees of freedom to use in your basis functions, etc. The defaults are set to the processing as done in _Haugne et al. (2018)_. 
2. **Populate Data Structure**: Put your data in the right locations and in the **right format** (see `?get.process.chunks` for information on correct file structures) - the data you want to project should be in `defaults$base.data.dir`, and the data you want to use to project should be in `defaults$mod.data.dir`. Make sure to not put `.nc` files with similar filename structures in those directories;... it could mess things up. 
3. **Generate Basis Functions:** Run `get.predictors()` to build the basis functions for the quantile estimates. You'll want to build one for every combination of degrees of freedom and number of files/runs in your analysis (all functions that need it will generate these by themselves if the files don't exist, but you get a major speed boost from pre-generating them). For the setup used in _Haugen et al. (2018)_, you will need 5: 
   - `get.predictors(n_files=40,dfs=c(14,6,3),save.predictors=TRUE)` (for model normalizing and 'bulk' fits for 40 runs)
   - `get.predictors(n_files=40,dfs=c(3,1,0),save.predictors=TRUE)` (for tail fits for 40 runs)
   - `get.predictors(n_files=1,dfs=c(14,6,3),save.predictors=TRUE)` (for model normalizing and 'bulk' fit evaluation)
   - `get.predictors(n_files=1,dfs=c(3,1,0),save.predictors=TRUE)` (for model tail fit evaluation) 
   - `get.predictors(n_files=1,dfs=c(10,1,0),save.predictors=TRUE)` (for reanalysis normalizing fit (1 run)) 
4. **Estimate quantiles**: Run `estimate.quantiles(defaults)`, to estimate the quantiles in the model data.
5. **Project climate**: Run `build.projection(defaults)`, to use those quantile estimates to project the weather/base data
6. **Merge and save files**: Run `combine.locs(defaults)`, to merge all the output files from `build.projection()` into a single NetCDF output file. Done! 


## Motivation
Researchers and policymakers studying the economic impacts of climate change need useful projections of the future climate. However, raw data from even state-of-the-art climate models (commonly referred to as GCMs - global climate models / general circulation models) cannot be used without processing since their base climates are unrealistic, esepecially in their higher moments (climate models are generally 'tuned' to existing data to a certain extent to improve parametrizations, but 'overtuning' is frowned upon, since it muddles the internal validity of the model). As a result, projections incorporate information from existing weather records using two methodologies. 

- 'Bias correction': future GCM output is scaled by the difference in some metrics between the historical GCM output and historical weather data
- 'Delta Method': historical weather data is scaled by the difference in some metrics between the future GCM output and historical GCM output

Both methods are limited by the projecting metrics they use. The simplest projection - a scaling of the mean (long used in early damage projections) - ignores changes in higher moments of the climate distribution, which have a greater impact on the frequency of policy-relevant extremes than a similar mean change by itself (see i.e. _Katz and Brown (2002)_). Increasingly sophisticated projection techniques have therefore been developed, scaling data by the standard deviation and other metrics of variabiltiy or even scaling by changes in specific quantiles. However, these methods are limited by data availability - due to computation limitations, most GCMs are only run once or a few times (the historical climate record is also only 'run' once), representing but one realization out of a possible distribution of climate outcomes. As a result, gaining a fine-scaled understanding of how models think the tails of distributions are going to change is difficult. 

**Enter Large Ensembles**: Due to recent advances in computing power, some modeling groups have begun running their GCMs many times, with only miniscule differences in the initial conditions between the runs. All resultant differences in the model are the result of chaotic internal model variability. A 40-member ensemble like the [CESM Large Ensemble](http://www.cesm.ucar.edu/projects/community-projects/LENS/) (CESM-LE/LENS) therefore has 40 times as much data for each day than a single run. This allows us to estimate even extreme quantiles of the model ensemble to a high degree of accuracy (within 'model-world' at least), and create a projection of the future climate that provides a much better idea of changes in extremes (again, at least in 'model-world') than scaling by less sophisticated methods would give us. 

## Method Description
For a more complete and rigorous description of the methodology, please consult the original paper by _Haugen et al_, [here](https://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-17-0782.1). We used ERA-INTERIM for the weather data and CESM-LE for the model, but you can use any data that follows the file structure conventions (for example, the new [ERA5](https://www.ecmwf.int/en/forecasts/datasets/archive-datasets/reanalysis-datasets/era5) product and the GFDL large ensemble).  

1. Both the model and weather data are normalized, by subtracting the estimated median and dividing by the estimated IQR
2. Around two dozen quantiles are estimated over the historical and future periods of the model, across 40 runs. Tail quantiles are estimated using thresholds beyond the max and min non-tail quantiles estimated.
3. Each day in the normalized weather data is 'mapped' to the corresponding normalized quantile in the model data for the same day in the historical period, using linear interpolation between estimated quantiles
4. Base periods for the projection are chosen (i.e. 1979-2010 can be projected day-by-day to, say, 2058-2099, or each future year can be projected from a resampled historic year, in a sequence 2011-2099, etc.) 
5. For each future day, the model change in the quantile mapped to the base day in the historical period is added to the base day's value
6. The future data is un-normalized, creating the projected future time series

Crucially, each quantile is assumed to be changing **slowly** and **smoothly** over time, ensured by using cubic splines basis functions. This allows each quantile function to be estimated using `[(number of runs) x (number of days)]` data points, generating pretty stable estimates. The basis functions depend on three parameters:
- a periodic component based only on the _day of year_ (modeling the seasonal cycle)
- a non-periodic component based only on the _year_ (modeling the long-term secular change due to global warming)
- an interaction term (modeling the long-term change in the seasonal cycle)
The degrees of freedom for each of these segments is set in the code through the `set.defaults()` function. 






