#' Hot and dry indexes from ERA5
#'
#' A subset of data from the ERA5 reanalysis : list containing the time series of the hot and the dry index, for one grid point located in France at -3.25°E/48.5°N.
#'
#' @format ## `dataHD_3_25_W_48_5_N`
#' A list containing two vectors of 219 values (JJA summer months from 1950 to 2022) :
#' \describe{
#'   \item{dataH}{hot index (Tmax) : the monthly maximum of daily maximum of temperature (in Celcius degree)}
#'   \item{dataD}{drought index (S) : S=-SPEI6 (SPEI : Standardized Precipitation Evapotranspiration Index, without unit)}
#' }
#' @source Hersbach et al., 2020
"dataHD_3_25_W_48_5_N"
