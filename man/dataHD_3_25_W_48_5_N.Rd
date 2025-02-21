% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataHD_3_25_W_48_5_N.R
\docType{data}
\name{dataHD_3_25_W_48_5_N}
\alias{dataHD_3_25_W_48_5_N}
\title{Hot and dry indexes from ERA5}
\format{
\subsection{\code{dataHD_3_25_W_48_5_N}}{

A list containing two vectors of 219 values (JJA summer months from 1950 to 2022) :
\describe{
\item{dataH}{hot index (Tmax) : the monthly maximum of daily maximum of temperature (in Celcius degree)}
\item{dataD}{drought index (S) : S=-SPEI6 (SPEI : Standardized Precipitation Evapotranspiration Index, without unit)}
}
}
}
\source{
Hersbach et al., 2020
}
\usage{
dataHD_3_25_W_48_5_N
}
\description{
A subset of data from the ERA5 reanalysis : list containing the time series of the hot and the dry index, for one grid point located in France at -3.25°E/48.5°N.
}
\keyword{datasets}
