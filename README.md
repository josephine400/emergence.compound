
<!-- README.md is generated from README.Rmd. Please edit that file -->

# emergence.compound

<!-- badges: start -->
<!-- badges: end -->

The goal of emergence.compound is to analyse emergence of compound event
probability (the signal). Periods of emergence (PoE) are the moments
during which the signal exceeds the natural variability. The probability
signal is modelled with copula functions in order to disentangle the
different contributions (both marginals and dependence). The package
gives the PoE features (like the frequency and the duration) and the
contribution of each statistical driver to the signal emergence. For
more details on the method and the metrics output by the functions,
refer to the article (<https://doi.org/10.5194/egusphere-2025-461>).

## Installation

You can install the development version of emergence.compound from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("josephine400/emergence.compound")
```

## Testing the package

You can test the package using the data provided in the `data` folder,
which includes two time series: the heat index and the drought index,
for a grid point located in France (-3.25°E/48.5°N). Both indices are
derived from the ERA5 reanalysis dataset.

The main function is `analyse_emergence()` and detailed information can
be assessed using `help(analyse_emergence)`.

``` r
library(emergence.compound)
help(analyse_emergence)
#> ℹ Rendering development documentation for "analyse_emergence"
```

## Contributors

This package was developed thanks to the contributions of:

- **Mathieu Vrac** – <mathieu.vrac@lsce.ipsl.fr>  
- **Bastien François** – <bastien.francois@knmi.nl>

## Aknowledgment

This software product is developed in the context the European Union
Horizon Europe Programme – Grant Agreement number 101058386
(’ÌnterTwin”), the National Research Agency under France 2030 bearing
the reference ANR-22-EXTR-0005 (TRACCS-PC4-EXTENDING project), the
European Union’s Horizon 2020 research and innovation programme under
grant agreement No. 101003469 (“XAIDA”), and the “COESION” project
funded by the French National program LEFE (Les Enveloppes Fluides et
l’Environnement).
