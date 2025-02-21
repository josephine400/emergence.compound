% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_5prob.R
\name{determine_prob_bivar}
\alias{determine_prob_bivar}
\title{Bivariate exceedance probability}
\usage{
determine_prob_bivar(name_cop, theta, prob_var1, prob_var2)
}
\arguments{
\item{name_cop}{character string specifying the copula name.}

\item{theta}{copula parameter}

\item{prob_var1}{probability for the variable 1 where the joint exceedance probability needs to be evaluated}

\item{prob_var2}{probability for the variable 1 where the joint exceedance probability needs to be evaluated}
}
\value{
the bivariate exceedance probability
}
\description{
Compute the bivariate exceedance probability with copula
}
