% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_5prob.R
\name{compute_5prob}
\alias{compute_5prob}
\title{CE probabilities}
\usage{
compute_5prob(
  threshold1,
  threshold2,
  y_start,
  y_end,
  length_sliding_window,
  step_sliding_window,
  length_smoothing_window,
  res_distr_max_minAIC
)
}
\arguments{
\item{threshold1}{threshold for series1, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).}

\item{threshold2}{threshold for series2, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).}

\item{y_start}{first year of the whole studied period.}

\item{y_end}{last year of the whole studied period.}

\item{length_sliding_window}{number of years of the sliding window.}

\item{step_sliding_window}{number of years separating two sliding windows.}

\item{length_smoothing_window}{width of the smoothing window. Defaults to 5.}

\item{res_distr_max_minAIC}{output of the function sel_distr_max_minAIC(), giving the fit features of each component (marginal 1, marginal 2 and copula).}
}
\value{
a list containing the following probability signals :
\itemize{
\item \code{prob_margdep} CE probability when the three statistical component evolve,
\item \code{prob_marg} CE probability when both marginals evolve (and the copula parameter is held constant),
\item \code{prob_marg1} CE probability when only the marginal 1 varies,
\item \code{prob_marg2} CE probability when only the marginal 2 varies.
}
}
\description{
Give five compound event (CE) probability signals, each based on different scenarios where one or more components (marginal 1, marginal 2, and dependence) varies while others are held constant.
The probability of a bivariate compound event relies on three components : the two marginal (i.e., univariate) distributions and the dependence structure coupling them, modeled here with copula.
It is possible to compute the probability of a specific CE over a given period by assuming that only one component (e.g, marginal 1) has changed since the reference period, while the two other components (e.g., marginal 2 and the copula) remain unchanged, as they were in the reference period.
}
