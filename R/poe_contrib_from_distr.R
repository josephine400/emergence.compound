#' @title Detect emergence (poe and toe) and quantify the contributions of each statistical component (marginals and dependence)
#'
#' @description Give periods of emergence (PoE) features (e.g., the number of PoEs, their duration). The signal can either emerge below the lower bound of the natural variability (PoE-low) or above the upper bound of the natural variability (PoE-up).
#' Give also the values of the statistical drivers associated with PoE-low and PoE-up.
#'
#' @param res_distr_max_minAIC outputs of the function fit_distr(), giving a list with the estimated parameters, the loglikelihood value, AIC and BIC for each tested distribution.
#' @param series1 numeric vector.
#' @param series2 numeric vector.
#' @param y_start first year of the whole studied period.
#' @param y_end last year of the whole studied period.
#' @param length_sliding_window number of years of the sliding window.
#' @param step_sliding_window number of years separating two sliding windows.
#' @param threshold1 threshold for series1, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2). Defaults to NaN means that the threshold is associated with the 95th quantile of the series 1 during the reference period.
#' @param threshold2 threshold for series2, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).  Defaults to NaN means that the threshold is associated with the 95th quantile of the series 2 during the reference period.
#' @param ci_level the confidence level required. Defaults to 0.68.
#' @param method character string, either "ml" (default) or "bootstrap.
#' @param B number of iterations for the bootstrap. Defaults to 500.
#' @param length_smoothing_window width of the smoothing window. Defaults to 5.
#'
#' @return a list containing the fitting features, CE probability signals and the emergence analysis:
#' \itemize{
#'    \item \code{biv_thresold} a list containing the two thresholds for both variables, used to computed CE probability
#'    \item \code{res_5prob} a list containing the following probability signals :
#'      \itemize{
#'        \item \code{prob_margdep} CE probability when the three statistical components evolve (marginal 1, marginal 2, copula)
#'        \item \code{prob_marg} CE probability when both marginals evolve (and the copula parameter is held constant)
#'        \item \code{prob_marg1} CE probability when only the marginal 1 varies
#'        \item \code{prob_marg2} CE probability when only the marginal 2 varies
#'        \item \code{prob_dep} CE probability when only the copula varies
#'       }
#'    \item \code{res_ci} vector of 3 values : the estimated CE probability during the reference period and the confidence interval associated to this value
#'    \item \code{res_poe_contrib} a list containing the emergence analysis.
#'      PoE features :
#'      \itemize{
#'        \item you can choose the signal for which you want to visualize the PoE (prob_margdep, prob_marg, prob_marg1, prob_marg2, prob_dep)
#'        \item you can choose to visualize either the case where the signal emerges above the natural variability (upper) or below it (lower).
#'        \item you can select the feature to visualize: the number of PoEs (number), the duration of each PoE (length), the starting year of each PoE (starting_year), the date of ToE if it exists (toe)
#'      }
#'      Contribution values during PoE:
#'      \itemize{
#'        \item you can choose the type of PoE on which the contribution values are calculated: the PoE that starts with a ToE (ipoe), PoEs not starting with a ToE (poe), all poe (poe_ipoe)
#'        \item you can choose to visualize either the case where the signal emerges above the natural variability (upper) or below it (lower).
#'        \item you can visualize the contribution of each statistical component change during an emergence (marginal1, marginal2, dep, marginal1&2)
#'      }
#'    }
#'
poe_contrib_from_distr <- function(res_distr_max_minAIC,
                                    series1,
                                    series2,
                                    y_start,
                                    y_end,
                                    length_sliding_window,
                                    step_sliding_window,
                                    threshold1=NaN,
                                    threshold2=NaN,
                                    ci_level=0.68,
                                    method="ml",
                                    B=500,
                                    length_smoothing_window=5) {

  n=length(res_distr_max_minAIC$var1$par)

  if (is.na(threshold1)){
    threshold1=q_marg(res_distr_max_minAIC$var1$name, res_distr_max_minAIC$var1$par[[1]], 0.9)
  }
  if (is.na(threshold2)){
    threshold2=q_marg(res_distr_max_minAIC$var2$name, res_distr_max_minAIC$var2$par[[1]], 0.9)
  }

  biv_threshold = list("threshold1"=threshold1,"threshold2"=threshold2)

  # compute CE probabilities, when one, two or the three statistical components evolve (marginal1, marginal2 and copula)
  res_5prob = compute_5prob(threshold1,threshold2,y_start,y_end,length_sliding_window,step_sliding_window,length_smoothing_window,res_distr_max_minAIC)

  # compute the probability associated to the natural variability (confidence interval of CE probability during the reference period)
  res_ci = determine_uncertainty_init(series1,series2,y_start,y_end,length_sliding_window,threshold1, threshold2,res_distr_max_minAIC, ci_level,method, B)

  # detect poe and compute the contribution metrics
  res_poe_contrib = poe_contrib_from_prob(res_5prob, res_ci, y_start, y_end,length_sliding_window,step_sliding_window)

  return(list(biv_threshold=biv_threshold,res_5prob=res_5prob,res_ci=res_ci,res_poe_contrib=res_poe_contrib))
}
