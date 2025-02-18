#' @title Detect emergence (poe and toe) and quantify the contributions of each statistical component (marginals and dependence)
#'
#' @description
#' Give periods of emergence (PoE) features (e.g., the number of PoEs, their duration). Compound event (CE) probability signal can either emerge below the lower bound of the natural variability (PoE-low) or above the upper bound of the natural variability (PoE-up).
#' Give also the values of the statistical drivers associated with PoE-low and PoE-up.
#'
#' @param series1 numeric vector.
#' @param series2 numeric vector.
#' @param y_start first year of the whole studied period.
#' @param y_end last year of the whole studied period.
#' @param length_sliding_window number of years of the sliding window.
#' @param step_sliding_window number of years separating two sliding windows.
#' @param threshold1 threshold for series1, used to compute CE probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2). Defaults to NaN means that the threshold is associated with the 95th quantile of the series 1 during the reference period.
#' @param threshold2 threshold for series2, used to compute CE probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).  Defaults to NaN means that the threshold is associated with the 95th quantile of the series 2 during the reference period.
#' @param list_marg1 character vector specifying the names of distribution family to be tested for series1. Defaults to c("normal","GEV","gamma","log-normal",GPD").
#' @param list_marg2 character vector specifying the names of distribution family to be tested for series2. Defaults to c("normal","GEV","gamma","log-normal",GPD").
#' @param list_cop character vector specifying the names of copula family to be tested. Defaults to c("Normal","Clayton","Gumbel","Joe","Frank").
#' @param thresh1_GPD threshold value to be used if GPD family is tested for series 1 (i.e. included in list_marg1). Defaults to NaN.
#' @param thresh2_GPD threshold value to be used if GPD family is tested for series 2 (i.e. included in list_marg2). Defaults to NaN.
#' @param ci_level the confidence level required. Defaults to 0.68.
#' @param method character string, either "ml" (default) or "bootstrap.
#' @param B number of iterations for the bootstrap. Defaults to 500.
#' @param length_smoothing_window width of the smoothing window. Defaults to 5.
#'
#' @return a list containing the fitting features, CE probability signals and the emergence analysis:
#' \itemize{
#'    \item \code{res_distr}  a list containing the fitting features of each tested distribution for each component (marginal 1, marginal 2, copula) :
#'      \itemize{
#'        \item \code{par} a list of estimated parameters,
#'        \item \code{logLik} a list of loglikelihood values,
#'        \item \code{AIC} a list of Akaike information criterion values,
#'        \item \code{BIC} a list of Bayesian information criterion values.
#'        }
#'    \item \code{res_distr_max_minAIC} a list containing the fitting features of the selected distribution for each component (marginal 1, marginal 2, copula). The selected distribution is the one that obtains the highest number of periods with the minimum AIC.
#'    \item \code{biv_thresold} a list containing the two thresholds for both variables, used to computed CE probability.
#'    \item \code{res_5prob} a list containing the following probability signals :
#'      \itemize{
#'        \item \code{prob_margdep} CE probability when the three statistical components evolve (marginal 1, marginal 2, copula),
#'        \item \code{prob_marg} CE probability when both marginals evolve (and the copula parameter is held constant),
#'        \item \code{prob_marg1} CE probability when only the marginal 1 varies,
#'        \item \code{prob_marg2} CE probability when only the marginal 2 varies,
#'        \item \code{prob_dep} CE probability when only the copula varies,
#'       }
#'    \item \code{res_ci} vector of 3 values : the estimated CE probability during the reference period and the confidence interval associated to this value.
#'    \item \code{res_poe_contrib} a list containing the emergence analysis.
#'      PoE features :
#'      \itemize{
#'        \item you can choose the signal for which you want to visualize the PoE (prob_margdep, prob_marg, prob_marg1, prob_marg2, prob_dep),
#'        \item you can choose to visualize either the case where the signal emerges above the natural variability (upper) or below it (lower),
#'        \item you can select the feature to visualize: the number of PoEs (number), the duration of each PoE (length), the starting year of each PoE (starting_year), the date of ToE if it exists (toe).
#'      }
#'      Contribution values during PoE:
#'      \itemize{
#'        \item you can choose the type of PoE on which the contribution values are calculated: the PoE that starts with a ToE (ipoe), PoEs not starting with a ToE (poe), all poe (poe_ipoe),
#'        \item you can choose to visualize either the case where the signal emerges above the natural variability (upper) or below it (lower),
#'        \item you can visualize the contribution of each statistical component change during an emergence (marginal1, marginal2, dep, marginal1&2).
#'      }
#'    }
#'
#' @export
#' @examples
#' series1 = dataHD_3_25_W_48_5_N$dataH
#' series2 = dataHD_3_25_W_48_5_N$dataD
#' result = analyse_emergence(series1,series2,y_start=1950,y_end=2022,
#'                           length_sliding_window=20,step_sliding_window=1)
#'
#' # Print the selected distribution family names
#' cat("variable 1:",result$res_distr_max_minAIC$var1$name," / variable 2:",
#' result$res_distr_max_minAIC$var2$name, " / copula:", result$res_distr_max_minAIC$copula$name)
#'
#' # Print the thresholds used to copute CE probabilities
#' cat("thresholds for variable 1:", result$biv_threshold$threshold1,
#' "and for variable 2:", result$biv_threshold$threshold2)
#'
#' # Plot CE probabilities
#' years = middle_sliding_window(names(result$res_5prob$prob_margdep))
#' y_max=0.05
#' plot(unlist(result$res_5prob$prob_margdep)~years,col="black",ylim=c(0,y_max),type="l",
#' lwd=3,xlab = "Years", ylab = "Probability",main=paste0("CE probabilities "))
#' abline(h=result$res_ci[2],col="red",lty=2)
#' abline(h=result$res_ci[3],col="red",lty=2)
#' lines(unlist(result$res_5prob$prob_marg)~years,type="l",col="blue",lwd=1.5,lty=1)
#' lines(unlist(result$res_5prob$prob_marg1)~years,type="l",col="orange",lwd=1.5,lty=1)
#' lines(unlist(result$res_5prob$prob_marg2)~years,type="l",col="magenta2",lwd=1.5,lty=1)
#' lines(unlist(result$res_5prob$prob_dep)~years,type="l",col="chartreuse3",lwd=1.5,lty=1)
#' legend(x="topleft",legend=c("natural variability","signal p","p_(Tmax,S)","p_Tmax","p_S","p_C"),
#' lty=c(2,1,1,1,1,1),col=c("red", "black","blue","orange","magenta2","green"),lwd = 2, cex=1)
#'
#' # Add upper-PoE to the graph
#' starting_years = result$res_poe_contrib$poe$margdep$upper$starting_year
#' abline(v=starting_years,col="black",lty=2)
#' abline(v=starting_years+result$res_poe_contrib$poe$margdep$upper$length-1,col="black",lty=2)
#'
#' # Print the year of the ToE
#' cat("upper-ToE:",result$res_poe_contrib$poe$margdep$upper$toe)
#'
#' # Print the value of the dependence contribution to the last PoE and to every PoE
#' cat(result$res_poe_contrib$contrib$ipoe$upper$dep,
#' result$res_poe_contrib$contrib$poe_ipoe$upper$dep)
#'
#' # Compare ToE when the dependence is considered or not
#' cat("ToE without considering dependence change:",
#' result$res_poe_contrib$poe$marg$upper$toe)
#' cat("ToE with dependence change:",
#' result$res_poe_contrib$poe$margdep$upper$toe)
#'
analyse_emergence <- function(series1,
                               series2,
                               y_start,
                               y_end,
                               length_sliding_window,
                               step_sliding_window,
                               threshold1=NaN,
                               threshold2=NaN,
                               list_marg1 = c("normal","GEV","gamma","log-normal","GPD"),
                               list_marg2 = c("normal","GEV","gamma","log-normal","GPD"),
                               list_cop = c("Normal","Clayton","Gumbel","Joe","Frank"),
                               thresh1_GPD = NaN,
                               thresh2_GPD = NaN,
                               ci_level=0.68,
                               method="ml",
                               B=500,
                               length_smoothing_window=5) {

  # test different distributions for series1, series2, and the copula
  res_distr = fit_distr(series1,series2,y_start,y_end,length_sliding_window,step_sliding_window,list_marg1,list_marg2,list_cop,thresh1_GPD,thresh2_GPD)

  # select the best distribution family for series1, series2, and the copula acccording to AIC criteria
  res_distr_max_minAIC = sel_distr_max_minAIC(res_distr)

  # if the thresholds are not given as input, we select the 95th percentile of the data during the first period
  if (is.na(threshold1)){
    threshold1=q_marg(res_distr_max_minAIC$var1$name, res_distr_max_minAIC$var1$par[[1]], 0.95)
  }
  if (is.na(threshold2)){
    threshold2=q_marg(res_distr_max_minAIC$var2$name, res_distr_max_minAIC$var2$par[[1]], 0.95)
  }

  biv_threshold = list("threshold1"=threshold1,"threshold2"=threshold2)

  # compute CE probabilities, when one, two or the three statistical components evolve (marginal1, marginal2 and copula)
  res_5prob = compute_5prob(threshold1,threshold2,y_start,y_end,length_sliding_window,step_sliding_window,length_smoothing_window,res_distr_max_minAIC)

  # compute the probability associated to the natural variability (confidence interval of CE probability during the reference period)
  res_ci = determine_uncertainty_init(series1,series2,y_start,y_end,length_sliding_window,threshold1, threshold2,res_distr_max_minAIC, ci_level,method, B)

  # detect poe and compute the contribution metrics
  res_poe_contrib = poe_contrib_from_prob(res_5prob, res_ci, y_start, y_end,length_sliding_window,step_sliding_window)

  return(list(res_distr=res_distr,res_distr_max_minAIC=res_distr_max_minAIC,biv_threshold=biv_threshold,res_5prob=res_5prob,res_ci=res_ci,res_poe_contrib=res_poe_contrib))
}

