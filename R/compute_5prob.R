#' @title CE probabilities
#'
#' @description
#' Give five compound event (CE) probability signals, each based on different scenarios where one or more components (marginal 1, marginal 2, and dependence) varies while others are held constant.
#' The probability of a bivariate compound event relies on three components : the two marginal (i.e., univariate) distributions and the dependence structure coupling them, modeled here with copula.
#' It is possible to compute the probability of a specific CE over a given period by assuming that only one component (e.g, marginal 1) has changed since the reference period, while the two other components (e.g., marginal 2 and the copula) remain unchanged, as they were in the reference period.
#'
#' @param threshold1 threshold for series1, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).
#' @param threshold2 threshold for series2, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).
#' @param y_start first year of the whole studied period.
#' @param y_end last year of the whole studied period.
#' @param length_sliding_window number of years of the sliding window.
#' @param step_sliding_window number of years separating two sliding windows.
#' @param length_smoothing_window width of the smoothing window. Defaults to 5.
#' @param res_distr_max_minAIC output of the function sel_distr_max_minAIC(), giving the fit features of each component (marginal 1, marginal 2 and copula).
#'
#' @return a list containing the following probability signals :
#' \itemize{
#'    \item \code{prob_margdep} CE probability when the three statistical component evolve,
#'    \item \code{prob_marg} CE probability when both marginals evolve (and the copula parameter is held constant),
#'    \item \code{prob_marg1} CE probability when only the marginal 1 varies,
#'    \item \code{prob_marg2} CE probability when only the marginal 2 varies.
#'    }
#'
compute_5prob <- function(threshold1,
                          threshold2,
                          y_start,
                          y_end,
                          length_sliding_window,
                          step_sliding_window,
                          length_smoothing_window,
                          res_distr_max_minAIC) {
  nb_y <- y_end - y_start + 1

  prob_margdep <- list()
  prob_marg <- list()
  prob_marg1 <- list()
  prob_marg2 <- list()
  prob_dep <- list()

  marg1_name=res_distr_max_minAIC$var1$name
  marg2_name=res_distr_max_minAIC$var2$name
  cop_name=res_distr_max_minAIC$copula$name

  for (k in seq(1, (nb_y))) {
    if (k == 1) {
      label_period <- paste0(y_start, "_", y_start + length_sliding_window - 1)
      prob_var1_init <- p_marg(marg1_name, res_distr_max_minAIC$var1$par[[k]], threshold1)
      prob_var2_init <- p_marg(marg2_name, res_distr_max_minAIC$var2$par[[k]], threshold2)
      par_cop_init = res_distr_max_minAIC$copula$par[[k]]
    } else {
      label_period <- paste0(y_start + (k - 1) * step_sliding_window, "_", y_start + (k - 1) * step_sliding_window + length_sliding_window - 1)
    }
    # only if we have the whole period
    if (as.numeric(substr(label_period, 6, 9)) < (y_start + nb_y)) {

      prob_var1 <- p_marg(marg1_name, res_distr_max_minAIC$var1$par[[k]], threshold1)
      prob_var2 <- p_marg(marg2_name, res_distr_max_minAIC$var2$par[[k]], threshold2)
      par_cop = res_distr_max_minAIC$copula$par[[k]]

      # compute bivariate probability
      prob_margdep[[label_period]] <- determine_prob_bivar(cop_name, par_cop, prob_var1, prob_var2)
      prob_marg[[label_period]] <- determine_prob_bivar(cop_name, par_cop_init, prob_var1, prob_var2)
      prob_marg1[[label_period]] <- determine_prob_bivar(cop_name, par_cop_init, prob_var1, prob_var2_init)
      prob_marg2[[label_period]] <- determine_prob_bivar(cop_name, par_cop_init, prob_var1_init, prob_var2)
      prob_dep[[label_period]] <- determine_prob_bivar(cop_name, par_cop, prob_var1_init, prob_var2_init)
    }
  }

  # Smooth the signals
  n_time=length(unlist(prob_margdep))
  ind_smooth <- floor(length_smoothing_window / 2)

  if (length_smoothing_window %% 2 == 0){
    stop("length_smoothing_window has to be odd")
  }

  if (length_smoothing_window >= n_time){
    stop(paste0("length_smoothing_window has to be lower than the signal length : ",length_smoothing_window," \u2265 ",n_time))
  }

  proba_margdep_sm <- apply(stats::embed(unlist(prob_margdep), length_smoothing_window), 1, mean)
  proba_marg_sm <- apply(stats::embed(unlist(prob_marg), length_smoothing_window), 1, mean)
  proba_marg1_sm <- apply(stats::embed(unlist(prob_marg1), length_smoothing_window), 1, mean)
  proba_marg2_sm <- apply(stats::embed(unlist(prob_marg2), length_smoothing_window), 1, mean)
  proba_dep_sm <- apply(stats::embed(unlist(prob_dep), length_smoothing_window), 1, mean)

  # add the first and the last values that can not be used for smoothing
  prob_margdep_SM=prob_margdep
  prob_marg_SM=prob_marg
  prob_marg1_SM=prob_marg1
  prob_marg2_SM=prob_marg2
  prob_dep_SM=prob_dep

  prob_margdep_SM[(ind_smooth+1):(n_time-ind_smooth)]=proba_margdep_sm
  prob_marg_SM[(ind_smooth+1):(n_time-ind_smooth)]=proba_marg_sm
  prob_marg1_SM[(ind_smooth+1):(n_time-ind_smooth)]=proba_marg1_sm
  prob_marg2_SM[(ind_smooth+1):(n_time-ind_smooth)]=proba_marg2_sm
  prob_dep_SM[(ind_smooth+1):(n_time-ind_smooth)]=proba_dep_sm

  return(list(prob_margdep = prob_margdep_SM, prob_marg = prob_marg_SM, prob_marg1 = prob_marg1_SM, prob_marg2 = prob_marg2_SM, prob_dep = prob_dep_SM))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Bivariate exceedance probability
#'
#' @description
#' Compute the bivariate exceedance probability with copula
#'
#' @param name_cop character string specifying the copula name.
#' @param theta copula parameter
#' @param prob_var1 probability for the variable 1 where the joint exceedance probability needs to be evaluated
#' @param prob_var2 probability for the variable 1 where the joint exceedance probability needs to be evaluated
#'
#'@return the bivariate exceedance probability
#'
determine_prob_bivar <- function(name_cop, theta, prob_var1, prob_var2) {
  d <- 2
  expand_proba <- c(prob_var1, prob_var2)

  if (name_cop == "Normal") {
    number_f <- 1
  }
  if (name_cop == "Clayton") {
    number_f <- 3
  }
  if (name_cop == "Gumbel") {
    number_f <- 4
  }
  if (name_cop == "Frank") {
    number_f <- 5
  }
  if (name_cop == "Joe") {
    number_f <- 6
  }
  if (name_cop == "indep") {
    number_f <- 0
  }

  test_fit_Cop <- tryCatch(
    {
      VineCopula::BiCop(number_f, theta)
    },
    error = function(e) {}
  )

  if (is.null(test_fit_Cop)) {
    fitted_cop <- copula::onacopulaL(name_cop, list(theta, 1:d))
    res_prob <- 1 - copula::pCopula(c(1, expand_proba[2]), fitted_cop) - copula::pCopula(c(expand_proba[1], 1), fitted_cop) + copula::pCopula(c(expand_proba[1], expand_proba[2]), fitted_cop)
  } else {
    fitted_cop <- test_fit_Cop
    res_prob <- 1 - VineCopula::BiCopCDF(1, expand_proba[2], fitted_cop) - VineCopula::BiCopCDF(expand_proba[1], 1, fitted_cop) + VineCopula::BiCopCDF(expand_proba[1], expand_proba[2], fitted_cop)
  }

  epsilon <- .Machine$double.eps
  res_prob[res_prob <= epsilon] <- 0
  return(res_prob)
}

