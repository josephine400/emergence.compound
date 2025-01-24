#' @title Fitting marginals and copula to data
#'
#' @description Give the fit features for several distributions in order to be able to choose the best marginal and copula family depending on the criteria (loglikelihood, AIC, BIC).
#'
#' @param series1 numeric vector.
#' @param series2 numeric vector.
#' @param y_start first year of the whole studied period.
#' @param y_end last year of the whole studied period.
#' @param length_sliding_window number of years of the sliding window.
#' @param step_sliding_window number of years separating two sliding windows.
#' @param list_marg1 character vector specifying the names of distribution family to be tested for series1.
#' Defaults to c("normal","GEV","gamma","log-normal",GPD").
#' @param list_marg2 character vector specifying the names of distribution family to be tested for series2.
#' Defaults to c("normal","GEV","gamma","log-normal",GPD").
#' @param list_cop character vector specifying the names of copula family to be tested.
#' Defaults to c("Normal","Clayton","Gumbel","Joe","Frank").
#' @param thresh1_GPD threshold value to be used if GPD family is tested for series 1 (i.e. included in list_marg1).
#' Defaults to NaN.
#' @param thresh2_GPD threshold value to be used if GPD family is tested for series 2 (i.e. included in list_marg2).
#' Defaults to NaN.
#'
#' @return a list with the estimated parameters, the loglikelihood value, AIC and BIC for each tested distribution family of each component (marginal 1, marginal 2, copula).
#'
fit_distr <- function(series1,
                      series2,
                      y_start,
                      y_end,
                      length_sliding_window,
                      step_sliding_window,
                      list_marg1 = c("normal","GEV","gamma","log-normal","GPD"),
                      list_marg2 = c("normal","GEV","gamma","log-normal","GPD"),
                      list_cop = c("Normal","Clayton","Gumbel","Joe","Frank"),
                      thresh1_GPD = NaN,
                      thresh2_GPD = NaN) {

  nb_y <- y_end - y_start + 1
  nb_data_per_y <- length(series1) / nb_y
  res_fit=list()

  if (length(series1)!=length(series2)){
    stop("both data series must have the same dimension")
  }
  if (length_sliding_window >=nb_y){
    stop("length_sliding_window should be lower than the number of years")
  }

  series1[is.nan(series2)]=NaN
  series2[is.nan(series1)]=NaN

  # Construct the label names associated with each sliding window
  for (k in seq(1, (nb_y))) {
    if (k == 1) {
      index_start <- k
      index_end <- k * nb_data_per_y + length_sliding_window * nb_data_per_y - 1 * nb_data_per_y
      label_period <- paste0(y_start, "_", y_start + index_end / nb_data_per_y - 1)
    } else {
      index_start <- (k - 1) * step_sliding_window * nb_data_per_y
      index_end <- (k - 1) * step_sliding_window * nb_data_per_y + length_sliding_window * nb_data_per_y
      label_period <- paste0(y_start + index_start / nb_data_per_y, "_", y_start + index_end / nb_data_per_y - 1)
    }

    series1_window <- stats::na.omit(series1[index_start:index_end])
    series2_window <- stats::na.omit(series2[index_start:index_end])

    # only if we have the whole period of a sliding window : for example,
    #if the sliding window is 10 and we have data drom 1950 to 1999, the last window will be 1989-1999
    if (as.numeric(substr(label_period, 6, 9)) < (y_start + nb_y)) {
      for (f in list_marg1) {
        res <- determine_par_marg(series1_window, f, thresh1_GPD)
        res_fit[["var1"]][[f]][["par"]][[label_period]] <- res$par
        res_fit[["var1"]][[f]][["logLik"]][[label_period]] <- res$logLik
        res_fit[["var1"]][[f]][["AIC"]][[label_period]] <- 2 * length(res$par) - 2 * res$logLik
        res_fit[["var1"]][[f]][["BIC"]][[label_period]] <- length(res$par) * log(length(series2_window)) - 2 * res$logLik
      }

      for (f in list_marg2) {
        res <- determine_par_marg(series2_window, f, thresh2_GPD)
        res_fit[["var2"]][[f]][["par"]][[label_period]] <- res$par
        res_fit[["var2"]][[f]][["logLik"]][[label_period]] <- res$logLik
        res_fit[["var2"]][[f]][["AIC"]][[label_period]] <- 2 * length(res$par) - 2 * res$logLik
        res_fit[["var2"]][[f]][["BIC"]][[label_period]] <- length(res$par) * log(length(series2_window)) - 2 * res$logLik
      }

      for (f in list_cop) {
        seriesU_window <- matrix(NaN, ncol = 2, nrow = length(series1_window))
        seriesU_window[, 1] <- VineCopula::pobs(series1_window)
        seriesU_window[, 2] <- VineCopula::pobs(series2_window)
        res <- determine_theta_copula(seriesU_window, f)
        res_fit[["copula"]][[f]][["par"]][[label_period]] <- res$par
        res_fit[["copula"]][[f]][["logLik"]][[label_period]] <- res$logLik
        res_fit[["copula"]][[f]][["AIC"]][[label_period]] <- 2 * length(res$par) - 2 * res$logLik
        res_fit[["copula"]][[f]][["BIC"]][[label_period]] <- length(res$par) * log(length(series2_window)) - 2 * res$logLik
      }
    }
  }

  return(res_fit)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Estimating copula parameter
#'
#' @description
#' Parameter estimation of copulas, i.e., fitting of a copula model to bivariate “pseudo” observations.
#'
#' @param data_U numeric vector of pseudo-observations.
#' @param name_family character string specifying the copula name.
#'
#' @return a list with the estimated parameters (par) and the loglikelihood (logLik)
#'
determine_theta_copula <- function(data_U, name_family) {
  res <- list()
  if (name_family == "Clayton") {
    cop_object3 <- copula::claytonCopula(dim = 2)
    number_f <- 3
    par_lim_inf <- -0.9999
  }
  if (name_family == "Gumbel") {
    cop_object3 <- copula::gumbelCopula(dim = 2)
    number_f <- 4
    par_lim_inf <- 1.0001
  }
  if (name_family == "Frank") {
    cop_object3 <- copula::frankCopula(dim = 2)
    number_f <- 5
    par_lim_inf <- -Inf
  }
  if (name_family == "Joe") {
    cop_object3 <- copula::joeCopula(dim = 2)
    number_f <- 6
    par_lim_inf <- 1.0001
  }

  if (name_family == "Normal") {
    cop_object3 <- copula::normalCopula(dim = 2)
    number_f <- 6
    par_lim_inf <- -0.9999
  }

  test_fit_BFGS <- tryCatch(
    {
      fit <- copula::fitCopula(cop_object3, data_U, method = "ml", optim.method = "BFGS")
    },
    error = function(e) {}
  )

  if (is.null(test_fit_BFGS)) {
    test_fit_NelderMead1 <- tryCatch(
      {
        fit <- copula::fitCopula(cop_object3, data_U, method = "ml", optim.method = "Nelder-Mead")
      },
      error = function(e) {}
    )

    if (is.null(test_fit_NelderMead1)) {
      test_fit_NelderMead2 <- tryCatch(
        {
          fit <- copula::fitCopula(cop_object3, data_U, method = "mpl", optim.method = "Nelder-Mead")
        },
        error = function(e) {}
      )

      if (is.null(test_fit_NelderMead2)) {
        est_BiCop <- VineCopula::BiCopEst(data_U[, 1], data_U[, 2], number_f, method = "mle")
        res[["logLik"]] <- est_BiCop$logLik
        res[["par"]] <- est_BiCop$par
      } else {
        fit <- test_fit_NelderMead2
        res[["par"]] <- fit@estimate
        res[["logLik"]] <- fit@loglik
      }
    }
    else{
    fit <- test_fit_NelderMead1
    res[["par"]] <- fit@estimate
    res[["logLik"]] <- fit@loglik
    }
  }
  else {
    fit <- test_fit_BFGS
    res[["par"]] <- fit@estimate
    res[["logLik"]] <- fit@loglik
  }

  if (res[["par"]] <= par_lim_inf) {
    res[["par"]] <- par_lim_inf
  }
  return(res)
}

#' @title Estimating marginal parameters
#'
#' @description
#' Parameter estimation of univariate distributions, i.e., fitting of a distribution to observations.
#'
#' @param data_x numeric vector.
#' @param name_family character string specifying the distribution family name.
#' @param thresh_GPD threshold value to be used for GPD distribution. Defaults to NaN.
#'
#' @return a list with the estimated parameters (par) and the loglikelihood (logLik)
#'
determine_par_marg <- function(data_x, name_family, thresh_GPD = NaN) {
  res <- list()
  family_distr <- c("exponential", "log-normal","logistic", "gamma","normal")
  if (name_family=="gamma" &  all(data_x>0)==FALSE) {
    res[["par"]]=NaN
    res[["logLik"]]=NaN
    return(res)
  }
  if (name_family=="log-normal" & all(data_x>0)==FALSE) {
    res[["par"]]=NaN
    res[["logLik"]]=NaN
    return(res)
  }

  if (name_family %in% family_distr) {
    fitted <- MASS::fitdistr(data_x, name_family)
  }

  if (name_family == "GPD") {
    if (is.na(thresh_GPD)){
      res[["par"]]=NaN
      res[["logLik"]]=NaN
      return(res)
    }
    else{
      # can have pb as: "observed information matrix is singular; use std.err = FALSE"
      test_fit <- tryCatch(
        {
          (fitted <- evd::fpot(data_x, thresh_GPD))
        },
        error = function(e) {}
      )
      if (is.null(test_fit)) {
        (fitted <- evd::fpot(data_x, thresh_GPD, std.err = FALSE))
      } else {
        fitted <- test_fit
      }
    }
  }
  if (name_family == "GEV") {
    test_fit <- tryCatch(
      {
        (fitted <- evd::fgev(data_x))
      },
      error = function(e) {}
    )
    if (is.null(test_fit)) {
      (fitted <- evd::fgev(data_x, std.err = FALSE))
    } else {
      fitted <- test_fit
    }
  }
  res[["par"]] <- fitted$estimate
  res[["logLik"]] <- stats::logLik(fitted)[1]
  return(res)
}
