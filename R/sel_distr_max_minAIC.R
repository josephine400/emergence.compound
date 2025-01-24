#' @title Select the best family according to AIC criteria
#'
#' @description Give the distribution family that gets the highest number of times the lowest AIC values.
#'
#' @param res_distr outputs of the function fit_distr(), giving a list with the estimated parameters, the loglikelihood value, AIC and BIC for each tested distribution.
#'
#' @return a list with the name of the selected distribution, the estimated parameters and the loglikelihood value for each component (marginal 1, marginal 2, copula).
#'
sel_distr_max_minAIC <- function(res_distr) {
  res <- list()
  for (var in c("var1", "var2", "copula")) {
    AIC_data <- lapply(res_distr[[var]], function(x) x[['AIC']])
    AIC_list <- lapply(AIC_data, function(x) unlist(x, use.names = FALSE))
    AIC_array <- do.call(cbind, AIC_list)
    family_minAIC <- apply(AIC_array, 1, function(data) names(which.min(data)))
    family_name <- names(table(family_minAIC))[which.max(table(family_minAIC))]
    if (sum(table(family_minAIC)==max(table(family_minAIC)))!=1){
      family_meanAIC=apply(AIC_array,2,mean)
      family_name=names(which.min(family_meanAIC))
    }
    res[[var]][["name"]]=family_name
    res[[var]][["par"]]=res_distr[[var]][[family_name]][["par"]]
    res[[var]][["logLik"]]=res_distr[[var]][[family_name]][["logLik"]]
  }
  return(res)
}

