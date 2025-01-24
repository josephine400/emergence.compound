#' @title Density distribution and quantile function
#'
#' @description Density distribution and quantile function for the following cases : "normal","log"normal","GEV","GPD","gamma".
#'
#' @param name_marg character string specifying the distribution family name. Examples : "normal","log-normal","GEV","GPD","gamma".
#' @param par list of parameters.
#' @param quantile quantile value.
#' @param thresh_GPD threshold value for GPD distribution. Default to NaN.
#'
#' @return p_marg gives the probability associated to the quantile, q_marg gives the quantile associated to the probability
#'
#'@aliases p_marg q_marg
#'
p_marg <- function(name_marg, par, quantile, thresh_GPD = NaN) {
  if (name_marg == "normal") {
    proba <- stats::pnorm(quantile, mean = par[[1]], sd = par[[2]])
  }
  if (name_marg=="log-normal"){
    proba <- stats::plnorm(quantile, meanlog = par[[1]], sdlog = par[[2]])
  }
  if (name_marg == "GEV") {
    proba <- evd::pgev(quantile, loc = par[[1]], scale = par[[2]], shape = par[[3]])
  }
  if (name_marg == "GPD") {
    proba <- evd::pgpd(quantile, loc = thresh_GPD, scale = par[[1]], shape = par[[2]])
  }
  if (name_marg == "gamma") {
    proba <- stats::pgamma(quantile, shape = par[[1]], rate = par[[2]])
  }
  return(proba)
}

q_marg <- function(name_marg, par, prob_to_estimate, thresh_GPD = NaN) {
  if (name_marg == "normal") {
    quantile <- stats::qnorm(prob_to_estimate, mean = par[[1]], sd = par[[2]])
  }
  if(name_marg == "log-normal"){
    quantile <- stats::qlnorm(prob_to_estimate, meanlog = par[[1]], sdlog = par[[2]])
  }
  if (name_marg == "GEV") {
    quantile <- evd::qgev(prob_to_estimate, loc = par[[1]], scale = par[[2]], shape = par[[3]])
  }
  if (name_marg == "GPD") {
    quantile <- evd::qgpd(prob_to_estimate, loc = thresh_GPD, scale = par[[1]], shape = par[[2]])
  }
  if (name_marg == "gamma") {
    quantile <- stats::qgamma(prob_to_estimate, shape = par[[1]], rate = par[[2]])
  }
  return(quantile)
}

#' @title Middle of each sliding window
#'
#' @description Give the middle year of each period.
#'
#' @param label vector of character string specifying each sliding period.
#'
#' @return numeric vector specifying the midpoint year of each period.
#'
#' @export
#'
#' @examples
#' label = c("1950_1969","1951_1970","1952_1971","1953_1972")
#' years = middle_sliding_window(label)
#'
middle_sliding_window <- function(label){
  years_list=c()
  for (k in seq(1,length(label))){
    y_1=as.numeric(substr(label[k],1,4))
    y_2=as.numeric(substr(label[k],6,9))
    y_int=ceiling((y_2+y_1)/2) #round up to the next integer
    years_list=c(years_list,y_int)
  }
  return(years_list)
}
