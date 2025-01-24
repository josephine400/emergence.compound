#' @title Uncertainty of the joint exceedance probability during the reference period
#'
#' @description
#' Compute uncertainty of CE probability with two methods ("ml" and "bootstrap").
#'
#' @param series1 numeric vector.
#' @param series2 numeric vector.
#' @param y_start first year of the whole studied period.
#' @param y_end last year of the whole studied period.
#' @param length_sliding_window number of years of the sliding window.
#' @param threshold1 threshold for series1, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).
#' @param threshold2 threshold for series2, used to compute the joint exceedance probability (the probability that variable1 exceeds threshold1 and variable2 exceeds threshold2).
#' @param res_distr_max_minAIC output of the function sel_distr_max_minAIC(), giving the fit features of each component (marginal 1, marginal 2 and copula).
#' @param ci_level the confidence level required. Defaults to 0.68.
#' @param method character string, either "ml" (default) or "bootstrap.
#' @param B number of iterations for the bootstrap. Defaults to 500.
#'
#' @return vector of 3 values : the estimated CE probability during the reference period and the confidence interval associated to this value.
#'
determine_uncertainty_init <- function(series1,
                                       series2,
                                       y_start,
                                       y_end,
                                       length_sliding_window,
                                       threshold1,
                                       threshold2,
                                       res_distr_max_minAIC,
                                       ci_level=0.68,
                                       method="ml",
                                       B=500){
  #T1<-Sys.time()
  nb_data_per_y <- length(series1) / (y_end - y_start + 1)
  p_ <- length_sliding_window * nb_data_per_y
  data1_init = stats::na.omit(series1[1:p_])
  data2_init = stats::na.omit(series2[1:p_])
  p=length(data1_init)
  dataU_init = array(NaN,dim=c(p,2))
  dataU_init[,1]=VineCopula::pobs(data1_init)
  dataU_init[,2]=VineCopula::pobs(data2_init)
  object_marg1_init=list("name"=res_distr_max_minAIC$var1$name,"par"=res_distr_max_minAIC$var1$par[[1]],"logLik"=res_distr_max_minAIC$var1$logLik[[1]])
  object_marg2_init=list("name"=res_distr_max_minAIC$var2$name,"par"=res_distr_max_minAIC$var2$par[[1]],"logLik"=res_distr_max_minAIC$var2$logLik[[1]])
  object_cop_init=list("name"=res_distr_max_minAIC$copula$name,"par"=res_distr_max_minAIC$copula$par[[1]],"logLik"=res_distr_max_minAIC$copula$logLik[[1]])
  th_prob1 <- p_marg(object_marg1_init$name, object_marg1_init$par, threshold1)
  th_prob2 <- p_marg(object_marg2_init$name, object_marg2_init$par, threshold2)

  lower_bound = (1-ci_level)/2
  upper_bound = 1-lower_bound

  ## estimation with theta uncertainty
  if (method=="bootstrap"){
    ci_theta <- determine_ci_theta_boot(object_cop_init$name, object_cop_init$par, n_sample=B, nb_pt=p, ci_level=ci_level)
  }

  if (method=="ml"){
    ci_theta <- determine_ci_theta_confint(object_cop_init, dataU_init, ci_level=ci_level)
  }

  ci_theta_prob <- c(
    determine_prob_bivar(object_cop_init$name, object_cop_init$par, th_prob1, th_prob2),
    determine_prob_bivar(object_cop_init$name, ci_theta[2], th_prob1, th_prob2),
    determine_prob_bivar(object_cop_init$name, ci_theta[3], th_prob1, th_prob2)
  )

  ## estimation with marginal parameters uncertainty
  if (method=="bootstrap"){
    ci_par_marg1=determine_ci_par_marg_boot(object_marg1_init$name, object_marg1_init$par, n_sample=B, nb_pt=p, ci_level=ci_level)
    ci_par_marg2=determine_ci_par_marg_boot(object_marg2_init$name, object_marg2_init$par, n_sample=B, nb_pt=p, ci_level=ci_level)

  }
  if (method=="ml"){
    ci_par_marg1=determine_ci_par_marg_confint(object_marg1_init$name, data1_init, ci_level=ci_level)
    ci_par_marg2=determine_ci_par_marg_confint(object_marg2_init$name, data2_init, ci_level=ci_level)
  }

  B_m=B

  b_m1_prob <- boot_prob_marg(object_marg1_init$name, ci_par_marg1, threshold1, B_m, ci_level=ci_level)
  b_m2_prob <- boot_prob_marg(object_marg2_init$name, ci_par_marg2, threshold2, B_m, ci_level=ci_level)

  ## combine both uncertainty
  b_prob_margdep <- c()
  for (b in 1:B_m) {
    b_prob_margdep[b] <- determine_prob_bivar(object_cop_init$name, object_cop_init$par, b_m1_prob[b], b_m2_prob[b])
  }

  coord_b_sup_margdep <- which(b_prob_margdep == sort(b_prob_margdep)[trunc(B_m * upper_bound)])[1]
  coord_b_inf_margdep <- which(b_prob_margdep == sort(b_prob_margdep)[trunc(B_m * lower_bound)])[1]

  ci_prob_margdep <- c()
  ci_prob_margdep[1] <- determine_prob_bivar(object_cop_init$name, object_cop_init$par, th_prob1, th_prob2)
  ci_prob_margdep[2] <- determine_prob_bivar(object_cop_init$name, ci_theta[2], b_m1_prob[coord_b_inf_margdep], b_m2_prob[coord_b_inf_margdep])
  ci_prob_margdep[3] <- determine_prob_bivar(object_cop_init$name, ci_theta[3], b_m1_prob[coord_b_sup_margdep], b_m2_prob[coord_b_sup_margdep])

  #T2<-Sys.time()
  #cat("baseline period uncertainty - Time difference : ",T2-T1,"secs \n")
  return(ci_prob_margdep)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Uncertainty of the copula parameter ("bootstrap")
#'
#' @description
#' Give the confidence interval of the copula parameter with a bootstrap approach.
#'
#' @param name_cop character string specifying the copula name.
#' @param theta copula parameter.
#' @param n_sample number of samples for the bootstrap.
#' @param nb_pt number of points to simulate.
#' @param ci_level the confidence level required.
#'
#' @return vector of 3 values : theta (input) and the confidence interval associated to this value estimated with a bootstrap approach.
#'
determine_ci_theta_boot <- function(name_cop, theta, n_sample, nb_pt, ci_level) {
  d <- 2
  if (name_cop == "Clayton") {
    par_lim_inf <- -0.9999
    cop_object2 <- copula::claytonCopula(dim = 2)
    cop_object3 <- copula::claytonCopula(theta, dim = 2)
  }
  if (name_cop == "Gumbel") {
    par_lim_inf <- 1.0001
    cop_object2 <- copula::gumbelCopula(dim = 2)
    cop_object3 <- copula::gumbelCopula(theta, dim = 2)
  }
  if (name_cop == "Frank") {
    par_lim_inf <- theta - 3
    cop_object2 <- copula::frankCopula(dim = 2)
    cop_object3 <- copula::frankCopula(theta, dim = 2)
  }
  if (name_cop == "Joe") {
    par_lim_inf <- 1.0001
    cop_object2 <- copula::joeCopula(dim = 2)
    cop_object3 <- copula::joeCopula(theta, dim = 2)
  }
  if (name_cop == "Normal") {
    par_lim_inf <- -0.9999
    cop_object2 <- copula::normalCopula(dim = 2)
    cop_object3 <- copula::normalCopula(theta, dim = 2)
  }

  sim_theta <- list()
  for (k in 1:n_sample) {
    sim_data <- copula::rCopula(nb_pt, cop_object3)
    sim_theta[[k]] <- determine_theta_copula(sim_data, name_cop)$par
  }

  lower_bound = (1-ci_level)/2
  upper_bound = 1-lower_bound

  ci <- c("MLE" = theta, c(stats::quantile(unlist(sim_theta), lower_bound), stats::quantile(unlist(sim_theta), upper_bound)))

  return(ci)
}

#' @title Uncertainty of the copula parameter ("ml")
#'
#' @description
#' Give the confidence interval of the copula parameter.
#'
#' @param cop_object copula object.
#' @param data_U numeric vector of pseudo-observations.
#' @param ci_level the confidence level required.
#'
#' @return vector of 3 values : theta (input) and the confidence interval associated to this value estimated with the confint() function from the package 'stats'.
#'
determine_ci_theta_confint<-function(cop_object, data_U, ci_level=ci_level){
  # print(paste0("CI ", ci_level, ":"))
  l_mle= cop_object[["logLik"]]
  theta_mle= cop_object[["par"]]
  d=2
  if(cop_object[["name"]]=="Normal"){
    name_family="Normal"
    par_lim_inf=-0.9999
    cop_object3 <- copula::normalCopula(theta_mle, dim=2)
  }
  if(cop_object[["name"]]=="Clayton"){
    name_family="Clayton"
    par_lim_inf=-0.9999
    cop_object3 <- copula::claytonCopula(theta_mle, dim=2)
  }
  if(cop_object[["name"]]=="Gumbel"){
    name_family="Gumbel"
    par_lim_inf=1.0001
    cop_object3 <- copula::gumbelCopula(theta_mle, dim=2)
  }
  if(cop_object[["name"]]=="Frank"){
    name_family="Frank"
    par_lim_inf=theta_mle-3
    cop_object3 <- copula::frankCopula(theta_mle, dim=2)
  }
  if(cop_object[["name"]]=="Joe"){
    name_family="Joe"
    par_lim_inf=1.0001
    cop_object3 <- copula::joeCopula(theta_mle, dim=2)
  }
  if (cop_object[["name"]]!="Normal"){
    cop_object2 <- copula::onacopulaL(name_family, list(theta_mle,1:d))

  }
  ### Version 1
  test_efm=tryCatch({
    fit.tau <- copula::fitCopula(cop_object3, data_U, method="ml",
                                 optim.method="BFGS")
    ci=stats::confint(fit.tau, level=ci_level)
  }, error=function(e){})
  if(is.null(test_efm)){
    ci=c( "MLE"=theta_mle, c(NaN,NaN))
  }else{
    fit.tau <- copula::fitCopula(cop_object3, data_U, method="ml",
                                 optim.method="BFGS")
    ci=stats::confint(fit.tau, level=ci_level)
    ci<-c("MLE"=theta_mle, ci)
  }
  # print(paste0("v1: ", round(ci[1],4), ", [", round(ci[2],4), ", ", round(ci[3],4), "]"))
  ### Version 2
  if(sum(is.na(ci))>0){
    test_efm=tryCatch({
      efm <- copula::emle(data_U, cop_object2)
      pfm <- stats::profile(efm)
    }, error=function(e){})
    if(is.null(test_efm)){
      ci=c( "MLE"=theta_mle, c(NaN,NaN))
    }else{
      efm <- copula::emle(data_U, cop_object2)
      pfm <- stats::profile(efm)
      ci  <- stats::confint(pfm, level=ci_level)
      if(length(ci)==0){
        ci<-c("MLE"=theta_mle, NaN, NaN)
      }else{
        ci<-c("MLE"=theta_mle, ci)
      }
    }
  }
  # print(paste0("v2: ", round(ci[1],4), ", [", round(ci[2],4), ", ",  round(ci[3],4), "]"))

  ### Version 3: if ci still has NaNs
  if(sum(is.na(ci))>0){
    if(is.na(ci[2])){
      print("ci_low a la main")
      th4=seq(theta_mle,par_lim_inf-0.0001,by=-0.0001) ### from theta mle to lower bound
      Chi=stats::qchisq(ci_level, df=1)
      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
      LR=2*(l_mle-Lt1)
      k=1
      while(LR<=Chi & th4[k]>=par_lim_inf){
        k=k+1
        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
        LR=2*(l_mle-Lt1)
      }
      if(k==length(th4)){
        ci[2]=par_lim_inf
      }else{
        ci[2]=th4[k]
      }
      # print(paste0("CI comp: ", ci[2]))
    }
    if(is.na(ci[3])){
      print("ci_high a la main")
      th4=seq(theta_mle,theta_mle+3,by=0.0001)
      Chi=stats::qchisq(ci_level, df=1)
      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
      LR=2*(l_mle-Lt1)
      k=1
      while(LR<=Chi){
        k=k+1
        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
        LR=2*(l_mle-Lt1)
      }
      ci[3]=th4[k]
      # print(paste0("CI comp: ", ci[3]))
    }
  }
  # print(paste0("v3: ", round(ci[1],4), ", [", round(ci[2],4), ", ", round(ci[3],4), "]"))
  ### to avoid cases where ci goes beyond lim inf
  if(ci[2]<par_lim_inf & name_family %in% c("Frank", "Joe", "Gumbel", "Clayton")){ci[2]<-par_lim_inf}
  # print(paste0("vf: ", round(ci[1],4), ", [", round(ci[2],4), ", ",  round(ci[3],4), "]"))
  return(ci)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Uncertainty of the marginal parameters ("bootstrap")
#'
#' @description
#' Give the confidence interval of the marginal parameters with a bootstrap approach.
#'
#' @param name_marg character string specifying the name of the distribution family.
#' @param par list of parameters.
#' @param n_sample number of samples for the bootstrap.
#' @param nb_pt number of points to simulate.
#' @param ci_level the confidence level required.
#' @param thresh_GPD threshold value to be used for GPD distribution. Defaults to NaN.
#'
#' @return vector of 3 values : initial parameters (input) and the confidence interval associated to this value estimated with a bootstrap approach.
#'
determine_ci_par_marg_boot <- function(name_marg, par, n_sample, nb_pt, ci_level,thresh_GPD=NaN){
  ci_par=list()
  sim_par=list()
  q1=(1-ci_level)/2
  q2=1-q1

  if (name_marg=="log-normal"){
    for (k in 1:n_sample) {
      sim_data <- stats::rlnorm(nb_pt, meanlog=par[1],sdlog=par[2])
      fitted = MASS::fitdistr(sim_data, "log-normal")
      sim_par[[k]] = fitted$estimate
    }
    ci_par[["min_meanlog"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q1)
    ci_par[["max_meanlog"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q2)
    ci_par[["min_sdlog"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q1)
    ci_par[["max_sdlog"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q2)
  }

  if (name_marg == "gamma") {
    for (k in 1:n_sample){
      sim_data <- stats::rgamma(nb_pt,shape=par[1],rate=par[2])
      fitted = MASS::fitdistr(sim_data, "gamma")
      sim_par[[k]] = fitted$estimate
    }
    ci_par[["min_shape"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q1)
    ci_par[["max_shape"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q2)
    ci_par[["min_rate"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q1)
    ci_par[["max_rate"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q2)
  }

  if (name_marg=="normal"){
    for (k in 1:n_sample) {
      sim_data <- stats::rnorm(nb_pt, mean=par[1],sd=par[2])
      fitted = MASS::fitdistr(sim_data, "normal")
      sim_par[[k]] = fitted$estimate
    }
    ci_par[["min_mean"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q1)
    ci_par[["max_mean"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q2)
    ci_par[["min_sd"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q1)
    ci_par[["max_sd"]]=stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q2)
  }

  if (name_marg=="GEV"){
    for (k in 1:n_sample){
      sim_data = evd::rgev(nb_pt,loc=par[1],scale=par[2],shape=par[3])
      test_fit <- tryCatch(
        {
          (fitted <- evd::fgev(sim_data))
        },
        error = function(e) {}
      )
      if (is.null(test_fit)) {
        (fitted <- evd::fgev(sim_data, std.err = FALSE))
      } else {
        fitted <- test_fit
      }
      sim_par[[k]] = fitted$estimate
    }
    ci_par[["min_loc"]] <- stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q1)
    ci_par[["max_loc"]] <- stats::quantile(unlist(lapply(sim_par,function(x)x[1]),use.names=FALSE), q2)
    ci_par[["min_scale"]] <- stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q1)
    ci_par[["max_scale"]] <- stats::quantile(unlist(lapply(sim_par,function(x)x[2]),use.names=FALSE), q2)
    ci_par[["min_shape"]] <- stats::quantile(unlist(lapply(sim_par,function(x)x[3]),use.names=FALSE), q1)
    ci_par[["max_shape"]] <- stats::quantile(unlist(lapply(sim_par,function(x)x[3]),use.names=FALSE), q2)

  }

  if (name_marg == "GPD") {
    for (k in 1:n_sample) {
      sim_data <- evd::rgpd(nb_pt, loc = thresh_GPD, scale = par[1], shape = par[2])
      test_fit <- tryCatch(
        {
          (fitted <- ismev::gpd.fit(sim_data, threshold = thresh_GPD))
        },
        error = function(e) {}
      )

      if (is.null(test_fit)) {
        (fitted <- ismev::gpd.fit(sim_data, threshold = thresh_GPD, std.err = FALSE, show = FALSE))
      } else {
        fitted <- test_fit
      }
      sim_par[[k]] = fitted$mle
    }
    # Calcul des intervalles de confiance pour les paramètres
    ci_par[["min_scale"]] <- stats::quantile(unlist(lapply(sim_par, function(x) x[1]), use.names = FALSE), q1)
    ci_par[["max_scale"]] <- stats::quantile(unlist(lapply(sim_par, function(x) x[1]), use.names = FALSE), q2)
    ci_par[["min_shape"]] <- stats::quantile(unlist(lapply(sim_par, function(x) x[2]), use.names = FALSE), q1)
    ci_par[["max_shape"]] <- stats::quantile(unlist(lapply(sim_par, function(x) x[2]), use.names = FALSE), q2)
  }

  return(ci_par)
}

#' @title Uncertainty of the marginal parameters ("ml")
#'
#' @description
#' Give the confidence interval of the marginal parameters.
#'
#' @param name_marg character string specifying the name of the distribution family.
#' @param data_x numeric vector.
#' @param ci_level the confidence level required.
#' @param thresh_GPD threshold value to be used for GPD distribution. Defaults to NaN.
#'
#' @return vector of 3 values : initial parameters (input) and the confidence interval associated to this value estimated with the confint() function from the package 'stats'.
#'
determine_ci_par_marg_confint <- function(name_marg,data_x,ci_level,thresh_GPD=NaN){
  ci_par=list()

  if (name_marg == "log-normal"){
    fitted <- MASS::fitdistr(data_x, "log-normal")
    meanlog_hat <- fitted$estimate[1]
    sdlog_hat <- fitted$estimate[2]
    confint_lognorm <- stats::confint(fitted, level = ci_level)
    ci_par[["min_meanlog"]] <- confint_lognorm[1]
    ci_par[["max_meanlog"]] <- confint_lognorm[3]
    ci_par[["min_sdlog"]] <- confint_lognorm[2]
    ci_par[["max_sdlog"]] <- confint_lognorm[4]
  }

  if (name_marg == "gamma"){
    fitted <- MASS::fitdistr(data_x,"gamma")
    shape_hat <- fitted$estimate[1]
    rate_hat <- fitted$estimate[2]
    confint_gamma <- stats::confint(fitted,level = ci_level)
    ci_par[["min_shape"]] <- confint_gamma[1]
    ci_par[["max_shape"]] <- confint_gamma[3]
    ci_par[["min_rate"]] <- confint_gamma[2]
    ci_par[["max_rate"]] <- confint_gamma[4]
  }

  if (name_marg == "normal") {
    fitted <- MASS::fitdistr(data_x, "normal")
    mean_hat <- fitted$estimate[1]
    sd_hat <- fitted$estimate[2]
    confint_gauss <- stats::confint(fitted, level = ci_level)
    ci_par[["min_mean"]] <- confint_gauss[1]
    ci_par[["max_mean"]] <- confint_gauss[3]
    ci_par[["min_sd"]] <- confint_gauss[2]
    ci_par[["max_sd"]] <- confint_gauss[4]
  }

  if (name_marg == "GEV") {
    test_fit <- tryCatch(
      {
        (fitted <- evd::fgev(data_x))
      },
      error = function(e) {
        NULL
      }
    )
    if (is.null(test_fit)) {
      gev_fit <- ismev::gev.fit(data_x)
      loc_hat <- gev_fit$mle[1]
      scale_hat <- gev_fit$mle[2]
      shape_hat <- gev_fit$mle[3]

      loc_se <- gev_fit$se[1]
      scale_se <- gev_fit$se[2]
      shape_se <- gev_fit$se[3]
      z_value <- stats::qnorm(1 - (1 - ci_level) / 2)  # 1.96 pour 95%

      loc_ci <- loc_hat + c(-1, 1) * z_value * loc_se
      scale_ci <- scale_hat + c(-1, 1) * z_value * scale_se
      shape_ci <- shape_hat + c(-1, 1) * z_value * shape_se

      ci_par[["min_loc"]] <- loc_ci[1]
      ci_par[["max_loc"]] <- loc_ci[2]
      ci_par[["min_scale"]] <- scale_ci[1]
      ci_par[["max_scale"]] <- scale_ci[2]
      ci_par[["min_shape"]] <- shape_ci[1]
      ci_par[["max_shape"]] <- shape_ci[2]

    } else {
      fitted <- test_fit
      confint_gev <- stats::confint(fitted, level = ci_level)

      ci_par[["min_loc"]] <- confint_gev[1, 1]
      ci_par[["max_loc"]] <- confint_gev[1, 2]
      ci_par[["min_scale"]] <- confint_gev[2, 1]
      ci_par[["max_scale"]] <- confint_gev[2, 2]
      ci_par[["min_shape"]] <- confint_gev[3, 1]
      ci_par[["max_shape"]] <- confint_gev[3, 2]
    }
  }

  if (name_marg == "GPD") {
    test_fit <- tryCatch(
      {
        (fitted <- evd::fpot(data_x, thresh_GPD))
      },
      error = function(e) {
        NULL }
    )
    if (is.null(test_fit)) {
      gpd_fit <- ismev::gpd.fit(data_x, threshold = thresh_GPD, std.err = TRUE, show = FALSE)
      scale_hat <- gpd_fit$mle[1]
      shape_hat <- gpd_fit$mle[2]

      scale_se <- gpd_fit$se[1]
      shape_se <- gpd_fit$se[2]

      z_value <- stats::qnorm(1 - (1 - ci_level) / 2)  # 1.96 pour 95%

      scale_ci <- scale_hat + c(-1, 1) * z_value * scale_se
      shape_ci <- shape_hat + c(-1, 1) * z_value * shape_se

      ci_par[["min_scale"]] <- scale_ci[1]
      ci_par[["max_scale"]] <- scale_ci[2]
      ci_par[["min_shape"]] <- shape_ci[1]
      ci_par[["max_shape"]] <- shape_ci[2]

    } else {
      fitted <- test_fit
      confint_gpd <- stats::confint(fitted, level = ci_level)

      ci_par[["min_scale"]] <- confint_gpd[1, 1]
      ci_par[["max_scale"]] <- confint_gpd[1, 2]
      ci_par[["min_shape"]] <- confint_gpd[2, 1]
      ci_par[["max_shape"]] <- confint_gpd[2, 2]
    }

  }

  return(ci_par)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Univariate probabilities estimated through bootstrap
#'
#' @description
#' Generate a bootstrap sample of probabilities
#'
#' @param name_marg character string specifying the name of the distribution family.
#' @param ci_par list of the estimated parameters and their confidence interval
#' @param quantile quantile value.
#' @param B_m number of samples
#' @param ci_level the confidence level required.
#' @param thresh_GPD threshold value used for the GPD distribution.
#'
#' @return vector of dimension B_m containing probabilities estimated through bootstrap sampling.
#'
boot_prob_marg <- function(name_marg, ci_par,quantile, B_m, ci_level, thresh_GPD = FALSE) {
  res_boot_p_of_x <- c()

  if (name_marg == "log-normal") {
    for (b in 1:B_m){
      b_meanlog <- stats::runif(1, ci_par[["min_meanlog"]], ci_par[["max_meanlog"]])
      b_sdlog <- stats::runif(1, ci_par[["min_sdlog"]], ci_par[["max_sdlog"]])
      res_boot_p_of_x[b] <- stats::plnorm(quantile, meanlog = b_meanlog, sdlog = b_sdlog)
    }
  }

  if (name_marg == "gamma") {
    for (b in 1:B_m){
      b_shape <- stats::runif(1, ci_par[["min_shape"]], ci_par[["max_shape"]])
      b_rate <- stats::runif(1, ci_par[["min_rate"]], ci_par[["max_rate"]])
      res_boot_p_of_x[b] <- stats::pgamma(quantile, shape = b_shape, rate = b_rate)
    }
  }

  if (name_marg == "normal") {
    for (b in 1:B_m) {
      b_mean <- stats::runif(1, ci_par[["min_mean"]], ci_par[["max_mean"]])
      b_sd <- stats::runif(1, ci_par[["min_sd"]], ci_par[["max_sd"]])
      res_boot_p_of_x[b] <- stats::pnorm(quantile, mean = b_mean, sd = b_sd)
    }
  }

  if (name_marg == "GEV") {
    for (b in 1:B_m) {
      b_loc <- stats::runif(1, ci_par[["min_loc"]], ci_par[["max_loc"]])
      b_scale <- stats::runif(1, ci_par[["min_scale"]], ci_par[["max_scale"]])
      b_shape <- stats::runif(1, ci_par[["min_shape"]], ci_par[["max_shape"]])
      res_boot_p_of_x[b] <- evd::pgev(quantile, loc = b_loc, scale = b_scale, shape = b_shape)
    }
  }

  if (name_marg == "GPD") {
    for (b in 1:B_m) {
      b_scale <- stats::runif(1, ci_par[["min_scale"]], ci_par[["max_scale"]])
      b_shape <- stats::runif(1, ci_par[["min_shape"]], ci_par[["max_shape"]])
      res_boot_p_of_x[b] <- evd::pgpd(quantile,
                                      loc = thresh_GPD,
                                      scale = b_scale, shape = b_shape
      )
    }
  }

  return(res_boot_p_of_x)
}

# from Hofert 2011:  Likelihood inference for Archimedean copulas
# theta : copula object.
# acop : an object of class acopula representing an Archimedean copula family.
# the loglikelihood value.
LogL <- function(theta, acop, u, n.MC=0, ...) {
  sum(acop@dacopula(u, theta, n.MC=n.MC, log=TRUE, ...))
}

