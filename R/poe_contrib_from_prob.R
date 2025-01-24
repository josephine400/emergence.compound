#' @title Detect emergence (poe and toe) and quantify the contributions of each statistical component (marginals and dependence)
#'
#' @description Give periods of emergence (PoE) features (e.g., the number of PoEs, their duration). The signal can either emerge below the lower bound of the natural variability (PoE-low) or above the upper bound of the natural variability (PoE-up).
#' Give also the values of the statistical drivers associated with PoE-low and PoE-up.
#'
#' @param res_5prob output of the function compute_5prob() giving CE probabilities when one, two or the three statistical components (marginal 1, marginal 2 and the dependence) vary while the others are held constant.
#' @param res_ci output of the function determine_uncertainty_init() giving the probability associated with the natural variability.
#' @param y_start first year of the whole studied period.
#' @param y_end last year of the whole studied period.
#' @param length_sliding_window number of years of the sliding window.
#' @param step_sliding_window number of years separating two sliding windows.
#'
#' @return a list containing the features associated with periods of emergence (poe) and the values of each statistical's component contribution.
#'
poe_contrib_from_prob <- function(res_5prob, res_ci, y_start, y_end, length_sliding_window, step_sliding_window) {
  res <- list()

  low_bound = res_ci[2]
  up_bound = res_ci[3]

  #### poe and toe ####
  poe_margdep <- detect_poe(res_5prob$prob_margdep, low_bound, up_bound,y_start, y_end, length_sliding_window,step_sliding_window)
  poe_marg <- detect_poe(res_5prob$prob_marg, low_bound, up_bound,y_start, y_end, length_sliding_window,step_sliding_window)
  poe_marg1 <- detect_poe(res_5prob$prob_marg1, low_bound, up_bound,y_start, y_end, length_sliding_window,step_sliding_window)
  poe_marg2 <- detect_poe(res_5prob$prob_marg2, low_bound, up_bound,y_start, y_end, length_sliding_window,step_sliding_window)
  poe_dep <- detect_poe(res_5prob$prob_dep, low_bound, up_bound,y_start, y_end, length_sliding_window,step_sliding_window)


  if (step_sliding_window>1){
    print("step_sliding_window>1 : the output gives information on the emergence of the signal (but does not refer to poe)")
    for (margdep_or_marg in c("margdep","marg","marg1","marg2","dep")){
      if(margdep_or_marg=="margdep"){
        poe=poe_margdep
      }
      if (margdep_or_marg=="marg"){
        poe=poe_marg
      }
      if (margdep_or_marg=="marg1"){
        poe=poe_marg1
      }
      if (margdep_or_marg=="marg2"){
        poe=poe_marg2
      }
      if (margdep_or_marg=="dep"){
        poe=poe_dep
      }
      res[["emergence"]][[margdep_or_marg]][["upper"]][["nb_group"]] <- poe$nb_poe_above
      res[["emergence"]][[margdep_or_marg]][["upper"]][["nb_periods_per_group"]] <- poe$length_poe_above
      res[["emergence"]][[margdep_or_marg]][["upper"]][["first_period_per_group"]] <- poe$year_start_poe_above
      res[["emergence"]][[margdep_or_marg]][["upper"]][["toe"]] <- poe$toe_above
      res[["emergence"]][[margdep_or_marg]][["lower"]][["nb_group_of_consec_periods"]] <- poe$nb_poe_below
      res[["emergence"]][[margdep_or_marg]][["lower"]][["nb_periods_per_group"]] <- poe$length_poe_below
      res[["emergence"]][[margdep_or_marg]][["lower"]][["first_period_per_group"]] <- poe$year_start_poe_below
      res[["emergence"]][[margdep_or_marg]][["lower"]][["toe"]] <- poe$toe_below
    }
  }

  if (step_sliding_window==1){
    print("step_sliding_window=1 : the output gives poe features")
    for (margdep_or_marg in c("margdep","marg","marg1","marg2","dep")){
      if(margdep_or_marg=="margdep"){
        poe=poe_margdep
      }
      if (margdep_or_marg=="marg"){
        poe=poe_marg
      }
      if (margdep_or_marg=="marg1"){
        poe=poe_marg1
      }
      if (margdep_or_marg=="marg2"){
        poe=poe_marg2
      }
      if (margdep_or_marg=="dep"){
        poe=poe_dep
      }
      res[["poe"]][[margdep_or_marg]][["upper"]][["number"]] <- poe$nb_poe_above
      res[["poe"]][[margdep_or_marg]][["upper"]][["length"]] <- poe$length_poe_above
      res[["poe"]][[margdep_or_marg]][["upper"]][["starting_year"]] <- poe$year_start_poe_above
      res[["poe"]][[margdep_or_marg]][["upper"]][["toe"]] <- poe$toe_above
      res[["poe"]][[margdep_or_marg]][["lower"]][["number"]] <- poe$nb_poe_below
      res[["poe"]][[margdep_or_marg]][["lower"]][["length"]] <- poe$length_poe_below
      res[["poe"]][[margdep_or_marg]][["lower"]][["starting_year"]] <- poe$year_start_poe_below
      res[["poe"]][[margdep_or_marg]][["lower"]][["toe"]] <- poe$toe_below
    }
  }

  #### contrib ####
  n_time=length(unlist(res_5prob$prob_margdep))
  p0 <- res_5prob$prob_margdep[[1]]

  contrib_marg <- contribution(unlist(res_5prob$prob_marg),p0,unlist(res_5prob$prob_margdep))
  contrib_marg1 <- contribution(unlist(res_5prob$prob_marg1),p0,unlist(res_5prob$prob_margdep))
  contrib_marg2 <- contribution(unlist(res_5prob$prob_marg2),p0,unlist(res_5prob$prob_margdep))
  contrib_dep <- contribution(unlist(res_5prob$prob_dep),p0,unlist(res_5prob$prob_margdep))

  for (upper_or_lower in c("upper", "lower")) {
    if (upper_or_lower == "upper") {
      ind_poe <- (poe_margdep$ind_start_poe_above)
      length_poe <- (poe_margdep$length_poe_above)
      y_toe <- poe_margdep$toe_above
    }
    if (upper_or_lower == "lower") {
      ind_poe <- (poe_margdep$ind_start_poe_below)
      length_poe <- (poe_margdep$length_poe_below)
      y_toe <- poe_margdep$toe_below
    }

    if (length(ind_poe) == 1) {
      contrib_marg_poe <- mean(contrib_marg[ind_poe:(ind_poe + length_poe - 1)])
      contrib_dep_poe <- mean(contrib_dep[ind_poe:(ind_poe + length_poe - 1)])
      contrib_marg1_poe <- mean(contrib_marg1[ind_poe:(ind_poe + length_poe - 1)])
      contrib_marg2_poe <- mean(contrib_marg2[ind_poe:(ind_poe + length_poe - 1)])
    } else if (length(ind_poe) > 1) {
      margdep_list <- c()
      marg_list <- c()
      dep_list <- c()
      int_list <- c()
      marg1_list <- c()
      marg2_list <- c()
      int_marg_list <- c()
      for (k in 1:length(ind_poe)) {
        marg_list <- c(marg_list, mean(contrib_marg[ind_poe[k]:(ind_poe[k] + length_poe[k] - 1)]))
        dep_list <- c(dep_list, mean(contrib_dep[ind_poe[k]:(ind_poe[k] + length_poe[k] - 1)]))
        marg1_list <- c(marg1_list, mean(contrib_marg1[ind_poe[k]:(ind_poe[k] + length_poe[k] - 1)]))
        marg2_list <- c(marg2_list, mean(contrib_marg2[ind_poe[k]:(ind_poe[k] + length_poe[k] - 1)]))
      }
      contrib_marg_poe <- mean(marg_list)
      contrib_dep_poe <- mean(dep_list)
      contrib_marg1_poe <- mean(marg1_list)
      contrib_marg2_poe <- mean(marg2_list)

    } else {
      contrib_marg_poe <- NA
      contrib_dep_poe <- NA
      contrib_marg1_poe <- NA
      contrib_marg2_poe <- NA
    }

    res[["contrib"]][["poe_ipoe"]][[upper_or_lower]][["marginal1&2"]] <- contrib_marg_poe
    res[["contrib"]][["poe_ipoe"]][[upper_or_lower]][["marginal1"]] <- contrib_marg1_poe
    res[["contrib"]][["poe_ipoe"]][[upper_or_lower]][["marginal2"]] <- contrib_marg2_poe
    res[["contrib"]][["poe_ipoe"]][[upper_or_lower]][["dep"]] <- contrib_dep_poe

    if ((!is.na(y_toe)) & (length(ind_poe) > 1)){
      # if ToE and several PoE -> we can compute contrib of ipoe and poe
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal1&2"]] <- marg_list[length(marg_list)]
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal1"]] <- marg1_list[length(marg1_list)]
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal2"]] <- marg2_list[length(marg2_list)]
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["dep"]] <- dep_list[length(dep_list)]

      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal1&2"]] <- mean(marg_list[1:(length(marg_list)-1)])
      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal1"]] <- mean(marg1_list[1:(length(marg1_list)-1)])
      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal2"]] <- mean(marg2_list[1:(length(marg2_list)-1)])
      res[["contrib"]][["poe"]][[upper_or_lower]][["dep"]] <- mean(dep_list[1:(length(dep_list)-1)])
    }

    else if ((!is.na(y_toe)) & (length(ind_poe) == 1)){
      # if ToE and one PoE -> there is no contrib of poe, only contrib of ipoe
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal1&2"]] <- contrib_marg_poe
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal1"]] <- contrib_marg1_poe
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal2"]] <- contrib_marg2_poe
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["dep"]] <- contrib_dep_poe

      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal1&2"]] <- NA
      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal1"]] <- NA
      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal2"]] <- NA
      res[["contrib"]][["poe"]][[upper_or_lower]][["dep"]] <- NA
    }
    else {
      # if there is no ToE -> there is no contrib of ipoe, only poe
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal1&2"]] <- NA
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal1"]] <- NA
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["marginal2"]] <- NA
      res[["contrib"]][["ipoe"]][[upper_or_lower]][["dep"]] <- NA

      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal1&2"]] <- contrib_marg_poe
      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal1"]] <- contrib_marg1_poe
      res[["contrib"]][["poe"]][[upper_or_lower]][["marginal2"]] <- contrib_marg2_poe
      res[["contrib"]][["poe"]][[upper_or_lower]][["dep"]] <- contrib_dep_poe
    }

  }

  return(res)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Contribution metric
#'
#' @description Give the percentage of each statistical component (marginal1, marginal2, dependence) change to the overall signal change.
#'
#' @param proba CE probability when one or two statistical components evolve.
#' @param p0 CE probability during the reference period.
#' @param prob_margdep CE probability when the three statistical components evolve.
#'
#' @return the contribution value.
#'
contribution <- function(proba,p0,prob_margdep){
  return ((proba - p0)/(prob_margdep - p0)*100)
}

#' @title Periods of emergence (PoE) detection
#'
#' @description Detect PoEs, periods during which the signal exceeds the natural variability, either below the lower bound or above the upper bound.
#'
#' @param proba CE probability when one, two or the three statistical components evolve.
#' @param low_bound lower bound of the natural variability (confidence interval of the probability estimated during the reference period).
#' @param up_bound upper bound of the natural variability (confidence interval of the probability estimated during the reference period).
#' @param y_start first year of the whole studied period.
#' @param y_end last year of the whole studied period.
#' @param length_sliding_window number of years of the sliding window.
#' @param step_sliding_window number of years separating two sliding windows.
#'
#' @return a list containing PoE features (number of PoEs, their duration, their starting year, the potential ToE).
#'
detect_poe <- function(proba, low_bound, up_bound, y_start, y_end, length_sliding_window, step_sliding_window) {

  if (length(proba)==1){
    stop("the signal is a single value (step_sliding_window is probably too high) - no emergence can be detected")
  }

  length_poe_below <- 0
  length_poe_above <- 0
  length_poe_below_all <- c()
  length_poe_above_all <- c()
  ind_poe_below <- c()
  ind_poe_above <- c()
  ind_start_poe_below <- c()
  ind_start_poe_above <- c()
  for (k in 1:length(proba)) {
    if (proba[k] < low_bound) {
      length_poe_below <- length_poe_below + 1
      ind_poe_below <- c(ind_poe_below, k)
    } else if (proba[k] > up_bound) {
      length_poe_above <- length_poe_above + 1
      ind_poe_above <- c(ind_poe_above, k)
    } else {
      if (length_poe_below != 0) {
        length_poe_below_all <- c(length_poe_below_all, length_poe_below)
        ind_start_poe_below <- c(ind_start_poe_below, ind_poe_below[1])
        length_poe_below <- 0
      }
      if (length_poe_above != 0) {
        length_poe_above_all <- c(length_poe_above_all, length_poe_above)
        ind_start_poe_above <- c(ind_start_poe_above, ind_poe_above[1])
        length_poe_above <- 0
      }
      ind_poe_below <- c()
      ind_poe_above <- c()
    }
  }

  if (length_poe_below != 0) {
    length_poe_below_all <- c(length_poe_below_all, length_poe_below)
    ind_start_poe_below <- c(ind_start_poe_below, ind_poe_below[1])
  }
  if (length_poe_above != 0) {
    length_poe_above_all <- c(length_poe_above_all, length_poe_above)
    ind_start_poe_above <- c(ind_start_poe_above, ind_poe_above[1])
  }

  year_start_poe_below = y_start + (ind_start_poe_below-1)*step_sliding_window + floor(length_sliding_window/2)
  year_start_poe_above = y_start + (ind_start_poe_above-1)*step_sliding_window + floor(length_sliding_window/2)


  #### toe ####
  last_start_y_below <- utils::tail(year_start_poe_below, 1)
  last_start_y_above <- utils::tail(year_start_poe_above, 1)
  last_length_below <- utils::tail(length_poe_below_all, 1)
  last_length_above <- utils::tail(length_poe_above_all, 1)
  last_year = middle_sliding_window(paste0(y_end-length_sliding_window+1,"_",y_end))

  if ((length(year_start_poe_below) != 0) && (last_start_y_below + last_length_below - 1) == last_year) {
    toe_below <- last_start_y_below
  } else {
    toe_below <- NA
  }

  if ((length(year_start_poe_above) != 0) && (last_start_y_above + (last_length_above-1)*step_sliding_window) == last_year) {
    toe_above <- last_start_y_above
  } else {
    toe_above <- NA
  }

  return(list(
    nb_poe_below = length(length_poe_below_all), nb_poe_above = length(length_poe_above_all),
    length_poe_below = length_poe_below_all, length_poe_above = length_poe_above_all,
    length_tot_poe_below = sum(length_poe_below_all), length_tot_poe_above = sum(length_poe_above_all),
    ind_start_poe_below = ind_start_poe_below, ind_start_poe_above = ind_start_poe_above,
    year_start_poe_below = year_start_poe_below, year_start_poe_above = year_start_poe_above,
    toe_below = toe_below, toe_above = toe_above
  ))
}
