MAX_YEARS <- 200

# Finds the empirical distribution survival distributions
.fit_dist_empirical <- function(last_test, ever_test, weights, time_at_risk, aids_dist){

  # Construct empirical distribution for TID | TESTER
  tm <- last_test
  event <- !is.na(tm) & tm < time_at_risk
  event[(ever_test | is.na(ever_test)) & is.na(tm)] <- NA
  tm[!event & !is.na(event)] <- time_at_risk[!event & !is.na(event)]
  tdist <- survival::survfit(survival::Surv(tm, event) ~ 1, weights=weights)
  times <- seq(from=0,to=MAX_YEARS*12, by=1)
  inds <- findInterval(times, tdist$time, all.inside = TRUE)
  test_surv <- tdist$surv[inds]
  test_surv_cond <- (test_surv - min(test_surv)) / (1-min(test_surv))

  # Calculate tau=P(TESTER) Taking into account censoring
  et <- ever_test
  et[tm >= time_at_risk & !is.na(tm >= time_at_risk)] <- FALSE
  po <- 1 - test_surv_cond[pmin(MAX_YEARS*12-1, ceiling(time_at_risk) + 1) ]
  tau_lik <- function(tau) {
    val <- (et * log(tau * po) + (!et) * log(1-tau * po)) * weights
    val[!is.finite(val)] <- NA
    sum(val, na.rm=TRUE)
  }
  tau <- optimise(tau_lik, interval = c(0,1), maximum = TRUE)$maximum

  # Calculate competing risk survival function
  aids_surv <- 1 - aids_dist(times)
  surv <- test_surv_cond * aids_surv * tau + aids_surv * (1 - tau)

  list(survival=surv,
       times=times,
       tau=tau,
       test_surv_cond=test_surv_cond,
       aids_surv=aids_surv)
}

# Fits a weibull distribution to possibly windowed data
.fit_dist_weibull <- function(last_test, last_upper, ever_test, weights, time_at_risk, aids_dist, initial_rate=1/14){

  exact_time <- all(na.omit(last_test == last_upper))

  miss_test <- is.na(last_test) | is.na(last_upper)
  time_at_risk[time_at_risk == 0] <- 1
  event <- ever_test
  censored <- last_test > time_at_risk & !is.na(last_test) | is.na(time_at_risk)
  event[censored] <- FALSE
  last_upper[censored] <- last_test[censored] <- NA

  if(!exact_time){
    last_upper <- ifelse(!is.na(time_at_risk) & !is.na(last_upper) & time_at_risk < last_upper, time_at_risk, last_upper)
    last_upper[!is.na(last_upper) & last_upper <= last_test] <- last_test[!is.na(last_upper) & last_upper <= last_test] + 1
  }

  # Fit a weibull distribtion to get P(TID | TESTER)
  weibull_lik <- function(scale, k, tau){
    p <- function(x) suppressWarnings(pweibull(x, scale=scale,shape=k))
    log_d <- function(x) suppressWarnings(dweibull(x, scale=scale,shape=k, log = TRUE))
    if(exact_time){
      lik <- sum((log_d(last_test[event])) * weights[event], na.rm=TRUE)
    }else{
      lik <- sum((log(p(last_upper[event]) - p(last_test[event]))) * weights[event], na.rm=TRUE)
    }
    nna <- is.na(last_test) & !is.na(event)
    lik <- lik + sum(log(p(time_at_risk[nna & event])) * weights[nna & event], na.rm=TRUE)
    lik <- lik + sum(log(tau) * weights[event], na.rm=TRUE)
    lik <- lik + sum( log((1 - p(time_at_risk[!event])) * tau + 1 - tau) * weights[!event], na.rm=TRUE )
    lik <- max(-.Machine$double.xmax/10, min(.Machine$double.xmax/10, lik))
    -lik
  }

  opt <- optim(function(x)weibull_lik(x[1],x[2], x[3]),
               par = c(1/initial_rate,1, .5),
               lower=c(0,0,0) + .0000000001,
               upper=c(Inf, Inf, 1 - .0000000001),
               method="L-BFGS-B")
  times <- seq(from=0,to=MAX_YEARS*12, by=1)
  test_surv_cond <- 1 - pweibull(times, scale=opt$par[1],shape=opt$par[2])
  tau <- opt$par[3]

  # Calculate competing risk survival function
  aids_surv <- 1 - aids_dist(times)
  surv <- test_surv_cond * aids_surv * tau + aids_surv * (1 - tau)

  list(survival=surv,
       times=times,
       tau=tau,
       test_surv_cond=test_surv_cond,
       aids_surv=aids_surv)
}

.assert <- function(..., msg=NULL){
  if(is.null(msg)){
    msg <- paste(deparse(substitute(...)),"is not true.")
  }
  tr <- try(stopifnot(...))
  if(inherits(tr, "try-error"))
    stop(msg)

}

#checks testing_incidence parameters and throws errors if they are invalid
.check_params <- function(report_pos, biomarker_art, low_viral, hiv,
                          ever_test, last_test, last_upper,
                          age, testing_debut_age,
                          weights,distribution, test_pop,
                          ptruth, ptreated,
                          aids_dist,
                          age_breaks, subset, uniform_missreport){
  n <- length(report_pos)
  if(n < 3)
    stop("Too few observations to perform estimation")
  n1 <- length(biomarker_art)
  .assert(n == n1, msg = paste0("biomarker_art is of length ", n1, " expecting ", n))
  n1 <- length(low_viral)
  .assert(n == n1, msg = paste0("low_viral is of length ", n1, " expecting ", n))
  n1 <- length(hiv)
  .assert(n == n1, msg = paste0("hiv is of length ", n1, " expecting ", n))
  n1 <- length(ever_test)
  .assert(n == n1, msg = paste0("ever_test is of length ", n1, " expecting ", n))
  n1 <- length(last_test)
  .assert(n == n1, msg = paste0("last_test is of length ", n1, " expecting ", n))
  n1 <- length(last_upper)
  .assert(n == n1, msg = paste0("last_upper is of length ", n1, " expecting ", n))
  if(!is.null(age)){
    n1 <- length(age)
    .assert(n == n1, msg = paste0("age is of length ", n1, " expecting ", n))
  }
  if(!is.null(subset)){
    .assert(is.logical(subset), length(subset) == n)
    if(uniform_missreport && (is.null(ptruth) || is.null(ptreated))){
      has_trt <- sum((low_viral | biomarker_art) * weights, na.rm=TRUE) > 0
      .assert(has_trt, msg="Unable to estimate missreporting as no treated cases are present")
    }
    subset[is.na(subset)] <- FALSE
    report_pos <- report_pos[subset]
    biomarker_art <- biomarker_art[subset]
    low_viral <- low_viral[subset]
    hiv <- hiv[subset]
    ever_test <- ever_test[subset]
    last_test <- last_test[subset]
    last_upper <- last_upper[subset]
    age <- age[subset]
    weights <- weights[subset]
    if(!uniform_missreport && (is.null(ptruth) || is.null(ptreated))){
      has_trt <- sum((low_viral | biomarker_art) * weights, na.rm=TRUE) > 0
      .assert(has_trt, msg="Unable to estimate missreporting as no treated cases are present. Try using uniform_missreport=TRUE.")
    }
  }else{
    if(is.null(ptruth) || is.null(ptreated)){
      has_trt <- sum((low_viral | biomarker_art) * weights, na.rm=TRUE) > 0
      .assert(has_trt, msg="Unable to estimate missreporting as no treated cases are present.")
    }
  }
  tbl <- wtd.table(report_pos, weights=weights)
  .assert(length(tbl) == 2 & all(tbl > 0), msg="report_pos should have two unique values.")

  tbl <- wtd.table(hiv, weights=weights)
  .assert(length(tbl) == 2 & all(tbl > 0), msg="hiv should have two unique values.")

  tbl <- wtd.table(ever_test, weights=weights)
  .assert(length(tbl) == 2 & all(tbl > 0), msg="ever_test should have two unique values.")

  .assert(all(na.omit(last_upper >= last_test)))
  .assert(all(na.omit(age >= testing_debut_age)))
  .assert(all(na.omit(age < 100)))
  .assert(distribution != "empirical" || testing_debut_age == 0, msg = "The empirical distribution may only be used when testing_debut_age=0")

  if(!is.null(age_breaks)){
    .assert(age_breaks > min(age, na.rm=TRUE), msg= "age_breaks: break points must all be greater than the smallest age value")
    .assert(age_breaks < max(age, na.rm=TRUE), msg= "age_breaks: break points must all be smaller than the largest age value")
  }
  TRUE
}

#' Incindence from testing history
#'
#' @param report_pos A logical vector indicating whether each subject reported a positive hiv status
#' @param biomarker_art A logical vector indicating whether ART antibodies are present. NA if test not done.
#' @param low_viral A logical vector indicating whether viral load is <=1000
#' @param hiv A logical vector indicating hiv status
#' @param ever_test A logical vector indicating whether the subject had eer had an hiv test.
#' @param last_test A numeric vector indicating the time since last hiv test in months. If testing times are binned into buckets, this is the lower bound of the months since last hiv test.
#' @param age A numeric vector indicating the age of each subject in years.
#' @param testing_debut_age The age at which individuals begin engaging in regular testing.
#' @param last_upper A numeric vector indicating the upper bound of the months since last hiv test for each individual.
#' @param weights Survey weights
#' @param distribution Either "empirical", or "weibull." This controls the family of distribution used to
#' model time since last test. "empirical" may not be used with binned testing times.
#' @param test_pop If "negative', the time since last negative is calculated amoung the HIV- population, otherwise it is calculated of those who report being undiagnosed.
#' @param ptruth The proportion of the diagnosed hiv positive population that would report being hiv positive.
#' If NULL, this is estimated using biomarker_art and low_viral.
#' @param ptreated The proportion of hiv positive individuals with postive biomarker_art or low_viral.
#' @param age_breaks age stratification break points.
#' If NULL, this is estimated from the data.
#' @param aids_dist The distribution function of time from infection to aids in months.
#' @param uniform_missreport If true, the rate of missreporting of undiagnosed status is conisuered uniform over the age strata.
#' @param subset An optional vector specifying a subset of observations on which to perform the analysis. If uniform_missreport is TRUE, ptruth and ptreated are calculated over the full data.
#' @details
#' When time since last test is grouped into bins, last_upper should always be greater than last_test,
#'  and may be infinite (e.g. test was  > 24 months ago). Those with missing data or who never had an hiv test
#' should be assigned NA for both their lower and upper bound.
#' @return A data.frame of class test_inc with elements:
#' 'incidence': the estimated incidence.
#' 'transmission rate': the estimated transmission rate.
#' 'pundiag': The estimated proportion of positive cases that are undiagnosed.
#' 'psay_undiag': The proportion of positive cases that report being undiagnosed.
#' 'pmiss_class': the proportion whose diagnosis status is incorrectly reported by the individual and are observed as treated due to viral load or art biomarkers.
#' 'phiv': The proportion with a positive diagnosis.
#' 'ptester': The proportion who have ever been tested.
#' 'mean_time_since_last_test': the mean tie since last test in years.
#' 'tid': mean time between infection and diagnosis in years.
#' 'ptruth' the proportion of positive individuals that correctly report their status.
#' 'ptreated': the proportion of positive indivudals who are identified as treated by viral load or biomarker.
#'
#' @examples
#'   data(tstdat)
#'   tstdat$age <- rep(15:64, 200)
#'
#'
#'   inc <- with(tstdat, testing_incidence(report_pos, biomarker_art, low_viral, hiv,
#'                                         ever_test, last_test))
#'   inc
#'
#'   # Using last test times divided into bins
#'   tstdat$last_test_lower <- 24
#'   tstdat$last_test_upper <- Inf
#'   tstdat$last_test_lower[tstdat$last_test < 6] <- 0
#'   tstdat$last_test_upper[tstdat$last_test < 6] <- 6
#'   tstdat$last_test_lower[tstdat$last_test >= 6 & tstdat$last_test < 12] <- 6
#'   tstdat$last_test_upper[tstdat$last_test >= 6 & tstdat$last_test < 12] <- 12
#'   tstdat$last_test_lower[tstdat$last_test >= 12 & tstdat$last_test < 24] <- 12
#'   tstdat$last_test_upper[tstdat$last_test >= 12 & tstdat$last_test < 24] <- 24
#'   inc <- with(tstdat, testing_incidence(report_pos, biomarker_art, low_viral, hiv,
#'                                         ever_test, last_test_lower, last_test_upper))
#'   inc
#'
#'   # HIV testing starts at age 13
#'   inc <- with(tstdat, testing_incidence(report_pos, biomarker_art, low_viral, hiv,
#'                                         ever_test, last_test,
#'                                         age=age,testing_debut_age=13))
#'   inc
#'
#'   # Stratify by age
#'   inc <- with(tstdat, testing_incidence(report_pos, biomarker_art, low_viral, hiv,
#'                                         ever_test, last_test,
#'                                         age=age,testing_debut_age=13,age_breaks=c(25,35,45,55)))
#'   inc
#'
#'   # Pooled estimate probability of miss-reporting diagnosis status across groups
#'   inc <- with(tstdat, testing_incidence(report_pos, biomarker_art, low_viral, hiv,
#'                                         ever_test, last_test,
#'                                         age=age,testing_debut_age=13,
#'                                         age_breaks=c(25,35,45,55),
#'                                         uniform_missreport=TRUE))
#'   inc
#' @export
testing_incidence <- function(report_pos, biomarker_art, low_viral, hiv,
                              ever_test, last_test, last_upper = last_test,
                              age=NULL, testing_debut_age=0,
                              weights=rep(1, length(report_pos)) / length(report_pos),
                              distribution=c("weibull","empirical"), test_pop=c("negative","undiagnosed"),
                              ptruth=NULL, ptreated=NULL,
                              aids_dist = function(x) pweibull(x / 12, scale=1/0.086, shape=2.516),
                              age_breaks=NULL,
                              subset=NULL,
                              uniform_missreport=FALSE){

  test_pop <- match.arg(test_pop)
  distribution <- match.arg(distribution)
  if(!is.null(age_breaks))
    age_breaks <- sort(age_breaks)

  .check_params(report_pos, biomarker_art, low_viral, hiv,
                ever_test, last_test, last_upper,
                age, testing_debut_age,
                weights,distribution, test_pop,
                ptruth, ptreated,
                aids_dist,
                age_breaks, subset, uniform_missreport)

  # Function Parameters
  params <- list(report_pos=report_pos, biomarker_art=biomarker_art, low_viral=low_viral, hiv=hiv,
                 ever_test=ever_test, last_test=last_test, last_upper=last_upper,
                 age=age, testing_debut_age=testing_debut_age,
                 weights=weights,
                 distribution=distribution, test_pop=test_pop,
                 ptruth=ptruth, ptreated=ptreated,
                 aids_dist=aids_dist,
                 age_breaks=age_breaks,
                 subset=subset,
                 uniform_missreport=uniform_missreport)

  treated <- low_viral | biomarker_art
  treated[is.na(treated)] <- FALSE

  if(!is.null(age_breaks)){
    tlie <- wtd.table(treated, !report_pos, weights = weights)
    if(is.null(ptruth) & uniform_missreport)
      ptruth <- tlie[2,1] / sum(tlie[2,])
    if(is.null(ptreated) & uniform_missreport)
      ptreated <- tlie[2,1] / sum(tlie[,1])
    if(!is.null(subset)){
      subset[is.na(subset)] <- FALSE
      report_pos <- report_pos[subset]
      biomarker_art <- biomarker_art[subset]
      low_viral <- low_viral[subset]
      hiv <- hiv[subset]
      ever_test <- ever_test[subset]
      last_test <- last_test[subset]
      last_upper <- last_upper[subset]
      age <- age[subset]
      weights <- weights[subset]
    }
    age_breaks <- sort(age_breaks)
    sub_estimates <- list()
    lower <- c()
    upper <- c()
    mu <- c()
    k <- length(age_breaks) + 1
    for(i in 1:k){
      ua <- if(i > length(age_breaks)) Inf else age_breaks[i]
      la <- if(i == 1) testing_debut_age else age_breaks[i-1]
      sub <- !is.na(age) & age < ua & age >= la
      lower[i] <- la
      upper[i] <- ua
      mu[i] <- sum(sub * weights, na.rm=TRUE) / sum((!is.na(age)) * weights, na.rm=TRUE)
      if(sum(sub) == 0){
        stop("No observations in age group")
      }
      ti <- testing_incidence(report_pos[sub], biomarker_art[sub], low_viral[sub], hiv[sub],
                                    ever_test[sub], last_test[sub], last_upper[sub],
                                    age[sub], testing_debut_age,
                                    weights[sub],
                                    distribution, test_pop,
                                    ptruth=ptruth, ptreated=ptreated,
                                    aids_dist,
                                    age_breaks=NULL)
      sub_estimates[[i]] <- ti
    }
    ph <- sapply(sub_estimates, function(x) x$phiv)
    pu <- sapply(sub_estimates, function(x) x$pundiag)
    lower <- lower * 12
    upper <- pmin(MAX_YEARS, upper) * 12
    Q <- matrix(0,nrow=k,ncol=k)

    for(i in 1:k){
      for(j in i:k){
        if(length(sub_estimates[[i]]) == 1)
          next
        tid <- attr(sub_estimates[[i]], "survival_distribution")$survival
        cumtid <- c(0, cumsum(tid))
        asub <- floor(age[!is.na(age) & (age < upper[i] / 12) & (age >= lower[i] / 12)] * 12)
        int <- cumtid[upper[j] - asub] - cumtid[pmax(1, lower[j] - asub)]
        Q[j,i] <- mean(int)
      }
    }
    #B <- sweep(Q, 1, (1 - ph) * mu, "*")
    c <- pu * ph * mu
    incidence <- solve(Q / 12) %*% c / ((1-ph)*mu)
    transmission <- incidence * (1 - ph) / ph
    result <- do.call(rbind, sub_estimates)
    result$incidence <- incidence
    result$transmission <- transmission
    pooled_inc <- sum(result$incidence * mu * (1 - result$phiv)) / sum(mu * (1 - result$phiv))
    pooled_phiv <- sum(result$phiv * mu)
    pooled_trans <- pooled_inc * (1 - pooled_phiv) / pooled_phiv
    result <- rbind(result, c(pooled_inc, pooled_trans, rep(NA, ncol(result) - 2)))
    attr(result, "Q") <- Q
    attr(result, "c") <- c
    attr(result, "lower") <- lower
    attr(result, "upper") <- upper
    attr(result, "survival_distribution") <- lapply(sub_estimates,
                                                    function(x) attr(x, "survival_distribution"))
    upper[upper == MAX_YEARS * 12] <- Inf
    row.names(result) <- c(paste0("[", lower / 12, ", ", upper / 12, ")"), "pooled")
    attr(result, "params") <- params
    class(result) <- c("test_inc","data.frame")
    return(result)
  }

  if(!is.null(subset)){
    subset[is.na(subset)] <- FALSE
    if(uniform_missreport){
      # Adjust for missreporting of HIV status over total population
      tlie <- wtd.table(treated, !report_pos, weights = weights)
      if(is.null(ptruth))
        ptruth <- tlie[2,1] / sum(tlie[2,])
      if(is.null(ptreated))
        ptreated <- tlie[2,1] / sum(tlie[,1])
    }
    report_pos <- report_pos[subset]
    biomarker_art <- biomarker_art[subset]
    low_viral <- low_viral[subset]
    hiv <- hiv[subset]
    ever_test <- ever_test[subset]
    last_test <- last_test[subset]
    last_upper <- last_upper[subset]
    age <- age[subset]
    weights <- weights[subset]
    treated <- treated[subset]
  }

  # Adjust for missreporting of HIV status
  tlie <- wtd.table(treated, !report_pos, weights = weights)
  if(is.null(ptruth))
    ptruth <- tlie[2,1] / sum(tlie[2,])
  if(is.null(ptreated))
    ptreated <- tlie[2,1] / sum(tlie[,1])
  pmiss_class <- 1 - (ptruth + ptreated * ( 1 - ptruth))

  if(!is.null(age)){
    time_at_risk <- (age - testing_debut_age) * 12
  }else{
    time_at_risk <- rep(MAX_YEARS*12, length(last_test))
  }
  time_at_risk[is.na(time_at_risk)] <- MAX_YEARS*12

  ptester <- as.vector(prop.table(wtd.table(ever_test[!hiv],weights=weights[!hiv]))[2])

  tab <- wtd.table(!report_pos, hiv, weights=weights)
  psay_undiag <- tab[2,2] / sum(tab[,2])

  # proportion that say undiagnosed and arn't observed to be treated
  tab <- wtd.table(!report_pos & !treated, hiv, weights=weights)
  ppos_undiag <- tab[2,2] / sum(tab[,2])

  phiv <- as.vector(prop.table(wtd.table(hiv, weights=weights))[2])

  # Calculate mean time since last test
  if(test_pop == "negative")
    ln_sub <- !hiv & !is.na(hiv)
  else
    ln_sub <- !report_pos & !is.na(report_pos) & hiv & !is.na(hiv) & !treated
  if(distribution == "empirical"){
    if(!all(na.omit(last_test == last_upper))){
      stop("The empirical distribution can only be used for non-binned testing times.")
    }
    surv_dist <- .fit_dist_empirical(last_test[ln_sub], ever_test[ln_sub],
                                     weights[ln_sub], time_at_risk[ln_sub], aids_dist)
  }else{
    surv_dist <- .fit_dist_weibull(last_test[ln_sub], last_upper[ln_sub], ever_test[ln_sub],
                                     weights[ln_sub], time_at_risk[ln_sub], aids_dist)
  }
  m2 <- sum(surv_dist$test_surv_cond)
  tid <- sum(surv_dist$survival)
  ptester <- surv_dist$tau


  #pundiag <- (ppos_undiag - pmiss_class) / (1 - pmiss_class)
  pundiag <- 1 - (1 - ppos_undiag) / (ptruth + ptreated * ( 1 - ptruth))
  # Calculate incidence
  tid <- tid / 12
  inc <- pundiag * phiv / (tid * (1-phiv))

  # Calculate transmission rate
  trans <- pundiag / tid

  # Format return object
  result <- data.frame(list(incidence=inc,
                 transmission=trans,
                 pundiag=pundiag,
                 psay_undiag=psay_undiag,
                 pmiss_class=pmiss_class,
                 phiv=phiv,
                 ptester=ptester,
                 mean_time_since_last_test=m2 / 12,
                 tid=tid,
                 ptruth=ptruth,
                 ptreated=ptreated,
                 n=length(report_pos)))
  row.names(result) <- "estimate"
  attr(result, "survival_distribution") <- surv_dist
  attr(result, "params") <- params
  class(result) <- c("test_inc","data.frame")
  result
}


#' Calculates bootstrap variability estimates
#' @param incidence a test_inc object returned by testing_incidence.
#' @param nrep The number of bootstrap replications (ignored if rep_weights is not null.
#' @param rep_weights A dataframe of replicate weights.
#' @param type The type of resampling weights. See svrepdesign.
#' @param combined_weights TRUE if the rep_weights already include the sampling weights. This is usually the case.
#' @param show_progress print bootstrap progress
#' @param ... additional parameters to svrepdesign.
#' @examples
#'   data(tstdat)
#'   tstdat$age <- rep(15:64, 200)
#'
#'
#'   inc <- with(tstdat, testing_incidence(report_pos, biomarker_art, low_viral, hiv,
#'                                         ever_test, last_test))
#'
#'   # Simple random sample bootstrap with 5 replicates.
#'   bootstrap_incidence(inc,nrep=5, show_progress=FALSE)
#'
#' @export
bootstrap_incidence <- function(incidence,
                                nrep=1000,
                                rep_weights=NULL,
                                type=c("BRR", "Fay", "JK1","JKn","bootstrap","other"),
                                combined_weights=TRUE,
                                show_progress=TRUE, ...){
  type <- match.arg(type)

  params <- attr(incidence, "params")
  testing_debut_age <- params$testing_debut_age
  distribution <- params$distribution
  test_pop <- params$test_pop
  ptruth <- params$ptruth
  ptreated <- params$ptreated
  aids_dist <- params$aids_dist
  age_breaks <- params$age_breaks
  subset <- params$subset
  uniform_missreport <- params$uniform_missreport


  dat <- params[c('report_pos', 'biomarker_art', 'low_viral', 'hiv',
                                'ever_test', 'last_test', 'last_upper',
                                'age','weights','subset')]
  fun <- function(data, weights=NULL){
    if(!is.null(weights))
      data$weights <- weights
    with(data, testing_incidence(report_pos=report_pos, biomarker_art=biomarker_art,
                                 low_viral=low_viral, hiv=hiv, ever_test=ever_test,
                                 last_test=last_test, last_upper=last_upper, weights=weights,
                                 age=age, testing_debut_age=testing_debut_age, distribution = distribution,
                                 age_breaks=age_breaks,
                                 subset=subset,
                                 test_pop=test_pop, uniform_missreport = uniform_missreport))
  }
  if(is.null(rep_weights)){
    result <- .srs_bootstrap(dat, fun, nrep=nrep, show_progress=show_progress)
  }else{
    result <- .replicate_bootstrap(dat, dat$weights, rep_weights, fun, type=type,
                                  combined_weights = combined_weights, show_progress=show_progress, ...)
  }
  result
}


# Performs a bootstrap assuming a simple random sample
# @param dat a data.frame
# @param fun a function taking a dataframe as a parameter
# @param nrep the number of replicates
# @param show_progress print progress
.srs_bootstrap <- function(dat, fun, nrep=1000, show_progress=TRUE){
  value <- as.matrix(fun(dat))
  n <- max(sapply(dat, length))
  nr <- nrow(value)
  nc <- ncol(value)

  errors <- list()
  estimates <- array(NA, dim=c(nr, nc, nrep))
  for(i in 1:nrep){
    if(show_progress)
      cat(".")
    samp <- sample.int(n, n, replace=TRUE)
    boot <- lapply(dat, function(x) x[samp])
    val <- try(as.matrix(fun(boot)))
    if(!inherits(val, "try-error"))
      estimates[,,i] <- as.matrix(fun(boot))
    else
      errors[[length(errors) + 1]] <- val
  }
  if(show_progress)
    cat("\n")
  vars <- matrix(NA,nrow=nr,ncol=nc)
  for(i in 1:nr){
    for(j in 1:nc){
      if(!all(is.na(estimates[i,j,])))
        vars[i,j] <- var(estimates[i,j,], na.rm=TRUE)
    }
  }
  bb <- list(value=value,
             var = vars,
             replicates=estimates,
             nrep=nrep,
             errors=errors)
  class(bb) <- c("inc_boot_est","list")
  bb
}

# Performs bootstrap variance estimation using replicate weights
# @param dat a dataframe
# @param weights the survey wieghts
# @param rep_weights a dataframe of replicate wieghts
# @param fun The function to be applied. takes a data.frame and a vector of wieghts as parameters, and returns a matrix
# @param type The type of resampling weights. See svrepdesign.
# @param combined_weights TRUE if the rep_weights already include the sampling weights. This is usually the case.
# @param show_progress print bootstrap progress
# @param ... additional parameters to svrepdesign.
.replicate_bootstrap <- function(dat, weights, rep_weights, fun,
                                type=c("BRR", "Fay", "JK1","JKn","bootstrap","other"),
                                combined_weights=TRUE,
                                show_progress=TRUE, ...){
  values <- as.matrix(fun(dat, weights))

  type <- match.arg(type)
  not_miss <- (rowSums(is.na(rep_weights)) + is.na(weights)) < 0.5
  dat <- lapply(dat,function(x) x[not_miss])#dat[not_miss,]
  weights <- weights[not_miss]
  rep_weights <- rep_weights[not_miss,]

  nr <- nrow(values)
  nc <- ncol(values)

  rep_design <- survey::svrepdesign(repweights=rep_weights,weights=weights,
                                    data=dat,
                                    combined.weights=combined_weights,
                                    type=type, ...)
  scale <- rep_design$scale
  rscales <- rep_design$rscales
  mse <- rep_design$mse

  rep_weights <- weights(rep_design)
  nrep <- ncol(rep_weights)

  errors <- list()
  estimates <- array(NA, dim=c(nr, nc, nrep))
  for(i in 1:nrep){
    if(show_progress)
      cat(".")
    val <- try(as.matrix(fun(dat, rep_weights[,i])))
    if(!inherits(val, "try-error"))
      estimates[,,i] <- as.matrix(fun(dat, rep_weights[,i]))
    else
      errors[[length(errors) + 1]] <- val
  }
  if(show_progress)
    cat("\n")
  vars <- matrix(NA,nrow=nr,ncol=nc)
  for(i in 1:nr){
    for(j in 1:nc){
      if(!all(is.na(estimates[i,j,])))
        vars[i,j] <- svrVar(estimates[i,j,], scale, rscales, mse=mse,coef=values[i,j])
    }
  }
  bb <- list(value=values,
             var = vars,
             replicates=estimates,
             nrep=nrep,
             errors=errors)
  class(bb) <- c("inc_boot_est","list")
  bb
}

#' print bootstrap
#' @param x a inc_boot_est object
#' @param ... additional parameters for summary.inc_boot_est
#' @export
print.inc_boot_est <- function(x, ...){
  res <- summary(x,...)
  print(res)
}

#' summary of incidence bootstrap
#' @param object a inc_boot_est object
#' @param conf_level confidence level for intervals
#' @param ... additional parameters for print.data.frame
#' @export
summary.inc_boot_est <- function(object, conf_level=.95, ...){
  cival <- qnorm(1 - (1 - conf_level) / 2)

  val <- object$value[,1]
  v <- object$var[,1]
  res <- data.frame(incidence=val,
                    se=sqrt(v),
                    ci_lower=val - cival * sqrt(v),
                    ci_upper=val + cival * sqrt(v))

  val <- object$value[,2]
  v <- object$var[,2]
  res2 <- data.frame(transmission=val,
                     se=sqrt(v),
                     ci_lower=val - cival * sqrt(v),
                     ci_upper=val + cival * sqrt(v))
  #res <- cbind(res,res2,x$value[,-c(1,2), drop=FALSE])
  res <- cbind(res,res2)
  res
}


