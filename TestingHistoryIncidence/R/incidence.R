MAX_YEARS <- 200


# Fits either a weibull or gamma distribution to possibly windowed data
.fit_dist <- function(last_test, last_upper, weights, distribution, initial_rate=1/14){
  #get counts for testing windows
  test_window <- paste(last_test,last_upper,sep="_")
  ln_count <- wtd.table(test_window, weights = weights)
  ln_windows <- strsplit(names(ln_count),"_")
  ln_lower <- as.numeric(sapply(ln_windows, function(x) x[1]))
  ln_upper <- as.numeric(sapply(ln_windows, function(x) x[2]))
  ln_count <- as.vector(ln_count)

  #The data likelihood for time since last negative test
  lik2 <- function(scale,k){

    if(distribution == "weibull"){
      p <- function(x) suppressWarnings(pweibull(x, scale=scale,shape=k))
      log_d <- function(x) suppressWarnings(dweibull(x, scale=scale,shape=k, log = TRUE))
    }else{
      p <- function(x) suppressWarnings(pgamma(x, scale=scale,shape=k))
      log_d <- function(x) suppressWarnings(dgamma(x, scale=scale,shape=k, log = TRUE))
    }

    result <- 0
    n <- length(ln_count)
    exact_time <- ln_upper == ln_lower
    for(i in 1:n){
      if(!exact_time[i]){
        result <- result - ln_count[i] * log(p(ln_upper[i]) - p(ln_lower[i]))
      }
    }
    result <- result - sum(ln_count[exact_time] * log_d(ln_upper[exact_time]) )
    result
  }
  opt <- optim(function(x)lik2(x[1],x[2]),par = c(1/initial_rate,1))
  opt$par
}

.mean_of_dist <- function(par, distribution){
  if(distribution == "weibull"){
    m2 <- par[1] * gamma(1 + 1/ par[2])
  }else{
    m2 <- par[1] * par[2]
  }
  m2
}

# Fits either a weibull or gamma distribution to possibly windowed data
.fit_dist_empirical <- function(last_test, ever_test, weights, time_at_risk, aids_dist){

  #construct empirical distribution for TID | TESTER
  tm <- last_test
  event <- !is.na(tm) & tm < time_at_risk
  event[(ever_test | is.na(ever_test)) & is.na(tm)] <- NA
  tm[!event & !is.na(event)] <- time_at_risk[!event & !is.na(event)]
  tdist <- survfit(Surv(tm, event) ~ 1, weights=weights)
  times <- seq(from=0,to=MAX_YEARS*12, by=iv)
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


.fit_dist_weibull <- function(last_test, last_upper, ever_test, weights, time_at_risk, aids_dist, initial_rate=1/14){

  exact_time <- all(na.omit(last_test == last_upper))

  #construct empirical distribution for TID | TESTER
  miss_test <- is.na(last_test) | is.na(last_upper)
  tm <- last_test
  event <- !miss_test & last_test < time_at_risk
  event[(ever_test | is.na(ever_test)) & miss_test] <- NA

  if(!exact_time){
    last_upper <- ifelse(!is.na(time_at_risk) & !is.na(last_upper) & time_at_risk < last_upper, time_at_risk, last_upper)
  }

  # Fit a weibull distribtion to get P(TSLT | TESTER)
  weibull_lik <- function(scale, k, tau){
    p <- function(x) suppressWarnings(pweibull(x, scale=scale,shape=k))
    log_d <- function(x) suppressWarnings(dweibull(x, scale=scale,shape=k, log = TRUE))
    if(exact_time){
      lik <- sum((log_d(last_test[event]) + log(tau)) * weights[event], na.rm=TRUE)
    }else{
      lik <- sum((log(p(last_upper[event]) - p(last_test[event])) + log(tau)) * weights[event], na.rm=TRUE)
    }
    lik <- lik + sum( log((1 - p(time_at_risk[!event])) * tau + 1 - tau) * weights[!event], na.rm=TRUE )
    -lik
  }
  opt <- optim(function(x)weibull_lik(x[1],x[2], x[3]),
               par = c(1/initial_rate,1, .5),
               lower=c(0,0,0),
               upper=c(Inf, Inf, 1),
               method="L-BFGS-B")
  times <- seq(from=0,to=MAX_YEARS*12, by=iv)
  test_surv_cond <- 1 - pweibull(times, scale=opt$par[1],shape=opt$par[2])
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

#' Incindence from testing history
#'
#' @param report_pos A logical vector indicating whether each subject reported a positive hiv status
#' @param biomarker_art A logical vector indicating whether ART antibodies are present. NA if test not done.
#' @param low_viral A logical vector indicating whether viral load is <=1000
#' @param hiv A logical vector indicating hiv status
#' @param ever_test A logical vector indicating whether the subject had eer had an hiv test.
#' @param last_test A numeric vector indicating the time since last hiv test in months. If testing times are binned into buckets, this is the lower bound of the months since last hiv test.
#' @param last_upper A numeric vector indicating the upper bound of the months since last hiv test for each individual.
#' @param weights Survey weights
#' @param distribution Either "empirical", "weibull" or "gamma." This controls the family of distribution used to
#' model time since last test. "empirical" may not be used with binned testing times.
#' @param test_pop If "negative', the time since last negative is calculated amoung the HIV- population, otherwise it is calculated of those who report being undiagnosed.
#' @param ptruth The proportion of the diagnosed hiv positive population that would report being hiv positive.
#' If NULL, this is estimated using biomarker_art and low_viral.
#' @param ptreated The proportion of hiv positive individuals with postive biomarker_art or low_viral.
#' If NULL, this is estimated from the data.
#' @param non_tester_tid mean nummber of months between infection and diagnosis for non_testers.
#' @details
#' When time since last test is grouped into bins, last_upper should always be greater than last_test,
#'  and may be infinite (e.g. test was  > 24 months ago). Those with missing data or who never had an hiv test
#' should be assigned NA for both their lower and upper bound.
#' @return A data.frame of class test_inc with elements:
#' 'incidence': the estimated incidence.
#' 'transmission rate': the estimated transmission rate.
#' 'pundiag': The estimated proportion of positive cases that are undiagnosed.
#' 'psay_undiag': The proportion of positive cases that report being undiagnosed.
#' 'pmiss_class': the proportion whose diagnosis status is incorrectly reported by the individual .
#' 'phiv': The proportion with a positive diagnosis.
#' 'ptester': The proportion who have ever been tested.
#' 'mean_time_since_last_test': the mean tie since last test in years.
#' 'tid': mean time between infection and diagnosis in years.
#' 'ptruth' the proportion of positive individuals that correctly report their status.
#' 'ptreated': the proportion of positive indivudals who are identified as treated by viral load or biomarker.
#' @export
testing_incidence <- function(report_pos, biomarker_art, low_viral, hiv,
                              ever_test, last_test, last_upper = last_test,
                              age=NULL, debut_age=0,
                              weights=rep(1, length(report_pos)) / length(report_pos),
                              distribution=c("weibull","gamma","empirical"), test_pop=c("negative","undiagnosed"),
                              ptruth=NULL, ptreated=NULL,
                              non_tester_tid= 123.824){
  test_pop <- match.arg(test_pop)
  distribution <- match.arg(distribution)

  treated <- low_viral | biomarker_art
  treated[is.na(treated)] <- FALSE

  if(!is.null(age)){
    time_at_risk <- (age - debut_age) * 12
  }else{
    time_at_risk <- rep(1000*12, length(last_test))
  }

  ptester <- as.vector(prop.table(wtd.table(ever_test[!hiv],weights=weights[!hiv]))[2])

  tab <- wtd.table(!report_pos, hiv, weights=weights)

  psay_undiag <- tab[2,2] / sum(tab[,2])

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

  # Adjust for missreporting of HIV status
  tlie <- wtd.table(treated, !report_pos, weights = weights)
  if(is.null(ptruth))
    ptruth <- tlie[2,1] / sum(tlie[2,])
  if(is.null(ptreated))
    ptreated <- tlie[2,2] / sum(tlie[,2])
  pmiss_class <- 1 - (ptruth + ptreated * ( 1 - ptruth))

  pundiag <- (psay_undiag - pmiss_class) / (1 - pmiss_class)

  # Calculate incidence
  tid <- tid / 12
  inc <- pundiag * phiv / (tid * (1-phiv))

  #Calculate transmission rate
  trans <- pundiag / tid

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
                 ptreated=ptreated))
  row.names(result) <- "estimate"
  class(result) <- c("test_inc","data.frame")
  result
}


#' Perorm a survey bootstrap
#' @param design an object of class svydesign
#' @param fun a function taking a dataframe as a parameter
#' @param show_progress print progress
#' @param ... additional parameters passed to as.srvrepdesign
#' @return
#' a list with class boot_est containing
#' 'value': the value of the function applied to dat.
#' 'var': the boostrap variance.
#' 'replicates': The bootstrap values of fun.
#' 'nrep': the number of replicates.
#' @examples
#' library(survey)
#' data(api,package="survey")
#' ## one-stage cluster sample
#' dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
#' ## convert to JK1 jackknife
#' rclus1<-as.svrepdesign(dclus1)
#' ## convert to bootstrap
#' set.seed(1)
#' bclus1<-as.svrepdesign(dclus1,type="bootstrap", replicates=100)
#' attr(svymean(~api00, bclus1),"var")
#' set.seed(1)
#' survey_bootstrap(dclus1, fun=function(dat, wts) {sum(wts*dat$api00) / sum(wts)},
#'   type="bootstrap", replicates=100)$var
#' @export
survey_bootstrap <- function(design, fun, show_progress=TRUE, ...){
  dat <- design$variables
  weights <- weights(design)

  coef <- fun(dat, weights)

  des1 <- as.svrepdesign(design, compress=FALSE, ...)
  scale <- des1$scale
  rscales <- des1$rscales
  mse <- des1$mse

  rep_weights <- weights(des1)
  nrep <- ncol(rep_weights)
  estimates <- matrix(NA, nrow=nrep, ncol=length(coef))
  for(i in 1:nrep){
    if(show_progress)
      cat(".")
    estimates[i,] <- fun(dat, weights * rep_weights[,i])
  }
  if(show_progress)
    cat("\n")
  vars <- apply(estimates, 2, function(x) svrVar(x, scale, rscales, mse=mse,coef=coef))
  bb <- list(value=coef,
       var = vars,
       replicates=estimates,
       nrep=nrep)
  class(bb) <- c("boot_est","list")
  bb
}

#' Bootstrap variance for incidence
#' @param frame a dataframe conatining the variables needed for id_formula, strata_formula and wieghts_formula
#' @param strata_formula a survey formula specifying the strata structure
#' @param weights_formula a survey formula specifying the weights
#' @param id_formula a survey formula specifying the unit of sampling
#' @param report_pos A logical vector indicating whether each subject reported a positive hiv status
#' @param biomarker_art A logical vector indicating whether ART antibodies are present. NA if test not done.
#' @param low_viral A logical vector indicating whether viral load is <=1000
#' @param hiv A logical vector indicating hiv status
#' @param last_test A numeric vector indicating the lower bound of the months since last hiv test bin in months.
#' @param last_upper A numeric vector indicating the upper bound of the months since last hiv test bin in months.
#' @param ever_test A logical vector indicating whether the subject had eer had an hiv test.
#' @param replicates The number of bootstrap replicated
#' @param type The type of survey bootstrap
#' @param show_progress print progress
#' @param ... additional parameters for testing_incidence
#' @export
bootstrap_incidence <- function(report_pos, biomarker_art, low_viral, hiv,
                                last_test, last_upper, ever_test, replicates=5000,
                                frame=NULL, id_formula=NULL, strata_formula=NULL, weights_formula=NULL, type="mrbbootstrap", show_progress=TRUE, ...){
  if(is.null(frame) && is.null(id_formula) && is.null(strata_formula) && is.null(weights_formula)){
    dat <- data.frame(report_pos,
      biomarker_art,
      low_viral,
      hiv,
      last_test,
      last_upper,
      ever_test,
      weights = rep(1,length(report_pos))
      )
    bb <- srs_bootstrap(dat,
                           fun=function(dat) unlist(testing_incidence(dat$report_pos, dat$biomarker_art, dat$low_viral, dat$hiv,
                                                                    dat$last_test, dat$last_upper, dat$ever_test, dat$weights, ...)),
                           replicates=replicates,
                        show_progress=show_progress)
  }else{
    include <- rowSums(is.na(frame)) == 0
    frame$report_pos <- report_pos
    frame$biomarker_art <- biomarker_art
    frame$low_viral <- low_viral
    frame$hiv <- hiv
    frame$last_test <- last_test
    frame$last_upper <- last_upper
    frame$ever_test <- ever_test

    frame <- frame[include,]
    des <- svydesign(ids = id_formula, strata=strata_formula, weights =weights_formula, data=frame)
    bb <- survey_bootstrap(des,
                           fun=function(dat, wts) {
                             inc <- testing_incidence(dat$report_pos, dat$biomarker_art,
                                                      dat$low_viral, dat$hiv, dat$last_test,
                                                      dat$last_upper, dat$ever_test,
                                                      wts, ...)
                             unlist(inc)
                           },
                           show_progress=show_progress,
                           type=type,
                           replicates=replicates)
  }
  bb
}

#' Performs a bootstrap assuming a simple random sample
#' @param dat a data.frame
#' @param fun a function taking a dataframe as a parameter
#' @param replicates the number of replicates
#' @param show_progress print progress
#' @return
#' a list with class boot_est containing
#' 'value': the value of the function applied to dat.
#' 'var': the boostrap variance.
#' 'replicates': The bootstrap values of fun.
#' 'nrep': the number of replicates.
#' @examples
#' df <- data.frame(a=rnorm(100))
#' srs_bootstrap(df, colMeans,replicates=100)$var
#' @export
srs_bootstrap <- function(dat, fun, replicates=1000, show_progress=TRUE){
  value <- fun(dat)
  n <- nrow(dat)
  estimates <- matrix(NA, nrow=replicates, ncol=length(value))
  for(i in 1:replicates){
    if(show_progress)
      cat(".")
    samp <- sample.int(n, n, replace=TRUE)
    boot <- dat[samp, , drop=FALSE]
    estimates[i,] <- fun(boot)
  }
  if(show_progress)
    cat("\n")
  bb <- list(value=value,
       var=diag(var(estimates)),
       replicates=estimates,
       nrep=replicates)
  class(bb) <- c("boot_est","list")
  bb
}

#' print bootstrap
#' @param x a boot_est object
#' @param probs A list of quantiles for the bootstrap
#' @param ... additional parameters for print.data.frame
#' @export
print.boot_est <- function(x,probs=c(.025,.975), ...){
  res <- data.frame(value=x$value,
                    se=sqrt(x$var))
  res <- cbind(res,t(apply(x$replicates, 2, quantile, probs=probs)))
  print(res, ...)
}

