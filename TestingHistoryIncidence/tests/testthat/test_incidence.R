library(testthat)
context("Incidence")

test_that(".fit_dist_weibull", {
  ss <- 50000
  last_test <- rweibull(ss, 1, 5) * 12
  mlt <- mean(last_test)
  ever_test <- runif(ss) > .2
  last_test[!ever_test] <- NA
  weights <- rep(1,ss)
  time_at_risk <- runif(ss,0,300)
  censored <- last_test > time_at_risk
  last_test[censored] <- NA
  ever_test[censored] <- FALSE
  aids_dist <- function(x) pweibull(x / 12, scale=1/0.086, shape=2.516)
  fit <- .fit_dist_weibull(last_test,last_test, ever_test, weights, time_at_risk, aids_dist)
  fit_mean <- sum(fit$test_surv_cond)
  fit_mean
  mlt
  fit$tau
  expect_true(abs(fit_mean - mlt) < 2)
  expect_true(abs(fit$tau - .8) < .02)
  expect_equal(1, 1)
})



test_that(".fit_dist_empirical", {
  ss <- 50000
  last_test <- rweibull(ss, 1, 5) * 12
  mlt <- mean(last_test)
  ever_test <- runif(ss) > .2
  last_test[!ever_test] <- NA
  weights <- rep(1,ss)
  time_at_risk <- rep(MAX_YEARS*12,ss)
  #censored <- last_test > time_at_risk
  #last_test[censored] <- NA
  #ever_test[censored] <- FALSE
  aids_dist <- function(x) pweibull(x / 12, scale=1/0.086, shape=2.516)
  fit <- .fit_dist_empirical(last_test, ever_test, weights, time_at_risk, aids_dist)
  fit_mean <- sum(fit$test_surv_cond)
  fit_mean
  mlt
  fit$tau
  expect_true(abs(fit_mean - mlt) < 1)
  expect_true(abs(fit$tau - .8) < .02)
})


test_that("testing_incidence", {
  data(tstdat)

  inc <- with(tstdat, testing_incidence(report_pos=report_pos, hiv=hiv,
                                        ever_test=ever_test, last_test=last_test, biomarker_art=biomarker_art,
                                        low_viral=low_viral))
  expect_equivalent(as.matrix(inc),
                    structure(c(0.00575807622336409, 0.0508601757271186, 0.291371713560825,
                                0.396263520157325, 0.0897422142330957, 0.1017, 0.707334433262165,
                                4.4364845698551, 5.72887744478369, 0.850271528316524, 0.400633605458755,
                                10000), .Dim = c(1L, 12L), .Dimnames = list("estimate", c("incidence",
                                                                                          "transmission", "pundiag", "psay_undiag", "pmiss_class", "phiv",
                                                                                          "ptester", "mean_time_since_last_test", "tid", "ptruth", "ptreated",
                                                                                          "n")))
  )

  bi <- bootstrap_incidence(inc,nrep=3, show_progress = FALSE)
  expect_equal(summary(bi)[1,1],inc[1,1])


  inc <- with(tstdat, testing_incidence(report_pos=report_pos, hiv=hiv,
                                        ever_test=ever_test, last_test=last_test))
})


