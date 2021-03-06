% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incidence.R
\name{testing_incidence}
\alias{testing_incidence}
\title{Incindence from testing history}
\usage{
testing_incidence(report_pos, hiv, ever_test, last_test,
  last_upper = last_test, biomarker_art = rep(FALSE, length(report_pos)),
  low_viral = rep(FALSE, length(report_pos)), age = NULL,
  testing_debut_age = 0, weights = rep(1,
  length(report_pos))/length(report_pos), distribution = c("weibull",
  "empirical"), test_pop = c("negative", "undiagnosed"), ptruth = NULL,
  ptreated = NULL, aids_dist = function(x) pweibull(x/12, scale = 1/0.086,
  shape = 2.516), age_breaks = NULL, subset = NULL,
  uniform_missreport = FALSE)
}
\arguments{
\item{report_pos}{A logical vector indicating whether each subject reported a positive hiv status}

\item{hiv}{A logical vector indicating hiv status}

\item{ever_test}{A logical vector indicating whether the subject had eer had an hiv test.}

\item{last_test}{A numeric vector indicating the time since last hiv test in months. If testing times are binned into buckets, this is the lower bound of the months since last hiv test.}

\item{last_upper}{A numeric vector indicating the upper bound of the months since last hiv test for each individual.}

\item{biomarker_art}{A logical vector indicating whether ART antibodies are present. NA if test not done.}

\item{low_viral}{A logical vector indicating whether viral load is <=1000}

\item{age}{A numeric vector indicating the age of each subject in years.}

\item{testing_debut_age}{The age at which individuals begin engaging in regular testing.}

\item{weights}{Survey weights}

\item{distribution}{Either "empirical", or "weibull." This controls the family of distribution used to
model time since last test. "empirical" may not be used with binned testing times.}

\item{test_pop}{If "negative', the time since last negative is calculated amoung the HIV- population, otherwise it is calculated of those who report being undiagnosed.}

\item{ptruth}{The proportion of the diagnosed hiv positive population that would report being hiv positive.
If NULL, this is estimated using biomarker_art and low_viral.}

\item{ptreated}{The proportion of hiv positive individuals with postive biomarker_art or low_viral.}

\item{aids_dist}{The distribution function of time from infection to aids in months.}

\item{age_breaks}{age stratification break points.
If NULL, this is estimated from the data.}

\item{subset}{An optional vector specifying a subset of observations on which to perform the analysis. If uniform_missreport is TRUE, ptruth and ptreated are calculated over the full data.}

\item{uniform_missreport}{If true, the rate of missreporting of undiagnosed status is conisuered uniform over the age strata.}
}
\value{
A data.frame of class test_inc with elements:
'incidence': the estimated incidence.
'transmission rate': the estimated transmission rate.
'pundiag': The estimated proportion of positive cases that are undiagnosed after removing treated cases and adjusting for miss-reporting.
'psay_undiag': The proportion of positive cases that report being undiagnosed.
'pmiss_class': the proportion whose diagnosis status is incorrectly reported by the individual and are observed as treated due to viral load or art biomarkers.
'phiv': The proportion with a positive diagnosis.
'ptester': The proportion who have ever been tested.
'mean_time_since_last_test': the mean tie since last test in years.
'tid': mean time between infection and diagnosis in years.
'ptruth' the proportion of positive individuals that correctly report their status.
'ptreated': the proportion of positive indivudals who are identified as treated by viral load or biomarker.
}
\description{
Incindence from testing history
}
\details{
When time since last test is grouped into bins, last_upper should always be greater than last_test,
 and may be infinite (e.g. test was  > 24 months ago). Those with missing data or who never had an hiv test
should be assigned NA for both their lower and upper bound.
}
\examples{
  data(tstdat)
  tstdat$age <- rep(15:64, 200)


  inc <- with(tstdat, testing_incidence(report_pos, hiv,
                                        ever_test, last_test,
                                        biomarker_art=biomarker_art, low_viral=low_viral))
  inc

  # Using last test times divided into bins
  tstdat$last_test_lower <- 24
  tstdat$last_test_upper <- Inf
  tstdat$last_test_lower[tstdat$last_test < 6] <- 0
  tstdat$last_test_upper[tstdat$last_test < 6] <- 6
  tstdat$last_test_lower[tstdat$last_test >= 6 & tstdat$last_test < 12] <- 6
  tstdat$last_test_upper[tstdat$last_test >= 6 & tstdat$last_test < 12] <- 12
  tstdat$last_test_lower[tstdat$last_test >= 12 & tstdat$last_test < 24] <- 12
  tstdat$last_test_upper[tstdat$last_test >= 12 & tstdat$last_test < 24] <- 24
  inc <- with(tstdat, testing_incidence(report_pos, hiv,
                                        ever_test, last_test_lower, last_test_upper,
                                        biomarker_art=biomarker_art, low_viral=low_viral))
  inc

  # HIV testing starts at age 13
  inc <- with(tstdat, testing_incidence(report_pos, hiv,
                                        ever_test, last_test,
                                        biomarker_art=biomarker_art, low_viral=low_viral,
                                        age=age,testing_debut_age=13))
  inc

  # Stratify by age
  inc <- with(tstdat, testing_incidence(report_pos, hiv,
                                        ever_test, last_test,
                                        biomarker_art=biomarker_art, low_viral=low_viral,
                                        age=age,testing_debut_age=13,age_breaks=c(25,35,45,55)))
  inc

  # Pooled estimate probability of miss-reporting diagnosis status across groups
  inc <- with(tstdat, testing_incidence(report_pos, hiv,
                                        ever_test, last_test,
                                        biomarker_art=biomarker_art, low_viral=low_viral,
                                        age=age,testing_debut_age=13,
                                        age_breaks=c(25,35,45,55),
                                        uniform_missreport=TRUE))
  inc
}
