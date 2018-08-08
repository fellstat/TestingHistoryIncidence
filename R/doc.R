
#' @author Ian E. Fellows \email{ian@fellstat.com}
#' @importFrom survey svrVar as.svrepdesign svydesign svymean
#' @importFrom stats weights optim pgamma pweibull var dgamma dweibull quantile na.omit optimise qnorm
#' @importFrom questionr wtd.table
NULL





#' A Synthetic Dataset for Incidence Estimation
#'
#' @usage
#' data(tstdat)
#' @docType data
#' @name tstdat
#' @section Variables:
#' 'hiv' - HIV status
#' 'ever_test' - Reported having ever had an HIV test.
#' 'last_test' - number of months since last HIV test.
#' 'undiag' - Actually has never been diagnosed with HIV in the past.
#' 'report_pos' - Reports having a previous positive test.
#' 'biomarker_art' - Tests positive for ART biomarkers.
#' 'low_viral' - Is classified as having low (undetectable) viral load.
#' @aliases tstdat
#' @keywords datasets
NULL
