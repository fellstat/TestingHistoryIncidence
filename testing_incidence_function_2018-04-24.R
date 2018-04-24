testing_incidence<-function (report_hiv_pos, biomarker_art, low_viral, hiv, last_hiv_lower, 
          last_hiv_upper, ever_hiv_test, weights = rep(1, length(report_hiv_pos))/length(report_hiv_pos), 
          distribution = "weibull", ptruth = NULL, ptreated = NULL, 
          non_tester_tid = 123.824) 
{
    treated <- low_viral | biomarker_art
    treated[is.na(treated)] <- FALSE
    tab <- wtd.table(!report_hiv_pos, hiv, weights = weights)
    #******************
    # error check
    if(prod(dim(tab))!=4)
    {
        result<-NA
        return(result)
    }
    #******************
    psay_undiag <- tab[2, 2]/sum(tab[, 2])
    phiv <- as.vector(prop.table(wtd.table(hiv, weights = weights))[2])
    ln_sub <- !hiv & !is.na(hiv) & !is.na(last_hiv_upper) & !is.na(last_hiv_lower) & 
        !is.na(weights)
    test_window <- paste(last_hiv_lower[ln_sub], last_hiv_upper[ln_sub], 
                         sep = "_")
    ln_count <- wtd.table(test_window, weights = weights[ln_sub])
    ln_windows <- strsplit(names(ln_count), "_")
    ln_lower <- as.numeric(sapply(ln_windows, function(x) x[1]))
    ln_upper <- as.numeric(sapply(ln_windows, function(x) x[2]))
    ln_count <- as.vector(ln_count)
    rate <- 1/14.03188
    lik2 <- function(scale, k) {
        if (distribution == "weibull") {
            p <- function(x) suppressWarnings(pweibull(x, scale = scale, 
                                                       shape = k))
            log_d <- function(x) suppressWarnings(dweibull(x, 
                                                           scale = scale, shape = k, log = TRUE))
        }
        else {
            p <- function(x) suppressWarnings(pgamma(x, scale = scale, 
                                                     shape = k))
            log_d <- function(x) suppressWarnings(dgamma(x, scale = scale, 
                                                         shape = k, log = TRUE))
        }
        result <- 0
        n <- length(ln_count)
        exact_time <- ln_upper == ln_lower
        for (i in 1:n) {
            if (!exact_time[i]) {
                result <- result - ln_count[i] * log(p(ln_upper[i]) - 
                                                         p(ln_lower[i]))
            }
        }
        result <- result - sum(ln_count[exact_time] * log_d(ln_upper[exact_time]))
        result
    }
    #***** end lik2 function
    opt <- optim(function(x) lik2(x[1], x[2]), par = c(1/rate,1))
    if (distribution == "weibull") {
        m2 <- opt$par[1] * gamma(1 + 1/opt$par[2])
    }    else {
        m2 <- opt$par[1] * opt$par[2]
    }
    ptester <- as.vector(prop.table(wtd.table(ever_hiv_test[!hiv],weights = weights[!hiv]))[2])
    tid <- m2 * ptester + non_tester_tid * (1 - ptester)
    tlie <- wtd.table(treated, !report_hiv_pos, weights = weights)
    
    #******************
    # error check
    
    if(prod(dim(tlie))!=4)
    {
        result<-NA
        return(result)
    }
    #******************
    
    if (is.null(ptruth)) 
        ptruth <- tlie[2, 1]/sum(tlie[2, ])
    if (is.null(ptreated)) 
        ptreated <- tlie[2, 2]/sum(tlie[, 2])
    pmiss_class <- 1 - (ptruth + ptreated * (1 - ptruth))
    pundiag <- (psay_undiag - pmiss_class)/(1 - pmiss_class)
    tid <- tid/12
    inc <- pundiag * phiv/(tid * (1 - phiv))
    result <- data.frame(list(incidence = inc, pundiag = pundiag, 
                              psay_undiag = psay_undiag, pmiss_class = pmiss_class, 
                              phiv = phiv, ptester = ptester, mean_time_since_last_test = m2/12, 
                              tid = tid, ptruth = ptruth, ptreated = ptreated))
    row.names(result) <- "estimate"
    class(result) <- c("test_inc", "data.frame")
    result
}