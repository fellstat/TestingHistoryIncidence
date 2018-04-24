survey.bootstrap<-function(design,fun,...)    
{   ##############################################
    # Function to compute survey replications
    # Parameters:
    # dat = design$variables from survey design object 
    # weights = weights from survey design object
    # est = estimate to return.
    #*********************************************
    # Return:
    # estimates
    # var=variance of replicates
    #*********************************************
    
    output<-fun(report_hiv_pos=design$variables$report_hiv_pos,biomarker_art=design$variables$biomarker_art,low_viral=design$variables$low_viral,hiv=design$variables$hiv,
                              last_hiv_lower=design$variables$last_hiv_lower,last_hiv_upper=design$variables$last_hiv_upper,ever_hiv_test=design$variables$ever_hiv_test,
                              weights=design$pweights,distribution="weibull",
                              ptruth=NULL,ptreated=NULL,non_tester_tid=123.824)
    if(anyNA(output))
    {
        return(list(value=rep(NA,7)))
    }
    result<-c(incidence=output$incidence,transmission=output$pundiag/output$tid,"P(U|H)"=output$pundiag,"P(H)"=output$phiv,"E(TID)"=output$tid)
    scale <- design$scale
    rscales <- design$rscales
    rep_weights <- design$repweights
    nrep <- ncol(rep_weights)
    k<-length(output)
    estimates<-matrix(NA,nrow=nrep,ncol=k+1)
    colnames(estimates)<-c(names(output),"transmission")
    for(i in 1:nrep)
    {
        if(i%%50==0) print(i);
        estimates[i,1:k]<-unlist(fun(report_hiv_pos=design$variables$report_hiv_pos,biomarker_art=design$variables$biomarker_art,low_viral=design$variables$low_viral,hiv=design$variables$hiv,
                                  last_hiv_lower=design$variables$last_hiv_lower,last_hiv_upper=design$variables$last_hiv_upper,ever_hiv_test=design$variables$ever_hiv_test,
                                  weights=rep_weights[,i],distribution="weibull",
                                  ptruth=NULL,ptreated=NULL,non_tester_tid=123.824))
        estimates[i,k+1]<-estimates[i,"pundiag"]/estimates[i,"tid"]

    }
    
    incidence.var.T<-svrVar(estimates[,1], scale, rscales, mse=T,coef=result["incidence"])
    incidence.var.F<-svrVar(estimates[,1], scale, rscales, mse=F,coef=result["incidence"])
    transmission.var.T<-svrVar(estimates[,k+1], scale, rscales, mse=T,coef=result["transmission"])
    transmission.var.F<-svrVar(estimates[,k+1], scale, rscales, mse=F,coef=result["transmission"])
    inc.se<-sqrt(incidence.var.T[1])
    trans.se<-sqrt(transmission.var.T[1])
    
    list(value=c("incidence"=result[[1]],"SE(lambda)"=inc.se,"transmission"=result[[2]],"SE(tau)"=trans.se,result[3:5]),inc.se=inc.se,trans.se=trans.se,
         output=output,estimates.rep=estimates,
         incidence.var.T=incidence.var.T,incidence.var.F=incidence.var.F,transmission.var.T=transmission.var.T,transmission.var.F=transmission.var.F,
         nrep=nrep)
}
