## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE---------------------------------------------------------
library(knitr)
library(pander)
panderOptions('table.alignment.default', function(df) rep("left",ncol(df)))

## ------------------------------------------------------------------------
library(TestingHistoryIncidence)
data(tstdat)

# Add age
tstdat$age <- rep(15:64, 200)

head(tstdat)
summary(tstdat)

## ------------------------------------------------------------------------
attach(tstdat)
inc <- testing_incidence(report_pos=report_pos, biomarker_art=biomarker_art, 
                         low_viral=low_viral, hiv=hiv,ever_test=ever_test,
                         last_test=last_test)
#print(inc)

## ---- echo=FALSE, results="asis"-----------------------------------------
#kable(inc)
pander(inc)

## ------------------------------------------------------------------------
inc <- testing_incidence(report_pos=report_pos, biomarker_art=biomarker_art, 
                         low_viral=low_viral, hiv=hiv, ever_test=ever_test, 
                         last_test=last_test, age=age, testing_debut_age=13)
#print(inc)

## ---- echo=FALSE, results="asis"-----------------------------------------
#kable(inc)
pander(inc)

## ------------------------------------------------------------------------
subset <- c(rep(TRUE,1000),rep(FALSE,9000))
inc1 <- testing_incidence(report_pos, hiv,ever_test, last_test,  
                          biomarker_art=biomarker_art, low_viral=low_viral,age=age, 
                          testing_debut_age=13, subset=subset)
inc2 <- testing_incidence(report_pos, hiv,
                                        ever_test, last_test, 
                          biomarker_art=biomarker_art, low_viral=low_viral, age=age, 
                                        testing_debut_age=13, subset=subset, uniform_missreport=TRUE)
inc <- rbind(inc1,inc2)
rownames(inc) <- c("restricted","uniform")
#print(inc[,c(1,2,3,10,11)])

## ---- echo=FALSE, results="asis"-----------------------------------------
#kable(inc)
pander(inc[,c(1,2,3,10,11)],caption=NULL, style = 'rmarkdown')

## ------------------------------------------------------------------------
inc <- testing_incidence(report_pos, hiv,
                                        ever_test, last_test, biomarker_art=biomarker_art, low_viral=low_viral, age=age, 
                                        testing_debut_age=13, age_breaks=c(25,35,45,55))

#print(inc[,c(1,2,3,6,9)])

## ---- echo=FALSE, results="asis"-----------------------------------------
#kable(inc)
pander(inc[,c(1,2,3,6,9)],caption=NULL, style = 'rmarkdown')

## ------------------------------------------------------------------------
inc <- testing_incidence(report_pos, hiv,
                                        ever_test, last_test, biomarker_art=biomarker_art, low_viral=low_viral, age=age, testing_debut_age=13)

boots <- bootstrap_incidence(inc, nrep=10)
#print(boots)

## ---- echo=FALSE, results="asis"-----------------------------------------
#kable(inc)
pander(summary(boots))

## ------------------------------------------------------------------------
weights <- runif(10000)
inc <- testing_incidence(report_pos, hiv,
                                        ever_test, last_test, biomarker_art=biomarker_art, low_viral=low_viral, age=age, testing_debut_age=13, weights=weights)


## ------------------------------------------------------------------------
library(survey)
stype <- sample(1:10, 10000, replace = TRUE)
# stratified sample
dstrat<-svydesign(id=~1,strata=~stype, weights=~weights, data=tstdat)

# generate replicate design
des1 <- as.svrepdesign(dstrat, type="bootstrap", compress=FALSE, replicates=25)
rep_weights <- des1$repweights * weights # create combined replication weights
head(rep_weights)

## ------------------------------------------------------------------------
boots <- bootstrap_incidence(inc, rep_weights = rep_weights, type="bootstrap")
#print(boots)

## ---- echo=FALSE, results="asis"-----------------------------------------
#kable(inc)
pander(summary(boots))

## ---- echo=FALSE---------------------------------------------------------
detach(tstdat)

