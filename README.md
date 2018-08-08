# Incidence from Testing History in Crosssectional Surveys
the `TestingHistoryIncidence` package utilizes crosssectional survey data containing information on partcipants' testing history and diagnosis to estimate incidence.



## Installation

To install the latest development version from the github repo run:
```
# If devtools is not installed:
# install.packages("devtools")

devtools::install_github("fellstat/ShinyAsyncTools")
devtools::install_github("fellstat/TestingHistoryIncidence")
```

## Resources


* For a more detailed description of what can be done with the ``ShinyAsyncTools`` package, **[see the introductory vignette](http://htmlpreview.github.io/?https://github.com/fellstat/TestingHistoryIncidence/blob/master/inst/doc/tst_hist_vig.html)**.

## Shiny User Interface

A user interface to the methods is provided in the package, and can be launched locally with
```
library(TestingHistoryIncidence)
shiny_testing_history()
```
