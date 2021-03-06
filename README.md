# Incidence from Testing History in Crosssectional Surveys
the `TestingHistoryIncidence` package utilizes crosssectional survey data containing information on partcipants' testing history and diagnosis to estimate incidence.



## Installation

To install the latest development version from the github repo run:
```
# If devtools is not installed:
# install.packages("devtools")

devtools::install_github("fellstat/ipc")
devtools::install_github("fellstat/TestingHistoryIncidence", dependencies="Suggests")
```

## Resources


* For a more detailed description of what can be done with the ``ShinyAsyncTools`` package, **[see the introductory vignette](http://htmlpreview.github.io/?https://github.com/fellstat/TestingHistoryIncidence/blob/master/inst/doc/tst_hist_vig.html)**.


## Shiny User Interface

A user interface to the methods is provided in the package, and can be launched locally with
```
library(TestingHistoryIncidence)
shiny_testing_history()
```
* A 15 minute video documenting the various options of the UI and how to use them is available **[here](https://www.youtube.com/watch?v=YVPcLLs9zxc&t=08s)**.
* See the wiki guide for the interface **[here](https://github.com/fellstat/TestingHistoryIncidence/wiki/Shiny-App-Documentation)**

