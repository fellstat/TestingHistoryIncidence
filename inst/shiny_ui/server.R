

library(shiny)
library(promises)
library(future)
library(ipc)
plan(multiprocess)
library(ggplot2)
library(TestingHistoryIncidence)
library(survey)
options(shiny.maxRequestSize=300*1024^2)

shinyServer(function(input, output, session) {

  get_raw_data <- reactive({
    inFile <- input$file1

    if (is.null(inFile))
      return(NULL)

    dat <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
    vars <- names(dat)

    for(name in c("hiv","biomarker_art","low_viral","ever_test",
                  "report_pos","last_test","last_test_lower",
                  "last_test_upper","age", "weights","strata",
                  "design_strata","design_clusters"))
      updateSelectizeInput(session, name, choices = c("Choose Variable"="",vars))
    updatePickerInput(session, "rep_weights", choices = vars)
    dat
  })

  get_numeric <- function(name){
    dat <- get_raw_data()
    if(is.null(dat) || input[[name]] == "")
      return(NULL)
    as.numeric(dat[[input[[name]]]])
  }

  get_categorical <- function(name){
    dat <- get_raw_data()
    if(is.null(dat) || input[[name]] == "")
      return(NULL)
    as.factor(dat[[input[[name]]]])
  }

  get_logical <- function(name){
    dat <- get_raw_data()
    if(is.null(dat) || input[[name]] == "")
      return(NULL)
    variable <- dat[[input[[name]]]]
    if(is.logical(variable))
      return(variable)
    if(is.numeric(variable))
      return(variable == max(variable, na.rm=TRUE))
    variable <- as.factor(variable)
    variable == max(levels(variable))
  }

  render_raw_table <- function(name){
    renderTable({
      dat <- get_raw_data()
      if(is.null(dat))
        return(NULL)
      v <- dat[[input[[name]]]]
      table(v, useNA = "always", dnn=name)
    })
  }

  render_table <- function(name){
    renderTable({
      v <- eval(parse(text=paste0("get_", name)))()
      table(v, useNA = "always", dnn=name)
    })
  }

  get_age <- reactive({
    get_numeric("age")
  })

  get_weights <- reactive({
    get_numeric("weights")
  })

  get_hiv <-  reactive({
    get_logical("hiv")
  })

  get_report_pos <-  reactive({
    get_logical("report_pos")
  })

  get_ever_test <-  reactive({
    get_logical("ever_test")
  })

  get_biomarker_art <-  reactive({
    get_logical("biomarker_art")
  })

  get_low_viral <-  reactive({
    get_logical("low_viral")
  })

  get_last_test <-  reactive({
    get_numeric("last_test")
  })

  get_last_test_lower <-  reactive({
    get_numeric("last_test_lower")
  })

  get_rep_weights <- reactive({
    dat <- get_raw_data()
    rw_names <- input$rep_weights
    if(is.null(dat) || length(rw_names) == 0)
      return(NULL)
    rw <- lapply(as.list(rw_names), function(name){
      as.numeric(dat[[name]])
    })
    as.data.frame(rw)
  })

  get_last_test_upper <-  reactive({
    low <- get_numeric("last_test_lower")
    up <- get_numeric("last_test_upper")
    up[is.na(up) & !is.na(low) & low == max(low, na.rm=TRUE)]
    up
  })

  get_strata <-  reactive({
    get_categorical("strata")
  })

  get_design_strata <- reactive({
    get_categorical("design_strata")
  })

  get_design_clusters <- reactive({
    get_categorical("design_clusters")
  })

  output$contents <- renderTable({
    if(is.null(get_raw_data()))
      return(NULL)

    head(get_raw_data())
  })

  output$table <- renderDataTable({
    if(is.null(get_raw_data()))
      return(NULL)

    head(get_raw_data())
  })

  output$age_desc_plot <- renderPlot({
    age <- get_age()
    if(is.null(age))
      return(NULL)
    print(qplot(age, bins = 30))
  })

  output$age_desc_errors <- renderText({
    age <- get_age()
    if(is.null(age))
      return("")
    txt <- ""
    if(sum(is.na(age)) > length(age) / 2)
      txt <- paste(txt, "Error: More than half of values are missing",sep="\n")
    age <- na.omit(age)
    if(any(age < 1))
      txt <- paste(txt, "Error: Age values < 1 detected",sep="\n")
    if(any(age > 100))
      txt <- paste(txt, "Error: Age values > 100 detected",sep="\n")
    txt
  })

  output$last_test_plot <- renderPlot({
    lt <- get_last_test()
    if(is.null(lt))
      return(NULL)
    print(qplot(lt, bins = 30) + xlab("Last HIV Test"))
  })

  output$last_test_errors <- renderText({
    lt <- get_last_test()
    if(is.null(lt))
      return("")
    txt <- ""
    if(sum(is.na(lt)) > length(lt) / 2)
      txt <- paste(txt, "Error: More than half of values are missing",sep="\n")
    lt <- na.omit(lt)
    if(any(lt <= 0))
      txt <- paste(txt, "Error: Last HIV test values <= 0 detected",sep="\n")
    if(any(lt > 35 * 12))
      txt <- paste(txt, "Error: Last HIV test values > 35 years ago",sep="\n")
    txt
  })

  output$hiv_desc <- render_table("hiv")
  output$hiv_desc_raw <- render_raw_table("hiv")

  output$report_pos_desc <- render_table("report_pos")
  output$report_pos_desc_raw <- render_raw_table("report_pos")

  output$ever_test <- render_table("ever_test")
  output$ever_test_desc_raw <- render_raw_table("ever_test")

  output$biomarker_art_desc <- render_table("biomarker_art")
  output$biomarker_art_desc_raw <- render_raw_table("biomarker_art")

  output$low_viral_desc <- render_table("low_viral")
  output$low_viral_desc_raw <- render_raw_table("low_viral")

  output$last_test_bound_desc <-     renderTable({
    low <- get_last_test_lower()
    up <- get_last_test_upper()
    hiv <- get_hiv()
    if(is.null(hiv))
      tbl <- table(low, up, useNA = "always",
                   dnn=c("Lower Bound","Upper Bound"))
    else
      tbl <- table(low, up, hiv, useNA = "always",
                   dnn=c("Lower Bound","Upper Bound", "HIV+"))
    tbl
  })

  incidence <- reactiveVal()
  nclicks <- reactiveVal(0)
  output$inc_results <- renderTable({
    nclicks(0)

    # Rerun on UI change
    get_design_clusters()
    get_design_strata()
    get_rep_weights()

    report_pos <- get_report_pos()
    biomarker_art <- get_biomarker_art()
    low_viral <- get_low_viral()
    hiv <- get_hiv()
    ever_test <- get_ever_test()
    last_test <- get_last_test()
    if(is.null(last_test)){
      last_test <- get_last_test_lower()
      last_upper <- get_last_test_upper()
      if(is.null(last_upper) || is.null(last_test))
        stop("Either an Exact Last Test, or Both a Lower and Upper Bound for Tests Must Be Specified.")
    }else
      last_upper <- last_test
    age <- get_age()
    strata <- get_strata()

    testing_debut_age <- as.numeric(input$testing_debut_age)
    distribution <- tolower(input$distribution)
    age_breaks <- as.numeric(input$age_breaks)
    if(length(age_breaks) ==0)
      age_breaks <- NULL
    uniform_missreport <- input$uniform_missreport == "True"

    weights <- get_weights()
    if(is.null(weights))
      weights <- rep(1, length(report_pos))

    subset=NULL
    if(is.null(get_raw_data()))
      stop("No Data Uploaded")
    if(is.null(report_pos))
      stop("Required Variable Not Specified: Report Positive")
    if(is.null(hiv))
      stop("Required Variable Not Specified: HIV Status")
    if(is.null(ever_test))
      stop("Required Variable Not Specified: Ever Tested")
    if(is.null(last_test))
      stop("Required Variable Not Specified: Last Test")
    if(is.null(report_pos))
      stop("Required Variable Not Specified: Report Positive")
    sm_strata <- min(table(strata))
    if(!is.null(strata) && sm_strata < 5){
      stop("Stratifying Variable Has Strata With Too Few Observations")
    }
    if(is.null(biomarker_art))
      biomarker_art <- rep(FALSE, length(report_pos))
    if(is.null(low_viral))
      low_viral <- rep(FALSE, length(report_pos))
    inc <- list()
    if(is.null(strata)){
      inc[["all"]] <- testing_incidence(report_pos=report_pos,
                                        biomarker_art=biomarker_art,
                                        low_viral=low_viral,
                                        hiv=hiv,
                                        ever_test=ever_test,
                                        last_test=last_test,
                                        last_upper=last_upper,
                                        age=age,
                                        testing_debut_age=testing_debut_age,
                                        weights=weights,
                                        distribution=distribution,
                                        test_pop="negative",
                                        age_breaks=age_breaks,
                                        subset=subset,
                                        uniform_missreport=uniform_missreport)
      result <- inc[[1]]
      if(!is.null(age_breaks)){
        result <- cbind(row.names(result), result)
        colnames(result)[1] <- "age_subgroup"
      }
    }else{
      lvls <- levels(strata)
      for(lv in lvls){
        inc[[lv]] <- testing_incidence(report_pos=report_pos,
                                       biomarker_art=biomarker_art,
                                       low_viral=low_viral,
                                       hiv=hiv,
                                       ever_test=ever_test,
                                       last_test=last_test,
                                       last_upper=last_upper,
                                       age=age,
                                       testing_debut_age=testing_debut_age,
                                       weights=weights,
                                       distribution=distribution,
                                       test_pop="negative",
                                       age_breaks=age_breaks,
                                       subset=strata == lv,
                                       uniform_missreport=uniform_missreport)
      }
      result <- inc
      for(i in 1:length(inc)){
        result[[i]] <- cbind(names(result)[[i]], result[[i]])
        if(!is.null(age_breaks)){
          result[[i]] <- cbind(row.names(result[[i]]), result[[i]])
        }else
          row.names(result[[i]]) <- NULL
      }
      names(result) <- NULL
      result <- do.call(rbind,result)
      if(!is.null(age_breaks))
        colnames(result)[1:2] <- c("age_subgroup", "strata")
      else
        colnames(result)[1] <- "strata"
    }
    incidence(inc)
    boot_result(NULL)
    as.data.frame(result)
  }, rownames=FALSE,
  width="400px", digits=4)

  interruptor <- AsyncInterruptor$new()

  boot_result <- reactiveVal()
  observeEvent(input$run,{
    if(nclicks() != 0){
      print("Already running")
      return(NULL)
    }
    nclicks(nclicks() + 1)
    boot_result(data.frame(Status="Running..."))
    if(is.null(incidence()))
      return(NULL)
    type <- input$type
    nrep <- as.numeric(input$nrep)
    incc <- incidence()
    rep_weights <- get_rep_weights()
    design_strata <- get_design_strata()
    design_clusters <- get_design_clusters()
    weights <- get_weights()
    if(is.null(rep_weights) & (!is.null(design_strata) || !is.null(design_clusters))){
      print("Generating survey replicate weights from design")
      if(is.null(weights)){
        showNotification("Error: Survey design, but no weights given.")
        nclicks(0)
        return()
      }
      if(nrep > 5000){
        showNotification("Error: number of replicates must be < 5000 for survey bootstraps without replicate weights")
        nclicks(0)
        return()
      }
      if(is.null(design_strata) & !is.null(design_clusters)){
        df <- data.frame(design_clusters,weights)
        des <- svydesign(id = ~design_clusters,
                         weights = ~weights, data=df)
      }else if(!is.null(design_strata) & is.null(design_clusters)){
        df <- data.frame(design_strata,weights)
        des <- svydesign(id = ~1,
                         strata = ~ design_strata,
                         weights = ~weights, data=df)
      }else{
        df <- data.frame(design_clusters,design_strata,weights)
        des <- svydesign(id = ~design_clusters,
                         strata = ~ design_strata,
                         weights = ~weights, data=df)
      }
      des1 <- as.svrepdesign(des, type="bootstrap", compress=FALSE, replicates=nrep)
      rep_weights <- des1$repweights * weights
      type <- "bootstrap"
    }
    if(!is.null(rep_weights))
      nrep <- ncol(rep_weights)

    progress <- AsyncProgress$new(message="Generating Boostrap Samples")
    prog <- function(i, strata, nstrata){
      interruptor$execInterrupts()
      progress$set(( (strata - 1) * nrep + i) / (nrep*nstrata))
    }
    stratified <- length(incc) > 1
    result <- finally(
      catch(
        future({
          blist <- list()
          nstrata <- length(incc)
          for(i in 1:nstrata){
            callback <- function(k) prog(k, i, nstrata)
            if(is.null(rep_weights))
              boot <- as.data.frame(summary(bootstrap_incidence(incc[[i]],
                                                                nrep=nrep,
                                                                show_progress=callback)))
            else
              boot <- as.data.frame(summary(bootstrap_incidence(incc[[i]],
                                                                rep_weights=rep_weights,
                                                                type=type,
                                                                show_progress=callback)))
            if(stratified){
              boot <- cbind(names(incc)[i], boot)
              names(boot)[1] <- "strata"
            }
            if(nrow(boot) > 1){
              boot <- cbind(row.names(boot), boot)
              names(boot)[1] <- "age_subgroup"
            }
            blist[[i]] <- boot
          }
          do.call(rbind, blist)
        })  %...>% boot_result,
        function(e) {
          boot_result(NULL)
          print(e$message)
          showNotification(e$message)
        }
      ),
      function(){
        progress$sequentialClose()
        nclicks(0)
      }
    )

    NULL
  })
  output$bootstrap <- renderTable({
    req(boot_result())
  },digits=4)

  observeEvent(input$cancel,{
    print("cancel")
    interruptor$interrupt("User Interrupt")
  })

})
