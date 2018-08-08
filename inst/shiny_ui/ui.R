library(shiny)
library(shinyWidgets)
# Define UI for application that draws a histogram
shinyUI(
   navbarPage('Testing History Incidence',
              #tag$script(HTML("updateFirstMultiInput = function(x){ var event = new Event('input'); $('.multi-wrapper .search-input').get(0).dispatchEvent(event);};Shiny.addCustomMessageHandler('updateFirstMultiInput', updateFirstMultiInput);")),
      tabPanel('Introduction',
      #includeScript("updateMulti.js"),
      h4('Welcome to The Testing History Incidence Tools'),
      br(),
      p("The purpose of this tool is the estimation of incidence from a single cross-sectional survey utilizing self-reported testing history."),
      br(),br(),
      p("You can get started with an example dataset available at:"),
      a("Download",href="https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/fellstat/TestingHistoryIncidence/blob/master/inst/shiny_ui/tstdat.csv"),
      br(),br(),p("Documentation for the application is available at:"),
      a("Documentation Wiki",href="https://github.com/fellstat/TestingHistoryIncidence/wiki/Shiny-App-Documentation"),
      br(),br(),p("and a youtube walkthrough can be found at"),
      a("Video",href=""),
      br(),br(),h3('Please proceed to the', em('Data'), 'tab')
    ),
    tabPanel('Data',
      fluidPage(
      # Application title
      titlePanel("Load Data"),

        # Sidebar with a slider input for number of bins
        sidebarLayout(
          sidebarPanel(
            fileInput("file1", "Choose CSV File",
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv")
            ),
            conditionalPanel("output.table != null",
            selectizeInput("hiv", "HIV+:",
                        c("")),
            selectizeInput("report_pos", "Report Diagnosed:",
                        c("")),
            selectizeInput("ever_test", "Report Having Had A Previous HIV Test:",
                        c("")),
            tabsetPanel(
              tabPanel("Exact",
                selectizeInput("last_test", "Time Since Last HIV Test (Months):",
                            c(""))
              ),
              tabPanel("Bounded",
                selectizeInput("last_test_lower", "Time Since Last HIV Test Lower Bound (Months):",
                             c("")),
                selectizeInput("last_test_upper", "Time Since Last HIV Test Upper Bound (Months):",
                               c(""))
              )
            ),tags$hr(),
            selectizeInput("biomarker_art", "ART Biomarker (Optional):",
                        c("")),
            selectizeInput("low_viral", "Low/Undetectable Viral Load (Optional):",
                        c("")),
            selectizeInput("age", "Age (Optional):",
                           c("")),
            selectizeInput("weights", "Weights (Optional):",
                           c("")),
            selectizeInput("strata", "Stratify By (Optional):",
                           c("")),
            pickerInput(inputId = "rep_weights",
                        label = "Replication Weights (Optional):",
                        choices = c(""),
                        multiple=TRUE,
                        options = list(`actions-box` = TRUE,
                                       `live-search`=TRUE,
                                       `none-selected-text`="Choose Variable"))
          ),
          width=5),
          mainPanel(
            conditionalPanel("input.hiv != \"\"",
                             h2("HIV Descriptives:"),
                             p("Raw Values:"),
                             tableOutput("hiv_desc_raw"),
                             p("Processed:"),
                             tableOutput("hiv_desc")
            ),
            conditionalPanel("input.report_pos != \"\"",
                             h2("Report Positive Descriptives:"),
                             p("Raw Values:"),
                             tableOutput("report_pos_desc_raw"),
                             p("Processed:"),
                             tableOutput("report_pos_desc")
            ),
            conditionalPanel("input.ever_test != \"\"",
                             h2("Ever Tested Descriptives:"),
                             p("Raw Values:"),
                             tableOutput("ever_test_desc_raw"),
                             p("Processed:"),
                             tableOutput("ever_test_desc")
            ),
            conditionalPanel("input.biomarker_art != \"\"",
                             h2("ART Biomarker Positive Descriptives:"),
                             p("Raw Values:"),
                             tableOutput("biomarker_art_desc_raw"),
                             p("Processed:"),
                             tableOutput("biomarker_art_desc")
            ),
            conditionalPanel("input.low_viral != \"\"",
                             h2("Low Viral Load Descriptives:"),
                             p("Raw Values:"),
                             tableOutput("low_viral_desc_raw"),
                             p("Processed:"),
                             tableOutput("low_viral_desc")
            ),
            conditionalPanel("input.last_test != \"\"",
                             h2("Last Test Descriptives:"),
                             plotOutput("last_test_plot", width = "100%", height = "200px"),
                             verbatimTextOutput("last_test_errors")
            ),
            conditionalPanel("input.last_test_lower != \"\" & input.last_test_upper != \"\"",
                             h2("Last Test (Bounded) Descriptives:"),
                             tableOutput("last_test_bound_desc")
            ),
            conditionalPanel("input.age != \"\"",
              h2("Age Descriptives:"),
              plotOutput("age_desc_plot", width = "100%", height = "200px"),
              verbatimTextOutput("age_desc_errors")
            ),
            dataTableOutput("table"),
            width = 7
          )
        )
      )
    ),
    tabPanel('Analysis',
             fluidPage(
               # Application title
               titlePanel("Analysis"),

               # Sidebar with a slider input for number of bins
               sidebarLayout(
                 sidebarPanel(
                   radioButtons(
                     inputId = "distribution",
                     label = "Distribution:",
                     choices = c("Weibull", "Empirical")
                   ),
                   conditionalPanel("input.age != \"\"",
                     numericInput("testing_debut_age",
                                  "Age of Testing Debut:",
                                  0, min=0),
                     selectizeInput("age_breaks", "Age Break Points:",
                                    as.character(1:100),
                                    multiple=TRUE)
                   ),
                   conditionalPanel("input.age_breaks != null || input.strata != \"\"",
                     radioButtons(
                       inputId = "uniform_missreport",
                       label = "Uniform Missreporting of Undiagnosed Status:",
                       choices = c("False", "True")
                     )
                   ),
                   wellPanel(
                    h3("Bootstrap Intervals"),
                    conditionalPanel("input.rep_weights == null",
                      numericInput("nrep",
                                "# of Bootstraps:",
                                50, min=2)
                    ),
                    conditionalPanel("input.rep_weights != null",
                      selectInput("type",
                        "Replicate Weight Type:",
                        choices = c("Bootstrap"="bootstrap","Jackknife"="JK1","BRR", "Fay"),
                        selected="bootstrap"
                      )
                    ),
                    actionButton('run', 'Run'),
                    actionButton('cancel', 'Cancel')
                   )


                 ),
                 mainPanel(
                   h3("Incidence Results"),
                   tableOutput("inc_results"),
                   h3("Bootstrap Intervals"),
                   tableOutput("bootstrap")
                 )
              )
             )
    )
  )
)
