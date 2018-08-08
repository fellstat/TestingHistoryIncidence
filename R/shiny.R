#' Run Shiny UI
#' @details
#' Runs a Shiny application
#' @export
shiny_testing_history <- function() {
  appDir <- system.file("shiny_ui", package = "TestingHistoryIncidence")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `TestingHistoryIncidence`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
