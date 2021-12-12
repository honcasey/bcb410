#' Run Shiny App for cellLineConsistency
#'
#' A function that launches the Shiny app for cellLineConsistency.
#'
#' @return No return value, opens up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' cellLineConsistency::runCellLineConsistency()
#' }
#'
#' @references
#' Wickham, H. (2020). Mastering Shiny. O'Reilly Media. \href{https://mastering-shiny.org/basic-app.html}{Link}
#'
#' @importFrom shiny runApp
#' @export
#'

runCellLineConsistency <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "cellLineConsistency")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

# [END]
