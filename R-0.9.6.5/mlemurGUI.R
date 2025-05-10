#' Run mlemur in a browser window
#' 
#' @description This function initializes mlemur in graphical mode.
#' @param options arguments passed to shiny::shinyApp(options=...)
#' @export
mlemur <- function(options=list("display.mode" = "normal", "launch.browser" = T)) {
  shiny::shinyApp(ui=mlemurUI, server=mlemurServer, options=options)
}
# mlemurGUI <- function() {
#   shiny::shinyApp(ui=mlemurUI, server=mlemurServer, options=list("display.mode" = "normal", "launch.browser" = T))
# }#' 