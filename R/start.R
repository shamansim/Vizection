#' Vizection
#'
#' This function launches the Vizection shiny application.
#' @param local T if you want host to be "0.0.0.0"
#' @import DT ade4 adegraphics colorspace data.table dendextend dplyr ggplot2
#'   magrittr plotly shiny shinydashboard smallCAGEqc
#' @export
#' @examples
#' vizection(genes, libs, F)
vizection <- function(genes = genes, libs = libs, local = F) {
  
  appDir <- system.file("shiny", "vizection", package = "vizection")
  if (appDir == "") {
    stop("Could not find vizection directory. Try re-installing `vizection`.", call. = FALSE)
  }
  
  ifelse(local, runApp(appDir, host="0.0.0.0"), runApp(appDir))

}
