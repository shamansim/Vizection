#' Vizection
#'
#' This function launches the Vizection shiny application.
#' @param local T if you want host to be "0.0.0.0"
#' @keywords
#' @export
#' @examples
#' vizection(genes, libs, F)
vizection <- function(genes = genes, libs = libs, local = F) {
  
  library(shiny)
  library(shinydashboard)
  library(ade4)
  library(adegraphics)
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(colorspace)
  library(smallCAGEqc)
  library(dendextend)
  library(DT)
  library(plotly)
  
  appDir <- system.file("shiny", "vizection", package = "vizection")
  if (appDir == "") {
    stop("Could not find vizection directory. Try re-installing `vizection`.", call. = FALSE)
  }
  
  ifelse(local, runApp(appDir, host="0.0.0.0"), runApp(appDir))

}
