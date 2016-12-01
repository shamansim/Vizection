vizection <- function(genes = genes, libs = libs) {
  
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
  
  app <- shinyApp(
    ui = source("ui.R")$value,
    server = source("server.R")$value
  )
  
  runApp(app)#, host="0.0.0.0")
  
}
