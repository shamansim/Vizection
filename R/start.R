#' Launch the Vizection shiny application.
#' 
#' By default, Vizection will look for two data frames called \dQuote{genes}
#' and \dQuote{libs}.  Alternatively, these data frames can be passed as
#' arguments to the \code{vizection()} function.  Their names will be used to
#' update the global options (\code{vizection.genes} and \code{vizection.libs})
#' that Vizection uses to find the objects in the global environment.
#' 
#' @param local T if you want host to be "0.0.0.0"
#' @import DT ade4 adegraphics colorspace data.table dendextend dplyr ggplot2
#'   magrittr shiny shinydashboard smallCAGEqc
#' @export
#' @examples
#' vizection(genes, libs, F)

vizection <- function(genes = genes, libs = libs, local = F) {
  
  options("vizection.genes" = deparse(substitute(genes)))
  options("vizection.libs"  = deparse(substitute(libs)))
  
  appDir <- system.file("shiny", "vizection", package = "vizection")
  if (appDir == "") {
    stop("Could not find vizection directory. Try re-installing `vizection`.", call. = FALSE)
  }
  
  ifelse(local, runApp(appDir, host="0.0.0.0"), runApp(appDir))

}

.onLoad <- function(libname, pkgname){
  options("vizection.genes" = "genes")
  options("vizection.libs"  = "libs")
}