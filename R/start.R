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
#' \dontrun{
#' library("airway")
#' data("airway")
#' se <- airway
#' vizection(se)
#' 
#' 
#' # Older way to run Vizection:
#' data(iris)
#' genes <- iris[, 1:4]
#' libs <- as.data.frame(iris[, "Species"])
#' colnames(libs) <- c('group')
#' libs$samplename <- as.character(1:nrow(genes))
#' rownames(libs) <- libs$samplename 
#' libs$counts <- rnorm(n = nrow(genes), mean = 1000, sd = 200)
#' genes <- as.data.frame(t(genes))
#' vizection(genes = genes, libs = libs, local = FALSE)
#' }

vizection <- function(se, genes, libs, local = F) {

  if ( !missing(se) & !(missing(genes) & missing(libs)))
    stop( "Too much input: ", dQuote("SummarizedExperiment"), " objects should be provided "
        , "without ", dQuote("genes"), " or ", dQuote("libs"), " data.frames.")
  
  options("vizection.dataIs" = "se") # Search for SummarizedExperiment by default.
  
  if (missing(se) & !missing(genes) & !missing(libs)) {
    warning( "Using ", dQuote("genes"), " and ", dQuote("libs"), " data.frames "
           , "is deprecated and will be removed in the future.")
    options("vizection.dataIs" = "genesAndLibs") # Search for genes and libs data frames.
  }
  
  options("vizection.genes" = deparse(substitute(genes)))
  options("vizection.libs"  = deparse(substitute(libs)))
  options("vizection.se"    = deparse(substitute(se)))
  
  appDir <- system.file("shiny", "vizection", package = "vizection")
  if (appDir == "") {
    stop("Could not find vizection directory. Try re-installing `vizection`.", call. = FALSE)
  }
  
  ifelse(local, runApp(appDir, host="0.0.0.0"), runApp(appDir))

}

.onLoad <- function(libname, pkgname){
  options("vizection.genes" = "genes")
  options("vizection.libs"  = "libs")
  options("vizection.se"    = "se")
}
