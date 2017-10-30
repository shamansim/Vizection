# Class definition -------------------------------------------------------------

setOldClass("dendrogram")
setOldClass("pca")
setOldClass("coa")

#' VizectionSummarizedExperiment, S4 class that extends \code{"SummarizedExperiment"} class.
#'
#' Extension is to add the results of the analyses runned in Vizection as well
#' as recording session version (new feature to come).
#'
#' @slot dendrogram The dendrogram computed in vizection.
#' @slot pca A principal components analysis runned with \code{ade4} package.
#' @slot coa A correspondance analysis runned with \code{ade4} package.
#' @slot vizectionParams List of parameters used as metadata for running Vizection sessions.
#'
#' @name VizectionSummarizedExperiment
#' @export
setClass(
  Class = "VizectionSummarizedExperiment",
  representation = representation(
    dendrogram = "dendrogram",
    pca = "pca",
    coa = "coa",
    vizectionParams = "list"),
  contains = "SummarizedExperiment"
)

# Constructor ------------------------------------------------------------------

vizectionSummarizedExperiment <- function(se){
  if(class(se) == "SummarizedExperiment"){
    initialDendrogram <- FALSE
    initialPca <- FALSE
    initialCoa <- FALSE
    class(initialDendrogram) <- "dendrogram"
    class(initialPca) <- "pca"
    class(initialCoa) <- "coa"
    new(Class = "VizectionSummarizedExperiment",
        dendrogram = initialDendrogram,
        pca = initialPca,
        coa = initialCoa,
        colData = colData(se),
        assays = se@assays,
        NAMES = names(se),
        elementMetadata = elementMetadata(se),
        metadata = metadata(se)
    )
  } else{
    stop("Data provided is not of Class 'SummarizedExperiment'.")
  }
}
