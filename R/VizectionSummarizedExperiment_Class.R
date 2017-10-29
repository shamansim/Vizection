# Class definition -------------------------------------------------------------

setOldClass("dendrogram")
setOldClass("pca")
setOldClass("coa")

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
