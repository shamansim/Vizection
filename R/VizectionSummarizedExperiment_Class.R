setOldClass("dendrogram")
setOldClass("pca")
setOldClass("coa")

setClass(
  Class = "VizectionSummarizedExperiment",
  representation = representation(
    dendrogram = "dendrogram",
    pca = "pca",
    ca = "coa"),
  contains = "SummarizedExperiment"
)
