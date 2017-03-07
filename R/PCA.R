#' contentgenesPCA
#' 
#' Vizection's Principal Component Analysis
#' 
#' First, the expression table, assumed to contain counts,
#' is normalised to tags per million (TPM) and transposed.
#' Then, the first line, containing the combined expression
#' of all the CTSS that were not contained in an annotated
#' feature (such as a gene), is removed.  Then the PCA is
#' run with the ade4 package's function dudi.pca().
#' 
#' @param df An expression table (see above)
#' 
#' @example 
#' vizection:::vizectionExampleGenes() %>% contentgenesPCA
#' 
#' @importFrom ade4 dudi.pca
#' @export contentgenesPCA

contentgenesPCA_1 <- function(subgenes)
  subgenes %>% smallCAGEqc::TPM() %>% t

contentgenesPCA_2 <- function (genesTpm)
  ade4::dudi.pca( genesTpm[, -1]
                , center = T
                , scale  = F
                , scannf = F
                , nf     = 3)
  
contentgenesPCA <- function(df)
  df %>% contentgenesPCA_1 %>% contentgenesPCA_2