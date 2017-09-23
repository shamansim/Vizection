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
#' @examples 
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


#' pcaCompGenesList
#' 
#' It's complicated...
#' 
#' @param components Principal components of a PCA object of class
#'        "dudi", extracted by subsetting it for its "co" member.
#' @param n Rank number of the principal component to select
#' 
#' @return The "components" data frame, sorted on the "n"th component.
#' 
#' @examples 
#' vizection:::vizectionExampleGenes() %>% contentgenesPCA %>%
#'   extract2("co") %>% pcaCompGenesList(1)
#' 
#' @importFrom dplyr mutate select

pcaCompGenesList <- function(components, n){
  stopifnot(ncol(components) == 3)
  
  genesCo <- components %>%
    dplyr::mutate(., geneNames = rownames(.)) %>% 
    dplyr::select(geneNames, Comp1, Comp2, Comp3)
  
  ifelse( abs(min(genesCo[n+1])) <= abs(max(genesCo[n+1]))
        , genesCo %<>% setorderv(., colnames(.)[n+1], order = -1)
        , genesCo %<>% setorderv(., colnames(.)[n+1], order = 1))
  genesCo
}


#' Plot Eigen Values
#' 
#' Plot the Eigen values of the "dudi" PCA object produced by contentgenesPCA.
#' 
#' @param orderedCompPca a The PCA object produced by contentgenesPCA.
#' 
#' @examples 
#' vizection:::vizectionExampleGenes() %>%
#'   contentgenesPCA %>%
#'   plotEigenValues
#'   

plotEigenValues <- function(orderedCompPca)
  barplot(orderedCompPca %$% eig, xlab = "Eigenvalues")


#' plot head and tail as barplot.
#' 
#' @param orderedCompPca A data frame where columns are principal component,
#'        the rows are genes or similar features, and the component to
#'        plot is already sorted.
#' @param comp The number of the principal component in that matrix.
#' @param nbDispGenes How many bars for the head and the tail (each).
#' 
#' @examples
#' vizection:::vizectionExampleGenes() %>% contentgenesPCA %>%
#'   extract2("co") %>% pcaCompGenesList(1) %>%
#'   plotHTB(1)
#' 
#' @export plotHTB

plotHTB <- function(orderedCompPca, comp, nbDispGenes = 25){
  par(mfrow=c(1, 2))
  
  bp1h <- orderedCompPca[, comp+1] %>% 
    head(nbDispGenes) %>%
    barplot(. ,
            ylim = c(min(orderedCompPca[, comp+1]), max(orderedCompPca[, comp+1])),
            axes = FALSE, axisnames = FALSE, main = paste0("comp ", comp, " head"))
  text(bp1h, par("usr")[3], labels = orderedCompPca$geneNames %>% head(nbDispGenes),
       srt = 90, adj = c(1.1,1.1), xpd = TRUE, cex=1)
  axis(2)
  
  bp1t <- orderedCompPca[, comp+1] %>%
    tail(nbDispGenes) %>%
    barplot(. ,
            ylim = c(min(orderedCompPca[, comp+1]), max(orderedCompPca[, comp+1])),
            axes = FALSE, axisnames = FALSE, main = paste0("comp ", comp, " tail"))
  text(bp1t, par("usr")[3], labels = orderedCompPca$geneNames %>% tail(nbDispGenes),
       srt = 90, adj = c(1.1,1.1), xpd = TRUE, cex=1)
  axis(4)
  
  par(mfrow=c(1, 1))
}
