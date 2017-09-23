#' Calculate a correlation matrix.
#' 
#' @param genes An expression table of discrete counts
#'              (tags, molecules, ...).
#'
#' The first row is removed, because it is expected to contain
#' the total count of the reads (or molecules, ...) that did not
#' match an annotation.
#' 
#' @seealso smallCAGEqc::TPM
#' 
#' @examples 
#' data.frame(1:3, 2:4, 6:4) %>% corMat()

corMat <- function(genes)
  corMat_1(genes) %>% corMat_2 %>% corMat_3

corMat_1 <- function(genes)
  genes %>% tail(-1) %>% smallCAGEqc::TPM()

corMat_2 <- log1p

corMat_3 <- cor


#' Convert correlation to distance
#' 
#' Transforms a correlation matrix into a Euclidian distance matrix.
#' 
#' The final conversion is done with the quasieuclid function of the
#' ade4 package.
#' 
#' @param m A corelation matrix.
#' 
#' @examples 
#' data.frame(1:3, 2:4, 6:4) %>% corMat %>% distCorMat

distCorMat <- function(m)
  m %>% distCorMat_1 %>% distCorMat_2

distCorMat_1 <- function(m) {
  m %>%
    subtract(1, .) %>%
    divide_by(., 2) %>% 
    as.dist
}
  
distCorMat_2 <- ade4::quasieuclid


#' Complete hierarchical clustering
#' 
#' Cluster a distance matrix with the complete method.
#' 
#' @param d a distance matrix
#' 
#' @examples
#' data.frame(1:3, 2:4, 6:4) %>% corMat %>% distCorMat %>% genesDend

genesDend <- function(d)
    d %>% hclust(method = "complete")


#' Cluster a distance matrix with the complete method.
#' 
#' @param d a distance matrix
#' @param x A list where the element "nbClusters" is the number of
#'          clusters to compute and "nbGroups" is the number of
#'          labeled groups.
#' 
#' @examples
#' dendr <- data.frame(1:3, 2:4, 6:4) %>% corMat %>% distCorMat %>% genesDend()
#' genesDend2(dendr, x = vizectionExampleEnv())
#' genesDend2(dendr, x = vizectionExampleEnv() %>% inset("showGroupsColor", FALSE))
#'   
#' @importFrom grDevices rainbow
#' @importFrom colorspace rainbow_hcl

genesDend2 <- function(d, x)
  genesDend2_4( d
                , x
                , nbGroups = genesDend2_1(x)
                , colGrps = genesDend2_2(genesDend2_1(x))
                , cols = genesDend2_3(x)) %>%
  genesDend2_5()

genesDend2_1 <- function(x)
  length(x$groupsCheck)
  
genesDend2_2 <- grDevices::rainbow

genesDend2_3 <- function(x)
  colorspace::rainbow_hcl(x$nbClusters, c=50, l=100)

genesDend2_4 <- function(d, x, nbGroups, colsGrps, cols) {
  # In order to pipe ifelse
  ife <- function(cond, x, y) {
    if(cond) return(x) 
    else return(y)
  }
  
  d %>%
    as.dendrogram %>%
    set("branches_k_color", k = x$nbClusters, with = cols) %>%
      { 
        ife(x$showGroupsColor ,
            set(., "labels_colors", k = nbGroups, with = colsGrps),
            set(., "labels_colors", k = x$nbClusters, with = genes)
        )
      } 
  
}

genesDend2_5 <- function(dend)
    dend %>% dendextend::ladderize(FALSE)


#' Plots Vizection's heatmap of correlations.
#' 
#' @param cormat A correlation matrix like the output of corMat().
#' @param dendr A dendrogram object like the ouptut of gendsDend2().
#' @param sublibs A "libs" table like the output of vizectionExampleLibs().
#' 
#' @examples
#' cormat <- vizectionExampleGenes() %>% corMat
#' dendr <- cormat %>% distCorMat %>% genesDend %>% genesDend2(x = vizectionExampleEnv())
#' contentheatmapGenes(cormat, dendr, vizectionExampleLibs())

contentheatmapGenes <- function(cormat, dendr, sublibs)
  NMF::aheatmap( cormat
               , annCol = list(Run=sublibs$Run, Group=sublibs$group)
               , Rowv   = dendr
               , Colv   = dendr)
