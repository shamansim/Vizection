
vizection <- function(genes, libs) {
  
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
  
  # In order to pipe ifelse
  ife <- function(cond, x, y) {
    if(cond) return(x) 
    else return(y)
  }
  
  showDendrColors <- function(dendro){
    dendrapply(dendro, function(X){
      if(is.leaf(X)){
        attr(X, "edgePar")[1]
      }
    }) %>% unlist
  }
  
  pcaCompOrientation <- function(compPca){
    ifelse(abs(range(compPca)[1]) > abs(range(compPca)[2]), FALSE, TRUE)
  }
  
  pcaCompGenesList <- function(pcaAde4co, comp){
    stopifnot(ncol(pcaAde4co) == 3)
    
    genesCo <- pcaAde4co %>%
      mutate(., geneNames = rownames(.)) %>% 
      select(geneNames, Comp1, Comp2, Comp3)
    
    ifelse(pcaCompOrientation(genesCo[comp+1]),
      genesCo %<>% setorderv(., colnames(.)[comp+1], order=-1),
      genesCo %<>% setorderv(., colnames(.)[comp+1], order=1))
    
    genesCo
  }
  
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
  
  app <- shinyApp(
    ui = source("ui.R")$value,
    server = source("server.R")$value
  )
  
  runApp(app, host="0.0.0.0")
  
}
