# loading required library
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

shinyServer(function(input, output, session) {
  
  # FILTERS
  # =======
  
  UNaddNumberOfSamplesOrGroup <- function(Checklist){
    result <- c()
    for(i in Checklist){
      result <- c(result, i %>% sapply(., function(x) gsub("\\s[:|:]\\s.*", "", x)))
    }
    return(result)
  }
  
  contentfilterSelectionBool <- reactive({
    withProgress(message = 'Updating pre-filter', {
      incProgress(1/2, detail = "updating")
      (libs$counts > input$nbFilterExtracted) &
        (rownames(libs) %in% (libs %>% select(samplename, group) %>% filter(group %in% UNaddNumberOfSamplesOrGroup(input$groupsCheck)) %$% samplename))
    })
  })
  filterSelectionBool <- reactive({
    contentfilterSelectionBool()
  })

  contentfilterSelectionBoolFinal <- reactive({
    withProgress(message = 'Updating filter', {
      incProgress(1/2, detail = "updating")
      filterSelectionBool() &
        (libs$samplename %in% UNaddNumberOfSamplesOrGroup(paste(input$samplesCheck)))
    })
  })
  filterSelectionBoolFinal <- reactive({
    contentfilterSelectionBoolFinal()
  })
  
  # SUBGENES SUBLIBS
  # ================
  
  subgenes <- reactive({
    withProgress(message = 'Updating subgenes', {
      incProgress(1/3, detail = "filtering")
      pre_subgenes <- genes[, filterSelectionBoolFinal()]
      incProgress(2/3, detail = "removing useless genes")
      pre_subgenes[apply(pre_subgenes, 1, sum) != 0, ] #removing useless genes
    })
  })

  sublibs <- reactive({
    withProgress(message = 'Updating sublibs', {
      incProgress(1/2, detail = "filtering")
      sublibs0 <- libs[filterSelectionBoolFinal(), ]
      sublibs0$group %<>% extract(drop = T)
      sublibs0
    })
  })
  
  # -> libsGroup
  addNumberOfSamples <- function(listOfGroups){
    result <- c()
    for(i in listOfGroups){
      result <- c(result, paste0(i, " | ", libs$samplename[libs$group==i] %>% length))
    }
    return(result)
  }
  contentlibsGroup <- reactive({
    withProgress(message = 'updating groups', {
      incProgress(1/3, detail = "extracting from filter")
      filterExtractedBool <- libs$counts > input$nbFilterExtracted
      incProgress(2/3, detail = "creating checkbox")
      myGroups <- addNumberOfSamples(paste(unique(libs$group[filterExtractedBool])))
      checkboxGroupInput(inputId = "groupsCheck", label = "",
        choices = myGroups,
        selected = myGroups
      )
    })
  })
  output$libsGroup <- renderUI({
    contentlibsGroup()
  })

  observe({
    filterExtractedBool <- libs$counts > input$nbFilterExtracted
    myGroups <- addNumberOfSamples(paste(unique(libs$group[filterExtractedBool])))
    updateCheckboxGroupInput(session,
      "groupsCheck",
      choices = myGroups,
      selected = if(input$bar) myGroups
    )
  })
  
  # -> libsSamplename
  addgroup <- function(listOfSamples){
    result = c()
    for(i in listOfSamples){
      result <- c(result, paste0(i, " | ", libs$group[libs$samplename==i]))
    }
    return(result)
  }
  contentlibsSamplename <- eventReactive(input$updateSamples, {
    withProgress(message = 'updating samples', {
      incProgress(1/3, detail = "extracting selection")
      filterSelectionNames <- rownames(libs)[filterSelectionBool()]
      incProgress(2/3, detail = "creating checkbox")
      mySamples <- addgroup(paste(filterSelectionNames))
      checkboxGroupInput(inputId = "samplesCheck", label = "",
        choices = mySamples,
        selected = mySamples
      )
    })
  })
  output$libsSamplename <- renderUI({
    contentlibsSamplename()
  })
  
  # SHARED
  # ======
  
  corMat <- eventReactive(input$updateCorMat, {
    withProgress(message = 'correlation matrice', {
      incProgress(1/4, detail = "TPM")
      a <- subgenes() %>% extract(-1, ) %>% TPM
      incProgress(2/4, detail = "log1p")
      b <- a %>% log1p
      incProgress(3/4, detail = "cor")
      b %>% cor
    })
  })
  
  distCorMat <- reactive({
    withProgress(message = 'distance matrice', value = 0, {
      incProgress(1/3, detail = "as.dist")
      a <- corMat() %>%
        subtract(1, .) %>%
        divide_by(., 2) %>% 
        as.dist 
      incProgress(2/3, detail = "quasieuclid")
      a %>% quasieuclid
    })
  })
  
  genesDend <- reactive({
    withProgress(message = 'cluster', value = 0, {
      incProgress(1/2, detail = "hclust")
      distCorMat() %>%
        hclust(method = "complete")
    })
  })
  
  genesDend2 <- reactive({
    withProgress(message = 'dendrogram', {
      incProgress(1/6, detail = "nbGroups")
      nbGroups <- length(input$groupsCheck)
      incProgress(2/6, detail = "colGroups")
      colsGrps <- rainbow(nbGroups)
      incProgress(3/6, detail = "colors")
      cols <- rainbow_hcl(input$nbClusters, c=50, l=100)
      incProgress(4/6, detail = "customization")
      a <- genesDend() %>% as.dendrogram %>%
        set("branches_k_color", k = input$nbClusters, with = cols) %>%
        { 
          ife(input$showGroupsColor ,
            set(., "labels_colors", k = nbGroups, with = colsGrps),
            set(., "labels_colors", k = input$nbClusters, with = cols)
          )
        } 
      incProgress(5/6, detail = "ladderize")
      a %>% ladderize(FALSE)
    })
  })
  
  colorsPcaLi <- reactive({
    withProgress(message = 'colors PCA', {
      incProgress(1/3, detail = "collecting nb clusters")
      ifelse(input$nbClusters!= 1, palette(rainbow_hcl(input$nbClusters, c=50, l=100)), palette(rainbow_hcl(2, c=50, l=100)))
      incProgress(2/3, detail = "generating colors")
      data.frame(colors = showDendrColors(genesDend2()), sampleIndex = order.dendrogram(genesDend2())) %>%
        setorder("sampleIndex") %$%
        return(colors)
    })
  })
  
  # HEADER
  # ======
  
  contentgeneral <- eventReactive(input$updateSelection, {
      withProgress(message = 'Updating selection information', {
        incProgress(1/5, detail = "collecting sublibs")
        sublibs <- sublibs()
        incProgress(2/5, detail = "collecting subgenes")
        subgenes <- subgenes()
        incProgress(3/5, detail = "generating dataframe")
        data <- data.frame(
          group = c("Samples", "Groups", "Genes"),
          value = c(sum(filterSelectionBoolFinal()) / nrow(libs) * 100,
            length(unique(sublibs$group)) / length(unique(libs$group)) * 100,
            nrow(subgenes[-1,]) / nrow(genes[-1,]) * 100),
          total = c(nrow(libs), length(unique(libs$group)), nrow(genes[-1,])),
          selection = c(sum(filterSelectionBoolFinal()), length(unique(sublibs$group)), nrow(subgenes[-1,]))
        )
        incProgress(4/5, detail = "final process")
        list(
          samples = data %>% filter(group == "Samples") %$% value %>% round(digits = 2),
          groups = data %>% filter(group == "Groups") %$% value %>% round(digits = 2),
          genes = data %>% filter(group == "Genes") %$% value %>% round(digits = 2),
          totalSamples = data %>% filter(group == "Samples") %$% total %>% round(digits = 2),
          totalGroups = data %>% filter(group == "Groups") %$% total %>% round(digits = 2),
          totalGenes = data %>% filter(group == "Genes") %$% total %>% round(digits = 2),
          selectionSamples = data %>% filter(group == "Samples") %$% selection %>% round(digits = 2),
          selectionGroups = data %>% filter(group == "Groups") %$% selection %>% round(digits = 2),
          selectionGenes = data %>% filter(group == "Genes") %$% selection %>% round(digits = 2)
        )#list
      })
    })
  contenttasksMenu <- reactive ({
    general <- contentgeneral()
    dropdownMenu(type = "tasks", badgeStatus = "success",
      taskItem(value = general$samples, color = "blue",
        paste("Samples: ", general$selectionSamples, "/", general$totalSamples)
      ),
      taskItem(value = general$groups, color = "green",
        paste("Groups: ", general$selectionGroups, "/", general$totalGroups)
      ),
      taskItem(value = general$genes, color = "red",
        paste("Genes: ", general$selectionGenes, "/", general$totalGenes)
      )
    )
  })
  output$tasksMenu <- renderMenu({
    contenttasksMenu()
  })
  
  # HOME
  # ====
  
  output$UIboxplotGroupsSub <- renderUI({
    withProgress(message = 'Updating boxplot list', {
      incProgress(1/2, detail = "parsing sublibs")
      selectInput("boxplotGroupsSub", "Groups (selection)", c("none", paste(unique(sublibs() %>% select(group) %>% extract(,1)))), selected = "none")
    })
  })
  
  # -> boxplotTotal
  output$boxplotTotal <- renderPlot({
    if(input$boxplotGroupsTotal != "none"){
      withProgress(message = 'Updating total boxplot', {
        incProgress(1/3, detail = "collecting sublibs")
        sublibs <- sublibs()
        incProgress(2/3, detail = "generating")
        ggplot(data = libs[libs$group == input$boxplotGroupsTotal, ], aes(input$boxplotGroupsTotal, counts)) +
          geom_boxplot() +
          xlab("") + ylab("") +
          ylim(c(0, max(libs[libs$group == input$boxplotGroupsTotal, "counts"], sublibs[sublibs$group == input$boxplotGroupsSub, "counts"]))) +
          theme_minimal()
      })
    }
  })
  
  output$boxplotSub <- renderPlot({
    if(input$boxplotGroupsSub != "none"){
      withProgress(message = 'Updating sub boxplot', {
        incProgress(1/3, detail = "collecting sublibs")
        sublibs <- sublibs()
        incProgress(2/3, detail = "generating")
        ggplot(data = sublibs[sublibs$group == input$boxplotGroupsSub, ], aes(input$boxplotGroupsSub, counts)) +
          geom_boxplot() +
          xlab("") + ylab("") +
          ylim(c(0, max(libs[libs$group == input$boxplotGroupsTotal, "counts"], sublibs[sublibs$group == input$boxplotGroupsSub, "counts"]))) +
          theme_minimal()
      })
    }
  })
  
  # DENDROGRAM
  # ==========
  
  contentdendrogram <- eventReactive(input$updateDendrogram, {
    withProgress(message = 'dendrogram plot', value = 0, {
      incProgress(1/3, detail = "modifying display parameters")
      par(mar = c(6,2,2,6))
      incProgress(2/3, detail = "generating plot")
      genesDend2() %>%
        dendextend::set("labels_cex", input$dendroSize) %>%
        plot(horiz = input$dendroHoriz)
    })
  })
  output$dendrogram <- renderPlot({
    contentdendrogram()
  })
  
  output$dendrogramPlot <- renderUI({
    plotOutput("dendrogram", height = paste0(input$heightDendro,"px"))
  })
  
  contentheight <- eventReactive(input$updateheight, {
    withProgress(message = 'dendrogram height', value = 0, {
      incProgress(1/3, detail = "collecting dendrogram")
      genesDend <- genesDend() 
      genesDendRev <- rev(genesDend$height)
      incProgress(2/3, detail = "drawing plot")
      plot(genesDendRev[1:input$heightlength], pch = 20, ylab = "Clusters height")
      abline(v = input$nbClusters + 0.5, col = "red", lty = 2)
      for(i in genesDendRev){abline(h = i, lty= 2, col = "grey")}
    })
  })
  output$height <- renderPlot({
    contentheight()
  })
  
  # HEATMAP
  # =======
  
  contentheatmapGenes <- eventReactive(input$updateHeatmap, {
    withProgress(message = 'heatmap', value = 0, {
      incProgress(1/2, detail = "construction")
      sublibs <- sublibs()
      
      NMF::aheatmap(corMat(),
                    annCol=list(Run=sublibs$Run, Group=sublibs$group),
                    Rowv = genesDend2(), Colv = genesDend2())
    })
  })
  output$heatmapGenes <- renderPlot({
    contentheatmapGenes()
  })
  
})
