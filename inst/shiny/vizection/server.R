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

genes <- get(getOption("vizection.genes"), .GlobalEnv)
libs  <- get(getOption("vizection.libs"),  .GlobalEnv)

vizectionValidate(genes = genes, libs = libs)

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
  
  filterSelectionBool <- reactive({
    withProgress(message = 'Updating pre-filter', {
      incProgress(1/2, detail = "updating")
      vizection:::filterSelectionBool(libs, input)
    })
  })

  filterSelectionBoolFinal <- reactive({
    withProgress(message = 'Updating filter', {
      incProgress(1/2, detail = "updating")
      vizection:::filterSelectionBoolFinal(libs, input)
    })
  })
  
  # SUBGENES SUBLIBS
  # ================
  
  subgenes <- reactive({
    withProgress(message = 'Updating subgenes', {
      incProgress(1/3, detail = "filtering")
      pre_subgenes <- vizection:::subgenes_1(genes, input)
      incProgress(2/3, detail = "removing useless genes")
      vizection:::subgenes_2(pre_subgenes) #removing useless genes
    })
  })

  sublibs <- reactive({
    withProgress(message = 'Updating sublibs', {
      incProgress(1/2, detail = "filtering")
      vizection:::sublibs(input)
    })
  })
  
  # -> libsGroup
  
  contentlibsGroup <- reactive({
    withProgress(message = 'updating groups', {
      incProgress(1/3, detail = "extracting from filter")
      filterExtractedBool <- vizection:::filterExtractedBool(libs, input)
      incProgress(2/3, detail = "creating checkbox")
      myGroups <- vizection:::addNumberOfSamples(libs, paste(unique(libs$group[filterExtractedBool])))
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
    filterExtractedBool <- vizection:::filterExtractedBool(libs, input)
    myGroups <- vizection:::addNumberOfSamples(libs, paste(unique(libs$group[filterExtractedBool])))
    updateCheckboxGroupInput(session,
      "groupsCheck",
      choices = myGroups,
      selected = if(input$bar) myGroups
    )
  })
  
  # -> libsSamplename
  
  contentlibsSamplename <- eventReactive(input$updateSamples, {
    withProgress(message = 'updating samples', {
      incProgress(1/3, detail = "extracting selection")
      filterSelectionNames <- rownames(libs)[filterSelectionBool()]
      incProgress(2/3, detail = "creating checkbox")
      mySamples <- vizection:::addGroupName(libs, paste(filterSelectionNames))
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
      a <- subgenes() %>% corMat_1
      incProgress(2/4, detail = "log1p")
      b <- a %>% corMat_2
      incProgress(3/4, detail = "cor")
      b %>% corMat_3
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
  
  # PCoA AND KMEANS
  # ===============
  
  contentgenesPCoA <- eventReactive(input$updatePCoA,{
    withProgress(message = 'PCoA', {
      incProgress(1/3, detail = "collecting data")
      distCorMat <- distCorMat()
      incProgress(2/3, detail = "calculating")
      dudi.pco(distCorMat, scannf = F, nf = 2)
    })
  })
  genesPCoA <- reactive({
    contentgenesPCoA()
  })
  output$pcoasummary <- renderPrint({
    summary(genesPCoA())
  })
  #
  rangespcoa12 <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$pcoa12dblclick, {
    brush <- input$pcoa12brush
    if(!is.null(brush)) {
      rangespcoa12$x <- c(brush$xmin, brush$xmax)
      rangespcoa12$y <- c(brush$ymin, brush$ymax)
    } else {
      rangespcoa12$x <- NULL
      rangespcoa12$y <- NULL
    }
  })
  pcoa12 <- eventReactive(input$updatePCoA, {
    genesPCoAli <- genesPCoA()$li
    ggplot(genesPCoAli, aes(x = A1, y = A2))
  })
  contentpcoagenes12 <- reactive({
    withProgress(message = 'plot PCoA', {
      incProgress(1/4, detail = "collecting PCoA")
      pcoa12 <- pcoa12()
      incProgress(2/4, detail = "collecting k-means colors")
      kmeansColor <- kmeansColor()
      incProgress(3/4, detail = "creating plot")
      pcoa12 +
        geom_point(color = kmeansColor) +
        coord_cartesian(xlim = rangespcoa12$x, ylim = rangespcoa12$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light() +
        theme(legend.position = "none")
    })
  })
  output$pcoagenes12 <- renderPlot({
    contentpcoagenes12()
  })
  contentdataPCoA <- reactive({
    withProgress(message = 'data PCoA', {
      incProgress(1/4, detail = "collecting sulibs")
      sublibs <- sublibs()
      incProgress(2/4, detail = "filtering")
      res0 <- brushedPoints(genesPCoA()$li, input$pcoa12brush, xvar = "A1", yvar = "A2")
      colour <- kmeansColor()
      resCol <- cbind(colour, sublibs[, -1])
      res <- resCol[rownames(resCol) %in% rownames(res0), ]
      colour2 <- res$colour
      incProgress(3/4, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE)) %>% formatStyle(
        "colour", target = 'row', backgroundColor = styleEqual(colour2, colour2)
      )
    })
  })
  output$dataPCoA <- renderDataTable({
    contentdataPCoA()
  })
  #
  contentSSE <- eventReactive(input$updateSSE, {
    SSE <- function(mydata, title = ""){
      wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
      for (i in 2:input$SSElength) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
      plot(1:input$SSElength, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares", main = title)
      for(i in wss){abline(h = i, lty = 2, col = "grey")}
    }
    SSE(genesPCoA()$li[, c(1, 2)])
  })
  output$SSE <- renderPlot({
    contentSSE()
  })
  #
  contentpcoakmeans <- eventReactive(input$updatekmeans, {
    withProgress(message = 'kmeans PCoA', {
      incProgress(1/4, detail = "performing kmeans")
      kmeanspco <- kmeans(genesPCoA()$li[, c(1, 2)], input$kmeansClusters)
      incProgress(2/4, detail = "extracting clusters")
      kmeanspcofitted <- fitted(kmeanspco)
      incProgress(3/4, detail = "grouping information")
      data.frame("sampleName" = rownames(genesPCoA()$li), "cluster" = kmeanspcofitted %>% rownames() %>% as.factor(), "centroidX" = kmeanspcofitted[, 1], "centroidY" = kmeanspcofitted[, 2])
    })
  })
  pcoakmeans <- reactive({
    contentpcoakmeans()
  })
  contentkmeansColor <- eventReactive(input$updatekmeans, {
    withProgress(message = 'colors PCoA', {
      incProgress(1/4, detail = "attribution")
      pcoakmeans <- pcoakmeans()
      colorvector <- rainbow(input$kmeansClusters) %>% substr(., 1, nchar(.)-2)
      colorvector[pcoakmeans$cluster]
    })
  })
  kmeansColor <- reactive({
    contentkmeansColor()
  })
  
  # PCA
  # ===
  
  contentgenesPCA <- eventReactive(input$updatePCASummary, {
    withProgress(message = 'PCA summary', {
      incProgress(1/3, detail = "TPM")
      genesTpm <- subgenes() %>% TPM %>% t
      incProgress(2/3, detail = "dudi.pca")
      dudi.pca(genesTpm[, -1], center = T, scale = F, scannf = F, nf = 3)
    })
  })
  genesPca <- reactive({
    contentgenesPCA()
  })
  
  output$pcasummary <- renderPrint({
    summary(genesPca())
  })
  
  output$eigenvalues <- renderPlot({
    barplot(genesPca() %$% eig, xlab = "Eigenvalues")
  })
  
  #
  contentcomponents1 <- eventReactive(input$updatePCAComponents, {
    withProgress(message = 'components1', {
      incProgress(1/4, detail = "collecting PCA")
      genesPca <- genesPca()
      incProgress(2/4, detail = 'generating list')
      genesCoComp1 <- pcaCompGenesList(genesPca$co, 1)
      incProgress(3/4, detail = 'generating plot')
      plotHTB(genesCoComp1, 1, input$nbDispGenes)
    })
  })
  output$components1 <- renderPlot({
    contentcomponents1()
  })
  
  contentcomponents2 <- eventReactive(input$updatePCAComponents, {
    withProgress(message = 'components2', {
      incProgress(1/4, detail = "collecting PCA")
      genesPca <- genesPca()
      incProgress(2/4, detail = 'generating list')
      genesCoComp2 <- pcaCompGenesList(genesPca$co, 2)
      incProgress(3/4, detail = 'generating plot')
      plotHTB(genesCoComp2, 2, input$nbDispGenes)
    })
  })
  output$components2 <- renderPlot({
    contentcomponents2()
  })
  
  contentcomponents3 <- eventReactive(input$updatePCAComponents, {
    withProgress(message = 'components3', {
      incProgress(1/4, detail = "collecting PCA")
      genesPca <- genesPca()
      incProgress(2/4, detail = 'generating list')
      genesCoComp3 <- pcaCompGenesList(genesPca$co, 3)
      incProgress(3/4, detail = 'generating plot')
      plotHTB(genesCoComp3, 3, input$nbDispGenes)
    })
  })
  output$components3 <- renderPlot({
    contentcomponents3()
  })
  
  pcaColor <- reactive({
    if(input$PCAcolor == 2){
      paste(colorsPcaLi())
    }
    else if(input$PCAcolor == 3){
      kmeansColor()
    }
    else{
      myColors <- sublibs()$group %>% levels %>% length %>% rainbow() %>% substr(., 1, nchar(.)-2)
      myColors[sublibs()$group]
    }
  })
  
  pcaGroup <- reactive({
    if(input$PCAcolor == 2){
      as.factor(colorsPcaLi())
    }
    else if(input$PCAcolor == 3){
      as.factor(kmeansColor())
    }
    else{
      sublibs()$group
    }
  })
  
  # ax 12
  ####
  ranges12li <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$PCA12lidblclick, {
    brush <- input$PCA12librush
    if (!is.null(brush)) {
      ranges12li$x <- c(brush$xmin, brush$xmax)
      ranges12li$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges12li$x <- NULL
      ranges12li$y <- NULL
    }
  })
  ####
  g12li <- eventReactive(input$updatePCAPlots, {
    genesPcali <- genesPca()$li
    if(input$showEllipse){
      ggplot(genesPcali, aes(x = Axis1, y = Axis2, group = pcaGroup(), color = pcaColor(), fill = pcaColor())) + stat_ellipse(aes(color = pcaColor(), fill = pcaColor()))}
    else {
      ggplot(genesPcali, aes(x = Axis1, y = Axis2, group = pcaGroup(), color = pcaColor()))
    }
  })
  contentinteractPCA12li <- reactive({
    withProgress(message = 'li axes 1-2', {
      incProgress(1/4, detail = "collecting PCA")
      g12li <- g12li()
      incProgress(2/4, detail = "collecting colors")
      pcaColor <- pcaColor()
      incProgress(3/4, detail = "creating plot")
      g12li +
        geom_point(color = pcaColor) +
        coord_cartesian(xlim = ranges12li$x, ylim = ranges12li$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light() +
        theme(legend.position = "none")
    })
  })
  output$interactPCA12li <- renderPlot({
    contentinteractPCA12li()
  })
  ####
  contentdataPCA12li <- reactive({
    withProgress(message = 'data 1-2', {
      incProgress(1/4, detail = "collecting sulibs")
      sublibs <- sublibs()
      incProgress(2/4, detail = "filtering")
      res0 <- brushedPoints(genesPca()$li, input$PCA12librush, xvar = "Axis1", yvar = "Axis2")
      colour <- pcaColor()
      resCol <- cbind(colour, sublibs[, -1])
      res <- resCol[rownames(resCol) %in% rownames(res0), ]
      colour2 <- res$colour
      incProgress(3/4, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE)) %>% formatStyle(
        "colour", target = 'row', backgroundColor = styleEqual(colour2, colour2)
      )
    })
  })
  output$dataPCA12li <- renderDataTable({
    contentdataPCA12li()
  })
  #####
  ranges12co <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$PCA12codblclick, {
    brush <- input$PCA12cobrush
    if (!is.null(brush)) {
      ranges12co$x <- c(brush$xmin, brush$xmax)
      ranges12co$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges12co$x <- NULL
      ranges12co$y <- NULL
    }
  })
  #####
  g12co <- eventReactive(input$updatePCAPlots, {
    genesPcaco <- genesPca()$co
    ggplot(genesPcaco, aes(x = Comp1, y = Comp2))
  })
  contentinteractPCA12co <- reactive({
    withProgress(message = 'co axes 1-2', {
      incProgress(1/3, detail = "collecting PCA")
      g12co <- g12co()
      incProgress(2/3, detail = "creating plot")
      g12co +
        geom_segment(aes(x=0, y=0, xend=Comp1, yend=Comp2)) +
        coord_cartesian(xlim = ranges12co$x, ylim = ranges12co$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light()
    })
  })
  output$interactPCA12co <- renderPlot({
    contentinteractPCA12co()
  })
  #####
  contentdataPCA12co <- reactive({
    withProgress(message = 'data 1-2', {
      incProgress(1/3, detail = "filtering")
      res <- brushedPoints(genesPca()$co, input$PCA12cobrush, xvar = "Comp1", yvar = "Comp2")
      incProgress(2/3, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE))
    })
  })
  output$dataPCA12co <- renderDataTable({
    contentdataPCA12co()
  })
  
  # ax 13
  ####
  ranges13li <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$PCA13lidblclick, {
    brush <- input$PCA13librush
    if (!is.null(brush)) {
      ranges13li$x <- c(brush$xmin, brush$xmax)
      ranges13li$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges13li$x <- NULL
      ranges13li$y <- NULL
    }
  })
  ####
  g13li <- eventReactive(input$updatePCAPlots, {
    genesPcali <- genesPca()$li
    ggplot(genesPcali, aes(x = Axis1, y = Axis3))
  })
  contentinteractPCA13li <- reactive({
    withProgress(message = 'li axes 1-3', {
      incProgress(1/4, detail = "collecting PCA")
      g13li <- g13li()
      incProgress(2/4, detail = "collecting colors")
      pcaColor <- pcaColor()
      incProgress(3/4, detail = "creating plot")
      g13li +
        geom_point(color = pcaColor) +
        coord_cartesian(xlim = ranges13li$x, ylim = ranges13li$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light()
    })
  })
  output$interactPCA13li <- renderPlot({
    contentinteractPCA13li()
  })
  ####
  contentdataPCA13li <- reactive({
    withProgress(message = 'data 1-3', {
      incProgress(1/4, detail = "collecting sulibs")
      sublibs <- sublibs()
      incProgress(2/4, detail = "filtering")
      res0 <- brushedPoints(genesPca()$li, input$PCA13librush, xvar = "Axis1", yvar = "Axis3")
      colour <- pcaColor()
      resCol <- cbind(colour, sublibs[, -1])
      res <- resCol[rownames(resCol) %in% rownames(res0), ]
      colour2 <- res$colour
      incProgress(3/4, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE)) %>% formatStyle(
        "colour", target = 'row', backgroundColor = styleEqual(colour2, colour2)
      )
    })
  })
  output$dataPCA13li <- renderDataTable({
    contentdataPCA13li()
  })
  #####
  ranges13co <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$PCA13codblclick, {
    brush <- input$PCA13cobrush
    if (!is.null(brush)) {
      ranges13co$x <- c(brush$xmin, brush$xmax)
      ranges13co$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges13co$x <- NULL
      ranges13co$y <- NULL
    }
  })
  #####
  g13co <- eventReactive(input$updatePCAPlots, {
    genesPcaco <- genesPca()$co
    ggplot(genesPcaco, aes(x = Comp1, y = Comp3))
  })
  contentinteractPCA13co <- reactive({
    withProgress(message = 'co axes 1-3', {
      incProgress(1/3, detail = "collecting PCA")
      g13co <- g13co()
      incProgress(2/3, detail = "creating plot")
      g13co +
        geom_segment(aes(x=0, y=0, xend=Comp1, yend=Comp3)) +
        coord_cartesian(xlim = ranges13co$x, ylim = ranges13co$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light()
    })
  })
  output$interactPCA13co <- renderPlot({
    contentinteractPCA13co()
  })
  #####
  contentdataPCA13co <- reactive({
    withProgress(message = 'data 1-3', {
      incProgress(1/3, detail = "filtering")
      res <- brushedPoints(genesPca()$co, input$PCA13cobrush, xvar = "Comp1", yvar = "Comp3")
      incProgress(2/3, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE))
    })
  })
  output$dataPCA13co <- renderDataTable({
    contentdataPCA13co()
  })
  
  # ax 32
  ####
  ranges32li <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$PCA32lidblclick, {
    brush <- input$PCA32librush
    if (!is.null(brush)) {
      ranges32li$x <- c(brush$xmin, brush$xmax)
      ranges32li$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges32li$x <- NULL
      ranges32li$y <- NULL
    }
  })
  ####
  g32li <- eventReactive(input$updatePCAPlots, {
    genesPcali <- genesPca()$li
    ggplot(genesPcali, aes(x = Axis3, y = Axis2))
  })
  contentinteractPCA32li <- reactive({
    withProgress(message = 'li axes 3-2', {
      incProgress(1/4, detail = "collecting PCA")
      g32li <- g32li()
      incProgress(2/4, detail = "collecting colors")
      pcaColor <- pcaColor()
      incProgress(3/4, detail = "creating plot")
      g32li +
        geom_point(color = pcaColor) +
        coord_cartesian(xlim = ranges32li$x, ylim = ranges32li$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light()
    })
  })
  output$interactPCA32li <- renderPlot({
    contentinteractPCA32li()
  })
  ####
  contentdataPCA32li <- reactive({
    withProgress(message = 'data 3-2', {
      incProgress(1/4, detail = "collecting sulibs")
      sublibs <- sublibs()
      incProgress(2/4, detail = "filtering")
      res0 <- brushedPoints(genesPca()$li, input$PCA32librush, xvar = "Axis3", yvar = "Axis2")
      colour <- pcaColor()
      resCol <- cbind(colour, sublibs[, -1])
      res <- resCol[rownames(resCol) %in% rownames(res0), ]
      colour2 <- res$colour
      incProgress(3/4, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE)) %>% formatStyle(
        "colour", target = 'row', backgroundColor = styleEqual(colour2, colour2)
      )
    })
  })
  output$dataPCA32li <- renderDataTable({
    contentdataPCA32li()
  })
  #####
  ranges32co <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$PCA32codblclick, {
    brush <- input$PCA32cobrush
    if (!is.null(brush)) {
      ranges32co$x <- c(brush$xmin, brush$xmax)
      ranges32co$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges32co$x <- NULL
      ranges32co$y <- NULL
    }
  })
  #####
  g32co <- eventReactive(input$updatePCAPlots, {
    genesPcaco <- genesPca()$co
    ggplot(genesPcaco, aes(x = Comp3, y = Comp2))
  })
  contentinteractPCA32co <- reactive({
    withProgress(message = 'co axes 3-2', {
      incProgress(1/3, detail = "collecting PCA")
      g32co <- g32co()
      incProgress(2/3, detail = "creating plot")
      g32co +
        geom_segment(aes(x=0, y=0, xend=Comp3, yend=Comp2)) +
        coord_cartesian(xlim = ranges32co$x, ylim = ranges32co$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light()
    })
  })
  output$interactPCA32co <- renderPlot({
    contentinteractPCA32co()
  })
  #####
  contentdataPCA32co <- reactive({
    withProgress(message = 'data 3-2', {
      incProgress(1/3, detail = "filtering")
      res <- brushedPoints(genesPca()$co, input$PCA32cobrush, xvar = "Comp3", yvar = "Comp2")
      incProgress(2/3, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE))
    })
  })
  output$dataPCA32co <- renderDataTable({
    contentdataPCA32co()
  })
  
  contentpca3D <- eventReactive(input$generatepca3d, {
    withProgress(message = 'pca 3D', {
      incProgress(1/4, detail = "collecting PCA")
      pcaGenesli <- genesPca()$li
      incProgress(2/4, detail = "collecting colors")
      pcaColor <- pcaColor()
      pcaGroup <- pcaGroup()
      incProgress(3/4, detail = "creating 3D plot")
      plot_ly(data = pcaGenesli, x = pcaGenesli$Axis1, y = pcaGenesli$Axis2, z = pcaGenesli$Axis3,
              type = "scatter3d", mode = "markers", marker = list(size = input$pca3ddotsize),
              color = pcaGroup, colors = pcaColor,
              text = sublibs()$samplename) %>% 
        layout(scene = list(
                 xaxis = list(title = "Axis1"), 
                 yaxis = list(title = "Axis2"), 
                 zaxis = list(title = "Axis3")))
    })
  })
  output$pca3D <- plotly::renderPlotly({
    contentpca3D()
  })
  
  #
  contentcheckplot <- eventReactive(input$updatecheckplot, {
    withProgress(message = 'checkplot', {
      incProgress(1/4, detail = "data collection")
      if(input$dataCheckplot == "total"){
        d <- genes[paste(input$geneNamescheckplot), ] %>% t %>% tbl_df() %T>% setnames("geneName")
        d$group <-as.factor(libs$group)
      }
      else {
        subgenes <- subgenes()
        sublibs <- sublibs()
        d <- subgenes[paste(input$geneNamescheckplot), ] %>% t %>% tbl_df() %T>% setnames("geneName")
        d$group <- as.factor(sublibs$group)
      }
      incProgress(2/4, detail = "generation")
      g <- d %>% ggplot(aes(geneName), group = group) +
        geom_histogram(binwidth = 1) +
        theme_light() +
        xlab(paste(input$geneNamescheckplot))
      incProgress(3/4, detail = "(faceting) and annotating")
      if(input$facetcheckplot){
        g <- g + facet_grid(group ~ .)
      }
      if(input$sampleNamescheckplot != "None"){
        if(input$dataCheckplot == "total"){
          xvalue <- genes[paste(input$geneNamescheckplot), paste(input$sampleNamescheckplot)]
          yvalue <- sum(genes[paste(input$geneNamescheckplot), ] == xvalue)
          g <- g + geom_vline(xintercept = xvalue, colour = "red", linetype = "dashed") +
            annotate("text", x = xvalue, y = yvalue + 1, label = paste(input$sampleNamescheckplot), colour = "red")
        } else {
          subgenes <- subgenes()
          xvalue <- subgenes[paste(input$geneNamescheckplot), paste(input$sampleNamescheckplot)]
          yvalue <- sum(subgenes[paste(input$geneNamescheckplot), ] == xvalue)
          g <- g + geom_vline(xintercept = xvalue, colour = "red", linetype = "dashed") +
            annotate("text", x = xvalue, y = yvalue + 1, label = paste(input$sampleNamescheckplot), colour = "red")
        }
      }
      g
    })
  })
  output$checkplot <- renderPlot({
    contentcheckplot()
  })
  
  contentgeneNamescheckplotUI <- eventReactive(input$updatelistcheckplot, {
    withProgress(message = 'checkplot genes', {
      incProgress(1/3, detail = "searching...")
      if(input$dataCheckplot == "total"){
        checkplotGrep <- rownames(genes) %>% grep(input$geneNameCheckplot, .) %>% rownames(genes)[.]
      } else {
        checkplotGrep <- rownames(subgenes()) %>% grep(input$geneNameCheckplot, .) %>% rownames(subgenes())[.]
      }
      incProgress(2/3, detail = "creating UI")
      selectInput("geneNamescheckplot", "Gene name:", c("None", checkplotGrep), selected = "None")
    })
  })
  output$geneNamescheckplotUI <- renderUI({
    contentgeneNamescheckplotUI()
  })
  
  contentsampleNamescheckplotUI <- eventReactive(input$updatelistcheckplot, {
    withProgress(message = 'checkplot samples', {
      incProgress(1/2, detail = "data collection")
      if(input$dataCheckplot == "total"){
        selectInput("sampleNamescheckplot", "Sample name:",
                    choices = c("None", paste(rownames(libs))), selected = "None"
        )
      } else {
        selectInput("sampleNamescheckplot", "Sample name:",
                    choices = c("None", paste(rownames(sublibs()))), selected = "None"
        )
      }
    })
  })
  output$sampleNamescheckplotUI <- renderUI({
    contentsampleNamescheckplotUI()
  })
  
  output$checkplotUI <- renderUI({
    plotOutput("checkplot", height = input$heightcheckplot)
  })
  
  # CA
  # ==
  
  contentcontribDataFrame <- reactive({
    withProgress(message = 'calculating contribution', {
      incProgress(1/5, detail = "collecting PCA")
      genesPca <- genesPca()
      selectedAxis <- as.numeric(input$selectAxisCoA)
      incProgress(2/5, detail = "calculating")
      contribution <- abs(genesPca$co[,selectedAxis])/sum(abs(genesPca$co[,selectedAxis])) * 100
      incProgress(3/5, detail = "checking results")
      stopifnot(all.equal(sum(contribution), 100))
      incProgress(4/5, detail = "creating data frame")
      data.frame("geneName" = rownames(genesPca$co), "contribution" = contribution)
    })
  })
  contribDataFrame <- reactive({
    contentcontribDataFrame()
  })
  
  contentcontributionBoxplot <- eventReactive(input$generateContributionBoxplot, {
    withProgress(message = 'contrib boxplot', {
      incProgress(1/3, detail = "collecting contrib")
      contribDataFrame <- contribDataFrame()
      incProgress(2/3, detail = "generating boxplot")
      boxplot(contribDataFrame$contribution, horizontal = T, main = paste0("Contribution on axis ", input$selectAxisCoA))
    })
  })
  output$contributionBoxplot <- renderPlot({
    contentcontributionBoxplot()
  })
  
  contentthresholded <- eventReactive(input$applyThreshold, {
    withProgress(message = 'applying threshold', {
      incProgress(1/4, detail = "collecting threshold")
      threshold <- input$nbGenesToKeep
      contribDataFrame <- contribDataFrame()
      incProgress(2/4, detail = "filtering contrib data frame")
      indexesThresholded <- which(contribDataFrame$contribution >= threshold)
      incProgress(3/4, detail = "creating list")
      list(
        names = contribDataFrame$geneName[indexesThresholded],
        values = contribDataFrame$contribution[indexesThresholded]
      )
    })
  })
  thresholded <- reactive({
    contentthresholded()
  })
  
  output$thresholdedPrint <- renderDataTable({
    as.data.frame(thresholded())
  })
  
  contentnumberGenesThresholded <- eventReactive(input$applyThreshold, {
    thresholded <- thresholded()
    paste("Selection of", length(thresholded$names), "genes.")
  })
  output$numberGenesThresholded <- renderPrint({
    contentnumberGenesThresholded()
  })
  contentthresholdBoxplot <- eventReactive(input$applyThreshold, {
    contribDataFrame <- contribDataFrame()
    boxplot(contribDataFrame$contribution, horizontal = T, main = paste0("Contribution on axis ", input$selectAxisCoA))
    abline(lty = 2, col = "red", v = input$nbGenesToKeep)
  })
  output$thresholdBoxplot <- renderPlot({
    contentthresholdBoxplot()
  })
  
  #
  contentcoaGenes <- eventReactive(input$updateCoA, {
    withProgress(message = 'CoA summary', {
      incProgress(1/3, detail = "creating thresholded genes and libs")
      subgenes <- subgenes()
      thresholded <- thresholded()
      mainGenes <- subgenes[rownames(subgenes) %in% thresholded$names, ]
      incProgress(2/3, detail = "dudi.coa")
      dudi.coa(mainGenes %>% t %>% as.data.frame, scannf = F, nf = 3)
    })
  })
  coaGenes <- reactive({
    contentcoaGenes()
  })
  
  output$coasummary <- renderPrint({
    summary(coaGenes())
  })
  
  output$coaeigenvalues <- renderPlot({
    barplot(coaGenes() %$% eig, xlab = "Eigenvalues")
  })
  
  
  coaColor <- reactive({
    if(input$COAcolor == 2){
      paste(colorsPcaLi())
    }
    else if(input$COAcolor == 3){
      kmeansColor()
    }
    else{
      myColors <- sublibs()$group %>% levels %>% length %>% rainbow() %>% substr(., 1, nchar(.)-2)
      myColors[sublibs()$group]
    }
  })
  
  coaGroup <- reactive({
    if(input$COAcolor == 2){
      as.factor(colorsPcaLi())
    }
    else if(input$COAcolor == 3){
      as.factor(kmeansColor())
    }
    else{
      sublibs()$group
    }
  })
  
  ####
  rangescoa12 <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$COA12dblclick, {
    brush <- input$COA12brush
    if (!is.null(brush)) {
      rangescoa12$x <- c(brush$xmin, brush$xmax)
      rangescoa12$y <- c(brush$ymin, brush$ymax)
    } else {
      rangescoa12$x <- NULL
      rangescoa12$y <- NULL
    }
  })
  ####
  coa12 <- eventReactive(input$updateCoAPlots, {
    coaGenesli <- coaGenes()$li
    ggplot(coaGenesli, aes(x = Axis1, y = Axis2))
  })
  contentinteractCOA12 <- reactive({
    withProgress(message = 'coa axes 1-2', {
      incProgress(1/4, detail = "collecting COA")
      coa12 <- coa12()
      incProgress(2/4, detail = "collecting colors")
      coaColor <- coaColor()
      incProgress(3/4, detail = "creating plot")
      coa12 +
        geom_point(color = coaColor) +
        coord_cartesian(xlim = rangescoa12$x, ylim = rangescoa12$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light() +
        geom_point(data = coaGenes()$co, aes(x = Comp1, y = Comp2)) +
        geom_text(data = coaGenes()$co, aes(x = Comp1, y = Comp2, label = rownames(coaGenes()$co)), hjust = 0, nudge_x = 0.05)
    })
  })
  output$interactCOA12 <- renderPlot({
    contentinteractCOA12()
  })
  ####
  contentdataCOA12 <- reactive({
    withProgress(message = 'data 1-2', {
      incProgress(1/4, detail = "collecting sulibs")
      sublibs <- sublibs()
      incProgress(2/4, detail = "filtering")
      res0 <- brushedPoints(coaGenes()$li, input$COA12brush, xvar = "Axis1", yvar = "Axis2")
      colour <- coaColor()
      resCol <- cbind(colour, sublibs[, -1])
      res <- resCol[rownames(resCol) %in% rownames(res0), ]
      colour2 <- res$colour
      incProgress(3/4, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE)) %>% formatStyle(
        "colour", target = 'row', backgroundColor = styleEqual(colour2, colour2)
      )
    })
  })
  output$dataCOA12 <- renderDataTable({
    contentdataCOA12()
  })
  
  ####
  rangescoa13 <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$COA13dblclick, {
    brush <- input$COA13brush
    if (!is.null(brush)) {
      rangescoa13$x <- c(brush$xmin, brush$xmax)
      rangescoa13$y <- c(brush$ymin, brush$ymax)
    } else {
      rangescoa13$x <- NULL
      rangescoa13$y <- NULL
    }
  })
  ####
  coa13 <- eventReactive(input$updateCoAPlots, {
    coaGenesli <- coaGenes()$li
    ggplot(coaGenesli, aes(x = Axis1, y = Axis3))
  })
  contentinteractCOA13 <- reactive({
    withProgress(message = 'coa axes 1-2=3', {
      incProgress(1/4, detail = "collecting COA")
      coa13 <- coa13()
      incProgress(2/4, detail = "collecting colors")
      coaColor <- coaColor()
      incProgress(3/4, detail = "creating plot")
      coa13 +
        geom_point(color = coaColor) +
        coord_cartesian(xlim = rangescoa13$x, ylim = rangescoa13$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light() +
        geom_point(data = coaGenes()$co, aes(x = Comp1, y = Comp3)) +
        geom_text(data = coaGenes()$co, aes(x = Comp1, y = Comp3, label = rownames(coaGenes()$co)), hjust = 0, nudge_x = 0.05)
    })
  })
  output$interactCOA13 <- renderPlot({
    contentinteractCOA13()
  })
  ####
  contentdataCOA13 <- reactive({
    withProgress(message = 'data 1-3', {
      incProgress(1/4, detail = "collecting sulibs")
      sublibs <- sublibs()
      incProgress(2/4, detail = "filtering")
      res0 <- brushedPoints(coaGenes()$li, input$COA13brush, xvar = "Axis1", yvar = "Axis3")
      colour <- coaColor()
      resCol <- cbind(colour, sublibs[, -1])
      res <- resCol[rownames(resCol) %in% rownames(res0), ]
      colour2 <- res$colour
      incProgress(3/4, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE)) %>% formatStyle(
        "colour", target = 'row', backgroundColor = styleEqual(colour2, colour2)
      )
    })
  })
  output$dataCOA13 <- renderDataTable({
    contentdataCOA13()
  })
  
  ####
  rangescoa32 <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$COA32dblclick, {
    brush <- input$COA32brush
    if (!is.null(brush)) {
      rangescoa32$x <- c(brush$xmin, brush$xmax)
      rangescoa32$y <- c(brush$ymin, brush$ymax)
    } else {
      rangescoa32$x <- NULL
      rangescoa32$y <- NULL
    }
  })
  ####
  coa32 <- eventReactive(input$updateCoAPlots, {
    coaGenesli <- coaGenes()$li
    ggplot(coaGenesli, aes(x = Axis3, y = Axis2))
  })
  contentinteractCOA32 <- reactive({
    withProgress(message = 'coa axes 3-2', {
      incProgress(1/4, detail = "collecting COA")
      coa32 <- coa32()
      incProgress(2/4, detail = "collecting colors")
      coaColor <- coaColor()
      incProgress(3/4, detail = "creating plot")
      coa32 +
        geom_point(color = coaColor) +
        coord_cartesian(xlim = rangescoa32$x, ylim = rangescoa32$y) +
        geom_vline(xintercept = 0, alpha = 0.2) +
        geom_hline(yintercept = 0, alpha = 0.2) +
        theme_light() +
        geom_point(data = coaGenes()$co, aes(x = Comp3, y = Comp2)) +
        geom_text(data = coaGenes()$co, aes(x = Comp3, y = Comp2, label = rownames(coaGenes()$co)), hjust = 0, nudge_x = 0.05)
    })
  })
  output$interactCOA32 <- renderPlot({
    contentinteractCOA32()
  })
  ####
  contentdataCOA32 <- reactive({
    withProgress(message = 'data 3-2', {
      incProgress(1/4, detail = "collecting sulibs")
      sublibs <- sublibs()
      incProgress(2/4, detail = "filtering")
      res0 <- brushedPoints(coaGenes()$li, input$COA32brush, xvar = "Axis3", yvar = "Axis2")
      colour <- coaColor()
      resCol <- cbind(colour, sublibs[, -1])
      res <- resCol[rownames(resCol) %in% rownames(res0), ]
      colour2 <- res$colour
      incProgress(3/4, detail = "creating datatable")
      datatable(res, options = list(scrollX = TRUE)) %>% formatStyle(
        "colour", target = 'row', backgroundColor = styleEqual(colour2, colour2)
      )
    })
  })
  output$dataCOA32 <- renderDataTable({
    contentdataCOA32()
  })
  
  contentcoa3D <- eventReactive(input$generatecoa3d, {
    withProgress(message = 'coa 3D', {
      incProgress(1/4, detail = "collecting")
      coaGenesli <- coaGenes()$li
      incProgress(2/4, detail = "collecting colors")
      coaColor <- coaColor()
      coaGroup <- coaGroup()
      incProgress(3/4, detail = "creating 3D plot")
      plot_ly(data = coaGenesli, x = coaGenesli$Axis1, y = coaGenesli$Axis2, z = coaGenesli$Axis3,
              type = "scatter3d", mode = "markers", marker = list(size = input$coa3ddotsize),
              color = coaGroup, colors = coaColor,
              text = sublibs()$samplename)  %>% 
        layout(scene = list(
          xaxis = list(title = "Axis1"), 
          yaxis = list(title = "Axis2"), 
          zaxis = list(title = "Axis3")))
    })
  })
  output$coa3D <- plotly::renderPlotly({
    contentcoa3D()
  })
  
  # EXPORT
  # ======
  
  observeEvent(input$exportGenes, {
    withProgress(message = "Exporting genes", {
      incProgress(1/2, detail = "processing")
      saveRDS(subgenes(), file = file.path(getwd(), paste(input$genesRDSName)))
    })
  })
  
  observeEvent(input$exportLibs, {
    withProgress(message = "Exporting libs", {
      incProgress(1/2, detail = "processing")
      saveRDS(sublibs(), file = paste(input$libsRDSName))
    })
  })
  
})
