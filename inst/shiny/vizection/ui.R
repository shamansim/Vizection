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

genes <- get(getOption("vizection.genes"), .GlobalEnv)
libs  <- get(getOption("vizection.libs"),  .GlobalEnv)

# BEGIN shiny app
dashboardPage(
  
  # HEADER
  # ------
  dashboardHeader(title = tags$img(src = "vizection_white.png", alt = "Vizection", height = "40"),
    dropdownMenuOutput("tasksMenu")
  ),#dashboardHeader
  
  # SIDEBAR
  # -------
  dashboardSidebar(#width = 450,
    
    sidebarMenu(id = "sidebarmenu",
      
      menuItem("Selection", tabName = "selection", icon = icon("check-square-o")),

      menuItem("General", tabName = "home", icon = icon("home")),

      menuItem("Clustering", tabName = "clustering", icon = icon("sort-amount-desc"),
        menuSubItem("Hierarchical clustering", tabName = "hierarchicalClustering", icon = icon("sitemap")),
        menuSubItem("K-means", tabName = "k-means", icon = icon("dot-circle-o"))
      ),#clustering

      menuItem("PCA", tabName = "PCA", icon = icon("area-chart")
      ),#PCA

      menuItem("CA", tabName = "CA", icon = icon("area-chart")
      ),#CA

      menuItem("Export", tabName = "export", icon = icon("share-square-o")
      ),#Export
      
      menuItem("About", tabName = "about", icon = icon("question-circle")
      )#About

    )#sidebarMenu
    
  ),#dashboardSidebar
  
  # BODY
  # ----
  dashboardBody(
    
    tags$head(
      tags$style(HTML('
        /* 
        body 
        */
        
        body {
        height: auto;
        overflow: auto;
        font-family: sans-serif;
        }
        
        /* 
        header
        */
        
        .skin-blue .main-header .logo {
        background-color: black;
        }
        
        .skin-blue .main-header .navbar {
        background-color: black;
        }
        
        /*
        sidebar
        */
        
        .skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper {
        background-color: black;
        }
        
        .skin-blue .sidebar-menu > li.active > a, .skin-blue .sidebar-menu > li:hover > a {
        background: black;
        }
        
        /*
        content
        */
        
        .content {
        background-color: white;
        }
        
        .content-wrapper, .right-side {
        background-color: white;
        }
        
        img {
        display: block;
        margin: auto;
        }
        
        .fa-5x {
        text-align: center;
        width: 100%;
        }
        
        .titre {
        text-align: center;
        width: 100%;
        font-family: sans-serif;
        font-variant: small-caps;
        }
        
        /* tabs */
        .nav > li > a:active, .nav > li > a:focus, .nav > li > a:hover {
        background: #ff9d00;
        }
        
        /* boxplotGeneral */
        #boxplotGeneral {
        text-align: center;
        }

        /* label above input (e.g. dendrogram nb clusters) */
        .centeredInput {
        text-align: center;
        display: block;
        margin: auto;
        width: 10%;
        }

        /* solid boxes primary and info */
        .box.box-solid.box-primary > .box-header {
        background-color: #111987;
        }

        .box.box-solid.box-info > .box-header {
        background-color: #078e53;
        }

      '
      ))
    ),#tags$head
    
    # SELECTION
    # =========
    conditionalPanel(condition = "input.sidebarmenu == 'selection'",
      tags$i(class = "fa fa-check-square-o fa-5x"),
      tags$h4(class = "titre", "Selection"),
      helpText("Please choose a filter on the 'count' column [libs table] (0 for none) then select the groups. You can then update the samples and make your selection. Do not forget to update the selection and validate."),
      tags$hr(),
      tags$br(),
      numericInput("nbFilterExtracted", "Filter (counts)", min = 0, value = 0, step = 1),
      tags$br(),
      fluidRow(
        box(title = "Groups", solidHeader = T, collapsible = T, status = "primary",
          checkboxInput('bar', 'All groups/None'),
          uiOutput("libsGroup")
          
        ),#box Groups
        box(title = "Samples", solidHeader = T, collapsible = T, status = "info",
          actionButton("updateSamples", "Update samples"),
          uiOutput("libsSamplename")
        )
      ),#fluidRow
      fluidRow(
        column(1, actionButton("updateSelection", "Update selection")),
        column(1, actionButton("updateCorMat", "Validate selection"))
      )#fluidRow
    ),#conditionalPanel displayGroup
    
    tabItems(

      # HOME
      # ====
      tabItem(tabName = "home",
        tags$i(class = "fa fa-home fa-5x"),
        tags$h4(class = "titre", "General information"),
        helpText("Need a description..."),
        tags$hr(),
        tags$div(id = "boxplotGeneral",
          fluidRow(
            column(2, selectInput("boxplotGroupsTotal", "Groups (total)", c("none", paste(unique(libs$group))), selected = "none")),
            column(1, plotOutput("boxplotTotal"), height = "50px"),
            column(1, plotOutput("boxplotSub"), height = "50px"),
            column(2, uiOutput("UIboxplotGroupsSub"))
          )#fluidRow
        )#tags$div boxplot
      ),#tabItem home

      # HIERARCHICAL CLUSTERING
      # =======================
      tabItem(tabName = "hierarchicalClustering",
        # tags$i(class = "fa fa-sort-amount-desc fa-5x"),
        # tags$h4(class = "titre", "Hierarchical clustering"),
        # tags$hr(),

        # DENDROGRAM
        tabsetPanel(
          tabPanel("Dendrogram",
            tags$i(class = "fa fa-sitemap fa-5x"),
            tags$h4(class = "titre", "Dendrogram"),
            helpText("Need a description..."),
            tags$hr(),
            tags$div(class = "centeredInput", 
              numericInput("nbClusters", "Number of clusters", min = 1, max = nrow(libs), value = 1, step = 1)
            ),#tags$div centeredInput
            tags$br(),
            # fluidRow(
              box(title = "Dendrogram", solidHeader = T, collapsible = T, status = "primary",
                fluidRow(
                  column(3, sliderInput("dendroSize", "size", min = 0.1, max = 3, step = .1, value = 1)),
                  column(3, numericInput("heightDendro", "Plot height:", min = 20, max = 5000, value = 250)),
                  column(3, checkboxInput("dendroHoriz", "horizontal", FALSE)),
                  column(5, checkboxInput("showGroupsColor", "Color labels by group", FALSE))
                ),#fluidRow
                actionButton("updateDendrogram", "Update dendrogram"),
                tags$br(),
                uiOutput("dendrogramPlot")
              ),#box Dendrogram
              box(title = "Nodes height", solidHeader = T, collapsible = T, status = "warning",
                sliderInput("heightlength", "Number of clusters to display", min = 2, max = 100, value = 10, step = 1),
                actionButton("updateheight", "Update plot"),
                plotOutput("height")
              )#box Nodes height
            # )#fluidRow
          ),#tabPanel Dendrogram
          
          # HEATMAP
          tabPanel("Heatmap",
            tags$i(class = "fa fa-th-large fa-5x"),
            tags$h4(class = "titre", "Heatmap"),
            helpText("Need a description..."),
            tags$hr(),
            box(title = "Heatmap", solidHeader = T, collapsible = T, status = "primary", width = "950px", height = "950px",
                fluidRow(column(1, actionButton("updateHeatmap", "Update heatmap"))),
                tags$br(),
                plotOutput("heatmapGenes", width = "950px", height = "750px")
            )#box Heatmap
          )#tabPanel Heatmap
        )#tabsetPanel

      ),#tabItem hierarchicalClustering

      # K-MEANS
      # =======
      tabItem(tabName = "k-means",
        tags$i(class = "fa fa-dot-circle-o fa-5x"),
        tags$h4(class = "titre", "K-means"),
        helpText("Need a description..."),
        tags$hr(),
        helpText("Do not forget to 'validate parameters'."),
        tabsetPanel(
          tabPanel("PCoA and k-means clustering",
                   helpText("1 - Perform PCoA first"),
                   helpText("2 - Check for optimal number of clusters"),
                   helpText("3 - Perform k-means"),
                   fluidRow(column(1, actionButton("updatePCoA", "Update PCoA"))),
                   fluidRow(
                     column(4, numericInput("kmeansClusters", h4("Number of clusters", style = "color: rgb(0, 123, 16)"), min = 1, max = nrow(libs), value = 1, step = 1)
                     )
                   ),
                   fluidRow(column(1, actionButton("updatekmeans", "Update k-means"))),
                   plotOutput("pcoagenes12", brush = brushOpts(id = "pcoa12brush", resetOnNew = T), width = 470, dblclick = "pcoa12dblclick"),
                   dataTableOutput("dataPCoA"),
                   verbatimTextOutput("pcoasummary")
          ),
          tabPanel("SSE",
                   helpText("Within sum of squared error plot."),
                   helpText("Needs to perform PCoA first."),
                   sliderInput("SSElength", "Number of clusters to display", min = 2, max = 100, value = 10, step = 1),
                   fluidRow(column(1, actionButton("updateSSE", "Update SSE plot"))),
                   plotOutput("SSE")
          )
        )
      ),#tabItem k-means

      # PCA
      # ===
      tabItem(tabName = "PCA",
        tags$i(class = "fa fa-area-chart fa-5x"),
        tags$h4(class = "titre", "Principal Components Analysis"),
        helpText("Need a description..."),
        tags$hr(),
        helpText("Do not forget to 'validate parameters'."),
        tabsetPanel(
          tabPanel("Summary",
                   fluidRow(column(1, actionButton("updatePCASummary", "Update PCA"))),
                   plotOutput("eigenvalues"),
                   verbatimTextOutput("pcasummary")
          ),
          tabPanel("Components",
                   fluidRow(column(1, actionButton("updatePCAComponents", "Update components"))),
                   fluidRow(
                     sliderInput("nbDispGenes", "How many genes:", min = 1, max = 100, value = 25, step = 1)
                   ),
                   fluidRow(
                     plotOutput("components1"),
                     plotOutput("components2"),
                     plotOutput("components3")
                   )
          ),
          tabPanel("Plots",
                   fluidRow(
                     column(2, radioButtons("PCAcolor", "Color points by:", choices = list("biological groups" = 1, "dendrogram clusters" = 2, "k-means clusters" = 3), selected = 1)),
                     column(2, checkboxInput("showEllipse", "Show ellipse(s)", value = F))
                   ),
                   tabsetPanel(
                     tabPanel("2D",
                              fluidRow(column(1, actionButton("updatePCAPlots", "Update plots"))),
                              #
                              fluidRow(
                                column(6, plotOutput("interactPCA12li", brush = brushOpts(id = "PCA12librush", resetOnNew = T), width = 470, dblclick = "PCA12lidblclick"), align="center"),
                                column(6, plotOutput("interactPCA12co", brush = brushOpts(id = "PCA12cobrush", resetOnNew = T), width = 470, dblclick = "PCA12codblclick"), align="center")
                              ),
                              fluidRow(
                                column(6, dataTableOutput("dataPCA12li")),
                                column(6, dataTableOutput("dataPCA12co"))
                              ),
                              #
                              fluidRow(
                                column(6, plotOutput("interactPCA13li", brush = brushOpts(id = "PCA13librush", resetOnNew = T), width = 470, dblclick = "PCA13lidblclick"), align="center"),
                                column(6, plotOutput("interactPCA13co", brush = brushOpts(id = "PCA13cobrush", resetOnNew = T), width = 470, dblclick = "PCA13codblclick"), align="center")
                              ),
                              fluidRow(
                                column(6, dataTableOutput("dataPCA13li")),
                                column(6, dataTableOutput("dataPCA13co"))
                              ),
                              #
                              fluidRow(
                                column(6, plotOutput("interactPCA32li", brush = brushOpts(id = "PCA32librush", resetOnNew = T), width = 470, dblclick = "PCA32lidblclick"), align="center"),
                                column(6, plotOutput("interactPCA32co", brush = brushOpts(id = "PCA32cobrush", resetOnNew = T), width = 470, dblclick = "PCA32codblclick"), align="center")
                              ),
                              fluidRow(
                                column(6, dataTableOutput("dataPCA32li")),
                                column(6, dataTableOutput("dataPCA32co"))
                              )
                     ),#2D
                     tabPanel("3D",
                              sliderInput("pca3ddotsize", "Dot size: ", min = 1, max = 20, step = 1, value = 2),
                              actionButton("generatepca3d", "Update plot"),
                              plotlyOutput("pca3D", width = "950px", height = "750px")
                     )#3D
                   )#tabsetPanel Plots
          ),#tabsetPanel Plots
          tabPanel("Checkplot",
                   fluidRow(column(2, radioButtons("dataCheckplot", "Dataset:", list("Total" = "total", "Selection" = "selection"), selected = "selection"))),
                   fluidRow(
                     column(2, textInput("geneNameCheckplot", label = "Gene name (or part of it):", value = "geneName...")),
                     column(1, actionButton("updatelistcheckplot", "Search gene"))
                   ),
                   fluidRow(
                     column(2, uiOutput("geneNamescheckplotUI")),
                     column(2, uiOutput("sampleNamescheckplotUI"))
                   ),
                   fluidRow(
                     column(2, numericInput("heightcheckplot", "Plot height", min = 250, max = 10000, value = 250)),
                     column(2, checkboxInput("facetcheckplot", label = "Separate by group", value = FALSE)),
                     column(1, actionButton("updatecheckplot", "Update checkplot"))
                   ),
                   uiOutput("checkplotUI")
          )
        )
      ),#tabItems PCA

      # CA
      # ==
      tabItem(tabName = "CA",
        tags$i(class = "fa fa-area-chart fa-5x"),
        tags$h4(class = "titre", "Correspondence Analysis"),
        helpText("Need a description..."),
        tags$hr(),
        helpText("Needs to have generated PCA first."),
        tabsetPanel(
          tabPanel("Settings",
                   fluidRow(
                     selectInput("selectAxisCoA", h4("PCA Axis of interest on which we find the interesting genes", style = "color: rgb(0, 123, 16)"), choices = list("Axis 1" = 1, "Axis 2" = 2, "Axis 3" = 3), selected = 1)
                   ),
                   fluidRow(column(1, actionButton("generateContributionBoxplot", "Generate contribution boxplot"))),
                   fluidRow(
                     plotOutput("contributionBoxplot", width = "950px", height = "200px")
                   ),
                   fluidRow(
                     numericInput("nbGenesToKeep", h4("Threshold", style = "color: rgb(0, 123, 16)"), min = 0, max = 100, value = 20, step = 1)
                   ),
                   fluidRow(column(1, actionButton("applyThreshold", "Apply threshold"))),
                   fluidRow(
                     column(3, verbatimTextOutput("numberGenesThresholded"))
                   ),
                   fluidRow(
                     plotOutput("thresholdBoxplot", width = "950px", height = "200px")
                   ),
                   fluidRow(
                     column(6, dataTableOutput("thresholdedPrint"))
                   )
          ),#Settings
          tabPanel("Summary",
                   fluidRow(column(1, actionButton("updateCoA", "Update CoA"))),
                   plotOutput("coaeigenvalues"),
                   verbatimTextOutput("coasummary")
          ),#Summary
          tabPanel("Plots",
                   radioButtons("COAcolor", "Color points by:", choices = list("biological groups" = 1, "dendrogram clusters" = 2, "k-means clusters" = 3), selected = 1),
                   tabsetPanel(
                     tabPanel("2D",
                              fluidRow(column(1, actionButton("updateCoAPlots", "Update plots"))),
                              tags$br(),
                              plotOutput("interactCOA12", brush = brushOpts(id = "COA12brush", resetOnNew = T), width = 470, dblclick = "COA12dblclick"),
                              dataTableOutput("dataCOA12"),
                              plotOutput("interactCOA13", brush = brushOpts(id = "COA13brush", resetOnNew = T), width = 470, dblclick = "COA13dblclick"),
                              dataTableOutput("dataCOA13"),
                              plotOutput("interactCOA32", brush = brushOpts(id = "COA32brush", resetOnNew = T), width = 470, dblclick = "COA32dblclick"),
                              dataTableOutput("dataCOA32")
                     ),#2D
                     tabPanel("3D",
                              sliderInput("coa3ddotsize", "Dot size: ", min = 1, max = 20, step = 1, value = 2),
                              actionButton("generatecoa3d", "Update plot"),
                              plotlyOutput("coa3D", width = "950px", height = "750px")
                     )#3D
                   )#tabsetPanel
          )#Plots
        )#TabsetPanel
      ),#tabItems CA

      # EXPORT
      # ======
      tabItem(tabName = "export",
        tags$i(class = "fa fa-share-square-o fa-5x"),
        tags$h4(class = "titre", "Export"),
        helpText("Need a description..."),
        tags$hr(),
        helpText("Do not forget to 'validate parameters'."),
        helpText("Exportation is done in RDS format."),
        fluidRow(
          column(3, textInput("genesRDSName", label = h4("Genes : choose a name", style = "color: rgb(0, 16, 148)", align = "center"), value = paste0("genes_", gsub("\\s", "_", gsub("\\d{2}[:::]\\d{2}[:::]\\d{2}\\s", "", date())), ".rds")))
        ),
        fluidRow(
          column(1, actionButton("exportGenes", "Export 'genes' with current selection"))
        ),
        fluidRow(
          column(3, textInput("libsRDSName", label = h4("Libs : choose a name", style = "color: rgb(0, 16, 148)", align = "center"), value = paste0("libs_", gsub("\\s", "_", gsub("\\d{2}[:::]\\d{2}[:::]\\d{2}\\s", "", date())), ".rds")))
        ),
        fluidRow(
          column(1, actionButton("exportLibs", "Export 'libs' with current selection"))
        )
      ),#tabItems export
      
      # ABOUT
      # =====
      
      tabItem(tabName = "about",
              tags$i(class = "fa fa-question-circle fa-5x"),
              tags$h4(class = "titre", "About"),
              helpText("Need a description..."),
              tags$hr(),
              helpText("Nothing for now...")
      )#tabItems about

    )#tabItems
    
  )#dashboardBody
  
)#dashboardPage
