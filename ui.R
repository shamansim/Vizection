dashboardPage(
  
  dashboardHeader(title = tags$img(src = "vizection_white.png", alt = "Vizection", height = "40"),
    dropdownMenuOutput("tasksMenu")
  ),#dashboardHeader
  
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
      )#Export

    )#sidebarMenu
    
  ),#dashboardSidebar
  
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

        tabsetPanel(
          tabPanel("Dendrogram",
            tags$i(class = "fa fa-sitemap fa-5x"),
            tags$h4(class = "titre", "Dendrogram"),
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
          tabPanel("Heatmap",
            tags$i(class = "fa fa-th-large fa-5x"),
            tags$h4(class = "titre", "Heatmap"),
            tags$hr()
          )#tabPanel Heatmap
        )#tabsetPanel

      ),#tabItem hierarchicalClustering

      # K-MEANS
      # =======
      tabItem(tabName = "k-means",
        tags$i(class = "fa fa-dot-circle-o fa-5x"),
        tags$h4(class = "titre", "K-means"),
        tags$hr(),
        helpText("Nothing more here !")
      ),#tabItem k-means

      # PCA
      # ===
      tabItem(tabName = "PCA",
        tags$i(class = "fa fa-area-chart fa-5x"),
        tags$h4(class = "titre", "Principal Components Analysis"),
        tags$hr(),
        helpText("Nothing for now...")
      ),#tabItems PCA

      # CA
      # ==
      tabItem(tabName = "CA",
        tags$i(class = "fa fa-area-chart fa-5x"),
        tags$h4(class = "titre", "Correspondence Analysis"),
        tags$hr(),
        helpText("Nothing more here !")
      ),#tabItems CA

      # EXPORT
      # ======
      tabItem(tabName = "export",
        tags$i(class = "fa fa-share-square-o fa-5x"),
        tags$h4(class = "titre", "Export"),
        tags$hr(),
        helpText("Nothing for now...")
      )#tabItems export

    )#tabItems
    
  )#dashboardBody
  
)#dashboardPage
