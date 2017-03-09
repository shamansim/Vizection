#' vizectionExampleLibs
#' 
#' Toy example for the "libs" table to be used in vizection.
#' 
#' This is for debugging purposes; the example does not reflect well
#' the kind of real data that is expected.
#' 
#' @example 
#' summary(vizectionExampleLibs())
#' head(vizectionExampleLibs())

vizectionExampleLibs <- function() {
  data(iris)
  genes <- iris[, 1:4]
  libs <- as.data.frame(iris[, "Species"])
  colnames(libs) <- c('group')
  libs$samplename <- as.character(1:nrow(genes))
  rownames(libs) <- libs$samplename 
  libs$counts <- rnorm(n = nrow(genes), mean = 1000, sd = 200)
  libs$Run <- "Run1" %>% factor
  libs
}

#' vizectionExampleGenes
#' 
#' Toy example for the "genes" table to be used in vizection.
#' 
#' This is for debugging purposes; the example does not reflect well
#' the kind of real data that is expected.
#' 
#' @example 
#' summary(vizectionExampleGenes())
#' head(vizectionExampleGenes())

vizectionExampleGenes <- function() {
  data(iris)
  genes <- iris[, 1:4]
  genes <- as.data.frame(t(genes))
}

#' vizectionExampleEnv
#' 
#' Toy example of "input" object passed by Shiny
#' 
#' This is for debugging purposes; the example does not reflect well
#' the kind of real data that is expected.

vizectionExampleEnv <- function()
  list( nbFilterExtracted = 0
      , nbClusters        = 3
      , nbDispGenes       = 10
      , showGroupsColor   = TRUE
      , groupsCheck       = c("setosa | 5", "virginica | 10")
      , samplesCheck      = c("1 | setosa", "150 | virginica"))

#' filterSelectionBool
#' 
#' @example 
#' filterSelectionBool( libs = vizectionExampleLibs()
#'                    , input = vizectionExampleEnv())

filterSelectionBool <- function(libs, input) {
  filterByCounts <- libs$counts > input$nbFilterExtracted
  filterByGroup  <- rownames(libs) %in% (
    libs %>%
      select(samplename, group) %>%
      filter(group %in% UNaddNumberOfSamplesOrGroup(input$groupsCheck)) %$%
      samplename)
  filterByCounts & filterByGroup
}

#' filterSelectionBoolFinal
#' 
#' @example 
#' filterSelectionBoolFinal( libs  = vizectionExampleLibs()
#'                         , input = vizectionExampleEnv())

filterSelectionBoolFinal <- function(libs, input) {
  filterSelectionBool(libs, input) &
    (libs$samplename %in%
       UNaddNumberOfSamplesOrGroup(input$samplesCheck))
}

#' filterExtractedBool

filterExtractedBool <- function(libs, input)
  libs$counts > input$nbFilterExtracted

#' subgenes

subgenes_1 <- function(libs, input, genes)
  genes[, filterSelectionBoolFinal(libs, input)]

subgenes_2 <- function(pre_subgenes)
  pre_subgenes[apply(pre_subgenes, 1, sum) != 0, ] # removing useless genes

subgenes <- function(libs, input, genes)
  subgenes_1 %>% subgenes_2

#' sublibs

sublibs <- function(libs, input) {
  sublibs0 <- libs[filterSelectionBoolFinal(libs, input), ]
  sublibs0$group %<>% extract(drop = T)
  sublibs0
}

#' addNumberOfSamples
#' 
#' Displays something like "groupname | n".
#' 
#' Takes a vector of group names name and constructs a vector of
#' strings made of the group name, a pipe separator and the number
#' of samples in the group.

addNumberOfSamples <- function(libs, groups){
  result <- c()
  for(i in groups){
    result <- c( result
               , paste0( i
                       , " | "
                       , libs$samplename[libs$group==i] %>% length))}
  result
}

#' addGroupName
#' 
#' Displays something like "samplename | groupname".
#' 
#' Takes a vector of sample names name and constructs a vector of
#' strings made of the sample name, a pipe separator and the name
#' of its group.

addGroupName <- function(libs, samples){
  result = c()
  for(i in samples){
    result <- c(result, paste0(i, " | ", libs$group[libs$samplename==i]))
  }
  return(result)
}

#' UNaddNumberOfSamplesOrGroup
#'
#' Same as above except that it does not keep name attributes.
#'
#' @param names Group or sample names to which other information
#'              have been added by the functions addGroupName or
#'              addNumberOfSamples.
#'
#' @example 
#' c("toto | 5", "H12 | toto") %>% UNaddNumberOfSamplesOrGroup

UNaddNumberOfSamplesOrGroup <- function(names)
  gsub("\\s[:|:]\\s.*", "", names)
