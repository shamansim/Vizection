UNaddNumberOfSamplesOrGroup <- function(Checklist){
  result <- c()
  for(i in Checklist){
    result <- c(result, i %>% sapply(., function(x) gsub("\\s[:|:]\\s.*", "", x)))
  }
  return(result)
}

#' filterSelectionBool

filterSelectionBool <- function(input) {
  filterByCounts <- libs$counts > input$nbFilterExtracted
  filterByGroup  <- rownames(libs) %in% (
    libs %>%
      select(samplename, group) %>%
      filter(group %in% UNaddNumberOfSamplesOrGroup(input$groupsCheck)) %$%
      samplename)
  filterByCounts & filterByGroup
}

#' filterSelectionBoolFinal

filterSelectionBoolFinal <- function(input) {
  filterSelectionBool(input) &
    (libs$samplename %in%
       UNaddNumberOfSamplesOrGroup(paste(input$samplesCheck)))
}

#' subgenes

subgenes_1 <- function(input, genes)
  genes[, filterSelectionBoolFinal(input)]

subgenes_2 <- function(pre_subgenes)
  pre_subgenes[apply(pre_subgenes, 1, sum) != 0, ] # removing useless genes

subgenes <- function(input, genes)
  subgenes_1 %>% subgenes_2

#' sublibs

sublibs <- function(input) {
  sublibs0 <- libs[filterSelectionBoolFinal(input), ]
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
