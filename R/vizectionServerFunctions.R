UNaddNumberOfSamplesOrGroup <- function(Checklist){
  result <- c()
  for(i in Checklist){
    result <- c(result, i %>% sapply(., function(x) gsub("\\s[:|:]\\s.*", "", x)))
  }
  return(result)
}

filterSelectionBool <- function(input) {
  filterByCounts <- libs$counts > input$nbFilterExtracted
  filterByGroup  <- rownames(libs) %in% (
    libs %>%
      select(samplename, group) %>%
      filter(group %in% UNaddNumberOfSamplesOrGroup(input$groupsCheck)) %$%
      samplename)
  filterByCounts & filterByGroup
}

filterSelectionBoolFinal <- function(input) {
  filterSelectionBool(input) &
    (libs$samplename %in%
       UNaddNumberOfSamplesOrGroup(paste(input$samplesCheck)))
}