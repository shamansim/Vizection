#' vizectionValidate
#' 
#' Validate objects passed to Vizection
#' 
#' Makes an error with objects that do not fit the requirements for being
#' used in Vizection
#' 
#' @param genes A \sQuote{genes} table.
#' @param libs A \sQuote{libs} table.
#' 
#' @importFrom assertthat assert_that
#' @importFrom magrittr "%>%"
#' 
#' @export
#' 
#' @examples
#' # Make arbitrary input data using the "iris table".
#' ilibs <- iris
#' ilibs$group <- ilibs$Species
#' ilibs$samplename <- paste("sample", 1:150)
#' rownames(ilibs) <- ilibs$samplename
#' ilibs$Run <- "Run1"
#' ilibs$counts <- 10
#' 
#' igenes <- t(iris[,1:4]) %>% data.frame
#' colnames(igenes) <- rownames(ilibs)
#' 
#' # Check that the format is valid (silently returns TRUE).
#' vizectionValidate(genes=igenes, libs=ilibs)

vizectionValidate<- function(genes, libs) {
  
  validateGenes <- function() {
    assert_that(all(genes >= 0))
  }
  
  validateLibs <- function() {
    assert_that("counts" %in% colnames(libs))
    assert_that("group" %in% colnames(libs))
    assert_that("samplename" %in% colnames(libs))
  }
  
  ( ifelse(missing(genes), TRUE, validateGenes()) &
    ifelse(missing(libs),  TRUE, validateLibs())    ) %>% invisible
}
