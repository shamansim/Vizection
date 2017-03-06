#' corMat
#' 
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
#' @example 
#' data.frame(1:3, 2:4, 6:4) %>% corMat()

corMat_1 <- function(genes)
  genes %>% tail(-1) %>% smallCAGEqc::TPM

corMat_2 <- log1p

corMat_3 <- cor
  
corMat <- function(genes)
  corMat_1(genes) %>% corMat_2 %>% corMat_3
