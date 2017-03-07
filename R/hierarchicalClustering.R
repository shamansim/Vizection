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
  genes %>% tail(-1) %>% smallCAGEqc::TPM()

corMat_2 <- log1p

corMat_3 <- cor
  
corMat <- function(genes)
  corMat_1(genes) %>% corMat_2 %>% corMat_3

#' distCorMat
#' 
#' Transforms a correlation matrix into a Euclidian distance matrix.
#' 
#' The final conversion is done with the quasieuclid function of the
#' ade4 package.
#' 
#' @param m A corelation matrix.
#' 
#' @example 
#' data.frame(1:3, 2:4, 6:4) %>% corMat %>% distCorMat

distCorMat_1 <- function(m) {
  m %>%
    subtract(1, .) %>%
    divide_by(., 2) %>% 
    as.dist
}
  
distCorMat_2 <- ade4::quasieuclid

distCorMat <- function(m)
  m %>% distCorMat_1 %>% distCorMat_2