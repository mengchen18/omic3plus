#' @title calculating distance
#' @description distance calculation
#' @param l a list of character or numeric vector 
#' @details the distance between two vectors are 
#'   d = size(union(A, B))/min(size(A), size(B))
#' @return a \code{dist} object
#' @example 
#' set.seed(5)
#' li <- lapply(4:7, function(x) sample(letters, x))
#' minjaccard(li)
#' li

minjaccard <- function(l) {
  
  if (!inherits(l, "list"))
    stop("l should be list.")
  
  a <- unique(unlist(l))
  am <- sapply(l, function(x) as.integer(a %in% x))
  inter <- crossprod(am)
  ms <- rep(colSums(am), ncol(am))
  ms <- matrix(ms, nrow = ncol(am))
  mst <- t(ms)
  ms <- pmin(ms, mst)
  d <- 1-inter/ms
  as.dist(d)
}

setGeneric("minjaccard")

setMethod(minjaccard, "GeneSetCollection", 
          definition = 
            function(l) {
              l <- lapply(l, geneIds)
              minjaccard(l)
            })



