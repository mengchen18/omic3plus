#' @title Normalized a vector to unit length
#' @author Chen Meng
#' @description Normalized a vector or matrix with 1 row/column to unit length
#' @param x a vector or matrix has one row/column
#' @param center whether x should be centered
#' @return a normalized vector/matrix has the same length as x
#' @keywords internal
#' @export

normvec <- function(x, center = FALSE) {
  # normalize a vector to unit length
  if (NROW(x) > 1 & NCOL(x) > 1)
    stop("x should be a vector or matrix with 1 row/column.")
  if (center)
    x <- x - mean(x, na.rm = TRUE)
  if (NCOL(x) > 1)
    prd <- tcrossprod(x) else
      prd <- crossprod(x)
    
  length <- sqrt(c(prd))
  v <- x/length
  attr(v, "norm") <- length
  v
}
