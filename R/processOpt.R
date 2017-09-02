#' @title Preprocess a list of matrices for concordance analysis
#' 
#' @description Internal function used to normalize a list of matrices passed into concordance analysis
#' 
#' @param x a list of matrices to be preprocessed. The rows are variables and columns are observations
#' @param center if the variables should be centered
#' @param scale if the variables should be scaled
#' @param option options for normalizing matrices, possible value are "lambda1", "inertia", "uniform"
#' @param value Only used if option is "lambda1" or inertia, the matrix are weighted by value/lambda1 
#' (or inertia)
#' 
#' @author Chen Meng
#' @return a list of preprocess matrices has the same dimension as input x
#' 
#' @keywords internal
#' @export
#' @examples 
#' data("NCI60_4arrays")
#' v <- processOpt(NCI60_4arrays, option = "lambda1")
#' 
#' processOpt(NCI60_4arrays, option = "inertia")

processOpt <- 
function(x, center=TRUE, scale=FALSE, option = c("lambda1", "inertia", "uniform", "nrow"), value = 1) {
  
  opt <- match.arg(option, c("lambda1", "inertia", "uniform", "nrow"))
  
  if (is.null(names(x)))
    names(x) <- paste("data", 1:length(x), sep = "_")

  x <- lapply(x, scale, center, scale)
  if (opt == "lambda1") {
    w <- sapply(x, function(xx) value/svd(xx)$d[1])
  } else if (opt == "inertia") {
    w <- sapply(x, function(xx) value/sqrt(sum(xx^2)))
  } else if (opt == "uniform") {
    w <- rep(1, length(x))
  } else if (opt == "uniform") {
    w <- sapply(x, function(xx) value/sqrt(nrow(xx)*length(x)))
  }
  mapply(SIMPLIFY = FALSE, function(xx, ww) xx*ww, xx=x, ww=w)
}



