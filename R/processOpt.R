#' @title Preprocess a list of matrices for concordance analysis
#' 
#' @description Internal function used to normalize a list of matrices passed into concordance analysis
#' 
#' @param x a list of matrices to be preprocessed. The columns are variables and rows are observations
#' @param center if the variables should be centered
#' @param scale if the variables should be scaled
#' @param option options for normalizing matrices, possible value are "lambda1", "inertia", "uniform", "nk"
#' @param value Only used if option is "lambda1" or inertia, the matrix are weighted by value/lambda1 
#' @param kx used if option = "nk", normalize according to the number of kept variables
#' (or inertia)
#' 
#' @author Chen Meng
#' @return a list of preprocess matrices has the same dimension as input x
#' 
#' @keywords internal
#' @export
#' @examples 
#' data("NCI60_4arrays")
#' dat <- lapply(NCI60_4arrays, t)
#' v <- processOpt(dat, option = "lambda1")
#' 

processOpt <- 
  function(x, center=TRUE, scale=FALSE, option = c("lambda1", "inertia", "uniform", "nk"), value = 1, kx = NULL) {
    
    opt <- match.arg(option, c("lambda1", "inertia", "uniform", "nk"))
    nd <- length(x)
    
    if (is.null(names(x)))
      names(x) <- paste("data", 1:length(x), sep = "_")
    
    x <- lapply(x, scale, center, scale)
    if (opt == "lambda1") {
      w <- sapply(x, function(xx) value/svd(xx)$d[1])
    } else if (opt == "inertia") {
      w <- sapply(x, function(xx) value/sqrt(sum(xx^2)))
    } else if (opt == "uniform") {
      w <- rep(1, nd)
    } else if (opt == "nk") {
      nx <- sum(sapply(x, ncol))
      if (!is.null(kx)) 
        nx <- value/sqrt(min(kx, nx))
      w <- rep(nx, nd)
    }
    mapply(SIMPLIFY = FALSE, function(xx, ww) xx*ww, xx=x, ww=w)
  }



