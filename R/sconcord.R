#' @title Sparse concordance analysis with cross-validation
#' 
#' @description Sparse concordance analysis, for each component, the optimal sparsity level is 
#'   optimized using cross-validation
#' 
#' @param x A list of predictive matrices
#' @param y the dependent matrix
#' @param opt.kx a numeric vector, candidate kx to be evaluated
#' @param opt.ky a numeric vector, candidate ky to be evaluated
#' @param fold the number of fold for cross validation
#' @param ncores the number of cores to be used, passed to \code{mclapply} in "parallel" package.
#' @param dmod see \code{\link{concord}}
#' @param center see \code{\link{concord}}
#' @param scale see \code{\link{concord}}
#' @param option see \code{\link{concord}}
#' @param ncomp see \code{\link{concord}}
#' @param kx see \code{\link{concord}}
#' @param ky see \code{\link{concord}}
#' @param wx see \code{\link{concord}}
#' @param wy see \code{\link{concord}}
#' @param pos see \code{\link{concord}}
#' 
#' @author Chen Meng
#' @return concordance with results of optimization of kx, ky params 
#' @importFrom parallel mclapply
#' @importFrom stats cor
#' @importFrom matrixStats rowSds rowMedians
#' @export
#' @examples
#' # example here
#' 

sconcord <- function(x, y, opt.kx, opt.ky, fold = 5, ncores = 1, ncomp = 1,
                     dmod = 1, center = TRUE, scale = FALSE, option = "uniform", 
                     wx = 1, wy = 1, pos = FALSE) {

  call <- match.call()
  
  "loading.x" <- "loading.y" <- "score.xsep" <- 
    "score.x" <- "score.y" <- "loading.x.index" <- 
    "score.x.index" <- "var" <- c()             
  # "normed" <- "deflated" <- list()
  opt <- list()
  for (i in 1:ncomp) {
    cat(paste("calculating component", i, "...\n"))
    cat("cross validation ...")
    ctv <- cv.concord(x = x, y = y, fold = fold, opt.kx = opt.kx, opt.ky = opt.ky, 
                      ncores = ncores, 
                      center = center, scale = scale, option = option, 
                      dmod = 0, pos = pos, wx = wx, wy = wy)
    
    okx <- opt.kx[which.max(rowMedians(ctv$cvx))]
    oky <- opt.ky[which.max(rowMedians(ctv$cvy))]
    cat("fit model ...")
    res <- concord(x, y, ncomp = 1, kx = okx, ky = oky, verbose = FALSE, 
                   center = center, scale = scale, option = option, 
                   dmod = dmod, pos = pos, wx = wx, wy = wy)
    enter <- FALSE
    scale <- FALSE
    option <- "uniform"
    
    y <- res$deflated$y
    x <- res$deflated$x
    
    loading.x <- cbind(loading.x, res$loading.x)
    loading.y <- cbind(loading.y, res$loading.y)
    score.xsep <- cbind(score.xsep, res$score.xsep)
    score.x <- cbind(score.x, res$score.x)
    score.y <- cbind(score.y, res$score.y)
    var <- c(var, res$var)
    # normed[[i]] <- res$normed
    # deflated[[i]] <- res$deflated
    opt[[i]] <- ctv
  }
  loading.x.index <- res$loading.x.index
  score.x.index <- res$score.x.index
  
  list(loading.x = loading.x, 
       loading.y = loading.y,
       score.xsep = score.xsep,
       score.x = score.x,
       score.y = score.y, 
       loading.x.index = loading.x.index, 
       score.x.index = score.x.index,
       var = var,
       cv = opt,
       # normed = normed, 
       # deflated = deflated, 
       call = call)
  
  
}