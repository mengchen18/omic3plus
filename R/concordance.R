#' concordance analysis
#' @param x a list of predictive matrices. The columns are observations, rows are varaibles. The columns (observations)
#' has to be matched. 
#' @param y a response matrix. Rows are variables, columns are observations. The columns should be matched with columns in x. 
#' @param ncomp the number of components want to retain
#' @param dmod the deflation mode, dmod = 2 is the original publication of 
#' @param center logical values, whether the variables should be centered
#' @param scale logical values, whether the variables should be scaled
#' @param option the option for normalizing matrix
#' @param kx the 
#' @param ky the number (if it is an integer > 1) or the proportion (if 0 < ky < 1) of 
#' @param wx weight for the rows of x
#' @param wy weight for the rows of y
#' @param pos logical value, whether only non-negative loadings retained
#' @param verbose if the process of calculation should be printed
#' 
#' @author Chen Meng
#' @export

concord <- function(x, y, ncomp=2, dmod = 1, center = TRUE, scale = FALSE, option = "uniform", 
                    kx = "all", ky = "all", wx = 1, wy = 1, pos = FALSE, verbose = TRUE) {

  option <- match.arg(option, c("uniform", "lambda1", "inertia"))
  call <- match.call()
  if (kx == "all") kx <- Inf
  if (ky == "all") ky <- Inf
  
  nmat <- length(x)
  if (is.null(names(x)))
    names(x) <- paste0("X", 1:nmat)
  
  nr <- sapply(x, nrow)
  nc <- unique(sapply(x, ncol))
  if (length(nc) > 1)
    stop("Number of columns in X need to be the same.")
  
  i.sample <- rep(names(x), each = nc)
  i.feature <- rep(names(x), nr)
  
  Ynorm <- scale(t(y), center = center, scale = scale)
  
  val <- switch(option,
                "uniform" = 1,
                "lambda1" = svd(Ynorm)$d[1],
                "inertia" = sqrt(sum(Ynorm^2)), 
                "nrow" = sqrt(nrow(y)))
  
  Xnorm <- processOpt(x, center = center, scale = scale, option = option, value = val)
  
  Xnorm <- lapply(Xnorm, t)
  Xcat <- do.call("cbind", Xnorm)
  Ynorm.o <- Ynorm
  Xnorm.o <- Xnorm
  
  ys <- yloading <- gls <- bls <- loading <- var <- c()
  
  for (f in 1:ncomp) {
    print(f)
    if (f == 1 || dmod != 1 )
      S <- t(Ynorm) %*% Xcat
    if (f == 1)
      S.o <- S

    decom <- softSVD(x = S, nf = 1, kv = kx, ku = ky, wv = wx, wu = wy, 
                     pos = pos, maxiter = 1000, verbose = verbose)
    
    xa <- Xcat %*% decom$v[, 1]
    yb <- Ynorm %*% decom$u[, 1]
    var <- c(var, decom$d[1]^2)
    ys <- cbind(ys, yb)
    gls <- cbind(gls, xa)
    yloading <- cbind(yloading, decom$u[, 1])
    loading <- cbind(loading, decom$v[, 1])
    
    xai <- lapply(names(x), function(x) {
      ii <- i.feature == x
      Xcat[, ii] %*% decom$v[ii, 1]
    })
    
    xai.var <- sapply(xai, crossprod)
    bls <- cbind(bls, unlist(xai))
    
    if (dmod == 1) {
      # deflation of S, crossprod matrix, the save with SVD directly
      # the classical concordance approach
      S <- S - tcrossprod(decom$u[, 1]) %*% S
      # or the same 
      # Ynorm <- Ynorm - Ynorm %*% tcrossprod(decom$u[, 1])
    } else if (dmod == 2) {
      # deflaltion X using loading of X, as Lafosse & Hanafi 1997
      # but not possible to incorporate with sparse factor
      Xcat <- Xcat - Xcat %*% tcrossprod(decom$v[, 1])
    } else if (dmod == 3) {
      # defaltion Y using its normed score
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(normvec(yb)))
    } else if (dmod == 4) {
      # defaltion X and Y using normed score of Y, not suggested
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(normvec(yb)))
      Xcat <- Xcat - t(t(Xcat) %*% tcrossprod(normvec(yb)))
    } else if (dmod == 5) {
      # defaltion X and Y using normed score of Y, not suggested
      Ynorm <- Ynorm - t(t(Ynorm) %*% tcrossprod(normvec(yb)))
      Xcat <- Xcat - t(t(Xcat) %*% tcrossprod(normvec(xa)))
    } else {
      stop("unknown deflation mode")
    }
    
  }
  
  rownames(loading) <- colnames(Xcat)
  
  list(loading.x = loading, 
       loading.y = yloading,
       score.xsep = bls,
       score.x = gls,
       score.y = ys, 
       loading.x.index = i.feature, 
       score.x.index = i.sample,
       var = var,
       norm.y = Ynorm.o,
       norm.x = Xnorm.o, 
       call = call)
}


