#' concordance analysis
#' @param x a list of predictive matrices. The columns are observations, rows are varaibles. The columns (observations)
#' has to be matched. 
#' @param y a response matrix. Rows are variables, columns are observations. The columns should be matched with columns in x. 
#' @param ncomp the number of components want to retain
#' @param dmod the deflation mode, dmod = 2 is the original publication of 
#' @param center logical values, whether the variables should be centered
#' @param scale logical values, whether the variables should be scaled
#' @param option the option for normalizing matrix
#' @param kx the number (if it is an integer > 1) or the proportion (if 0 < ky < 1) of kept variables in x. It should be
#'   a numeric value.
#' @param ky the number (if it is an integer > 1) or the proportion (if 0 < ky < 1) of kept variables in y. It should be 
#'   a numeric value. 
#' @param wx weight for the rows of x
#' @param wy weight for the rows of y
#' @param pos logical value, whether only non-negative loadings retained
#' @param verbose if the process of calculation should be printed
#' @param ncores number of cores to be used, passed to \code{\link{mclapply}}
#' @param fold the number of fold to be used in cross-validation, only used if kx or ky is a vector
#' @param nstart how many time the k-fold cross validation should be done
#' @param seed set seed for random number generation
#' @param loorss if the Leave-one-out procedure should be used in matrix reconstruction
#' @param scan If the PRESS plot should be shown and used to determine the optimal k in CV
#' @param nsd the the n*sd for selecting k automatically
#' @param init how to initialize the algorithm. if no sparsity, svd is fast.
#' 
#' @author Chen Meng
#' @export
#' @examples 
#' library(omic3plus)
#' data("NCI60_4arrays")
#' 
#' y <- as.matrix(NCI60_4arrays$agilent)
#' x <- lapply(NCI60_4arrays[2:4], as.matrix)
#' 
#' # no sparsity
#' con1 <- concord(x, y, ncomp = 3)
#' 
#' # sparsity on rows of x, select 10% genes
#' con <- concord(x, y, ncomp = 3, kx = 0.1)
#' 
#' # sparsity on rows of both x and y, select 10% genes
#' con <- concord(x, y, ncomp = 3, kx = 0.1, ky = 0.1, option = "nk")
#' plot(con$score.x[, 1], con$score.y[, 1])
#' abline(a = 0, b = 1)


concord <- function(x, y, ncomp=2, dmod = 1, center = TRUE, scale = FALSE, option = "uniform", 
                    kx = "all", ky = "all", wx = 1, wy = 1, pos = FALSE, verbose = TRUE, 
                    init = c("svd", "average")[2],
                    # for cv
                    ncores = 1, fold = 5, nstart = 1, seed = NULL, loorss = FALSE, 
                    scan = TRUE, nsd = 2) {
  
  option <- match.arg(option, c("uniform", "lambda1", "inertia", "nk"))
  call <- match.call()
  #
  if (kx[1] == "all") kx <- Inf
  if (ky[1] == "all") ky <- Inf
  if (kx[1] > 0 & kx[1] < 1) 
    kx <- ceiling(sum(sapply(x, nrow)) * kx)
  if (ky[1] > 0 & ky[1] < 1)
    ky <- ceiling(nrow(y) * ky)
  
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
  
  ##
  val <- switch(option,
                "uniform" = 1,
                "lambda1" = svd(Ynorm)$d[1],
                "inertia" = sqrt(sum(Ynorm^2)), 
                "nk" = sqrt(min(nrow(y), ky)))
  
  Xnorm <- lapply(x, t)  
  Xnorm <- processOpt(Xnorm, center = center, scale = scale, option = option, value = val, kx = kx)
  Xcat <- do.call("cbind", Xnorm)
  Ynorm.o <- Ynorm
  Xnorm.o <- Xnorm
  
  ys <- yloading <- gls <- bls <- loading <- var <- c()
  
  for (f in 1:ncomp) {
    if (verbose)
      cat(paste("calculating component", f, "...\n"))
    
    S <- t(Ynorm) %*% Xcat
    ok <- cv.softSVD(S, nf = 1, kv.opt = kx, ku.opt = ky, wv = wx, wu = wy, pos = pos, 
                     maxiter = 1000, verbose = TRUE, ncores = ncores, fold = fold, init = init,
                     nstart = nstart, seed = seed, loorss = loorss, scan = scan, nsd = nsd)
    if (verbose)
      cat(paste0("optimal kx = ", ok$sel.v, "; optimal ky = ", ok$sel.u, ".\n"))
    decom <- softSVD(x = S, nf = 1, kv = ok$sel.v, ku = ok$sel.u, wv = wx, wu = wy,init = init, 
                     pos = pos, maxiter = 1000, verbose = FALSE)
    
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
      # S <- S - tcrossprod(decom$u[, 1]) %*% S
      # or the same 
      Ynorm <- Ynorm - Ynorm %*% tcrossprod(decom$u[, 1])
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
    } else if (dmod == 0) {
      cat("no deflation\n")
    } else {
      stop("unknown deflation mode")
    }
    
    # prepare output
    xd <- lapply(names(x), function(it) {
      t(Xcat[, i.feature == it])
    })
    names(xd) <- names(x)
    yd <- t(Ynorm)
    
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
       normed = list(y = Ynorm.o, x = Xnorm.o), 
       deflated = list(y = yd, x = xd), 
       call = call)
}

