#' @title cross-validation for softSVD
#' @description This function use k-fold cross-valiation method to optimize the sparsity
#'   of right singular values
#' @param x input matrix
#' @param nf number of component
#' @param kv.opt optional value for sparsity on right singular value
#' @param wv weight for columns
#' @param wu weight for rows
#' @param pos whether retein non-negative results
#' @param maiter maximum number of iteration
#' @param tol convergence tolerance
#' @param verbose if print the progress
#' @param ncores the number of cores used, passed to mclapply
#' @param fold fold number in cross validation
#' @param nstart how many time the k-fold cross validation should be done
#' @param seed set seed
#' @param loorss if the Leave-one-out procedure should be used in matrix reconstruction
#' 
#' @importFrom parallel mclapply
#' @seealso \code{\link{cv.softSVD}}


cv.softKv <- function(x, nf = 1, kv.opt = c(0.3, 0.5, 0.8), wv = 1, wu = 1, pos = FALSE,
                      maxiter = 50, tol = sqrt(.Machine$double.eps), verbose = FALSE,
                      ncores = 1, fold = 5, nstart = 1, seed = NULL, loorss = FALSE) {
  
  set.seed(NULL)
  
  if (any(ikv <- kv.opt > nrow(x))) {
    cat("Some optional k exceed the dimenson of x, removed. ")
    kv.opt <- kv.opt[!ikv]
  }
  
  nr <- nrow(x)
  lab <- 1:fold
  
  val <- replicate(nstart, simplify = FALSE, {
    
    ivec <- sample(rep(lab, ceiling(nr/fold))[1:nr])
    rl <- mclapply(lab, mc.cores = ncores, function(ii) {
      
      if (verbose) 
        cat(paste("fold", ii, "...\n"))
      
      i <- ii == ivec
      m.train <- x[!i, ]
      m.test <- x[i, ]
      
      sapply(kv.opt, function(k) {
        s0 <- softSVD(m.train, nf = nf, kv = k, maxiter = maxiter, wv = wv, 
                      wu = wu, pos = pos, tol = tol, verbose = verbose)
        if (!loorss) {
          sum((m.test - tcrossprod(m.test %*% s0$v, s0$v))^2)
        } else {
          q <- sapply(1:ncol(m.test), function(r) 
            sum((m.test[, r] - s0$v[r] * (m.test[, -r] %*% s0$v[-r, ])) ^2)
            )
          sum(q)
        }
      })
    })
    do.call("rbind", rl)
  })
  do.call("rbind", val)
}



#' @title cross-validation for softSVD
#' @description This function use k-fold cross-valiation method to optimize the sparsity
#'   of right singular values
#' @param x input matrix
#' @param nf number of component
#' @param kv.opt optional value for sparsity on right singular value
#' @param ku.opt optional value for sparsity on left singular value
#' @param wv weight for columns
#' @param wu weight for rows
#' @param pos whether retein non-negative results
#' @param maiter maximum number of iteration
#' @param tol convergence tolerance
#' @param verbose if print the progress
#' @param ncores the number of cores used, passed to mclapply
#' @param fold fold number in cross validation
#' @param nstart how many time the k-fold cross validation should be done
#' @param seed set seed
#' @param loorss if the Leave-one-out procedure should be used in matrix reconstruction
#' 
#' @importFrom parallel mclapply
#' @export
#' @keywords internal
#' @return list consist of two matrix. 
#'  - cvv - The PRESS for right singular vector
#'  - cvu - the PRESS for left singular vector
#' @examples 
#' ###
#' seed <- rnorm(8)
#' noise1 <- matrix(rnorm(8*15, sd = 0.3), 15, 8)
#' noise2 <- matrix(rnorm(8*20, sd = 0.3), 8, 20)
#' noise1[1:4, ] <- noise1[1:4, ] + rbind(seed, seed, seed, seed)
#' noise2[, 1:6] <- noise2[, 1:6] + seed
#' mat <- noise1 %*% noise2
#' dim(mat)
#' 
#' cv <- cv.softSVD(x = mat, nf = 1, kv.opt = 1:10, ku.opt = 1:8)
#' boxplot(cv$cvu)
#' boxplot(cv$cvv)

cv.softSVD <- function(x, nf = 1, kv.opt = c(0.3, 0.5, 0.8), ku.opt = c(0.3, 0.5, 0.8), 
                       wv = 1, wu = 1, pos = FALSE, maxiter = 50, 
                       tol = sqrt(.Machine$double.eps), verbose = FALSE,
                       ncores = 1, fold = 5, nstart = 1, seed = NULL, loorss = FALSE) {
  
  ov <- ou <- "all" 
  if (length(kv.opt > 1)) {
    ov <- cv.softKv(x, nf = nf, kv.opt = kv.opt, wv = wv, wu = wu, pos = pos,
              maxiter = maxiter, tol = tol, verbose = FALSE,
              ncores = ncores, fold = fold, nstart = nstart, 
              seed = seed, loorss = loorss)
  }
  if (length(ku.opt > 1)) {
    xt <- t(x)
    wvt <- wu
    wut <- wv
    ou <- cv.softKv(xt, nf = nf, kv.opt = ku.opt, wv = wvt, wu = wut, pos = pos,
                    maxiter = maxiter, tol = tol, verbose = FALSE,
                    ncores = ncores, fold = fold, nstart = nstart, 
                    seed = seed, loorss = loorss)
  }
    
  list(cvu = ou, cvv = ov)
}


