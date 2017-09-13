#' @title resolving target space of drugs
#' @description resolving target space and identify biomarkers by integrating drug target, cell line sensitivity and 
#'   cell line molecular data
#' @param t target matrix in the form target (rows) x drug (columns), values should be the higher the more potent, e.g. pKD.
#' @param s the sensitivity matrix, in the form drug (rows) x cell line (columns), the higher the more sensitive (e.g. -log10(GI50))
#' @param m a list of matrices, abundance of molecular features in the form of variables (rows) x cell lines (columns)
#' @param ncomp the number of components want to retain
#' @param center logical values, whether the variables should be centered
#' @param scale logical values, whether the variables should be scaled
#' @param option the option for normalizing matrix
#' @param km the number (if it is an integer >= 1) or the proportion (if 0 < km < 1) of kept variables in m. It should be
#'   a numeric value.
#' @param kt the number (if it is an integer >= 1) or the proportion (if 0 < kt < 1) of kept variables in t. It should be 
#'   a numeric value. 
#' @param wm weight for the rows of m
#' @param wt weight for the rows of t
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

targetResolve <- function(t, s, m, ncomp = 1, center = TRUE, scale = TRUE, option = "uniform",
         kt = "all", km = "all", wt = 1, wm = 1, pos = FALSE, ncores = 1, fold = 5, init = "average", 
         nstart = 1, seed = NULL, loorss = FALSE, scan = FALSE, nsd = 1, verbose = TRUE) {
  
  call <- match.call()
  #
  if (km[1] == "all") km <- Inf
  if (kt[1] == "all") kt <- Inf
  if (km[1] > 0 & km[1] < 1) 
    km <- ceiling(sum(sapply(m, nrow)) * km)
  if (kt[1] > 0 & kt[1] < 1)
    kt <- ceiling(nrow(t) * kt)
  
  #
  nmat <- length(m)
  if (is.null(names(m)))
    names(m) <- paste0("X", 1:nmat)
  
  nr <- sapply(m, nrow)
  nc <- unique(sapply(m, ncol))
  if (length(nc) > 1)
    stop("Number of columns in X need to be the same.")
  
  i.sample <- rep(names(m), each = nc)
  i.feature <- rep(names(m), nr)
  
  ##
  val <- switch(option,
                "uniform" = 1,
                "lambda1" = svd(Ynorm)$d[1],
                "inertia" = sqrt(sum(Ynorm^2)), 
                "nk" = sqrt(min(nrow(y), kt)))
  
  Xnorm <- lapply(m, t)  
  Xnorm <- processOpt(Xnorm, center = center, scale = scale, option = option, value = val, kx = km)
  Xcat <- do.call("cbind", Xnorm)
  
  mat <- t %*% s %*% Xcat
  
  
  ll <- list(
    loading.m = matrix(NA, sum(nr), ncomp),
    loading.t = matrix(NA, nrow(t), ncomp),
    score.t = matrix(NA, ncol(t), ncomp),
    score.m = matrix(NA, nc, ncomp),
    score.ts = matrix(NA, ncol(t), ncomp),
    score.ms = matrix(NA, nc, ncomp),
    sdev = rep(NA, ncomp), 
    mindex = i.feature, 
    call = call
    )
  
  for (i in 1:ncomp) {
    ok <- cv.softSVD(mat, nf = 1, kv.opt = km, ku.opt = kt, wv = wm, wu = wt, pos = pos, 
                     maxiter = 1000, verbose = TRUE, ncores = ncores, fold = fold, init = init,
                     nstart = nstart, seed = seed, loorss = loorss, scan = scan, nsd = nsd)
    
    if (verbose)
      cat(paste0("Selected km = ", ok$sel.v, "; Selected kt = ", ok$sel.u, ".\n"))
    decom <- softSVD(x = mat, nf = 1, kv = ok$sel.v, ku = ok$sel.u, wv = wm, wu = wt,init = init, 
                     pos = pos, maxiter = 1000, verbose = FALSE)
    
    if (i < ncomp)
      mat <- mat - tcrossprod(decom$u[, 1]) %*% mat
    
    loading.m <- decom$v[, 1, drop = FALSE]
    loading.t <- decom$u[, 1, drop = FALSE]
    score.t <- crossprod(t, decom$u[, 1, drop = FALSE])
    score.m <- Xcat %*% decom$v[, 1, drop = FALSE]
    score.ts <- s %*% normvec(score.m)
    score.ms <- crossprod(s, normvec(score.t))
    
    ll$loading.m[, i] <- loading.m
    ll$loading.t[, i] <- loading.t
    ll$score.t[, i] <- score.t
    ll$score.m[, i] <- score.m
    ll$score.ts[, i] <- score.ts
    ll$score.ms[, i] <- score.ms
    ll$sdev[i] <- decom$d[1]
  }
  
  ll
}







