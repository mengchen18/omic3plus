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

targetResolve2 <- function(t, s, m, ncomp = 1, center = TRUE, scale = TRUE, option = "uniform",
         kt = "all", km = "all", wt = 1, wm = 1, pos = FALSE, ncores = 1, fold = 5, init = "average", 
         nstart = 1, seed = NULL, loorss = FALSE, scan = FALSE, nsd = 1, verbose = TRUE) {
  
  call <- match.call()

  y <- t %*% s 

  concord(
    m, y, ncomp=ncomp, dmod = 1, center = center, scale = scale, 
    center.y = FALSE, scale.y = FALSE, option = option, 
    kx = km, ky = kt, wx = wm, wy = wt, pos = pos, verbose = verbose, 
    init = init, ncores = ncores, fold = fold, nstart = nstart, 
    seed = seed, loorss = loorss, scan = scan, nsd = nsd)
  
  ###
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
  
  ll
}







