#' @title Cross-validation for concordance analysis
#' 
#' @description Cross-validation for concordance analysis to find optimal sparsity level of x and y
#' 
#' @param x A list of predictive matrices
#' @param y the dependent matrix
#' @param opt.kx a numeric vector, candidate kx to be evaluated
#' @param opt.ky a numeric vector, candidate ky to be evaluated
#' @param fold the number of fold for cross validation
#' @param ... other arguments to be passed to concord
#' 
#' @author Chen Meng
#' @return A list of two data.frames. One for x and one for y. Each data.frame has three columns:
#'   kx/ky (input values of opt.kx/opt.ky), mean - the cv (cross-validation) mean, sd - the sd
#'   of cv-error 
#' 
#' @export
#' @importFrom matrixStats rowSds
#' @examples
#' 

cv.concord <- function(x, y, fold = 7, opt.kx = seq(0.1, 0.6, length.out = 10), 
           opt.ky = seq(0.1, 0.9, length.out = 10), ...) {
  
  r <- concord(x, y, ncomp = 1, ...)
  nobs <- ncol(y)
  # 
  v <- sort(r$score.y[, 1])
  m <- replicate(ceiling(nobs/fold), 1:fold)
  m <- c(apply(m, 2, sample))[1:nobs]
  #
  cv <- lapply(unique(m), function(i) {
    ii <- m == i
    xn <- lapply(r$normed$x, t)
    yn <- t(r$normed$y)
    xx <- lapply(xn, function(mat) mat[, ii])
    yy <- yn[, ii]
    xt <- lapply(xn, function(mat) mat[, !ii])
    yt <- yn[, !ii]
    cvx <- sapply(opt.kx, function(kk) {
      r0 <- concord(xx, yy, ncomp = 1, dmod = 1, option = "uniform", center = FALSE, scale = FALSE, kx = kk)
      sum((yt - predictConcordance(r0, xt))^2)
    })
    cvy <- sapply(opt.ky, function(kk) {
      r0 <- concord(xx, yy, ncomp = 1, dmod = 1, option = "uniform", center = FALSE, scale = FALSE, ky = kk)
      sum((yt - predictConcordance(r0, xt))^2)
    })
    list(cvx = cvx, cvy = cvy)
  })
  
  cvx <- sapply(cv, "[[", "cvx")
  cvy <- sapply(cv, "[[", "cvy")
  
  mx <- rowMeans(cvx)
  sdx <- rowSds(cvx)
  
  my <- rowMeans(cvy)
  sdy <- rowSds(cvy)
  
  list(x = data.frame(kx = opt.kx, mean = mx, sd = sdx), 
       y = data.frame(ky = opt.ky, mean = my, sd = sdy))
}



