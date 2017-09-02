

predict.concordance <- function(res, x) {
  Xnorm <- lapply(x, function(x) scale(t(x), center = TRUE, scale = TRUE))
  
  
  Ynorm <- scale(t(y), center = center, scale = scale)
  
  val <- switch(option, 
                "uniform" = 1, 
                "lambda1" = svd(Ynorm)$d[1],
                "inertia" = sum(Ynorm^2))
  
  Xnorm <- processOpt(x, center = center, scale = scale, option = option, value = val)
  
  res <- r1
  
  xcat <- do.call(rbind, x)
  preds <- crossprod(xcat, res$loading.x)
  l <- tcrossprod(res$loading.y , preds)
}

