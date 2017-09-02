

predict.concordance <- function(res, x) {
  Xnorm <- lapply(x, function(x) scale(t(x), center = TRUE, scale = TRUE))
  xcat <- do.call(rbind, x)
  preds <- crossprod(xcat, res$loading.x^2)
  tcrossprod(res$loading.y^2 , preds)
}
