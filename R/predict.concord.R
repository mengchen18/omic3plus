#' @title predict y using result of concordance analysis
#' @description predict y using the result of concordance analysis
#' @param object an object returned by concord
#' @param newdata new x
#' @return the predicted y matrix
#' @author Chen Meng
#' @export


predictConcordance <- function(object, newdata) {
  xcat <- do.call(rbind, newdata)
  preds <- crossprod(xcat, object$loading.x)
  tcrossprod(object$loading.y, preds)
}


