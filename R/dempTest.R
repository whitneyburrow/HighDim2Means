#' Dempster Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' 
#' @importFrom highD2pop ChenQin.test
#'
#' @return Dempster statistic for different mean vectors.
#' @export

dempTest <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  dbar <- colMeans(x) - colMeans(y)
  sx <- cov(x)
  sy <- cov(y)
  s <- ((n1 - 1) * sx + (n2 - 1) * sy) / (n1 + n2)
  trace <- sum(diag(s))
  (1 / n1 + 1 / n2) ^ (-1) * t(dbar) * dbar / trace
}