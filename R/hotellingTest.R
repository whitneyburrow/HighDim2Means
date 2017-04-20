#' Performs Hotelling T2 Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @return
#' @export
hotellingT2 <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  p <- ncol(x)
  dbar <- colMeans(x) - colMeans(y)
  sPool <- ((n1 - 1) * cov(x) + (n2 - 1) * cov(y)) / n
  t(dbar) %*% solve((1 / n1 + 1 / n2) * sPool) %*% dbar
}
