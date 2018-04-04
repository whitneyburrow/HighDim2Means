#' Srivastava T+2 Test
#'
#' Calculates Srivastava test statstic base on Moore-Penrose
#' inverse for two high dimensional data sets.
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @export

smpTest <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  dbar <- colMeans(x) - colMeans(y)
  sx <- cov(x)
  sy <- cov(y)
  s <- ((n1 - 1) * sx + (n2 - 1) * sy) / (n1 + n2)
  splus <- MASS::ginv(s)
  (1 / n1 + 1 / n2) ^ (-1) * t(dbar) %*% splus %*% dbar
}

#' @rdname smpTest
#' @export

smpPval <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  p <- ncol(x)
  tplus <- smpTest(x, y)
  sx <- cov(x)
  sy <- cov(y)
  s <- ((n1 - 1) * sx + (n2 - 1) * sy) / (n1 + n2)
  ahat1 <- sum(diag(s)) / p
  ahat2 <- n^2 / ((n - 1) * (n + 2)) * 
    (sum(diag(s %*% s)) / p - p / n * (sum(diag(s)) / p)^2)
  bhat <- ahat1^2 / ahat2
  tplusS <- ((bhat * p / n) * tplus - n) / sqrt(2 * n)
  1 - pnorm(tplusS)
}