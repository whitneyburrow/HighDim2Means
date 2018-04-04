#' Srivastava-Du Test
#'
#' Calculates Srivastava-Du test statistic for two
#' high-dimensional data sets.
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @export
#'
sdTest <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  dbar <- colMeans(x) - colMeans(y)
  sx <- cov(x)
  sy <- cov(y)
  s <- ((n1 - 1) * sx + (n2 - 1) * sy) / (n1 + n2)
  diagonals <- diag(s)
  ds <- diag(diagonals)
  rhat <- solve(ds^(1/2)) %*% s %*% solve(ds^(1/2))
  cpn <- 1 + sum(diag(rhat))^2 / p^(3/2)
  num <- (1 / n1 + 1 / n2) ^ (-1) * t(dbar) %*% solve(ds) %*% dbar - p
  den <- (2 * (sum(diag(rhat))^2 - p^2 / (n1 + n2 - 2)) * cpn) ^ (1 / 2)
  num / den
}

#' @rdname sdTest
#' @export

sdPval <- function(x, y) {
  t <- sdTest(x, y)
  1 - pnorm(t)
}
