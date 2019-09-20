#' Bai Sarandasa Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @return Bai Sarandasa statistic for different mean vectors.
#' @export

bsTest <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  dbar <- colMeans(x) - colMeans(y)
  sx <- cov(x)
  sy <- cov(y)
  s <- ((n1 - 1) * sx + (n2 - 1) * sy) / (n1 + n2)
  trace <- sum(diag(s))
  trace2 <- sum(diag(s %*% s))
  num <- (1 / n1 + 1 / n2) ^ (-1) * t(dbar) %*% dbar - trace
  den <- sqrt(2 * (trace2 - 1/(n1 + n2 - 2) * trace^2))
  as.numeric(num / den)
}

#' @rdname bsTest
#' @export

bsPval <- function(x, y, B2 = 500) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  tobs <- bsTest(x, y)
  1 - pnorm(tobs)
}

#' @rdname bsTest
#' @export

bs.test <- function(x, y) {
  t <- bsTest(x, y)
  pval <- 1 - pnorm(t)
  c(tbs = t, pvalue = pval)
}