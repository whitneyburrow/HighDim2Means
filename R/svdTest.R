#' SVD Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @export

svdTest <- function(x, y, method = "svd") {
  reducedData <- svdReducedPair(x, y, method = "svd")
  as.numeric(hotellingT2(reducedData$x, reducedData$y))
}



#' @rdname svdTest
#' @export
hotellingT2 <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  p <- ncol(x)
  dbar <- .Internal(colMeans(x, n1, p, na.rm = FALSE)) -
    .Internal(colMeans(y, n2, p, na.rm = FALSE))
  sPool <- ((n1 - 1) * cov(x) + (n2 - 1) * cov(y)) / n
  t(dbar) %*% solve((1 / n1 + 1 / n2) * sPool) %*% dbar
}
