#' Thulin's Random Subspaces
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param k Number of dimensions to select. Defaults to floor((n1 + n2 - 2) / 2).
#' @param B1 Number of times to randomly subset full data. Defaults to 100.
#' @param B2 Number of permutations for calculating pvalue. Defaults to 100.
#'
#' @export

rsTest <- function(x, y, k, B1 = 100) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  p <- ncol(x)
  if(missing(k)) k <- floor(n / 2)
  rs <- sapply(1:B1, function(i) {
    cols <- sample(p, k)
    xSub <- x[, cols]
    ySub <- y[, cols]
    hotellingT2(xSub, ySub)
  })
  mean(rs)
}

#' @rdname rsTest
#' @export

rsPval <- function(x, y, k, B1 = 100, B2 = 100) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  p <- ncol(x)
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  z <- rbind(x, y)
  rsObs <- rsTest(x, y, k, B1)
  rsZ <- sapply(1:B2, function(i) {
      xRows <- sample(n1 + n2, n1)
      xNew <- z[xRows, ]
      yNew <- z[-xRows, ]
      rsTest(xNew, yNew, k, B1)
  })
  sum(rsZ >= rsObs) / B2
}

#' @rdname rsTest
#' @export

rs.test <- function(x, y, k, B1 = 100, B2 = 100) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  t <- rsTest(x, y, k) 
  pval <- rsPval(x, y, k) 
  c(trs = t, pvalue = pval)
}