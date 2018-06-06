#' SVD Reduction for High Dimensional Data Pairs.
#' 
#' The following performs dimension reduction on two samples of
#' high dimensional data.
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @export
#'

svdReducedPair2 <- function(x, y) {
  fullData <- rbind(x, y)
  sing <- svd(fullData)
  rank <- Matrix::rankMatrix(fullData)
  u <- sing$u[, 1:rank]
  d <- diag(sing$d[1:rank])
  v <- sing$v[, 1:rank]
  reducedFull <- u %*% d %*% v
  reducedX <- reducedFull[1:nrow(x), ]
  reducedY <- reducedFull[(nrow(x) + 1):nrow(reducedFull), ]
  list(x = reducedX, y = reducedY)
}
