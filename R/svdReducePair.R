#' SVD Reduction for High Dimensional Data Pairs.
#' 
#' The following performs dimension reduction on two samples of
#' high dimensional data.
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param method Method used for Singular Value Decompostion.
#'
#' @export
#'

svdReducedPair <- function(x, y, method = "svd") {
  fullData <- rbind(x, y)
  if(method == "svd") {
    sing <- svd(fullData)
    rank <- min(Matrix::rankMatrix(x), Matrix::rankMatrix(y))
    v <- sing$v[, 1:rank]
    reducedFull <- as.matrix(fullData) %*% v
    reducedX <- reducedFull[1:nrow(x), ]
    reducedY <- reducedFull[(nrow(x) + 1):nrow(reducedFull), ]
    list(x = reducedX, y = reducedY)
  }
  else if(method == "scatter") {
    scatter <- scatterMat(fullData)
    svdOfScatter <- svd(scatter)
    rank <- Matrix::rankMatrix(scatter)
    u1 <- svdOfScatter$u[, 1:rank]
    reducedX <- projectData(x, u1)
    reducedY <- projectData(y, u1)
    list(x = reducedX, y = reducedY)
  } else
    print("Error. Method of choice is not valid.")
}

centerMat <- function(n) {
  identity <- diag(1, n)
  centerTerm <- 1 / n * matrix(rep(1, n * n), nrow = n)
  identity - centerTerm
}

scatterMat <- function(data) {
  n <- nrow(data)
  t(data) %*% centerMat(n) %*% data
}

projectData <- function(data, projection) {
  data %*% projection
}
