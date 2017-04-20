#' Correlated Block Covariance Matrix
#'
#' @param p Number of dimensions.
#' @param inblock Covariance value inside of block.
#' @param outblock Covariance value outside of block.
#' @param numBlocks Total number of blocks. Defaults to p/25.
#'
#' @export

corrCovariance <- function(p, inblock, outblock, numBlocks) {
  if(missing(numBlocks)) numBlocks <- p / 25
  singleBlock <- matrix(inblock, p / numBlocks, p / numBlocks)
  diag(singleBlock) <- 1
  blockList <- lapply(1:numBlocks, function(i) singleBlock)
  diagonalMat <- Matrix::bdiag(blockList)
  diagonalMat <- as.matrix(diagonalMat)
  diagonalMat[which(diagonalMat == 0)] <- outblock
  diagonalMat
}
