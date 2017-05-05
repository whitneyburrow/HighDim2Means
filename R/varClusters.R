#' Clustering Variables
#'
#' @param data Data set.
#'
#' @return Clustered variables.
#' @export


varClusters <- function(data) {
  p <- ncol(data)
  n <- nrow(data)
  kc <- floor(2 * (n - 2) / 3)
  tree <- ClustOfVar::hclustvar(data)
  index <- 2
  cuts <- ClustOfVar::cutreevar(tree, k = index)
  while(min(table(cuts$cluster)) > kc) {
    index <- index + 1
    cuts <- ClustOfVar::cutreevar(tree, k = index)
  }
  cuts
}