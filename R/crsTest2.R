#' Burrow Clustered Random Subspaces
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param B1 Number of times to randomly subset. Defaults to 100.
#'
#' @return Clustered random subspaces statistic.
#' @export


crsTest2 <- function(x, y, k, B1 = 100) {
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  clusters <- burrowClusters2(x, y)
  cols <- lapply(unique(clusters$cluster), function(i) {
    which(clusters$cluster == i)
  })
  res <- sapply(seq_along(cols), function(i) {
    clusterCols <- cols[[i]]
    xSub <- x[, clusterCols]
    ySub <- y[, clusterCols]
    k <- min(floor((n1 + n2 - 2) / 2), ncol(xSub))
    rsTest(xSub, ySub, B1, k = k)
  })
  mean(res)
}

#' @rdname crsTest
#' @export
burrowClusters2 <- function(x, y) {
  p <- ncol(x)
  n <- nrow(x) + nrow(y) - 2
  distances <- dist(t(x), method = "euclidean")
  clusterStart <- flashClust::hclust(distances, method = "complete")
  allCuts <- cutree(clusterStart, k = 1:ncol(x))
  k <- 4
  cuts <- allCuts[, k]
  cuts <- cutree(clusterStart, k = 4)
  clusters <- data.frame(variable = names(cuts),
                         cluster = as.character(cuts),
                         stringsAsFactors = FALSE)
  clusters
}

