#' Burrow Clustered Random Subspaces
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param B1 Number of times to randomly subset. Defaults to 100.
#'
#' @return Clustered random subspaces statistic.
#' @export


crsTest2 <- function(x, y, B1 = 100) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  p <- ncol(x)
  k <- floor(n / 2)
  clusters <- burrowClusters(x, y)
  cols <- lapply(unique(clusters$cluster), function(i) {
    which(clusters$cluster == i)
  })
  
  res <- sapply(seq_along(cols), function(i) {
    clusterCols <- cols[[i]]
    xSub <- x[, clusterCols]
    ySub <- y[, clusterCols]
    if(ncol(xSub) < k) {
      t <- hotellingT2(xSub, ySub)
    } else {
      t <- rsTest(xSub, ySub, B1, k = k)
    }
    t
  })
  mean(res)
}


#' @rdname crsTest
#' @export
burrowClusters2 <- function(x, y) {
  df <- Reduce(f = rbind, list(x = x, y = y))
  p <- ncol(x)
  n <- nrow(df) - 2
  kc <- floor(n / 2)
  distances <- pearsonDistance(df)
  clusterStart <- hclust(distances)
  allCuts <- cutree(clusterStart, k = 1:ncol(x))
  means <- sapply(1:p, function(i) mean(table(allCuts[, i])))
  k <- max(which(means > kc))
  cuts <- allCuts[, k]
  clusters <- data.frame(variable = names(cuts),
                         cluster = as.character(cuts),
                         stringsAsFactors = FALSE)
  clusters
}
