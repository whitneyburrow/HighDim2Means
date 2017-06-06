#' Burrow Clustered Random Subspaces
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param B1 Number of times to randomly subset. Defaults to 1000.
#'
#' @return Clustered random subspaces statistic.
#' @export


crsTest2 <- function(x, y, k, B1 = 1000) {
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  p <- ncol(x)
  kc <- floor((n1 + n2 - 2) / 2)
  clusters <- burrowClusters2(x, y)
  cols <- lapply(unique(clusters$cluster), function(i) {
    which(clusters$cluster == i)
  })
  res <- sapply(seq_along(cols), function(i) {
    clusterCols <- cols[[i]]
    xSub <- x[, clusterCols]
    ySub <- y[, clusterCols]
    if(length(cols[[i]]) > kc) {
      rsTest(xSub, ySub, k = k, B1 = B1)
    } else {
      hotellingT2(xSub, ySub)
    }
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
  dc <- pearson(nrow(df) - 2, ncol(df))
  distances <- pearsonDistance(df)
  clusterStart <- flashClust::hclust(distances, method = "complete")
  cuts <- cutree(clusterStart, h = dc)
  clusters <- data.frame(variable = names(cuts),
                         cluster = as.character(cuts),
                         stringsAsFactors = FALSE)
  clusters
}

