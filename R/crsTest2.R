#' Burrow Clustered Random Subspaces
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @return
#' @export
#'
#' @examples

crsTest2 <- function(x, y, B1 = 100) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  clusters <- burrowClusters2(x, y)
  cols <- lapply(unique(clusters$cluster), function(i) {
    which(clusters$cluster == i)
  })
  res <- sapply(seq_along(cols), function(i) {
    clusterCols <- cols[[i]]
    xSub <- as.data.frame(x[, clusterCols])
    ySub <- as.data.frame(y[, clusterCols])
    k <- min(floor((n1 + n2 - 2)/2), ncol(xSub))
    rsTest(xSub, ySub, B1, k = k)
  })
  mean(res)
}

#' @rdname crs2Test
#' @export

crs2Pval <- function(x, y, k, B1 = 100, B2 = 100,
                   parallel = TRUE, mc.cores = 2) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  x <- as.data.frame(x)
  y <- as.data.frame(y)
  p <- ncol(x)
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  z <- rbind(x, y)
  crs2Obs <- crsTest2(x, y, B1)
  if(parallel == TRUE) {
    crs2Z <- mclapply(1:B2, function(i) {
      xRows <- sample(n1 + n2, n1)
      xNew <- z[xRows, ]
      yNew <- z[-xRows, ]
      crsTest2(xNew, yNew, B1)
    }, mc.cores = mc.cores)
    crs2Z <- unlist(crs2Z)
  } else {
    crs2Z <- sapply(1:B2, function(i) {
      xRows <- sample(n1 + n2, n1)
      xNew <- z[xRows, ]
      yNew <- z[-xRows, ]
      crsTest2(xNew, yNew, B1)
    })
  }
  sum(crs2Z >= crs2Obs) / B2
}


#' @rdname crsTest2
#' @export
burrowClusters2 <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  k <- floor((n1 + n2 - 2) / 2)
  df <- Reduce(f = rbind, list(x = x, y = y))
  p <- ncol(x)
  distances <- pearsonDistance(df)
  clusterStart <- cluster::diana(distances)
  allCuts <- stats::cutree(stats::as.hclust(clusterStart), k = 1:ncol(x))
  mins <- sapply(1:p, function(i) min(table(allCuts[, i])))
  cutSelect <- max(which(mins >= k))
  cuts <- allCuts[, cutSelect]
  clusters <- data.frame(variable = names(cuts),
                         cluster = as.character(cuts),
                         stringsAsFactors = FALSE)
  clusters
}
