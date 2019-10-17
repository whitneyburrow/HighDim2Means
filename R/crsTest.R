#' Worley Clustered Random Subspaces
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param B1 Number of times to randomly subset. Defaults to 100.
#'
#' @return Clustered random subspaces statistic.
#' @export


crsTest <- function(x, y, k, B1 = 100) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  clusters <- worleyClusters(x, y)
  cols <- lapply(unique(clusters$cluster), function(i) {
    which(clusters$cluster == i)
  })

  res <- sapply(seq_along(cols), function(i) {
    clusterCols <- cols[[i]]
    xSub <- x[, clusterCols]
    ySub <- y[, clusterCols]
    knew <- min(k, ncol(xSub))
    rsTest(xSub, ySub, B1, k = knew)
  })
  sum(res)
}

#' @rdname crsTest
#' @export
crs.test <- function(x, y, k, B1 = 100, B2 = 100) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  t <- crsTest(x, y, k, B1)
  z <- rbind(x, y)
  crsZ <- sapply(1:B2, function(i) {
    xRows <- sample(n1 + n2, n1)
    xNew <- z[xRows, ]
    yNew <- z[-xRows, ]
    crsTest(xNew, yNew, k, B1)
  })
  pval <- sum(crsZ >= t) / B2
  c(tcrs = t, pvalue = pval)
}


#' @rdname crsTest
#' @export
worleyClusters <- function(x, y) {
  df <- Reduce(f = rbind, list(x = x, y = y))
  p <- ncol(x)
  n <- nrow(df) - 2
  kc <- floor(n / 3)
  distances <- pearsonDistance(df)
  clusterStart <- flashClust::hclust(distances)
  allCuts <- cutree(clusterStart, k = 1:ncol(x))
  maxes <- sapply(1:p, function(i) max(table(allCuts[, i])))
  k <- min(which(maxes < kc))
  cuts <- allCuts[, k]
  clusters <- data.frame(variable = names(cuts),
                         cluster = as.character(cuts),
                         stringsAsFactors = FALSE)
  clusters
}
