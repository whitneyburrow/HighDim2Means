#' Zhang and Pan Cluster Subspaces Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param k k Number of dimensions to select. Defaults to floor((n1 + n2 - 2) / 2).
#' @param B1 Number of times to randomly subset full data. Defaults to 100.
#' @param B2 Number of permutations for calculating pvalue. Defaults to 100.
#' @param data Data frame on which to calculate Pearson distances between variables.
#' @param n Sample size.
#' @param p Dimension.
#'
#' @return Zhang and Pan cluster subspaces statistic.
#' @export
#'

csTest <- function(x, y, k, B1 = 100) {
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  clusters <- pearsonClusters(x, y)
  cols <- lapply(unique(clusters$cluster), function(i) {
    which(clusters$cluster == i)
  })
  res <- sapply(seq_along(cols), function(i) {
    clusterCols <- cols[[i]]
    xSub <- as.data.frame(x[, clusterCols])
    ySub <- as.data.frame(y[, clusterCols])
    hotellingT2(xSub, ySub)
  })
  sum(res)
}

#' @rdname csTest
#' @export

csPval <- function(x, y, k, B1 = 100, B2 = 100) {
  sum(rsZ >= rsObs) / B2
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2 - 2
  x <- as.matrix(x)
  y <- as.matrix(y)
  p <- ncol(x)
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  z <- rbind(x, y)
  csObs <- csTest(x, y, k, B1)
  csZ <- sapply(1:B2, function(i) {
      xRows <- sample(n1 + n2, n1)
      xNew <- z[xRows, ]
      yNew <- z[-xRows, ]
      csTest(xNew, yNew, k, B1)
  })
  sum(csZ >= csObs) / B2
}

#' @rdname csTest
#' @export

cs.test <- function(x, y, k, B1 = 100, B2 = 100) {
  if(missing(k)) k <- floor((n1 + n2 - 2) / 2)
  t <- csTest(x, y, k)
  pval <- csPval(x, y, k)
  c(tcs = t, pvalue = pval)
}

#' @rdname csTest
#' @importFrom stats cov
#' @export
hotellingCluster <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  dbar <- colMeans(x) - colMeans(y)
  sPool <- ((n1 - 1) * stats::cov(x) + (n2 - 1) * stats::cov(y)) / (n1 + n2 - 2)
  weight <- (n1 * n2) / (n1 + n2)
  as.numeric(weight %*% t(dbar) %*%
               solve((1 / n1 + 1 / n2) * sPool) %*% dbar)
}

#' @rdname csTest
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @export
pearsonClusters <- function(x, y) {
  df <- list(x = x, y = y)
  df <- Reduce(f = rbind, x = df)
  distances <- pearsonDistance(df)
  dc <- pearson(nrow(df) - 2, ncol(df))
  kc <- floor(2 * (nrow(df) - 2)/3)
  clusterStart <- flashClust::hclust(distances, method = "average")
  cuts <- stats::cutree(clusterStart, h = dc)

  clusters <- data.frame(variable = names(cuts),
                         cluster = as.character(cuts),
                         stringsAsFactors = FALSE)

  while(max(table(clusters$cluster)) > kc) {
    toSub <- names(which(table(clusters$cluster) > kc))
    for(i in seq_along(toSub)) {
      vars <- which(clusters$cluster == toSub[[i]])
      varNames <- clusters$variable[which(clusters$cluster == toSub[[i]])]
      newDF <- df[, vars]
      newDist <- pearsonDistance(newDF)
      subclusterStart <- flashClust::hclust(newDist, method = "average")
      cuts <- stats::cutree(subclusterStart, k = 2)
      cuts <- paste0(toSub[[i]], ".", cuts)
      clusters$cluster[vars] <- cuts
    }
  }

  return(clusters)
}

#' @rdname csTest
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @export
pearsonDistance <- function(data) {
  p <- ncol(data)
  data <- as.matrix(data)
  colPairs <- utils::combn(1:ncol(data), 2)
  distances <- sapply(1:ncol(colPairs), function(i) {
    col1 <- colPairs[1, i]
    col2 <- colPairs[2, i]
    1 - stats::cor(data[, col1], data[, col2])
  })
  results <- as.data.frame(t(rbind(colPairs, distances)))
  colnames(results) <- c("col1", "col2", "distance")

  distMat <- matrix(NA, nrow = p, ncol = p)
  rownames(distMat) <- colnames(distMat) <- paste0("V", 1:p)
  distMat[lower.tri(distMat)] <- results$distance
  return(as.dist(distMat))
}

#' @rdname csTest
#' @export
pearson <- function(n, p) {
  tc <- stats::pnorm(1 - 2 / (p * (p - 1)))
  zprimec <- tc / sqrt(n - 1)
  rc <- (exp(2 * zprimec) - 1) / (exp(2 * zprimec) + 1)
  dc <- 1 - rc
  dc
}
