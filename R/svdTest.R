#' SVD Test
#' 
#' The following performs a Hotelling T2 test on a high dimensional
#' dataset that has first been reduced via the singular value
#' decomposition.
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' @param method Method for reduction. Can use traditional svd or svd of scatter matrix.
#'
#' @export

svdTest <- function(x, y, method = "svd") {
  reducedData <- svdReducedPair(x, y, method = "svd")
  as.numeric(hotellingT2(reducedData$x, reducedData$y))
}
