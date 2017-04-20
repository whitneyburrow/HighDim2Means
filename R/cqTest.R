#' Chen-Qin Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' 
#' @importFrom highD2pop ChenQin.test
#'
#' @return Chen-Qin statistic for different mean vectors.
#' @export

cqTest <- function(x, y) {
  highD2pop::ChenQin.test(as.matrix(x), as.matrix(y))[[1]]
}

