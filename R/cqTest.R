#' Chen-Qin Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#'
#' @return
#' @export

cqTest <- function(x, y) {
  highD2pop::ChenQin.test(as.matrix(x), as.matrix(y))[[1]]
}

