#' Srivastava-Kubokawa Test
#'
#' @param x Data set 1.
#' @param y Data set 2.
#' 
#' @importFrom highD2pop SK.test
#'
#' @return
#' @export

skTest <- function(x, y) {
  highD2pop::SK.test(x, y)[[1]]
}

