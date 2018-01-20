#' Shifted Mean Vector for Data Generation
#'
#' This function shifts a given mean vector. The function shifts 4/5 of each
#' block in m out of the total number of blocks. (Thulin)
#'
#' @param mean Original mean vector.
#' @param shift Amount to shift by.
#' @param m Number of blocks to experience shift.
#'
#' @return
#' @export
#'
#' @examples
shiftedMean <- function(mean, shift, m) {
  browser()
  p <- length(mean)
  shiftVec <- c(rep(shift, 20), rep(0, 5))
  shiftVecXm <- rep(shiftVec, m)
  totalShiftVec <- c(shiftVecXm, rep(0, p - length(shiftVecXm)))
  mean + totalShiftVec
}
