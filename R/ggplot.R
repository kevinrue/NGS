
#' Oscillate Y nudge above and below data points
#'
#' @param y Amount of Y nudge
#' @param length Count of data points to nudge
#'
#' @return Vector of Y nudge alternated and balanced across 0.
#' @export
#'
#' @examples
nudge_balanced <- function(y, length){
  direction <- rep_len(c(-1,1), length.out=length)
  direction * y
}
