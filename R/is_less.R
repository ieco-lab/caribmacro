#' Identify Values that are Less than Previous Value
#'
#' This function will check to see if the values in a vector is less than the one before
#' it in the vector
#'
#' @param x Vector of length >0 to be evaluated
#'
#' @return A vector logical values (i.e. \code{TRUE} or \code{FALSE}) that denote if the
#' value is less than the one before it in the vector
#'
#' @examples
#' \dontrun{
#'
#' x <- c(1, 2, 1, 4)
#' y <- c(5, 4, 6, 5)
#'
#' is_less(x)
#' # [1] FALSE FALSE TRUE FALSE
#'
#' is_less(y)
#' # [1] FALSE TRUE FALSE TRUE
#'
#' # Return values that are less than previous in vector
#' x[is_less(x)]
#' # [1] 1
#'
#' y[is_less(y)]
#' # [1] 4 5
#'
#' # Return position of lesser values
#' which(is_less(x))
#' # [1] 3
#'
#' which(is_less(y))
#' # [1] 2 4
#'
#' }
#'
#' @export

is_less <- function (x) {
  out <- NULL
  for (i in 2:length(x)) {
    out <- c(out, x[i] < x[i - 1])
  }
  return(out)
}
