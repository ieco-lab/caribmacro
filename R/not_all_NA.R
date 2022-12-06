#' Evaluate Vector Emptiness
#'
#' This function will check to see if a vector has at least 1 non- \code{NA} elements
#'
#' @param x Vector of length >0 to be evaluated
#'
#' @return A logical value (i.e. \code{TRUE} or \code{FALSE}) that answers the question
#' that any element of vector \code{x} is not \code{NA}
#'
#' @examples
#' \dontrun{
#'
#' x <- c(NA, NA, 1, NA)
#' y <- c(NA, NA, NA, NA)
#'
#' not_all_NA(x)
#' # TRUE
#'
#' not_all_NA(y)
#' # FALSE
#'
#' }
#'
#' @export
not_all_NA <- function (x) {
  any(!is.na(x))
}
