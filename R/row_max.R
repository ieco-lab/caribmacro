#' Find the Maximum Value of a Row
#'
#' This function will determine the maximum value in a each row of a data frame
#'
#' @param data A data frame
#'
#' @param na,rm Logical. If \code{TRUE} the function will ignore \code{NA} values. If \code{FLASE}
#' then the function will return a \code{NA}.
#'
#' @return A vector of the maximum values of each row of \code{data}. If a row has no values (i.e.
#' all values are \code{NA} then the function will return a \code{NA}.
#'
#' @examples
#' \dontrun{
#'
#' set.seed(1)
#' df <- data.frame(X = rnorm(10, mean = 3, sd = 2),
#'                  Y = c(NA, rnorm(9, mean = 5, sd = 0.5)))
#'
#' row_max(df, na.rm = TRUE)
#'
#' row_max(df, na.rm = FALSE)
#'
#' }
#'
#' @export

row_max <- function(data, na.rm = FALSE) {
  out <- NULL
  for (i in 1:nrow(data)) {
    if (any(!is.na(data[i, ]))) {
      out <- c(out, max(data[i, ], na.rm = na.rm))
    } else {
      out <- c(out, NA)
    }
  }
  return(out)
}

