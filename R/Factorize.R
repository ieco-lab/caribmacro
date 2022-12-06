#' Turn Categorical Columns of a Data Frame to Factors
#'
#' This function is used to set columns of a data frame of class character or numerical to
#' class factor
#'
#' @param data A data frame
#'
#' @param columns A character or vector equal to the names of the column(s) of the variables to
#' be set as factors
#'
#' @return A the same data frame as \code{data} but with the class of the \code{columns} column(s)
#' set to factor
#'
#' @examples
#'
#' df <- data.frame(X = c(A, B, C),
#'                  Y = c(a, b, c))
#'
#' df <- Factorize(data = df, columns = 'X')
#'
#' class(df$X)
#' class(df$Y)
#'
#' @export


Factorize <- function(data, columns) {
  for (i in 1:length(columns)) {
    data[, columns[i]] <- as.factor(data[, columns[i]])
  }
  return(data)
}
