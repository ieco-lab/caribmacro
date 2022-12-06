#' Standardize Columns of a Data Frame
#'
#' This function is used to standardize data by mean centering and dividing by the standard deviation.
#'
#' @param vars A character or vector equal to the names of the column(s) of the variables to be standardized
#'
#' @param data A data frame that holds the variables to be standardized
#'
#' @param center A character equal to \code{'mean'}, \code{'min'}, \code{'max'} or a
#' conditional (e.g. \code{'Year == 2015'} or \code{'Sample == "D-4"'})
#' that sets the value on which the data should be centered. Default is \code{center = 'mean'}
#'
#' @return A the same data frame as \code{data} but with the data in the \code{vars} column(s) standardized
#'
#' @examples
#' \dontrun{
#'
#' set.seed(1)
#' df <- data.frame(X = rnorm(10, mean = 3, sd = 2),
#'                  Y = rnorm(10, mean = 5, sd = 0.5))
#'
#' s.df <- stndrd(vars = c('X', 'Y'), data = df, center = 'mean')
#'
#' }
#'
#' @export

stndrd <- function(vars, data, center = 'mean') {
  for (i in 1:length(vars)) {
    if (center == 'mean') {
      data[, vars[i]] <- (data[, vars[i]] - mean(data[, vars[i]], na.rm = TRUE)) / sd(data[, vars[i]], na.rm = TRUE)
    } else if (center == 'min') {
      data[, vars[i]] <- (data[, vars[i]] - min(data[, vars[i]], na.rm = TRUE)) / sd(data[, vars[i]], na.rm = TRUE)
    } else if(center == 'max') {
      data[, vars[i]] <- (data[, vars[i]] - max(data[, vars[i]], na.rm = TRUE)) / sd(data[, vars[i]], na.rm = TRUE)
    } else if (grepl('==', center)) {
      center <- paste0('data$', center)
      data[, vars[i]] <- (data[, vars[i]] - data[which(eval(parse(text = center))), vars[i]]) / sd(data[, vars[i]],
                                                                                                   na.rm = TRUE)
    }
  }
  return(data)
}
