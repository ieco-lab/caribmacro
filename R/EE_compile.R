#' Compile data from different files exported by Earth Engine
#'
#' This function is used to compile data exported from Google's Earth Engine as CSV files into a single tibble.
#' The function requires the package 'readr'
#'
#' @param filepath The file path to the folder that holds the exported CSV files from Earth Engine.
#' Should be in the format produced by the file.path() function (i.e. with '/' instead of '\\')
#'
#' @param group A character equal to the column name in which the feature names are stored (e.g. island names)
#'
#' @param columns A character or vector of the column names that have the data of interest
#'
#' @return A data frame of the selected data from the files exported from Google Earth Engine
#'
#' @examples
#' \dontrun{
#'
#' library(readr)
#'
#' filepath <- file.path("G:", "Shared drives", "Caribbean Macrosystem", "data", "islands", "island_covariates", "landcover")
#'
#' dat <- EE_compile(filepath,
#'                   group = "bank",
#'                   columns = c("max", "mean", "median", "min", "sd", "lower", "upper"))
#'
#' ## OR ##
#' dat <- EE_compile("G:/Shared drives/Caribbean Macrosystem/data/islands/island_covariates/landcover",
#'                   group = "bank",
#'                   columns = c("max", "mean", "median", "min", "sd", "lower", "upper"))
#'
#' }
#' @export
EE_compile <- function(filepath, group, columns) {

  files <- list.files(filepath, all.files = FALSE)
  out <- suppressMessages(as.data.frame(readr::read_csv(paste(filepath, files[1], sep = "/"))))
  out <- out[, c(group, columns)]
  names(out) <- c(group, paste(columns, gsub(".csv", "", files[1]), sep = "_"))

  for (i in 2:length(files)) {
    tmp <- suppressMessages(as.data.frame(readr::read_csv(paste(filepath, files[i], sep = "/"))))
    tmp <- tmp[, c(group, columns)]
    names(tmp) <- c(group, paste(columns, gsub(".csv", "", files[i]), sep = "_"))

    out <- merge(out, tmp, by = group, all = TRUE)
    rm(tmp)
  }

  return(out)
}
