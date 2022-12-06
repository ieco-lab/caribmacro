#' Calculate SR for geographic features
#'
#' This function is used to calculate the species richness (SR) for the geographic
#' features of interest from a list of community matrices.
#'
#' @param data A list of community matrices that have column with the geographic names. Each entry in
#' the list should be a different species group and that list entry should be named with the group name.
#' The resulting object from the \code{com_matrix()} is already formatted for this function.
#'
#' @param geo_name A character equal to the column name of the data frame(s) in 'data' in which the
#' geographic feature names are stored
#'
#' @return A data frame that holds the SR values for each geographic feature and species group.
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(reshape)
#'
#' dat <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_v6.csv'), header=TRUE)
#'
#' coms <- com_matrix(species = "binomial",
#'                    geo_group = "bank",
#'                    taxa_group = "class",
#'                    status = "bnk_status",
#'                    stat_levels = c('N', "E"),
#'                    total = TRUE,
#'                    data = dat)
#'
#' SR <- SR_geo(data = coms, geo_group = "bank")
#'
#' }
#'
#' @export
SR_geo <- function(data, geo_group) {
  out <- data.frame(geo_group = data[[1]][, geo_group])
  names(out)[1] <- geo_group

  for (i in 1:length(data)) {
    if (length(2:ncol(data[[i]]))>1) {
      tmp <- data.frame(geo_group = data[[i]][, geo_group],
                        SR = rowSums(data[[i]][, 2:ncol(data[[i]])], na.rm=TRUE))
      names(tmp)[1] <- geo_group
      names(tmp)[2] <- names(data)[i]

      out <- merge(out, tmp, by = geo_group, all=TRUE)
      rm(tmp)
    } else {
      tmp <- data.frame(geo_group = data[[i]][, geo_group],
                        SR = data[[i]][, 2])
      names(tmp)[1] <- geo_group
      names(tmp)[2] <- names(data)[i]

      out <- merge(out, tmp, by = geo_group, all=TRUE)
      rm(tmp)
    }
  }

  for (i in 2:ncol(out)) {
    out[is.na(out[, i]), i] <- 0
  }

  out
}
