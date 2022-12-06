#' Calculate the PSV for geographic features
#'
#' This function is used to calculate the phylogenetic species variability (PSV) for the geographic
#' features of interest from a list of community matrices such as produced by \code{com_matrix()}.
#' This function requires the package 'phyr'
#'
#' @param data A list of community matrices that have column with the geographic names. Each entry in
#' the list should be a different species group and that list entry should be named with the group name.
#' The resulting object from the \code{com_matrix()} is already formatted for this function.
#'
#' @param geo_name A character equal to the column name of the data frame(s) in 'data' in which the
#' geographic feature names are stored
#'
#' @param tree A phylogenetic tree that includes as many species as possible that are in the data. Should be
#' in the Newick format.
#'
#' @return A list of data frames. The first data frame holds the PSV values for each geographic feature
#' and species group. If there is less than two species of a community in  \code{tree} then the function
#' returns an \code{NA} for the PSV value for that community.
#'
#' The second data frame holds the number of taxa (e.g. species) used to calculate
#' the PSV value for each geographic feature and species group.
#'
#' The third data frame holds the 'Vars'
#' values for each geographic feature and species group.
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(reshape)
#' library(phyr)
#'
#' tetra <- read.tree(file.path(here(), 'data_raw', 'Tetrapoda_species.nwk'))
#' dat <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_final.csv'), header=TRUE)
#'
#' coms <- com_matrix(species = "binomial",
#'                    geo_group = "bank",
#'                    taxa_group = "class",
#'                    status = "bnk_status",
#'                    stat_levels = c('N', "E"),
#'                    total = TRUE,
#'                    data = dat)
#'
#' PSV <- PSV_geo(tree = tetra, data = coms, geo_group = 'bank')
#'
#' }
#'
#' @export
PSV_geo <- function(tree, data, geo_group) {

  out<-list(data.frame(geo_group=data[[1]][, geo_group]),
            data.frame(geo_group=data[[1]][, geo_group]),
            data.frame(geo_group=data[[1]][, geo_group]))
  for (i in 1:length(out)) {
    names(out[[i]])<-geo_group
  }
  names(out)<-c('PSV', 'Ntaxa', 'Vars')

  for (i in 1:length(data)) {
    if (length(intersect(names(data[[i]]), tree$tip.label)) > 1) {
      com <- as.matrix(as.data.frame(data[[i]][, 2:ncol(data[[i]])]))
      row.names(com) <- data[[i]][, geo_group]
      tmp <- phyr::psv(com,tree)
      names(tmp) <- paste(names(tmp), names(data)[i], sep=".")
      tmp[, geo_group] <- row.names(tmp)
      out[['PSV']] <- merge(out[['PSV']], tmp[, c(geo_group, names(tmp)[1])], by = geo_group, all=TRUE)
      out[['Ntaxa']] <- merge(out[['Ntaxa']], tmp[, c(geo_group, names(tmp)[2])], by = geo_group, all=TRUE)
      out[['Vars']] <- merge(out[['Vars']], tmp[, c(geo_group, names(tmp)[3])], by = geo_group, all=TRUE)

      rm(tmp, com)
    } else {
      tmp <- data.frame(geo_group = data[[i]][, geo_group],
                        PSVs = rep(NA, length(data[[i]][, geo_group])),
                        SR = rep(NA, length(data[[i]][, geo_group])),
                        vars = rep(NA, length(data[[i]][, geo_group])))

      names(tmp) <- c(geo_group, paste(names(tmp)[2:4], names(data)[i], sep="."))

      out[['PSV']] <- merge(out[['PSV']], tmp[, c(geo_group, names(tmp)[2])], by = geo_group, all=TRUE)
      out[['Ntaxa']] <- merge(out[['Ntaxa']], tmp[, c(geo_group, names(tmp)[3])], by = geo_group, all=TRUE)
      out[['Vars']] <- merge(out[['Vars']], tmp[, c(geo_group, names(tmp)[4])], by = geo_group, all=TRUE)
    }
  }
  names(out[['PSV']]) <- sub('PSVs.', '', names(out[['PSV']]))
  names(out[['Ntaxa']]) <- sub('SR.', '', names(out[['Ntaxa']]))
  names(out[['Vars']]) <- sub('vars.', '', names(out[['PSV']]))

  return(out)
}
