#' Retrieve COL and other database IDs
#'
#' This function is used to retrieve the Catelogue of Life id and any additional database ids for each species
#' in a data frame.
#' This function requires the package 'taxize'
#'
#' NOTE: Catelogue of Life is not supported by 'taxize' anymore... check back later for its incorporation
#'
#' @param data A data frame
#'
#' @param species A character equal to the column name of the data frame that holds the species names
#'
#' @param db Optional. A character or vector eaqual to the data base names for which you want the ids in addition
#' to Catelogue of Life (COL). By default only COL ids are returned. The possible data bases
#' are Integrated Taxonomic Information System ('itis'), Global Biodiversity Information
#' Facility ('gbif'), and The National Center for Biotechnology Information ('ncbi')
#'
#' @return A data frame with the COL IDs and the IDs of any other biodiversity database give in the \code{db}
#' parameter
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(taxize)
#'
#' herp <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_final.csv'), header=TRUE)
#'
#' ids <- div_db_ids(data = herp,
#'                   species = 'binomial',
#'                   db = c('itis', 'gbif', 'ncbi'))
#'
#' }
#'
#' @export
div_db_ids <- function(data, species, db = NULL) {

  stop('Catelogue of Life currently not supported:\n      div_db_ids function deprecated')

  sp <- levels(as.factor(data[, species]))
  tmp <- taxize::get_ids(names = sp[i], db = 'col')
  out <- data.frame(Species = names(c(tmp[[1]])), COL.id = c(tmp[[1]]))

  for (i in 2:length(sp)) {
    tmp <- taxize::get_ids(names = sp[i], db = 'col')
    tmp <- data.frame(Species = names(c(tmp[[1]])), COL.id = c(tmp[[1]]))
    out <- rbind(out, tmp)
    Sys.sleep(15)
  }

  if (length(db) == 1) {
    out$Species<-as.character(out$Species)
    out$COL_id<-as.character(out$COL_id)

    col <- paste(toupper(db), 'id', sep = "_")
    out[, col] <- NA

    for(i in 1:nrow(out)) {
      tmp <- taxize::get_ids(out[i, 'COL.id'], db = db)

      out[i, col]<-tmp[[1]][[1]]
    }
  } else if (length(db) == 2) {
    out$Species<-as.character(out$Species)
    out$COL_id<-as.character(out$COL_id)

    col <- paste(toupper(db), 'id', sep = "_")
    out[, col[1]] <- NA
    out[, col[2]] <- NA

    for(i in 1:nrow(out)) {
      tmp <- taxize::get_ids(out[i, 'COL.id'], db = db)

      out[i, col[1]]<-tmp[[1]][[1]]
      out[i, col[2]]<-tmp[[2]][[1]]
    }
  } else if (length(db) == 3) {
    out$Species<-as.character(out$Species)
    out$COL_id<-as.character(out$COL_id)

    col <- paste(toupper(db), 'id', sep = "_")
    out[, col[1]] <- NA
    out[, col[2]] <- NA
    out[, col[3]] <- NA

    for(i in 1:nrow(out)) {
      tmp <- taxize::get_ids(out[i, 'COL.id'], db = db)

      out[i, col[1]] <- tmp[[1]][[1]]
      out[i, col[2]] <- tmp[[2]][[1]]
      out[i, col[3]] <- tmp[[3]][[1]]
    }
  } else {NA}

  return(out)
}
