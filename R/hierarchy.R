#' Determine species hierarchy
#'
#' This function is used to retrieve the taxonomic hierarchy from an online database for each species
#' in a data frame.
#' This function requires the package 'taxize'
#'
#'
#' @param data A data frame
#'
#' @param species A character equal to the column name of the data frame that holds the species names
#'
#' @param db A character equal to the data base name for which you want to use to search for the taxonomic
#' hierarchy. The possible data bases are Integrated Taxonomic Information System (\code{'itis'}),
#' Global Biodiversity Information Facility (\code{'gbif'}), The National Center for Biotechnology
#' Information (\code{'ncbi'}), the Encyclopedia of Life (\code{'eol'}), National Biodiversity
#' Network (\code{'nbn'}), and Tropicos (\code{'tropicos'}). The default is \code{'gbif'}. The function
#' uses the first ID number returned by the \code{taxize::get_id} function.
#'
#' @param ranks A character or vector of rank names. The default is 'kingdom', 'phylum', 'class', 'order',
#' 'family', 'genus', and 'species'.
#'
#' @param error.ret Logical. If \code{TRUE} then a column is added to the returned data frame that indicates
#' if there was a some error in the retrieval of that species' hierarchy. \code{error.ret = TRUE} checks to see
#' if the species name matches the species name returned in the hierarchy from the database.
#'
#' @return A data frame that has a column for each taxonomic level of a species taxonomic hierarchy
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(taxize)
#'
#' herp <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_v6.csv'), header=TRUE)
#'
#' taxonomy <- hierarchy(data = herp, species = 'binomial', db = 'gbif')
#'
#' }
#'
#' @export
hierarchy <- function(data, species,
                      db = 'gbif',
                      ranks = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
                      ID = TRUE, error.ret = TRUE) {
  sp <- levels(as.factor(data[, species]))
  out <- data.frame(NULL)

  ranks <- tolower(ranks)

  for (i in 1:length(sp)) {
    if (db == 'gbif') {
      id <- c(taxize::get_gbifid_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])
    } else if (db == 'ncbi') {
      #id <- c(taxize::get_uid_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])
      stop('ID Function needs API for NCBI')
    } else if (db == 'itis') {
      id <- as.data.frame(taxize::get_tsn_(sp[i], messages = FALSE, rows = 1)[[1]][1, 1])[1, 1]
    } else if (db == 'eol') {
      id <- c(taxize::get_eolid_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])
    } else if (db == 'col') {
      stop(paste(db,"is currently not supported by taxize"))
    } else if (db == 'nbn') {
      id <- c(taxize::get_nbnid_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])
    } else if (db == 'tropicos') {
      #id <- c(taxize::get_tpsid_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])
      stop('ID Function not working with Tropicos')
    } else if (db == 'worms') {
      id <- c(taxize::get_wormsid_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])
    } else if (db == 'natserv') {
      id <- taxize::get_natservid_(sp[i], messages = FALSE, rows = 1)
      if (nrow(id[[1]]) > 0) {
        id <- c(id[[1]][1,1]$id)
      } else {
        id <- NULL
      }
    } else if (db == 'bold') {
      id <- c(taxize::get_boldid_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])
    } else if (db == 'wiki') {
      id <- as.data.frame(taxize::get_wiki_(sp[i], messages = FALSE, rows = 1)[[1]][1,1])[1,1]
    } else if (db == 'pow') {
      id <- taxize::get_pow_(sp[i], messages = FALSE, rows = 1)[[1]]$fqId
    } else {
      stop(paste(db,"is not an acceptable database name"))
    }

    if (is.null(id)) {
      id <- NA
    } else if (is.na(id)) {
      id <- NA
    }

    if (length(id) == 1 & !is.na(id)) {
      tmp <- taxize::classification(id, db = db)

      if (all(is.element(ranks, tmp[[1]]$rank))) {
        tmp2 <- data.frame(sp_in_data = sp[i])
        for (j in 1:length(ranks)) {
          tmp2[, ranks[j]] <- tmp[[1]][which(tmp[[1]]$rank == ranks[j]), 'name']
        }

      } else if (any(is.element(ranks, tmp[[1]]$rank))) {
        tmp2 <- data.frame(sp_in_data = sp[i])
        new.rnk <- intersect(ranks, tmp[[1]]$rank)
        for (j in 1:length(new.rnk)) {
          tmp2[, new.rnk[j]] <- tmp[[1]][which(tmp[[1]]$rank == new.rnk[j]), 'name']
        }
        bad.rnk <- setdiff(ranks, tmp[[1]]$rank)
        for (j in 1:length(bad.rnk)) {
          tmp2[, bad.rnk[j]] <- NA
        }
      }
    } else {
      tmp2 <- data.frame(sp_in_data = sp[i])

      for (j in 1:length(ranks)) {
        tmp2[, ranks[j]] <- NA
      }
    }

    tmp2 <- tmp2[, c('sp_in_data', ranks)]

    if (ID) {
      tmp2$ID <- id
      tmp2 <- tmp2[, c('sp_in_data', 'ID', ranks)]
    } else {NA}

    out <- rbind(out, tmp2)
  }


  if (error.ret) {
    out$Error<-'NO'

    for (i in as.numeric(rownames(out[!is.na(out$species), ]))) {
      if (out[i, 'sp_in_data'] != out[i, ranks[length(ranks)]]) {
        out[i, 'Error'] <- 'ERROR'
      } else {NA}
    }

    #MESSAGES DO NOT WORK ... NEED TO FIX
    #message('Could not find species: ')
    #mess<-levels(droplevels(out[is.na(out$Species), 'sp_in_data']))
    #for (i in 1:length(mess)) {
    #  message(mess[i])
    #}

    #message('')
    #message('--------------------------------------------')
    #message('')

    #message('Species names that do not match: ')
    #mess<-levels(droplevels(out[which(out$Error=='ERROR'), 'sp_in_data']))
    #for (i in 1:length(mess)) {
    #  message(mess[i])
    #}

    out[is.na(out$species), 'Error'] <- 'ERROR'

  } else {NA}

  return(out)
}
