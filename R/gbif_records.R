#' Download Records of Species from GBIF
#'
#' This function is used to download occurrence records from GBIF within a year range for
#' given taxonomic group(s) for a given area.
#' The function requires the packages 'rgbif', 'rgdal', and 'rgeos'
#'
#' @param groups A character or vector of taxonomic group names of the same taxonomic rank
#'
#' @param rank A character equal to the taxonomic level of the group(s) (e.g. 'class')
#'
#' @param bounds A shapefile of a simplified polygon that bounds the study area
#' (e.g. outline of the Caribbean bioregion)
#'
#' @param years A number equal to the number of years before the present year (default) or
#' before \code{yr_from} for which you want occurrences collected. See Details for more information.
#'
#' @param yr_from A number equal to the year of the most current records you want to download. The
#' default is the current year as determined by \code{Sys.Date()}. See Details for more information.
#'
#' @param limit A numerical value equal to the maximum number of records wanted. See Details for more
#' information.
#'
#' @param genus_level Logical. If \code{TRUE} then the records that were only identified to
#' the genus level are returned. If \code{FALSE}, then these records are omitted
#'
#' @return A data frame with the occurrence observations of species within the taxonomic group(s).
#' The function only returns the columns: "order", "family", "genus", "species", "eventDate", "year",
#' "decimalLatitude", "decimalLongitude", "basisOfRecord", "occurrenceID", and "references" from the
#' \code{occ_search()} function of the 'rgbif' package
#'
#' @details
#' In order to download records for a certain year range, you have to specify both \code{years}
#' and \code{yr_from}. The \code{years} argument tells the function to download all records from a
#' specified number of years prior to \code{yr_from}. For example, if you want records from 2016 through
#' 2018 then you would set \code{years = 2} and \code{yr_from = 2018}. However, the default for \code{yr_from}
#' is the year of the computer system (i.e. for most the current year). Therefore, if you only want records
#' from two years prior through the current year, you only need to set \code{years = 2} and nothing for \code{years}.
#'
#' The default number of records returned is 50,000, and therefore, if you want all records within a certain time
#' frame, you should check that the number of records returned does not equal the 50,000 If it does then you
#' should set \code{limit} to a larger number. For example, if the returned data frame of records for only one
#' taxonomic groups has \code{nrow(df) = limit}, which the default is 50,000, increase the limit to ensure all
#' of the records are downloaded for that group. However, the maximum number of records that can be downloaded
#' is 200,000 at a time.
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(rgbif)
#' library(rgdal)
#' library(rgeos)
#'
#' outline<-readOGR(file.path(here(), 'data_raw', 'gis', 'Carib_Poly.shp'))
#'
#' h.occ<-gbif_records(groups = c('Amphibia','Reptilia'),
#'                     rank = 'class',
#'                     years = 2,
#'                     bounds = outline)
#'
#' }
#'
#' @export

gbif_records <- function (groups, rank, bounds, years,
                          yr_from = as.numeric(substr(Sys.Date(), 1, 4)),
                          limit = 50000,
                          genus_level = FALSE) {
  bound <- rgeos::writeWKT(bounds, byid = FALSE)
  yr.range <- paste(yr_from - years, yr_from, sep = ",")

  out <- data.frame(NULL)

  pb <- winProgressBar(title = "Progress Bar", min = 0, max = length(groups), width = 300)
  for (i in 1:length(groups)) {
    key <- rgbif::name_suggest(q = groups[i], rank = rank)$key[1]
    temp <- rgbif::occ_search(taxonKey = key, return = 'data', limit = limit,
                              fields = c("order", "family", "genus", "species", "eventDate", "year",
                                         "decimalLatitude", "decimalLongitude", "basisOfRecord",
                                         "occurrenceID", "references"),
                              geometry = bound,
                              year = yr.range)

    temp$Group <- groups[i]

    out <- rbind(out, temp)
    rm(key, temp)

    setWinProgressBar(pb, i, title=paste( round(i/length(groups)*100, 0), "% Done", sep = ''))
    if (i == length(groups)) {
      close(pb)
      cat('DONE! \n')
    } else {NA}
  }

  out <- out[, c("Group", "order", "family", "genus", "species", "eventDate", "year",
                 "decimalLatitude", "decimalLongitude", "basisOfRecord",
                 "occurrenceID", "references")]

  if (genus_level) {
    out <- out
  } else {
    out <- out[!is.na(out$species),]
  }

  return(out)
}
