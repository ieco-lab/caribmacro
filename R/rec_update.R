#' Filter Out New GBIF Occurrence Records
#'
#' This function is used to relate the occurrence records from GBIF returned by the function \code{gbif_records()}
#' with geographic features of interest and select the species occurrences that are not in the data to be
#' updated for those geographic features.
#' The function requires the package 'sf'
#'
#' @param gbif_occ The resulting data frame from the \code{gbif_records()} funtion in the 'caribmacro' package
#'
#' @param geography A shapefile of the geographic features of interest. Must be read in using the \code{st_read()}
#'          function of the 'sf' package (e.g. multipolygon shapefile of Caribbean island banks)
#'
#' @param geo_name A character equal to the column name of the shapefile attribute table in which the geographic
#'            feature names are stored (e.g. bank names)
#'
#' @param data A data frame of the data to be updated
#'
#' @param sp_name A character equal to the column name of the 'data' that has the species' binomial names
#'
#' @param num_rec A number equal to the number of records needed to denote a credible species record on a
#' particular bank. The default is \code{num_rec = 5}
#'
#' @return A data frame of species records filtered from the data frame returned by \code{gbif_records}
#' for geographic features in a shapefile and matched to a database of records.
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(sf)
#'
#' banks <- st_read(file.path(here(), 'data_raw', 'gis', 'Carib_Banks.shp'))
#' dat <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_final.csv'), header=TRUE)
#'
#' ## records is a data frame returned by gbif_records()
#'
#' new_recs <- rec_update(gbif_occ = records,
#'                        geography = banks,
#'                        geo_name = 'bank',
#'                        data = dat,
#'                        sp_name = 'binomial.cleaned',
#'                        num_rec = 5)
#'
#' }
#'
#' @export

rec_update <- function(gbif_occ, geography, geo_name, data, sp_name, num_rec = 5) {
  temp <- sf::st_as_sf(gbif_occ, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(geography)$proj4string)
  att.tab <- sf::st_drop_geometry(geography)
  tmp1 <- sf::st_nearest_feature(temp, geography)
  temp[, geo_name] <- NA


  temp<-as.data.frame(temp)
  for (i in 1:nrow(temp)) {
    temp[i, geo_name] <- as.character(att.tab[tmp1[i], geo_name])
  }

  out<-data.frame(NULL)
  geo_group <- levels(as.factor(temp[, geo_name]))

  data[, sp_name] <- as.factor(data[, sp_name])

  for (i in 1:length(geo_group)){
    tmp2 <- setdiff(unique(temp[which(temp[, geo_name] == geo_group[i]), 'species']),
                    levels(droplevels(data[which(data[, geo_name] == geo_group[i]), sp_name])))

    for (j in 1:length(tmp2)){
      out <- rbind(out, temp[which(temp[, geo_name] == geo_group[i] & temp$species == tmp2[j]), ])
    }
  }

  coords <- as.character(out$geometry)
  coords <- sub('c(', '', coords, fixed = TRUE)
  coords <- sub(')', '', coords, fixed = TRUE)
  coords <- strsplit(coords, ',')
  coords <- matrix(unlist(coords), ncol=2, byrow=TRUE)
  coords <- as.data.frame(coords)
  names(coords) <- c("decimalLongitude", "decimalLatitude")

  out <- cbind(out, coords)
  out <- as.data.frame(out[, which(names(out)!='geometry')])


  out$Multiple_Records<-NA
  out$ref<-paste(out[, geo_name], out$species, sep="_")
  refs<-unique(out$ref)

  for (i in 1:length(refs)) {
    if (nrow(out[which(out$ref == refs[i]),]) >= num_rec) {
      out[which(out$ref == refs[i]),'Multiple_Records']<-'YES'
    } else {
      out[which(out$ref == refs[i]),'Multiple_Records']<-'NO'
    }
  }

  return(out[, which(names(out) != 'ref')])
}
