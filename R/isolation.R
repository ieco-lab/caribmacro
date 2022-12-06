#' Calculate the Isolation Metrics
#'
#' This function is used to calculate the isolation metrics based on a PCA of the isolation distances for the
#' geographic features of interest.
#' This function requires the package 'FactoMineR'
#'
#' @param dist a distance matrix in a data frame format that has the pairwise distance between all geographic
#' features of interest and the nearby continents. The row names of the data from should be the
#' names of the geographic features and should also be equal to the column names.
#'
#' @param source A character or vector of the column name(s) of the distance to the evolutionary source(s) of the
#' group of species that will be used for the calculation of diversity indexes.
#'
#' @param main A character or vector of the column name(s) of the distance to the nearby continent(s)
#'
#' @param geo.rm Optional. A character or vector of the column name(s) of the distance
#'
#' @param min Logical. If \code{TRUE} then the minimum distance to another bank is used in the PCA
#'
#' @param sum Logical. If \code{TRUE} then the sum of the all the distances to the other banks
#' is used in the PCA
#'
#' @param larger Logical. If \code{TRUE} then the minimum distance to a geographic feature of equal or larger size
#' is used in the PCA
#'
#' @param area Optional. A data frame with two columns. The first has to be the names of the geographic features used
#' in \code{dist} and the second the area for those geographic features. This argument is only needed if \code{larger = TRUE}
#'
#' @param sq.rt Logical. If \code{TRUE} then the distances will be square root transformed before the PCA is ran
#'
#' @param scale Logical. If \code{TRUE} then the data will be scaled in the PCA
#'
#' @param stats Logical. If \code{TRUE} then the the function will return a list of length 2 with the first element
#' holding the isolation metrics (PC1 - PC4), the second holding the eigenvalues and variance explained by the PCs,
#' and the third holding the loadings for the PCs.
#'
#' @param d.main Logical. If \code{TRUE} then the minimum distance to the mainland is used in the PCA
#'
#' @param export.dist Logical. If \code{TRUE} then all of the isolation distances are returned as a data frame in
#' the last element of the exported list.
#'
#' @return A data frame holding the PC1, PC2, and PC3 scores for each geographic feature of interest calculated
#' from a Principle Components Analysis of the isolation distances. If \code{stats = TRUE} then a list of
#' length 3 is returned with the second slot holding the fit statistics of the PCA and the third slot holding
#' the PC loadings for the isolation distances.
#'
#' Isolation distances are the summed distances between a geographic feature and all others (including the mainland),
#' the distance to nearest neighbor, minimum distance to mainland, and the minimum distance to a source geographic
#' feature.
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(sf)
#' library(FactoMineR)
#'
#' cent <- st_read(file.path(here(), 'gis', 'Centroids_Bank_Main_WGS_Merc.shp'))
#' cent.dist <- st_distance(x = cent, y = cent, by_element = FALSE)
#'
#' cent.dist <- data.frame(cent.dist)
#' names(cent.dist) <- gsub(' ', '.', cent$bank)
#' cent.dist$bank <- cent$bank
#' cent.dist < -cent.dist[order(cent.dist$bank),
#'                         c('bank', sort(names(cent.dist)[which(names(cent.dist) != 'bank')]))]
#'
#' for (i in 2:ncol(cent.dist)) {
#'   cent.dist[, i] <- as.numeric(cent.dist[, i])
#' }
#'
#' main <- c('central.america', 'north.america', 'south.america')
#' iso <- isolation(dist = cent.dist,
#'                  source = c('cuba', 'south.america'),
#'                  main = main)
#'
#' }
#'
#' @export

isolation <- function (dist, source, main,
                       geo.rm = NULL,
                       min = TRUE,
                       sum = TRUE,
                       larger = FALSE,
                       area = NULL,
                       sq.rt = TRUE,
                       scale = TRUE,
                       stats = TRUE,
                       d.main = TRUE,
                       export.dist = FALSE) {
  if (!min & !sum & !larger) {
    stop('At least "sum", "min", or "larger" has to be TRUE')
  } else {NA}

  if (larger & is.null(area)) {
    stop('Area not supplied for the calculation of distance to larger feature')
  } else {NA}

  tmp1 <- dist[setdiff(row.names(dist), geo.rm), setdiff(colnames(dist), geo.rm)]
  tmp1 <- tmp1[setdiff(row.names(tmp1), main), ]

  banks <- as.character(row.names(tmp1))

  tmp2 <- data.frame(bank = banks,
                     sum = rep(NA, length(banks)),
                     min = rep(NA, length(banks)),
                     main = rep(NA, length(banks)),
                     source = rep(NA, length(banks)))

  for (i in 1:length(banks)) {
    other <- names(tmp1)[which(tmp1[banks[i], ] == min(tmp1[banks[i], main]))]
    tmp2[i, 'sum'] <- sum(tmp1[banks[i], unique(c(banks, other))])
    tmp2[i, 'min'] <- min(tmp1[banks[i], c(banks[which(banks != banks[i])], main)])
    tmp2[i, 'main'] <- min(tmp1[banks[i], main])
    tmp2[i, 'source'] <- min(tmp1[banks[i], source])
  }

  if (larger) {
    tmp2$larger <- rep(NA, length(banks))

    for (i in 1:length(banks)) {
      lvl <- area[which(area[, 1] == banks[i]), 2]
      tmp2[i, 'larger'] <- min(tmp1[banks[i],
                                    c(main, unique(area[which(area[, 1] != banks[i] & area[, 2] >= lvl), 1]))])
    }
  } else {NA}

  if (sq.rt) {
    for (i in 2:ncol(tmp2)) {
      tmp2[, i] <- sqrt(tmp2[, i])
    }
  } else {NA}

  rownames(tmp2) <- tmp2$bank

  if (min & sum & larger) {
    mat <- tmp2[, c('sum', 'min', 'larger', 'main', 'source')]
  } else if (sum & min & !larger) {
    mat <- tmp2[, c('sum', 'min', 'main', 'source')]
  } else if (!sum & min & larger) {
    mat <- tmp2[, c('min', 'larger', 'main', 'source')]
  } else if (sum & !min & larger) {
    mat <- tmp2[, c('sum', 'larger', 'main', 'source')]
  } else if (sum & !min & !larger) {
    mat <- tmp2[, c('sum', 'main', 'source')]
  } else if (!sum & min & !larger) {
    mat <- tmp2[, c('min', 'main', 'source')]
  } else if (!sum & !min & larger) {
    mat <- tmp2[, c('larger', 'main', 'source')]
  } else {
    stop('Need at least min, sum, or distance to larger')
  }

  if (!d.main) {
    mat <- mat[, names(mat)[which(names(mat) != 'main')]]
  } else {NA}

  if (ncol(mat) < 3) {
    stop('Need at least 3 Isolation Distances')
  }


  temp <- FactoMineR::PCA(mat, scale = scale, graph = FALSE)

  iso <- data.frame(bank=row.names(as.data.frame(temp$ind$coord)),
                    PC1=as.data.frame(temp$ind$coord)[,1],
                    PC2=as.data.frame(temp$ind$coord)[,2],
                    PC3=as.data.frame(temp$ind$coord)[,3])

  var <- as.data.frame(temp[["eig"]])
  var$PC <- paste('PC', seq(1, nrow(temp[["eig"]]), 1), sep = ".")
  row.names(var) <- NULL
  names(var) <- c('Eigenvalue', 'X', 'Cum_Var', 'PC')
  var <- var[, c('PC', 'Eigenvalue', 'Cum_Var')]

  loadings <- as.data.frame(temp[["var"]]$cor)
  names(loadings) <- paste('PC', seq(1, nrow(temp[["var"]]$cor), 1), sep = ".")
  loadings$Variable <- row.names(loadings)
  row.names(loadings) <- NULL
  loadings <- loadings[, c('Variable', paste('PC', seq(1, nrow(temp[["var"]]$cor), 1), sep = "."))]

  if (stats) {
    out <- list(iso, var, loadings)
    names(out) <- c('Components', 'Variance', 'Loadings')
  } else {
    out <- iso
  }

  if (export.dist) {
    if (is.list(out)) {
      out <- c(out, list(tmp2))
      names(out)[length(out)] <- 'Distances'
    } else {
      out <- list(out, tmp2)
      names(out) <- c('Components', 'Distances')
    }
  } else {NA}

  return(out)
}
