#' Create a list of community matrices
#'
#' This function is used to create a list of i x j matrices where i is the number of geographical groups and j is the
#' total number of species. The matrix is filled in with the presence (1) and absence (0) of each species in the
#' data for a particular taxonomic group.
#' The function requires the package 'reshape'
#'
#' @param data A data frame that holds the species occurrence data with at least 4 columns of factors that holds the
#' geographic groups, taxonomic groups, species statuses, and species binomial names.
#'
#' @param species A character equal to the column name in which the species binomial names are stored
#'
#' @param geo_group A character equal to the column name in which the geographic groups are stored (e.g. bank names)
#'
#' @param taxa_group A character or vector equal to the column name(s) in which the names of the taxonomic groups
#' are stored (e.g. genus and/or family names)
#'
#' @param status A character equal to the column name in which the species' statuses are stored
#'
#' @param stat_levels A character or vector equal to the levels of the interest of the species' status column. The
#' default is to use all of the levels in the species' status column of the data.
#'
#' @param total Logical, If \code{TRUE} then the total community matrix for all species across all taxonomic groups is
#' included in the returned list. If \code{FLASE} then this complete community matrix is omitted.
#'
#' @return A list of length k where k is the number of taxonomic groups and statuses given to the function.
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(reshape)
#'
#' herp <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_final.csv'), header=TRUE)
#'
#' coms <- com_matrix(data = herp,
#'                    species = "binomial",
#'                    geo_group = "bank",
#'                    taxa_group = c("genus","taxon"),
#'                    status = "bnk_status",
#'                    stat_levels = c('N', "E"),
#'                    total = FALSE)
#'
#' }
#'
#' @export
com_matrix <- function(data,
                       species,
                       geo_group,
                       taxa_group,
                       status,
                       stat_levels = unique(data[, status]),
                       total = TRUE) {

  out<-list()

  df <- data.frame(NULL)
  for (m in 1:length(stat_levels)) {
    tmp <- droplevels(data[which(data[, status] == stat_levels[m]), ])
    df <- rbind(df, tmp)
  }

  data <- df
  rm(df)

  data[, species] <- gsub(' ', '_', data[, species], fixed = TRUE)
  data[, species] <- as.factor(data[, species])
  data[, geo_group] <- as.factor(data[, geo_group])
  data[, status] <- as.factor(data[, status])

  for (i in 1:length(taxa_group)) {
    data[, taxa_group[i]] <- as.factor(data[, taxa_group[i]])
  }

  data$num <- 1

  for (i in 1:length(taxa_group)) {
    nams <- NULL
    for (l in 1:length(stat_levels)) {
      grp <- data.frame(Group = levels(droplevels(data[which(data[, status] == stat_levels[l]),
                                                       taxa_group[i]])),
                        Status = rep(stat_levels[l],
                                     length(levels(droplevels(data[which(data[, status] == stat_levels[l]),
                                                                   taxa_group[i]])))))
      nams <- rbind(nams, grp)
    }

    grp <- levels(as.factor(nams$Group))
    for (k in 1:length(grp)) {
      for (l in 1:length(stat_levels)) {
        if (!is.element(stat_levels[l], nams[which(nams$Group == grp[k]), 'Status'])) {
          tmp <- data.frame(Group = grp[k], Status = stat_levels[l])
          nams <- rbind(nams, tmp)
        } else {NA}
      }
    }

    nams[, 1] <- as.character(nams[, 1])
    nams[, 2] <- as.character(nams[, 2])

    tmp <- list()
    for (j in 1:nrow(nams)) {
      if (is.element(nams[j, 2], data[which(data[, taxa_group[i]] == nams[j, 1]), status])) {
        tmp[[j]] <-  suppressMessages(reshape::cast(droplevels(data[which(data[, taxa_group[i]] == nams[j, 1] &
                                                                            data[, status] == nams[j, 2]), ]),
                                                    as.formula(paste(geo_group, "~", species)), min, fill = 0, value = 'num'))
      } else {
        rows <- levels(data[which(data[, taxa_group[i]] == nams[j, 1]), geo_group])
        cols <- levels(droplevels(data[which(data[, taxa_group[i]] == nams[j, 1]), species]))
        tmp[[j]] <- data.frame(matrix(0, nrow = length(rows), ncol = 1+length(cols)))
        names(tmp[[j]]) <- c(geo_group, cols)
        tmp[[j]][, geo_group] <- rows
        rm(rows, cols)
      }

    }
    names(tmp) <- paste(nams[, 1], nams[, 2], sep=".")
    out <- c(out, tmp)
    rm(tmp, nams)
  }


  if (total) {
    for (i in 1:length(taxa_group)) {
      nams <- data.frame(V1 = levels(droplevels(data[, taxa_group[i]])),
                         V2 = rep('T', length(levels(droplevels(data[, taxa_group[i]])))))

      nams[, 1] <- as.character(nams[, 1])
      nams[, 2] <- as.character(nams[, 2])

      tmp <- list()
      for (j in 1:nrow(nams)) {
        tmp[[j]] <- suppressMessages(reshape::cast(droplevels(data[which(data[, taxa_group[i]] == nams[j, 1]),]),
                                                   as.formula(paste(geo_group,"~",species)), min, fill = 0, value = 'num'))
      }
      names(tmp) <- paste(nams[, 1], nams[, 2], sep=".")
      out <- c(out, tmp)
      rm(tmp, nams)
    }

    nams <- expand.grid("All", stat_levels)
    nams[, 1] <- as.character(nams[, 1])
    nams[, 2] <- as.character(nams[, 2])

    tmp <- list()
    for (j in 1:nrow(nams)) {
      tmp[[j]] <- suppressMessages(reshape::cast(droplevels(data[which(data[, status] == nams[j, 2]),]),
                                                 as.formula(paste(geo_group,"~",species)), min, fill = 0, value = 'num'))
    }
    names(tmp) <- paste(nams[, 1], nams[, 2], sep=".")
    out <- c(out, tmp)
    rm(tmp, nams)

    tmp <- list()
    tmp[[1]] <- suppressMessages(reshape::cast(data, as.formula(paste(geo_group,"~",species)), min, fill = 0, value = 'num'))
    names(tmp)<-"All.T"

    out <- c(out, tmp)
    rm(tmp)
  } else {NA}

  out <- out[sort(names(out))]
  return(out)
}
