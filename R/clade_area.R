#' Area of Geographic Feature Occupancy
#'
#' This function calculates the summed land area of all of the geographic feature a taxonomic clade(s) occupies.
#'
#' @param data A data frame that holds the species occurrence data with at least 4 columns of factors that holds the
#' geographic groups, taxonomic groups, species statuses, and species binomial names.
#'
#' @param area A data frame with two columns. The first has to be the names of the geographic features used
#' in \code{data} and the second the area for those geographic features.
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
#' @return A data frame with two columns. The first column will hold the clade names and the second column the
#' area of bank occupancy in the units of the area data frame.
#'
#' @examples
#' \dontrun{
#' clade.area <- clade_area(data = herp,
#'                          area = bnk.area,
#'                          geo_group = 'bank',
#'                          taxa_group = c('Class', 'Order', 'Suborder', 'Family', 'Genus'),
#'                          status = 'bnk_status',
#'                          stat_levels = c("N", "E"))
#'
#' }
#'
#' @export
clade_area <- function(data,
                       area,
                       geo_group,
                       taxa_group,
                       status,
                       stat_levels = unique(data[, status])) {

  if (any(!complete.cases(area))) {
    stop('area cannot have missing values')
  } else if (!is.numeric(area[, 2])) {
    stop('area[, 2] is not numerical')
  } else {NA}

  out <- data.frame(NULL)

  df <- data.frame(NULL)
  for (m in 1:length(stat_levels)) {
    tmp <- droplevels(data[which(data[, status] == stat_levels[m]), ])
    df <- rbind(df, tmp)
  }

  data <- df
  rm(df)

  data <- caribmacro::Factorize(data = data, columns = c(geo_group, status, taxa_group))

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

    temp <- data.frame(NULL)
    for (j in 1:nrow(nams)) {
      if (is.element(nams[j, 2], data[which(data[, taxa_group[i]] == nams[j, 1]), status])) {
        bnks <- levels(droplevels(data[which(data[, taxa_group[i]] == nams[j, 1] &
                                               data[, status] == nams[j, 2]), geo_group]))
        areas <- NULL
        for (p in 1:length(bnks)) {
          areas <- c(areas, area[which(area[, 1] == bnks[p]), 2])
        }

        tmp <- data.frame(Clade = nams[j, 1],
                          Group = nams[j, 2],
                          Area = sum(areas))
        temp <- rbind(temp,tmp)
      } else {
        tmp <- data.frame(Clade = nams[j, 1],
                          Group = nams[j, 2],
                          Area = 0)
        temp <- rbind(temp,tmp)
      }
      rm(tmp)

      if (!is.element(nams[j, 1], temp[which(temp$Group == 'T'), 'Clade'])) {
        bnks <- levels(droplevels(data[which(data[, taxa_group[i]] == nams[j, 1]), geo_group]))

        areas <- NULL
        for (l in 1:length(bnks)) {
          areas <- c(areas, area[which(area[, 1] == bnks[l]), 2])
        }

        tmp <- data.frame(Clade = nams[j, 1],
                          Group = 'T',
                          Area = sum(areas))
        temp <- rbind(temp,tmp)
      } else {NA}

    }
    out <- rbind(out, temp)
    rm(temp, nams)
  }

  return(out)
}
