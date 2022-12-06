#' Pairwise Community Dissimilarity between Geogrphical Units
#'
#' This function is used to calculate the Sorensen or Jaccard dissimilairty indices that have also been partitioned for community
#' turnover and nestedness. This function requires the package \code{betapart} and is a wrapper for the function \code{beta.pair()}.
#'
#' @param coms A data frame that holds the species occurrence data with at least 4 columns of factors that holds the
#' geographic groups, taxonomic groups, species statuses, and species binomial names.
#'
#' @param geos A character equal to the column name in which the geographic groups are stored (e.g. island names or bank names)
#'
#' @param index.family A character that states the family of dissimilarity indices to be used. The options are "sorensen" or "jaccard".
#' Default is \code{index.family = "sorensen"}
#'
#' @param abund Logical. If \code{TRUE} then the Bray-Curtis Index will be returned that accounts for differences in abundance. Default is
#' \code{FALSE} so that presence/absence based dissimilarity measures (i.e., \code{family = 'sorensen'} or \code{family = 'jaccard'})
#' will be returned.
#'
#' @param return.mat Logical, If \code{TRUE} then the dissimilarity matrix or matrices is included in the returned list.
#' If \code{FLASE} then just the data frame(s) of the pairwise beta diversity measures are returned. Default is \code{FALSE}.
#'
#'
#' @return A list of length k where k is number of elements in \code{coms} if \code{coms} is a list (i.e. \code{length(coms)})
#' or if  \code{coms} is a data frame the list will be of the various dissimilarity metrics.
#'
#' If  \code{return.mat = TRUE}, then result is a list of length 2 with the first element holding a data frame(s) of the pairwise
#' beta diversity measures and the second element a  j x j matrix or matrices where j is the number of geographic groups you have.
#'
#' In the returned list(s) there are three measures of community dissimilarity. Therefore, in each element of the returned list there
#' are three elements (i.e. one for each kind dissimilarity index). The three idices correspond to (i) turnover (replacement), (ii)
#' nestedness-resultant component, and (iii) total dissimilarity (i.e. the sum of both components).
#'
#' For Sorensen indices (\code{index.family = "sorensen"} each element is named "beta.sim", "beta.sne", and "beta.sor".
#'
#' \code{beta.sim} is the dissimilarity due to species turnover (i.e. Simpson pair-wise dissimilarity)
#'
#' \code{beta.sim} is the dissimilarity due to community nestedness measured as the nestedness-fraction of Sorensen pair-wise
#' dissimilarity
#'
#' \code{beta.sim} is the total dissimilarity between communities (i.e. Sorensen pair-wise dissimilarity)
#'
#' ***
#'
#' For Jaccard indices (\code{index.family = "jaccard"} each element is named "beta.jtu", "beta.jne", and "beta.jac".
#'
#' \code{beta.jtu} is the dissimilarity due to species turnover measured as the turnover-fraction of Jaccard pair-wise dissimilarity
#'
#' \code{beta.jne} is the dissimilarity due to community nestedness measured as the nestedness-fraction of Jaccard pair-wise
#' dissimilarity
#'
#' \code{beta.jac} is the total dissimilarity between communities (i.e. Jaccard pair-wise dissimilarity)
#'
#' ***
#'
#' For abundance-based dissimilarity (\code{abund = TRUE}), each element is named "beta.bray.bal", "beta.bray.gra", and "beta.bray"
#'
#' \code{beta.bray.bal} is the dissimilarity accounting for the dissimilarity derived from balanced variation in abundance
#' between sites
#'
#' \code{beta.bray.gra} is the dissimilarity accounting for the dissimilarity derived from unidirectional abundance gradients
#'
#' \code{beta.bray} is the total dissimilarity accounting for total abundance-based dissimilarity between sites, measured as the
#' Bray-Curtis index
#'
#'
#' @examples
#' \dontrun{
#'
#' library(betapart)
#' library(here)
#'
#' recs <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_final.csv'), header = TRUE)
#'
#' bnk_coms <- com_matrix(data = recs,
#'                        species = "binomial",
#'                        geo_group = "bank",
#'                        taxa_group = "class",
#'                        status = "stat_new",
#'                        stat_levels = c("N", "E"),
#'                        total = TRUE)
#'
#' beta_div <- beta.pair_geo(coms = bnk_coms, geos = 'bank', index.family = "sorensen", return.mat = FALSE)
#'
#' }
#'
#' @export

beta.pair_geo <- function(coms, geos = NULL, index.family = "sorensen", abund = FALSE, return.mat = FALSE) {

  if (inherits(coms, "list")) {
    grps <- names(coms)
    com.df <- FALSE
  } else {
    coms <- list('A' = coms)
    grps <- names(coms)
    com.df <- TRUE
  }

  beta_mat <- vector('list', length(coms))
  names(beta_mat) <- grps

  beta_tabs <- vector('list', length(coms))
  names(beta_tabs) <- grps

  for (i in 1:length(coms)) {
    tmp <- coms[[grps[i]]]

    if (!is.null(geos)) {
      row.names(tmp) <- tmp[, geos]
      tmp <- tmp[, names(tmp)[which(names(tmp) != geos)]]
    }

    if (!abund) {
      beta_mat[[grps[i]]] <- betapart::beta.pair(tmp, index.family = index.family)
    } else {
      beta_mat[[grps[i]]] <- betapart::beta.pair.abund(tmp, index.family = 'bray')
    }

    beta_tabs[[grps[i]]] <- vector('list', length(names(beta_mat[[grps[i]]])))
    names(beta_tabs[[grps[i]]]) <- names(beta_mat[[grps[i]]])

    for (j in 1:3) {
      beta_mat[[grps[i]]][[j]] <- as.matrix(beta_mat[[grps[i]]][[j]])
      colnames(beta_mat[[grps[i]]][[j]]) <- gsub(' ', '.', row.names(tmp), fixed = TRUE)
      row.names(beta_mat[[grps[i]]][[j]]) <- gsub(' ', '.', row.names(tmp), fixed = TRUE)

      geos.tmp <- colnames(beta_mat[[grps[i]]][[j]])
      temp <- NULL
      for (k in 1:length(geos.tmp)) {
        if (k != length(geos.tmp)) {
          for (l in (k + 1):length(geos.tmp)) {
            tmp2 <- data.frame(comp = paste0(geos.tmp[k], '-', geos.tmp[l]),
                               beta.div = beta_mat[[grps[i]]][[j]][geos.tmp[k], geos.tmp[l]])
            temp <- rbind(temp, tmp2)
            rm(tmp2)
          }
        }
      }

      beta_tabs[[grps[i]]][[j]] <- temp
      rm(temp, geos.tmp)
    }
    rm(tmp)
  }

  if (com.df) {
    beta_tabs <- beta_tabs[['A']]
    beta_mat <- beta_mat[['A']]
  }

  if (return.mat) {
    return(list(Tables = beta_tabs, Matrices = beta_mat))
  } else {
    return(beta_tabs)
  }

}
