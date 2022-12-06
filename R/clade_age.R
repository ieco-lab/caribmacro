#' Calculate the crown age of clades
#'
#' This function is used to calculate the crown ages for a set of clades within a phylogenetic tree.
#' This function requires the package 'ape'
#'
#'
#' @param tree A time-resolved newick formatted tree with tip labels equal to the species binomial
#' names with an "_" such as the trees downloaded from timetree.org.
#'
#' @param taxonomy A data frame with each row holding the species' taxonomic hierarchy, such as what
#' is returned from \code{hierarchy()}, or the clade names in which each species belongs.
#'
#' @param species A character equal to the name of the column that holds the species binomial names.
#' This column can also to be a different taxonomic level if needed.
#'
#' @param clades A character or vector equal to the name(s) of the columns that holds the clade names
#'
#' @param total Logical. If \code{TRUE}, the crown age for the entire group of species is returned
#'
#' @param sp_sep A character equal to the separator of the binomial name in the species name column.
#' For example: for 'Apis mellifera' \code{sp_sep = ' '} or for 'Apis.mellifera' \code{sp_sep = '.'}
#'
#' @return A data frame with three columns. One for the clade name, one for the clade rank (i.e.
#' column name in the input data frame), and one for the crown age of the clade.
#'
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(ape)
#'
#' tree <- read.tree(file.path(here(), 'data_raw', 'Tetrapoda_species.nwk'))
#' herp <- read.csv(file.path(here(), 'data_out', 'supp_info', 'IBT_Carib_Taxonomy.csv'), header=TRUE)
#'
#' ages <- clade_age(tree = tree,
#'                   taxonomy = taxa,
#'                   species = 'sp_in_Data',
#'                   clades = c("Class", "Order", "Family", "Genus"),
#'                   total = TRUE,
#'                   sp_sep = " ")
#' }
#'
#' @export
clade_age <- function(tree, taxonomy, species, clades, total = TRUE, sp_sep = ' ') {

  row.names(taxonomy) <- sub(sp_sep, '_', taxonomy[, species])

  out <- data.frame(NULL)

  for (i in 1:length(clades)) {
    groups <- as.character(sort(unique(taxonomy[, clades[i]])))

    tmp1 <- data.frame(NULL)

    for (j in 1:length(groups)) {
      sub.tree <- ape::drop.tip(tree, setdiff(tree$tip.label,
                                              row.names(taxonomy[which(taxonomy[, clades[i]] == groups[j]), ])))

      if (!is.null(sub.tree)) {
        tmp2 <- data.frame(Clade = groups[j],
                           Rank = clades[i],
                           Age = max(ape::node.depth.edgelength(sub.tree)))
      } else {
        tmp2 <- data.frame(Clade = groups[j],
                           Rank = clades[i],
                           Age = NA)
      }

      tmp1 <- rbind(tmp1, tmp2)
    }

    out <- rbind(out, tmp1)
  }

  if (total) {
  tot.tree <- drop.tip(tree, setdiff(tree$tip.label, row.names(taxonomy)))

  out <- rbind(out, data.frame(Clade = 'All',
                               Rank = NA,
                               Age = max(ape::node.depth.edgelength(tot.tree))))
  } else {NA}

  return(out)
}










