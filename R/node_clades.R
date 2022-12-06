#' Descendant Species of Nodes
#'
#' This function is used to determine the descendant species of each node in a phylogenetic tree.
#' This function requires the package 'ape' and 'phangorn'
#'
#'
#' @param tree A time-resolved newick formatted tree
#'
#' @param data A data frame of occurrence or taxonomic data
#'
#' @param tip.lvl A character equal to the name of the column in \code{data} that holds the tip labels of \code{tree}.
#'
#' @param species A character equal to the name of the column in \code{data} that holds the species' names.
#'
#' @param sep A character equal to the separator of the tip labels and the separator to be used in the column headers
#' of the resultant data frame.
#' For example: for 'Apis_mellifera' \code{sep = '_'} (the default) or for 'Apis.mellifera' \code{sep = '.'}
#'
#' @param report Logical. If \code{TRUE}, the species or groups (i.e. tip labels) that are in \code{data} but not in \code{tree}
#' are returned.
#'
#' @param age Logical. If \code{TRUE}, the clade crown ages are calculated and returned in the matrix.
#'
#' @return A i x j data frame of presence/absence of j species in i clades determined by the phylogenetic nodes.
#'
#' @examples
#' tree <- read.tree(file.path(here(), 'data_raw', 'trees', 'Tetrapoda_genus.nwk'))
#' herp <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_final.csv'), header = TRUE)
#' taxa <- read.csv(file.path(here(), 'data_raw', 'IBT_Sp_Taxonomy_v1.csv'), header = TRUE)
#'
#' names(taxa)[which(names(taxa) == 'sp_in_Data')] <- 'binomial'
#' herp <- merge(herp, taxa[, c('binomial', 'Genus')], by = 'binomial', all = TRUE)
#'
#' clades <- node_clades(tree = tree, data = herp, tip.lvl = 'Genus', species = 'binomial', sep = '_', report = TRUE)
#'
#' @export

node_clades <- function(tree, data, tip.lvl, species, sep = '_', report = TRUE, age = FALSE) {

  tips <- gsub(' ', sep, levels(as.factor(data[, tip.lvl])), fixed = TRUE)

  new.tree <- ape::drop.tip(tree, setdiff(tree$tip.label, tips))

  if (tip.lvl != species) {
    sp <- levels(droplevels(as.factor(data[is.element(tips, data[, tip.lvl]), species])))
  } else {
    sp <- tips
  }

  sp <- gsub(' ', sep, sp, fixed = TRUE)

  a <- length(new.tree$tip.label) + 1
  b <- new.tree$Nnode + length(new.tree$tip.label)

  out <- as.data.frame(matrix(NA, new.tree$Nnode, length(sp)))
  names(out) <- sp
  row.names(out) <- a:b

  if (age) {
    out$Age <- NA
    out <- out[, c('Age', sort(names(out)[which(names(out) != 'Age')]))]
  } else {
    out <- out[, sort(names(out)[which(names(out) != 'Age')])]
  }

  des <- phangorn::Descendants(new.tree, a:b)
  names(des) <- a:b

  for (i in 1:length(des)) {
    des[[i]] <- new.tree$tip.label[des[[i]]]

    spp <- levels(droplevels(as.factor(data[is.element(data[, tip.lvl], des[[i]]), species])))
    spp <- gsub(' ', sep, spp, fixed = TRUE)

    out[names(des)[i], spp] <- 1
    out[names(des)[i], is.na(out[names(des)[i], ])] <- 0

    sub.tree <- ape::drop.tip(tree, setdiff(tree$tip.label, des[[i]]))
    if (age) {
      out[names(des)[i], 'Age'] <- max(ape::node.depth.edgelength(sub.tree))
    } else {NA}

    rm(spp, sub.tree)
  }

  if (report) {
    cat(paste0('Species or Groups not in Tree', ' (n = ',
               length(setdiff(tips, tree$tip.label)), ')', ': ', '\n',
               paste(setdiff(tips, tree$tip.label), collapse = " | ")))

  } else {NA}

  return(out)
}
