#' Run Constant Age or Size Bootstrapping
#'
#' This function is used bootstrap the curvature - age or curvature - size analyses while holding the clade size
#' or clade age constant, respectively. This functions requires the packages 'ape', 'ggplot2', and 'wesanderson.' In addition,
#' the packages required for the 'caribmacro' functions \code{com_matrix()}, \code{SR_geo()}, and \code{sr_LM()} are also
#' required for this function.
#'
#' @param occurrences A data frame that holds the species occurrence data with at least 4 columns of factors that holds the
#' geographic groups, taxonomic groups, species statuses, and species binomial names.
#'
#' @param real.dat A data frame that has the observed curvatures for each phylogenetic clade. The columns should be 'Clade'
#' for the clade names, 'Age' for the clade ages, 'SR' for the clade sizes, 'ABC.t' for the curvature of the total assemblage,
#' and 'ABC.n' for the curvature of the native assemblage.
#'
#' @param node.mat A data frame that has the presence/absence of each phylogenetic clade. The first column of this data
#' frame should hold the clade names and be named 'Clade'
#'
#' @param tree.gen A phylogenetic tree with genera as the tips.
#'
#' @param lineage A vector of characters that equal the clade names of the evolutionary lineage of interest.
#'
#' @param clade A character equal to the name of the smallest clade in the lineage.
#'
#' @param species_names A character equal to the column name in 'occurrences' in which the species binomial names are stored
#'
#' @param con_size Logical. If true constant size bootstrapping is performed. Default is \code{con_size = FALSE}.
#'
#' @param con_age Logical. If true constant age bootstrapping is performed. Default is \code{con_age = FALSE}.
#'
#' @param sp.num Numeric. If \code{con_age = TRUE}, 'sp.num' should be a single numeric equal to the clade size to be held
#' constant for the constant age bootstrapping.\cr
#' If \code{con_age = TRUE}, 'sp.num' should be a vector of varying clade sizes
#' (e.g., \code{seq(<size of smallest clade>, <size of clade at age>, (<size of clade at age> - <size of smallest clade>)/5)}.\cr
#' Default is \code{sp.num = NULL}.
#'
#' @param age Numeric. The age of a clade along the lineage that will be held constant for the constant age bootstrapping.
#' Required if \code{con_age = TRUE}. Default is \code{age = NULL}.
#'
#' @param runs Numeric. The number of random clades to be made for each age or clade size value. Default is \code{runs = 1000}.
#'
#' @param keep.gen.sp Logical. If true a species records that are only identified to their genus are retained in the
#' species-area curve analysis. Default is \code{keep.gen.sp = TRUE}
#'
#' @param make_plot Logical. If true a list is returned with the ggplot object saved in the third element named 'Plot'.
#' Default is \code{make_plot = FALSE}
#'
#' @param just.clades Logical. If true a just a data frame is returned with the presence/absence of species in each
#' of the randomly made clades. Default is \code{just.clades = FALSE}.
#'
#' @param bank.data A data frame that holds the explanatory variables for each geographic
#' feature of interest, specifically their area. Default is \code{bank.data = NULL}.
#'
#' @param status A character equal to the column name in 'occurrences' in which the species' statuses are stored.
#' Default is \code{status = NULL}.
#'
#' @param Area A character equal to the column name in 'bank.data' in which the area of the geographic feature is stored.
#' Default is \code{Area = NULL}.
#'
#' @param geo_group A character equal to the column name in 'bank.data' _AND_ 'occurrences' in
#' which the geographic feature names are stored. Default is \code{geo_group = NULL}.
#'
#' @param stat_levels A character or vector equal to the levels of the interest of the species' status column in
#' 'occurrences'. The default is to use all of the levels in the species' status column of the data.
#' Default is \code{stat_levels = NULL}.
#'
#' @return If \code{just.clades = TRUE}, a data frame with the with the presence/absence of species in each
#' of the randomly made clades. If \code{just.clades = TRUE}, a list of length 2 if \code{make_plot = FALSE} or
#' length 3 if \code{make_plot = TRUE}. The elements of this list are 'Matrix' (the presence/absence of species
#' in each clade), 'Curvature' (a data frame with the curvature of the species-area curves for each clade and
#' the clade ages and sizes), and if \code{make_plot = TRUE}, 'Plot' (a ggplot oblect of the results).
#'
#'
#' @examples
#' \dontrun{
#' if (FALSE) {
#'
#' node.dat <- read.csv(file.path(here(), 'data_out', 'results', 'sar_lin', 'supp_info', 'Node_ABC-Age_Data.csv'), header = TRUE)
#' nod.mat <- read.csv(file.path(here(), 'data_out', 'supp_info', 'IBT_Node_Clade_Data.csv'), header = TRUE)
#' tree.gen <- read.tree(file.path(here(), 'data_raw', 'trees', 'Tetrapoda_genus.nwk'))
#' herp <- read.csv(file.path(here::here(), 'data_raw', 'IBT_Herp_Records_final.csv'), header=TRUE)
#'
#' herp$stat_new <- NA
#' for (i in 1:nrow(herp)) {
#'   if (herp[i, 'bnk_status'] == 'E') {
#'     herp[i, 'stat_new'] <- 'E'
#'   } else if (herp[i, 'bnk_status'] == 'FE') {
#'     herp[i, 'stat_new'] <- 'N'
#'   } else if (herp[i, 'bnk_status'] == 'N') {
#'     herp[i, 'stat_new'] <- 'N'
#'   } else if (herp[i, 'bnk_status'] == 'PX') {
#'     herp[i, 'stat_new'] <- 'N'
#'   } else if (herp[i, 'bnk_status'] == 'U') {
#'     herp[i, 'stat_new'] <- 'E'
#'   } else {
#'     herp[i, 'stat_new'] <- 'X'
#'   }
#' }
#'
#' anol <- node.dat[which(node.dat$Lineage == clade), ]
#' lin <- anol$Clade[order(anol$Clade, decreasing = TRUE)]
#'
#' res <- curv_boot(occurrences = herp,
#'                  real.dat = node.dat,
#'                  node.mat = nod.mat,
#'                  tree.gen = tree.gen,
#'                  lineage = lin,
#'                  clade = 'Anolis',
#'                  species_names = 'binomial',
#'                  con_age = TRUE,
#'                  sp.num = seq(240, 600, 40),
#'                  age = 201,
#'                  runs = 1000,
#'                  make_plot = TRUE,
#'                  bank.data = bank_dat,
#'                  status = 'stat_new',
#'                  Area = 'Area',
#'                  geo_group = 'bank',
#'                  stat_levels = c("N", "E"))
#'
#' }
#' }
#'
#' @export

curv_boot <- function(occurrences, real.dat, node.mat, tree.gen, lineage, clade, species_names, con_age = FALSE, con_size = FALSE,
                     sp.num = NULL, age = NULL, runs = 1000, keep.gen.sp = TRUE, make_plot = FALSE, just.clades = FALSE,
                     bank.data = NULL, status = NULL, Area = NULL, geo_group = NULL, stat_levels = NULL) {
  # Determine Lineage
  n.dat <- real.dat[is.element(real.dat[, 1], lineage), ]

  # Determine Initial Clade Size
  cld.size <- data.frame(Clade = node.mat[, 1], SR = rowSums(node.mat[, 3:ncol(node.mat)]))
  cld.size <- max(cld.size[which(cld.size$Clade == clade), 'SR'])

  n.sr <- n.dat[which(n.dat$Clade == clade), 'SR.t']

  if (cld.size != n.sr) {
    stop("Disagreement between node.mat and real.dat -> \n SR in real.dat does not equal SR from node.mat for initial clade in lineage")
  }

  if (!any(con_age, con_size)) {
    stop("Need to set 'con_age' or 'con_size' as TRUE for constant age or size, respectively, bootstrapping")
  } else if (all(con_age, con_size)) {
    stop("'curv_boot' currently cannot run both constant age and constant size bootstrapping at the same time")
  }

  if (con_age & is.null(age)) {
    stop("Need to specify 'age' for constant age bootstrapping")
  }

  if (con_size & is.null(sp.num)) {
    stop("Need to specify 'sp.num' for constant size bootstrapping")
  }

  sp.list <- unique(occurrences[, species_names])

  gen.sp <- sp.list[which(substr(sp.list, nchar(sp.list)-2, nchar(sp.list)) == ' sp')]
  gen.gen <- substr(gen.sp, 1, nchar(gen.sp)-3)

  sp.list <- sp.list[which(substr(sp.list, nchar(sp.list)-2, nchar(sp.list)) != ' sp')]

  clad.mat <- NULL
  if (con_size) {
    pb <- winProgressBar(title = "Random Clades", min = 0, max = length(lineage), width = 300)
    for (i in 1:length(lineage)) {
      tmp.mat <- as.data.frame(matrix(0, runs, length(sp.list)))
      names(tmp.mat) <- gsub(' ', '_', sp.list)
      row.names(tmp.mat) <- paste(lineage[i], 1:runs, sep = '_')
      tmp.mat$Age <- NA

      for (j in 1:runs) {

        if (lineage[i] == clade) {
          sp.tmp <- names(nod.mat[, which(nod.mat[which(nod.mat$Clade == clade), ] == 1)])
          sp.tmp <- sp.tmp[which(substr(sp.tmp, nchar(sp.tmp)-2, nchar(sp.tmp)) != '_sp')]
          rand.sp <- gsub(' ', '_', sample(sp.tmp, sp.num, replace = FALSE))
          tmp.mat[j, rand.sp] <- 1

          # Clade age
          tips <- gsub(' ', '_', sp.tmp)
          sub.tree <- ape::drop.tip(tree.spp, setdiff(tree.spp$tip.label, tips))
          tmp.mat[j, 'Age'] <- max(ape::node.depth.edgelength(sub.tree))

        } else {
          tmp1 <- names(nod.mat[, which(nod.mat[which(nod.mat$Clade == lineage[i-1]), ] == 1)])
          tmp1 <- tmp1[which(substr(tmp1, nchar(tmp1)-2, nchar(tmp1)) != '_sp')]
          tmp2 <- names(nod.mat[, which(nod.mat[which(nod.mat$Clade == lineage[i]), ] == 1)])
          tmp2 <- tmp2[which(substr(tmp2, nchar(tmp2)-2, nchar(tmp2)) != '_sp')]
          tmp3 <- setdiff(tmp2, tmp1)
          tmp4 <- gsub(' ', '_', sample(tmp3, 1, replace = FALSE))

          rand.sp <- c(tmp4, gsub(' ', '_', sample(tmp2[!is.element(tmp2, tmp4)], (sp.num - 1), replace = FALSE)))
          tmp.mat[j, rand.sp] <- 1

          # Clade age
          tips <- unique(data.frame(matrix(unlist(strsplit(rand.sp, '_')), nrow = length(strsplit(rand.sp, '_')), byrow = TRUE))[, 1])
          sub.tree <- ape::drop.tip(tree.gen, setdiff(tree.gen$tip.label, tips))
          tmp.mat[j, 'Age'] <- max(ape::node.depth.edgelength(sub.tree))
        }
      }

      tmp.mat <- as.data.frame(tmp.mat)

      old <- length(names(as.data.frame(tmp.mat))[!is.element(names(as.data.frame(tmp.mat)), names(clad.mat))]) > 0
      new <- length(names(clad.mat)[!is.element(names(clad.mat), names(as.data.frame(tmp.mat)))]) > 0

      if (i == 1) {
        clad.mat <- rbind(clad.mat, as.data.frame(tmp.mat))
      } else if (old | new) {
        new_col_old <- c(names(as.data.frame(tmp.mat))[!is.element(names(as.data.frame(tmp.mat)), names(clad.mat))])
        new_col_new <- c(names(clad.mat)[!is.element(names(clad.mat), names(as.data.frame(tmp.mat)))])

        if (old) {
          clad.mat[, new_col_old] <- 0
        }

        if (new) {
          tmp.mat[, new_col_new] <- 0
        }

        clad.mat <- rbind(clad.mat, as.data.frame(tmp.mat))
      } else {
        clad.mat <- rbind(clad.mat, as.data.frame(tmp.mat))
      }

      setWinProgressBar(pb, i, title = paste0("Random Clades (", round(i/length(lineage)*100, 0), "% Done)"))
      if (i == length(lineage)) {
        close(pb)
      } else {NA}
    }

    clad.mat$Clade <- row.names(clad.mat)
    clad.mat <- clad.mat[, c('Clade', 'Age', sort(names(clad.mat[which(names(clad.mat) != 'Clade' & names(clad.mat) != 'Age')])))]

    err.age <- names(summary(as.factor(clad.mat$Age)))[which(summary(as.factor(clad.mat$Age)) != runs)]
    if (length(err.age != 0)) {
      age.correct <- NULL
      for (i in 1:length(err.age)) {
        tmp <- abs(as.numeric(err.age[i]) - as.numeric(err.age))
        age.correct <- rbind(age.correct, data.frame(err = err.age[i],
                                                     correct = err.age[which(tmp == min(tmp[-i]))],
                                                     keep = NA))
        age.correct[i, 'keep'] <- age.correct$correct[i] == min(age.correct$err[i], age.correct$correct[i])
      }
      age.correct <- age.correct[age.correct$keep, ]
    }

    for (i in 1:nrow(age.correct)) {
      clad.mat[which(clad.mat$Age == age.correct$err[i]), 'Age'] <- age.correct$correct[i]
    }

    x_lab <- 'Clade Age (mya)'

  }



  if (con_age & !length(sp.num) > 1) {
    stop("'sp.num' needs to be a vector of clade sizes for constant age bootstrapping (i.e. 'con_age = TRUE')")
  }

  if (con_age) {
    if (is.na(nchar(strsplit(as.character(age), "\\.")[[1]][2]))) {
      rnd <- 0
    } else {
      rnd <- nchar(strsplit(as.character(age), "\\.")[[1]][2])
    }
    node.at.age <- node.mat[which(round(node.mat$Age, rnd) == age), 'Clade']

    pb <- winProgressBar(title = "Random Cladses", min = 0, max = length(sp.num), width = 300)
    for (i in 1:length(sp.num)) {

      tmp.mat <- as.data.frame(matrix(0, runs, length(sp.list)))
      names(tmp.mat) <- gsub(' ', '_', sp.list)
      row.names(tmp.mat) <- paste(clade, sp.num[i], 1:runs, sep = '_')
      tmp.mat$Age <- NA

      for (j in 1:runs) {
        ps <- which(lineage == node.at.age)
        tmp1 <- names(node.mat[, which(node.mat[which(node.mat$Clade == lineage[ps-1]), ] == 1)])
        tmp1 <- tmp1[which(substr(tmp1, nchar(tmp1)-2, nchar(tmp1)) != '_sp')]
        tmp2 <- names(node.mat[, which(node.mat[which(node.mat$Clade == node.at.age), ] == 1)])
        tmp2 <- tmp2[which(substr(tmp2, nchar(tmp2)-2, nchar(tmp2)) != '_sp')]
        tmp3 <- setdiff(tmp2, tmp1)
        tmp4 <- gsub(' ', '_', sample(tmp3, 1, replace = FALSE))
        tmp5 <- gsub(' ', '_', sample(tmp1, 1, replace = FALSE))

        rand.sp <- c(tmp4, tmp5, gsub(' ', '_', sample(tmp2[!is.element(tmp2, tmp4) & !is.element(tmp2, tmp5)],
                                                       (sp.num[i] - 2),
                                                       replace = FALSE)))
        tmp.mat[j, rand.sp] <- 1

        # Clade age
        tips <- unique(data.frame(matrix(unlist(strsplit(rand.sp, '_')), nrow = length(strsplit(rand.sp, '_')), byrow = TRUE))[, 1])
        sub.tree <- ape::drop.tip(tree.gen, setdiff(tree.gen$tip.label, tips))
        tmp.mat[j, 'Age'] <- max(ape::node.depth.edgelength(sub.tree))
      }

      old <- length(names(as.data.frame(tmp.mat))[!is.element(names(as.data.frame(tmp.mat)), names(clad.mat))]) > 0
      new <- length(names(clad.mat)[!is.element(names(clad.mat), names(as.data.frame(tmp.mat)))]) > 0

      if (i == 1) {
        clad.mat <- rbind(clad.mat, as.data.frame(tmp.mat))
      } else if (old | new) {
        new_col <- c(names(as.data.frame(tmp.mat))[!is.element(names(as.data.frame(tmp.mat)), names(clad.mat))],
                     names(clad.mat)[!is.element(names(clad.mat), names(as.data.frame(tmp.mat)))])
        clad.mat[, new_col] <- 0
        clad.mat <- rbind(clad.mat, as.data.frame(tmp.mat))
      } else {
        clad.mat <- rbind(clad.mat, as.data.frame(tmp.mat))
      }
      setWinProgressBar(pb, i, title = paste0("Random Clades (", round(i/length(sp.num)*100, 0), "% Done)"))
      if (i == length(sp.num)) {
        close(pb)
      } else {NA}
    }

    clad.mat$Clade <- row.names(clad.mat)
    clad.mat <- clad.mat[, c('Clade', 'Age', sort(names(clad.mat[which(names(clad.mat) != 'Clade' & names(clad.mat) != 'Age')])))]

    err.age <- names(summary(as.factor(clad.mat$Age)))[which(summary(as.factor(clad.mat$Age)) != runs*length(sp.num))]
    if (length(err.age) != 0) {
      age.correct <- NULL
      for (i in 1:length(err.age)) {
        tmp <- abs(as.numeric(err.age[i]) - as.numeric(err.age))
        age.correct <- rbind(age.correct, data.frame(err = err.age[i],
                                                     correct = err.age[which(tmp == min(tmp[-i]))],
                                                     keep = NA))
        age.correct[i, 'keep'] <- age.correct$correct[i] == min(age.correct$err[i], age.correct$correct[i])
      }
      age.correct <- age.correct[age.correct$keep, ]

      for (i in 1:nrow(age.correct)) {
        clad.mat[which(clad.mat$Age == age.correct$err[i]), 'Age'] <- age.correct$correct[i]
      }

      if (length(unique(clad.mat$Age)) != 1) {
        stop('Non-constant resulting clade age')
      }
    }

    x_lab <- 'Clade Size (Num. of Sp.)'
  }

  # Remove replicate samples
  if(!(nrow(clad.mat) == nrow(unique(clad.mat[, 3:ncol(clad.mat)])))) {
    clad.mat <- unique(clad.mat[, 3:ncol(clad.mat)])
  }

  for (i in 1:ncol(clad.mat)) {
    clad.mat[is.na(clad.mat[, i]), i] <- 0
  }

  nams <- sort(names(clad.mat[which(names(clad.mat) != 'Clade' & names(clad.mat) != 'Age')]))
  clad.mat$SR <- rowSums(clad.mat[, !is.element(names(clad.mat), c('Clade', 'Age', 'SR'))])

  clad.mat <- clad.mat[, c('Clade', 'Age', 'SR', nams)]

  if (!just.clades) {
    if (!just.clades & any(is.null(bank.data), is.null(status), is.null(Area), is.null(geo_group), is.null(stat_levels))) {
      stop("To run analyses the options for the 'com_matrix' and 'SR_geo' and 'sr_LM' functions need to be included")
    }

    rand.clades <- clad.mat

    coms <- list()
    pb <- winProgressBar(title = "Community Matrices", min = 0, max = nrow(rand.clades), width = 300)
    for (i in 1:nrow(rand.clades)) {
      sp <- names(rand.clades[i, which(rand.clades[i, ] == 1)])
      sp <- gsub('_', ' ', sp, fixed = TRUE)

      if (keep.gen.sp) {
        for (m in 1:length(gen.gen)) {
          if (any(grepl(gen.gen[m], sp, fixed = TRUE))) {
            sp <- c(sp, paste0(gen.gen[m], ' sp'))
          }
        }
      }

      if (length(sp) > 0) {
        tmp <- droplevels(occurrences[is.element(occurrences$binomial, sp), ])
        tmp$Clade <- rand.clades[i, 'Clade']

        if (all(is.element(stat_levels, tmp[, status]))) {
          temp <- com_matrix(data = tmp,
                             species = species_names,
                             geo_group = geo_group,
                             taxa_group = 'Clade',
                             status = status,
                             stat_levels = stat_levels,
                             total = TRUE)

          bad.nams <- paste('All', c(stat_levels, 'T'), sep = '.')
          coms <- c(coms, temp[!is.element(names(temp), bad.nams)])
          suppressMessages(rm(temp))
        } else {NA}#stop('At least 1 status level not in data set for one of the clades')}
        suppressMessages(rm(tmp))
      } else {NA}#stop('At least 1 of the clades has no species')}

      suppressMessages(rm(sp))
      setWinProgressBar(pb, i, title = paste0("Community Matrices (", round(i/nrow(rand.clades)*100, 0), "% Done)"))
      if (i == nrow(rand.clades)) {
        close(pb)
      } else {NA}


    }

    coms <- coms[!grepl('.E', names(coms), fixed = TRUE)]


    sp.rich <- SR_geo(data = coms, geo_group = geo_group)

    # Rename column in SR for lm formula
    names(sp.rich)[grep('.T', names(sp.rich), fixed = TRUE)] <- paste0("T_",
                                                             names(sp.rich)[grep('.T', names(sp.rich),
                                                                            fixed = TRUE)])
    names(sp.rich)[grep('T_', names(sp.rich), fixed = TRUE)] <- gsub(".T", "",
                                                           names(sp.rich)[grep('T_', names(sp.rich),
                                                                          fixed = TRUE)], fixed = TRUE)
    names(sp.rich)[grep('.N', names(sp.rich), fixed = TRUE)] <- paste0("N_",
                                                             names(sp.rich)[grep('.N', names(sp.rich),
                                                                            fixed = TRUE)])
    names(sp.rich)[grep('N_', names(sp.rich), fixed = TRUE)] <- gsub(".N", "",
                                                           names(sp.rich)[grep('N_', names(sp.rich),
                                                                          fixed = TRUE)], fixed = TRUE)

    x_vars <- Area
    bnk.area <- bank.data[, c("bank", x_vars)]

    ## Curvilinear (squared Area)
    lm2.sr <- sr_LM(SR = sp.rich,
                    data = bnk.area,
                    x_vars = x_vars,
                    geo_group = geo_group,
                    area = Area,
                    log_area = TRUE,
                    sq_area = TRUE,
                    complete.case = TRUE)

    ## Linear
    lm.sr <- sr_LM(SR = sp.rich,
                   data = bnk.area,
                   x_vars = x_vars,
                   geo_group = geo_group,
                   area = Area,
                   log_area = TRUE,
                   complete.case = TRUE)

    mods <- names(lm.sr[['Models']])

    # Calculate the ABC values for each Model using the definite integral method
    abc <- data.frame(Model = mods, ABC = rep(NA, length(mods)))
    pb <- winProgressBar(title = "Curvature Calculation", min = 0, max = length(mods), width = 300)
    for (i in 1:length(mods)) {
      f_a <- lm.sr[["Models"]][[mods[i]]][["coefficients"]]
      g_a <- lm2.sr[["Models"]][[mods[i]]][["coefficients"]]

      z_eq <- function(x) {abs((g_a[1] + g_a[2]*x + g_a[3]*(x^2)) - (f_a[1] + f_a[2]*x))}

      tmp <- integrate(z_eq,
                       lower = min(lm.sr[["Data"]][[mods[i]]]$Area),
                       upper = max(lm.sr[["Data"]][[mods[i]]]$Area))

      abc[which(abc$Model == mods[i]), 'ABC'] <- tmp[[1]]
      setWinProgressBar(pb, i, title = paste0("Curvature Calculation (", round(i/length(mods)*100, 0), "% Done)"))
      if (i == length(mods)) {
        close(pb)
      } else {NA}
    }

    tmp <- data.frame(matrix(unlist(strsplit(abc$Model, '_')), nrow = length(strsplit(abc$Model, '_')), byrow = TRUE))

    if (ncol(tmp) == 3) {
      abc$Clade <- paste(tmp[, 2], tmp[, 3], sep = '_')
    } else if (ncol(tmp) == 4) {
      abc$Clade <- paste(tmp[, 2], tmp[, 3], tmp[, 4], sep = '_')
    } else {
      stop('Error in clade names; Suggest renaming clades in lineage')
    }

    abc$Group <- tmp[, 1]

    abc <- merge(abc, rand.clades[, c('Clade', 'Age', 'SR')], by = 'Clade', all.x = TRUE)
  }

  abc$Age <- as.numeric(abc$Age)

  if (make_plot) {
    # Plot with real data
    ages <- data.frame(Clade = lineage, Age = as.integer(unique(clad.mat$Age)))

    real.t.dat <- real.dat[which(real.dat$Lineage == clade), c("Clade", "ABC.t", 'SR.t')]
    real.n.dat <- real.dat[which(real.dat$Lineage == clade), c("Clade", "ABC.n", 'SR.t')]

    names(real.t.dat)[2:3] <- c('ABC', 'SR')
    real.t.dat$Group <- 'T'

    names(real.n.dat)[2:3] <- c('ABC', 'SR')
    real.n.dat$Group <- 'N'

    plot.dat <- rbind(real.n.dat, real.t.dat)
    plot.dat <- merge(plot.dat, ages, by = 'Clade', all.x = TRUE)

    plot.dat$Assemblage <- 'Native'
    plot.dat[which(plot.dat$Group == 'T'), 'Assemblage'] <- 'Total'

    abc$Assemblage <- 'Native'
    abc[which(abc$Group == 'T'), 'Assemblage'] <- 'Total'

    cols <- wesanderson::wes_palettes[['Cavalcanti1']][c(5, 1, 3)]
    cols <- c('Introduced' = cols[1], 'Native' = cols[2], 'Total' = cols[3])

    plot_grps <- c('Native', 'Total')

    if (con_size) {
      p <- ggplot(abc, aes(x = Age, y = ABC, color = Assemblage, fill = Assemblage)) +
        geom_smooth(linetype = 'dashed', method = lm, formula = y ~ log(x), se = FALSE) +
        geom_smooth(data = plot.dat, linetype = 'solid', method = lm, formula = y ~ log(x), se = FALSE) +
        geom_point(data = plot.dat, size = 2) +
        scale_color_manual(values = cols[plot_grps]) +
        scale_fill_manual(values = alpha(cols[plot_grps])) +
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y = c(0, max(plot.dat$ABC) + 0.25),
                      x = c(min(plot.dat$Age) - 10, max(plot.dat$Age) + 10)) +
        ylab("SAR Curvature (\U01B6)") + xlab(x_lab) +
        theme_classic() +
        theme(legend.position = c(0.8, 0.9),
              legend.text = element_text(size = 12, color = 'black'),
              legend.title = element_text(size = 12, color = 'black'),
              axis.text = element_text(size = 12, color = 'black'),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5),
              panel.grid.major.x = element_blank())
    }

    if (con_age) {
      p <- ggplot(abc, aes(x = SR, y = ABC, color = Assemblage, fill = Assemblage)) +
        geom_smooth(linetype = 'dashed', method = lm, formula = y ~ log(x), se = FALSE) +
        geom_smooth(data = plot.dat, linetype = 'solid', method = lm, formula = y ~ log(x), se = FALSE) +
        geom_point(data = plot.dat, size = 2) +
        scale_color_manual(values = cols[plot_grps]) +
        scale_fill_manual(values = alpha(cols[plot_grps])) +
        scale_y_continuous(expand = c(0, 0)) +
        expand_limits(y = c(0, max(plot.dat$ABC) + 0.25),
                      x = c(min(sp.num) - 10, max(plot.dat$SR) + 10)) +
        ylab("SAR Curvature (\U01B6)") + xlab(x_lab) +
        theme_classic() +
        theme(legend.position = c(0.8, 0.9),
              legend.text = element_text(size = 12, color = 'black'),
              legend.title = element_text(size = 12, color = 'black'),
              axis.text = element_text(size = 12, color = 'black'),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5),
              panel.grid.major.x = element_blank())
    }
  }

  if (just.clades) {
    out <- clad.mat
  }

  if (!just.clades) {
    out <- list(Matrix = clad.mat, Curvature = abc)
  }

  if (!just.clades & make_plot) {
    out <- list(Matrix = clad.mat, Curvature = abc, Plot = p)
  }

  return(out)
  cat('DONE! \n')
}
