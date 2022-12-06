#' Create Model Estimates ~ Clade Age Plots
#'
#' This function is used to create a scatter plot of model estimates against clade age.
#' Requires package 'ggplot2'
#'
#' @param data A data frame holding at least the beta estimates and age for each clade.
#'
#' @param X A character equal to the name of the column in \code{data} that holds the clade ages
#'
#' @param Y A character equal to the name of the column in \code{data} that holds the estimates
#'
#' @param Group Optional. A character equal to the name of the column in \code{data} that
#' holds the groups.
#'
#' @param abs.value Logical. If \code{TRUE}, then the absolute value of the estimates is plotted.
#'
#' @param var Optional. A character equal to the name of the variable the beta estimates are for
#'
#' @param Title Optional. A character to be displayed as the plot's title
#'
#' @param palette Optional. A vector with length equal to the number of groups to be displayed.
#' This vector is passed on to \code{scale_color_manual()} to dictate the color representation
#' of each group
#'
#' @param shapes Optional. A list of length 2 with the first element holding a character equal to
#' the name of the column in \code{data} that holds the groups. The second element should be a
#' vector with length equal to the number of groups to be displayed. This vector is passed on to
#' \code{scale_shape_manual()} to dictate the shape representation of each group
#'
#' @param fill Optional. A vector that holds the color(s) to fill the point shapes with
#'
#' @param plot_theme Optional. An object returned by the \code{theme()} function.
#' See \code{help(theme)} for more details
#'
#' @return A scatter plot of estimates against clade ages in the form of a ggplot object
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(ggplot2)
#'
#' ests <- read.csv(file.path(here(), 'data_out', 'Avg_Est_SR.csv'), header = TRUE)
#' ages <- read.csv(file.path(here(), 'data_out', 'Herp_Taxon_Ages.csv'), header = TRUE)
#' names(ages)[1] <- 'Taxon'
#'
#' dat <- merge(ests, ages, by = 'Taxon', all.x = TRUE)
#' dat <- droplevels(dat[which(dat$Group == 'Total' | dat$Group == 'Native'), ])
#' tmp <- droplevels(dat[which(dat$Variable == 'Area'), ])
#'
#' (plot_est_age(data = tmp,
#'               X = "Age",
#'               Y = "Estimate",
#'               Group = "Group",
#'               var = "Area",
#'               Title = "Effect of Area by Clade Age",
#'               palette = c("red3", "blue3"),
#'               shapes = list("Group", c(16, 23)),
#'               fill = 'white'))
#'
#' }
#'
#' @export

plot_est_age <- function(data, X, Y, Group = NULL, abs.value = FALSE, var = NULL, Title = NULL,
                         palette = NULL, shapes = NULL, fill = NULL, plot_theme = NULL) {

  if (abs.value) {
    data[, Y] <- abs(data[, Y])
    ltyp <- 0
  } else {
    ltyp <- 2
    }

  aes_plot <- aes_string(X, Y, colour = Group)

  if (!is.null(shapes)) {
    aes_point <- aes_string(shape = shapes[[1]])
    p <- ggplot(data, aes_plot) +
      geom_abline(intercept = 0, slope = 0, linetype = ltyp, size = 0.5) +
      geom_point(aes_point, fill = fill, size = 2) +
      scale_shape_manual(values = shapes[[2]])
  } else {
    p <- ggplot(data, aes_plot) +
      geom_abline(intercept = 0, slope = 0, linetype = 2, size = 0.5) +
      geom_point(size = 2)
  }

  if (!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  } else {NA}


  p <- p + xlab("Clade Age (myr)") +
            ylab(substitute('Standardized '~beta[v], list(v = var))) +
            ggtitle(Title)

  if (is.null(plot_theme)) {
    p <- p + theme(axis.line.y = element_line(colour = "black"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.text.y = element_text(colour = "black", size = 12),
                   axis.text.x = element_text(colour = "black", size = 12),
                   axis.title.y=element_text(size=14),
                   axis.title.x=element_text(size=14),
                   plot.title = element_text(hjust = 0.5, colour = "black", size = 16),
                   legend.title = element_text(colour="black", size = 14),
                   legend.text = element_text(colour="black", size = 12),
                   legend.key=element_blank())
  } else {
    p <- p + plot_theme
  }

  if (abs.value) {
    p <- p + scale_y_continuous(expand = c(0, 0),
                                limits = c(0, (max(data[, Y]) + 0.1*max(data[, Y])))) +
               theme(axis.line.x = element_line(colour = "black"),
                     axis.ticks.x = element_line(colour = "black"))
  } else {
    p <- p + theme(axis.ticks.x = element_blank())
  }

  return(p)
}
