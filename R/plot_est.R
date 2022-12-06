#' Create Model Estimates Plot
#'
#' This function is used to create a dot and whisker plot of model estimates.
#' requires package 'ggplot2'
#'
#' @param data A data frame holding the beta estimate and their confidence intervals
#'
#' @param aes_plot An aesthetic mapping object returned by the \code{aes()} function and passed on to
#' \code{ggplot()}.
#'
#' @param aes_err An aesthetic mapping object returned by the \code{aes()} function and passed on to
#' \code{geom_errorbar()}.
#'
#' @param Title Optional. A character to be displayed as the plot's title
#'
#' @param r_squared Optional. A number equal to the R squared for the model to be displayed as
#' the plot's subtitle
#'
#' @param xlabs Optional. A vector of labels for each value along the x axis
#'
#' @param palette Optional. A vector with length equal to the number of groups to be displayed.
#' This vector is passed on to \code{scale_color_manual()} to dictate the color representation
#' of each group
#'
#' @param shapes Optional. A list of length 2 with the first element holding \code{aes(shapes = group)}
#' to be passed on to \code{geom_points} where 'group' is the same as in \code{aes_plot}.
#' The second element should be a vector with length equal to the number of groups to be displayed.
#' This vector is passed on to \code{scale_shape_manual()} to dictate the shape representation
#' of each group
#'
#' @param fill Optional. A vector that holds the color(s) to fill the point shapes with
#'
#' @param plot_theme Optional. An object returned by the \code{theme()} function.
#' See \code{help(theme)} for more details
#'
#' @return A dot and whisker plot in the form of a ggplot object
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(ggplot2)
#'
#' ests <- read.csv(file.path(here(), 'data_out', 'Avg_Est_SR.csv'), header = TRUE)
#' tmp <- droplevels(ests[which(ests$Taxon == 'Anolis' & ests$Variable != '(Intercept)'), ])
#'
#' plot_est(data = tmp,
#'          aes_plot = aes(x = Variable, y = Estimate, colour = Group),
#'          aes_err = aes(ymin = LCL, ymax = UCL, col = Group),
#'          Title = "Anolis Lizards",
#'          r_squared = data.frame(V1 = c('Exotic', 'Native', 'Total'),
#'                                 R2 = c(0.23, 0.1, 0.123456)),
#'          palette = c("red3", "blue3", "chocolate3"),
#'          shapes = list(aes(shape = Group), c(16, 23, 15)),
#'          fill = 'white')
#'
#' }
#'
#' @export

plot_est <- function(data, aes_plot, aes_err,
                     Title = NULL, r_squared = NULL, xlabs = NULL,
                     palette = NULL, shapes = NULL, fill = NULL, plot_theme = NULL) {
  if (!is.null(shapes)){
    p <- ggplot(data, aes_plot) +
      geom_abline(intercept = 0, slope = 0, linetype = 2, size = 0.5) +
      geom_errorbar(aes_err, position = position_dodge(width = 0.5), lwd = 1, width = 0) +
      geom_point(shapes[[1]], position = position_dodge(width = 0.5), fill = fill, size = 2)
  } else {
    p <- ggplot(data, aes_plot) +
      geom_abline(intercept = 0, slope = 0, linetype = 2, size = 0.5) +
      geom_point(position = position_dodge(width = 0.5), size = 2) +
      geom_errorbar(aes_err, position = position_dodge(width = 0.5), lwd = 1, width = 0)
  }

  if (!is.null(xlabs)) {
    p <- p + scale_x_discrete(labels = xlabs)
  } else {NA}

  if (!is.null(r_squared)) {
    if (is.data.frame(r_squared)) {
      l <- c(NULL)
      for (i in 1:nrow(r_squared)) {
        r_squared <- r_squared[order(r_squared[, 1]), ]
        nam <- r_squared[i,1]
        rs <- format(round(r_squared[i,2], 3), nsmall = 3)
        l <- c(l, paste(nam, '\n', expression(R^2), ' = ', rs))
      }
      if (!is.null(shapes)) {
        p <- p + scale_shape_manual(values = shapes[[2]], labels = l)
        p <- p + ggtitle(Title)
      } else {
        p <- p + scale_shape_manual(labels = l)
        p <- p + ggtitle(Title)
      }

      if (!is.null(palette)){
        p <- p + scale_color_manual(values = palette, labels = l)
        p <- p + ggtitle(Title)
      } else {
        p <- p + scale_color_manual(labels = l)
        p <- p + ggtitle(Title)
      }

    } else if (is.numeric(r_squared)){
      if (!is.null(shapes)) {
        p <- p + scale_shape_manual(values = shapes[[2]])
        p <- p + ggtitle(Title, subtitle = substitute((R^{2}~~'='~~r), list(r = r_squared)))
      } else {
        p <- p + ggtitle(Title, subtitle = substitute((R^{2}~~'='~~r), list(r = r_squared)))
      }

      if (!is.null(palette)){
        p <- p + scale_color_manual(values = palette)
        p <- p + ggtitle(Title, subtitle = substitute((R^{2}~~'='~~r), list(r = r_squared)))
      } else {
        p <- p + ggtitle(Title, subtitle = substitute((R^{2}~~'='~~r), list(r = r_squared)))
      }

    } else {
      message('R Squared is not in a readable format: Ignoring r_squared argument')
      if (!is.null(shapes)) {
        p <- p + scale_shape_manual(values = shapes[[2]])
        p <- p + ggtitle(Title)
      } else {
        p <- p + ggtitle(Title)
      }

      if (!is.null(palette)){
        p <- p + scale_color_manual(values = palette)
        p <- p + ggtitle(Title)
      } else {
        p <- p + ggtitle(Title)
      }

    }
  } else {
    if (!is.null(shapes)) {
      p <- p + scale_shape_manual(values = shapes[[2]])
      p <- p + ggtitle(Title)
    } else {
      p <- p + ggtitle(Title)
    }

    if (!is.null(palette)){
      p <- p + scale_color_manual(values = palette)
      p <- p + ggtitle(Title)
    } else {
      p <- p + ggtitle(Title)
    }
  }

  p <- p + xlab("Explanatory Variable") + ylab("Standardized Beta")

  if (is.null(plot_theme)) {
    p <- p + theme(axis.line.y = element_line(colour = "black"),
                   axis.ticks.x = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.text.y = element_text(colour = "black", size = 12),
                   axis.text.x = element_text(colour = "black", size = 12,
                                              angle = 90, hjust = 1, vjust = 0.5),
                   axis.title.y=element_text(size=14),
                   axis.title.x=element_text(size=14),
                   plot.title = element_text(hjust = 0.5, colour = "black", size = 16),
                   plot.subtitle = element_text(hjust = 0.5, colour = "black", size = 12),
                   legend.title = element_text(colour="black", size = 14),
                   legend.text = element_text(colour="black", size = 12),
                   legend.key=element_blank(),
                   legend.key.height = unit(2.5, "line"),
                   legend.key.width = unit(2, "line"))
  } else {
    p <- p + plot_theme
  }

  return(p)
}
