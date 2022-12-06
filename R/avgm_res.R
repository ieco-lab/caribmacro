#' Diagnostic Plots and Residual Determination
#'
#' This function is used to create Diagnostic Plots and calculate the residuals for a set of model estimates.
#'
#' @param estimates A data frame with 2 columns. The first column must hold the variable names and if there is an
#' intercept, the intercept should also have a name in this column. The second column should hold the estimates
#'
#' @param data The data frame used in the model
#'
#' @param response A character equal to the name of the column in \code{data} that holds the response variable
#'
#' @param x_vars A character or vector of the explanatory variable names. Should also equal the column names
#' in \code{data} that holds the explanatory variable data. The name of the intercept should NOT be included
#' in this vector.
#'
#' @param plots Logical. If \code{TRUE}, then diagnostic plots are created and saved
#'
#' @param res.export Logical. If \code{TRUE}, then the residuals are saved in a data frame with the observed
#' and fitted values
#'
#' @param norm Logical. If \code{TRUE}, a Shapiro-Wilk normality test is performed on the residuals.
#'
#' @param show Logical. If \code{TRUE}, the plots will be show in a 1 x 3 paneled figure. If \code{FALSE},
#' the figures will still be saved but not shown so that the user has to plot them. Additionally, if
#' \code{norm = TRUE}, the results of a Shapiro-Wilk normality test on the residuals will be printed if
#' \code{show = TRUE}.
#'
#' @return A list or data frame with the information saved as designated by the \code{plots}, \code{res.export},
#' and \code{norm} arguments.
#'
#' @examples
#' \dontrun{
#'
#' library(MuMIn)
#' library(AICcmodavg)
#' library(ggplot2)
#' library(gridExtra)
#'
#' ## lm.sr is an object returned from sr.LM()
#' dat <- lm.sr[['Data']][['Anolis.N']]
#' out <- AIC_avg(lm.sr[['Models']][['Anolis.N']],
#'                data = dat,
#'                response = 'Anolis.N',
#'                x_vars = c(vars, 'sq_Area'),
#'                groups = data.frame(Genus = 'Anolis',
#'                                    Status = 'N'),
#'                cum.weight = 0.95,
#'                table = FALSE)
#'
#' vars <- attr(lm.sr[['Models']][['Anolis.N']]$terms,'term.labels')
#'
#' tmp <- avgm_res(out[, c('Variable', 'Estimate')],
#'                   data = dat,
#'                   response = 'Anolis.N',
#'                   x_vars = vars)
#'
#' }
#'
#' @export
avgm_res <- function(estimates, data, response, x_vars,
                     plots = TRUE, res.export = TRUE, norm = FALSE, show = FALSE) {
  if (nrow(data) > 0) {
    dat <- data[complete.cases(data), ]
    x_dat <- dat[, x_vars]

    intercept <- setdiff(estimates[, 1], x_vars)
    if (length(intercept) == 0) {
      int <- 0
    } else if (length(intercept) > 1) {
      stop('length(setdiff(estimates[, 1], x_vars)) > 1: cannot determine Y Intercept')
    } else {
      int <- estimates[which(estimates[, 1] == intercept), 2]
    }

    if (length(setdiff(c(intercept, x_vars), estimates[, 1])) > 0) {
      add <- data.frame(V1 = setdiff(c(intercept, x_vars), estimates[, 1]),
                        V2 = rep(0, length(setdiff(c(intercept, x_vars), estimates[, 1]))))
      names(add) <- names(estimates)
      estimates <- rbind(estimates, add)
    } else {NA}

    tmp <- data.frame(Y = dat[, response],
                      Y.hat = rep(NA, length(dat[, response])))

    for (i in 1:nrow(tmp)) {
      ests <- c(NULL)
      X <- c(NULL)
      for (j in 1:length(x_vars)) {
        ests <- c(ests, estimates[which(estimates[, 1] == x_vars[j]), 2])
        X <- c(X, x_dat[i, x_vars[j]])
      }

      tmp[i, 'Y.hat'] <- int + sum(ests*X)
    }

    tmp$res <- (tmp$Y-tmp$Y.hat)
    tmp$s.res <- (tmp$res - mean(tmp$res))/sd(tmp$res)
    tmp$sq.res <- sqrt(abs(tmp$s.res))

    if (plots) {
      plot_theme <- ggplot2::theme(axis.line = element_line(colour = "black"),
                                   axis.text = element_text(colour = "black"),
                                   axis.title = element_text(colour = "black"),
                                   plot.title = element_text(hjust = 0.5, colour = "black"),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   panel.background = element_blank(),
                                   legend.position = 'none')

      res.v.fit <- ggplot2::ggplot(tmp, aes(x = Y.hat, y = res)) +
        ggplot2::geom_abline(intercept = 0, slope = 0, linetype = 2, size = 0.5) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = 'loess', formula = y ~ x, color = 'red', se = FALSE) +
        ggplot2::labs(x = 'Fitted Values', y = 'Residuals', title = 'Residuals vs Fitted') +
        plot_theme

      qq_plot <- ggplot2::ggplot(tmp) +
        ggplot2::geom_qq_line(aes(sample = s.res), color = 'red', linetype = 2, size = 1) +
        ggplot2::geom_qq(aes(sample = s.res)) +
        ggplot2::labs(x = 'Theoretical Quantiles', y = 'Standardized Residuals', title = 'Normal Q-Q') +
        plot_theme

      scale <- ggplot2::ggplot(tmp, aes(x = Y.hat, y = sq.res)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = 'loess', formula = y ~ x, color = 'red', se = FALSE) +
        ggplot2::labs(x = 'Fitted Values', y = expression(sqrt('|Standardized Residuals|')),
                      title = 'Spread-Location') +
        plot_theme

      if (show) {
        gridExtra::grid.arrange(grobs = list(res.v.fit, qq_plot, scale),
                                nrow = 3, heights = c(3,3,3), widths = 4)
      } else {NA}
    } else {NA}

    if (norm) {
      tst <- shapiro.test(tmp$res)
      tst[["data.name"]] <- 'Residuals'
      if (show) {
        cat(paste0('Test of Residual Normality\n\n',
                   tst[["method"]], ': ','\n',
                   'W = ', round(tst[["statistic"]][["W"]], 4), '\n',
                   'P-value = ', round(tst[["p.value"]], 4), '\n'))
      } else {NA}
    } else {NA}

    tmp <- tmp[, c('Y', 'Y.hat', 'res')]
    names(tmp) <- c('Observed', 'Fitted', 'Residuals')

    if (plots) {
      if (res.export & norm) {
        out <- list(Residuals = tmp,
                    Plots = list(Residuals = res.v.fit, QQ = qq_plot, Spread = scale),
                    Test = tst)
      } else if (res.export & !norm) {
        out <- list(Residuals = tmp,
                    Plots = list(Residuals = res.v.fit, QQ = qq_plot, Spread = scale))
      } else if (!res.export & norm) {
        out <- list(Plots = list(Residuals = res.v.fit, QQ = qq_plot, Spread = scale),
                    Test = tst)
      } else {
        out <- list(Residuals = res.v.fit, QQ = qq_plot, Spread = scale)
      }
    } else if(res.export & norm) {
      out <- list(Residuals = tmp,
                  Test = tst)
    } else if (res.export & !norm) {
      out <- tmp
    } else if (!res.export & norm) {
      out <- tst
    } else {
      message(paste0('No Information to Export:', '\n',
                     'At least one of "plots", "res.export", and "norm" needs to = TRUE'))
    }
  } else {
    if (plots) {
      if (res.export & norm) {
        out <- list(Residuals = c(NA),
                    Plots = list(Residuals = c(NA), QQ = c(NA), Spread = c(NA)),
                    Test = c(NA))
      } else if (res.export & !norm) {
        out <- list(Residuals = c(NA),
                    Plots = list(Residuals = c(NA), QQ = c(NA), Spread = c(NA)))
      } else if (!res.export & norm) {
        out <- list(Plots = list(Residuals = c(NA), QQ = c(NA), Spread = c(NA)),
                    Test = c(NA))
      } else {
        out <- list(Residuals = c(NA), QQ = c(NA), Spread = c(NA))
      }
    } else if(res.export & norm) {
      out <- list(Residuals = c(NA),
                  Test = c(NA))
    } else if (res.export & !norm) {
      out <- c(NA)
    } else if (!res.export & norm) {
      out <- c(NA)
    } else {
      message(paste0('No Information to Export:', '\n',
                     'At least one of "plots", "res.export", and "norm" needs to = TRUE'))
    }

    message(paste0('data is empty:', '\n',
                   '    returning a list of NAs'))
  }

  return(out)
}
