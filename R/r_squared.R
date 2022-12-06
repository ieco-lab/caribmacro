#' Calculate R Squared for Averaged Model(s)
#'
#' This function is used to calculate the multiple and adjusted R squared values for a set of model estimates.
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
#' @param display Logical. If \code{TRUE}, then the R Squared values are printed.
#'
#' @return A data frame with 1 row and 2 columns. The columns hold the Multiple R Squared and the Adjusted R
#' Squared, respectively.
#'
#' @examples
#' \dontrun{
#'
#' library(MuMIn)
#' library(AICcmodavg)
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
#' tmp <- r_squared(out[, c('Variable', 'Estimate')],
#'                   data = dat,
#'                   response = 'Anolis.N',
#'                   x_vars = vars,
#'                   display = FALSE)
#'
#' }
#'
#' @export
r_squared <- function(estimates, data, response, x_vars, display = TRUE) {
  dat <- data[complete.cases(data), ]
  x_dat <- dat[, x_vars]
  n <- nrow(dat)
  p <- length(x_vars)

  intercept <- setdiff(estimates[, 1], x_vars)
  if (length(intercept) == 0) {
    int <- 0
  } else if (length(intercept) > 1) {
    stop('length(setdiff(estimates[, 1], x_vars)) > 1: cannot determine Y Intercept')
  } else {
    int <- estimates[which(estimates[, 1] == intercept), 2]
  }

  if (length(setdiff(x_vars, estimates[, 1])) > 0) {
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

  rss <- sum((tmp$Y-tmp$Y.hat)^2)
  tss <- sum((tmp$Y-mean(tmp$Y, na.rm = TRUE))^2)

  R.sq <- 1 - (rss/tss)
  adj.R.sq <- 1 - ((1 - R.sq)*((n-1)/(n-p-1)))

  out <- data.frame(R.sq = R.sq,
                    adj.R.sq = adj.R.sq)
  if (display) {
  cat(paste('Multiple R Squared =', round(R.sq, 4), '\n'))
  cat(paste('Adjusted R Squared =', round(adj.R.sq, 4), '\n'))
  cat('\n')
  } else {NA}

  return(out)
}
