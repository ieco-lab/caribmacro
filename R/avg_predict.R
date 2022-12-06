#' Predict for an Averaged Model
#'
#' This function is used to predict Y values from averaged estimates for specified X values.
#'
#' @param estimates A data frame with 2 columns. The first column must hold the variable names and if there is an
#' intercept, the intercept should also have a name in this column. The second column should hold the estimates
#'
#' @param data The data frame used in the model
#'
#' @param new.data A data frame that holds the new data to predict from. Should have the same number of columns
#' as \code{length(x_vars)}
#'
#' @param x_vars A character or vector of the explanatory variable names. Should also equal the column names
#' in \code{data} that holds the explanatory variable data. The name of the intercept should NOT be included
#' in this vector.
#'
#' @param response A character equal to the name of the column in \code{data} that holds the response variable
#'
#' @return A data frame with \code{length(x_vars)} + 1 columns and \code{nrow(new.data)} rows. The data frame has
#' the \code{new.data} and the estimated values (in a column named 'Y.hat').
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
#' preds <- avg_predict(out[, c('Variable', 'Estimate')],
#'                      new.data = dat,
#'                      x_vars = vars)
#'
#' }
#'
#' @export
avg_predict <- function(estimates, new.data, x_vars, se = FALSE, data = NULL, response = NULL) {
  if (se & is.null(data)) {
    stop('Need "data" and/or "response" to calculate SE')
  } else if (se & is.null(response)) {
    stop('Need "data" and/or "response" to calculate SE')
  } else {NA}

  x_dat <- new.data[, x_vars]

  intercept <- setdiff(estimates[, 1], x_vars)
  if (length(intercept) == 0) {
    int <- 0
  } else if (length(intercept) > 1) {
    stop('length(setdiff(estimates[, 1], x_vars)) > 1: cannot determine Y Intercept')
  } else {
    int <- estimates[which(estimates[, 1] == intercept), 2]
  }

  Yhat <- rep(NA, nrow(new.data))

  for (i in 1:nrow(new.data)) {
    ests <- c(NULL)
    X <- c(NULL)
    for (j in 1:length(x_vars)) {
      ests <- c(ests, estimates[which(estimates[, 1] == x_vars[j]), 2])
      X <- c(X, x_dat[i, x_vars[j]])
    }

    Yhat[i] <- int + sum(ests*X)
  }

  out <- data.frame(Y.hat = Yhat)

  if (se) {
    res <- avgm_res(estimates = estimates,
                    data = data,
                    response = response,
                    x_vars = x_vars,
                    plots = FALSE)

    mse <- sum(res$Residuals^2)/length(res$Residuals)

    se <- NULL

    X_n <- as.matrix(data[, x_vars])

    if (length(intercept) == 1) {
      X_n <- cbind(rep(1, nrow(X_n)), X_n)
    } else {NA}

    for (i in 1:length(Yhat)) {
      X_h <- X_n[i, ]
      err <- sqrt(mse + ((mse * (t(X_h) %*% solve(t(X_n) %*% X_n) %*% X_h)))^2)
      se <- c(se, err)
    }

    out <- data.frame(Y.hat = Yhat,
                      SE = se)
  }


  out <- cbind(out, new.data)

  return(out)
}
