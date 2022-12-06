#' Run AIC Model Averaging on a Model
#'
#' This function is used run all combinations of explanatory variables and rank them based on their AIC.
#' Then the function averages the coefficients of the models that hold a certain amount of the model weight
#' This function requires the packages 'AICcmodavg' and 'MuMIn'
#'
#' @param model A model object from \code{lm()} or \code{glm()}.
#'
#' @param data The data frame that holds the data used in \code{model}.
#'
#' @param groups A data frame with 1 row that holds the groups names. For example, for introduced Anolis lizards
#' \code{groups = data.frame(Genus = 'Anolis', Status = 'I')}
#'
#' @param cum.weight A number equal to the cumulative weight of the models to be averaged.
#'
#' @param table Logical. If true a list is returned with the averaged coefficients in the first element and the
#' AIC tables for each full model in the second element.
#'
#' @param std Value passed to \code{dredge} and \code{model.avg} to indicate whether and how the coefficients are standardized,
#' and must be one of \code{"none"}, \code{"sd"} or \code{"partial.sd"}. \code{"sd"} standardizes the coefficients by the standard
#' deviation, and \code{"partial.sd"} standardizes the coefficients by the partial standard deviation (recommended if multicollinearity
#' is present among the predictors).
#'
#' @return A data frame with the estimates and their 95 percent confidence interval of the averaged models for
#' each model in \code{x}. If \code{table = TRUE} then a list is returned with the second slot holding the
#' complete AIC tables for each model in \code{x}.
#'
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
#' }
#'
#' @export

AIC_avg <- function(model, data = NULL, groups = NULL, response = NULL, x_vars = NULL,
                    cum.weight = 0.95, table = FALSE, std = 'none') {
  options(na.action = "na.fail")

  if (!is.null(data)) {
    if (!is.null(response) & !is.null(x_vars)) {
      res <- response
      vars <- x_vars

      dat <- data[, c(response, x_vars)]
      dat <- droplevels(dat[complete.cases(dat), ])

      RHS <- vars[1]
      for (j in 2:length(vars)) {
        RHS <- paste(RHS, vars[j], sep = '+')
      }

      mod <- update(model, as.formula(paste(res, "~", RHS)), data = dat)
    } else {
      stop('response and/or x_vars not Provided')
    }
  } else {
    mod <- model
  }

  if (table) {
    table <- list()

    dd <- dredge(mod, beta = std)

    top <- get.models(dd, subset = cumsum(weight) <= cum.weight)

    if (length(top) > 0) {
      avgm <- model.avg(top, beta = std)
      ci <- confint(avgm)
      coeffs <- as.data.frame(coefTable(avgm, full = TRUE))[, c('Estimate', 'Std. Error')]

      temp <- data.frame(coeffs, ci)
      names(temp) <- c('Estimate', 'SE', 'LCL', 'UCL')
      temp$Variable <- row.names(temp)
      row.names(temp) <- NULL

      if (!is.null(groups)) {
        for (k in 1:ncol(groups)) {
          temp[, names(groups)[k]] <- groups[1, k]
        }
        temp <- temp[, c(names(groups), "Variable", "Estimate", 'SE', "LCL", "UCL")]
      } else {NA}

    } else {
      temp <- data.frame(Estimate = NA, SE = NA, LCL = NA, UCL = NA, Variable = NA)

      if (!is.null(groups)) {
        for (k in 1:ncol(groups)) {
          temp[, names(groups)[k]] <- groups[1, k]
        }
        temp <- temp[, c(names(groups), "Variable", "Estimate", 'SE', "LCL", "UCL")]
      } else {NA}
    }

    out <- list(temp, dd)
    names(out) <- c('Coefficients', 'AIC Table')
  } else {
    dd <- dredge(mod, beta = std)

    top <- get.models(dd, subset = cumsum(weight) <= cum.weight)
    if (length(top) > 0) {
      avgm <- model.avg(top, beta = std)
      ci <- confint(avgm)
      coeffs <- as.data.frame(coefTable(avgm, full = TRUE))[, c('Estimate', 'Std. Error')]

      temp <- data.frame(coeffs, ci)
      names(temp) <- c('Estimate', 'SE', 'LCL', 'UCL')
      temp$Variable <- row.names(temp)
      row.names(temp) <- NULL

      if (!is.null(groups)) {
        for (k in 1:ncol(groups)) {
          temp[, names(groups)[k]] <- groups[1, k]
        }
        temp<-temp[, c(names(groups), "Variable", "Estimate", 'SE', "LCL", "UCL")]
      } else {NA}

    } else {
      temp <- data.frame(Estimate = NA, SE = NA, LCL = NA, UCL = NA, Variable = NA)

      if (!is.null(groups)) {
        for (k in 1:ncol(groups)) {
          temp[, names(groups)[k]] <- groups[1, k]
        }
        temp <- temp[, c(names(groups), "Variable", "Estimate", 'SE', "LCL", "UCL")]
      } else {NA}
    }

    out <- temp
  }

  return(out)

  options(na.action = "na.omit")
}
