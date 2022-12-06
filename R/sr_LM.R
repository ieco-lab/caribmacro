#' Create a list of SR linear models
#'
#' This function is used to create a list of  linear models using the \code{lm()} function
#' with species richness (SR) as the response variable.
#' Requires the function \code{stndrd()} from the 'caribmacro' package
#'
#' @param SR A data frame in the format that is returned by \code{SR_geo()} that holds the
#' species richness values for each species group and geographic feature of interest.
#'
#' @param data A data frame that holds the explanatory variables for each geographic
#' feature of interest.
#'
#' @param x_vars A character or vector equal to the column name(s) in 'data' that hold
#' the explanatory variables of interest.
#'
#' @param geo_group A character equal to the column name in 'SR' and 'data' in which
#' the geographic feature names are stored. This needs to be equal for both the 'SR'
#' and 'data' data frames.
#'
#' @param area Optional. A character or vector equal to the column name in 'data' that
#' holds the area of each geographic feature. This argument is required if the
#' \code{log_area = TRUE} and/or \code{sq_area = TRUE}.
#'
#' @param standardize Logical. If \code{TRUE} then the explanatory variables are standardized by
#' mean centering and dividing by the standard deviation.
#'
#' @param log_sr Logical. If \code{TRUE} then the species richness is transformed by \code{log(SR+1)}
#'
#' @param log_area Logical. If \code{TRUE} then the geographic area is log transformed by \code{log(area)}
#'
#' @param sq_area Logical. If \code{TRUE} then the geographic area is squared (i.e. area^2).
#' NOTE: if \code{log_area = TRUE} then this option will square log(area).
#'
#' @param intercept Logical. If \code{TRUE} (default), then the y-intercept is estimated. If \code{FALSE},
#' then the y-intercept is set to 0.
#'
#' @param complete.case Logical. If \code{TRUE} only complete cases will be included in the models.
#' This argument should be set to TRUE if models will be passed to \code{AIC_analysis()}. Also, if
#' \code{TRUE} then the data used for each model will be exported in the list for use in \code{AIC_analysis()}.
#'
#' @examples
#' \dontrun{
#'
#' library(here)
#' library(reshape)
#'
#' dat <- read.csv(file.path(here(), 'data_raw', 'IBT_Herp_Records_v6.csv'), header=TRUE)
#' geo <- read.csv(file.path(here(), 'data_raw', 'bank_data.csv'), header=TRUE)
#'
#' coms <- com_matrix(species = "binomial",
#'                    geo_group = "bank",
#'                    taxa_group = "class",
#'                    status = "bnk_status",
#'                    stat_levels = c("N", "E"),
#'                    total = TRUE,
#'                    data = dat)
#'
#' SR <- SR_geo(data = coms, geo_group = "bank")
#'
#' ## Requires the function stndrd()
#' lm.sr <- sr_LM(SR = SR,
#'                data = geo,
#'                x_vars = c('area', 'isoPC1', 'isoPC2'),
#'                geo_group = 'bank',
#'                area = 'area',
#'                log_area = TRUE,
#'                sq_area = TRUE,
#'                complete.case = TRUE)
#'
#' }
#'
#' @export

sr_LM <- function(SR, data, x_vars, geo_group,
                  area = NULL,
                  standardize = TRUE,
                  log_sr = TRUE,
                  log_area = FALSE,
                  sq_area = FALSE,
                  intercept = TRUE,
                  complete.case = FALSE) {

  x_data <- data[, c(geo_group, x_vars)]

  if (log_sr) {
    SR <- cbind(SR[, geo_group], log((SR[, names(SR)[which(names(SR) != geo_group)]] + 1)))
    names(SR)[1] <- geo_group
  } else {NA}

  if (log_area) {
    x_data[, area] <- log(x_data[, area])
  } else {NA}

  if (standardize) {
    x_data <- stndrd(vars = names(x_data)[2:ncol(x_data)], data = x_data)
    SR <-  stndrd(vars = names(SR)[which(names(SR) != geo_group)], data = SR)
  } else {NA}

  if (sq_area) {
    nam <- paste('sq', area, sep = "_")
    x_data[, nam] <- x_data[, area]^2
  } else {NA}

  temp <- merge(SR, x_data, by = geo_group, all = TRUE)

  res <- names(SR)[which(names(SR) != geo_group)]
  vars <- names(x_data)[which(names(x_data) != geo_group)]

  out<-vector("list", length(res))
  mod.data <- vector("list", length(res))

  names(out) <- res
  names(mod.data) <- res

  for (i in 1:length(res)) {
    tmp <- temp[, c(res[i], vars)]

    if (complete.case) {
      tmp <- tmp[complete.cases(tmp), ]
    } else {NA}

    if (intercept) {
      RHS <- vars[1]
    } else {
      RHS <- paste('0', vars[1], sep = '+')
    }

    if (length(vars) > 1) {
      for (j in 2:length(vars)) {
        RHS <- paste(RHS, vars[j], sep = '+')
      }
    } else {NA}

    out[[res[i]]]<-lm(as.formula(paste(res[i], "~", RHS)), data = tmp)
    mod.data[[res[i]]] <- tmp

    rm(tmp)
  }

  if (complete.case) {
    out <- list(out, mod.data)
    names(out) <- c('Models', 'Data')
  } else {NA}

  return(out)
}
