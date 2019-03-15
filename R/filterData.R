#' Filter Data
#'
#' Discard rows or columns whose NA value is greater than a certain percentage.
#'
#' Normally, alpha1 and alpha2 only need to use one of them.
#'
#' @param data The data needs to be filter.
#' @param alpha1 The ratio of NA values in a column is greater than alpha1 and the column will be discarded.
#' @param alpha2 The ratio of NA values in a row is greater than alpha1 and the row  will be discarded.
#' @param filter_zero_as_NA If this parameter is true, the 0 in the data will be treated as NA.
#'
#' @return
#' Filtered data.
#'
#' @examples
#' data(COAD_Methy)
#' COAD_Methy <- filterData(COAD_Methy, alpha1 = 0.8)
#'
#' @export filterData
filterData <- function(data,alpha1 = 1,alpha2 = 1, filter_zero_as_NA = FALSE) {

  dt <- data

  max_na_col <- floor(nrow(dt) * alpha1)
  for(i in c(ncol(dt):1)) {
    if (filter_zero_as_NA) nZero <- sum(dt[,i] == 0)
    else nZero <- 0
    if(sum(is.na(dt[,i]), nZero) >= max_na_col) dt <- dt[,-i]
  }

  max_na_row <- floor(ncol(dt) * alpha2)
  for(i in c(nrow(dt):1)) {
    if (filter_zero_as_NA) nZero <- sum(dt[i,] == 0)
    else nZero <- 0
    if(sum(is.na(dt[i,]), nZero) >= max_na_row) dt <- dt[-i,]
  }

  return(dt)
}
