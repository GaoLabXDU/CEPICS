#' Normalize Data
#'
#' Using different methods to normalize the data. CEPICS provides four normalization methods.
#'
#' @param data A data matrix which needs to be normalized.
#' @param method A string means the method of normalization, choose from "log", "z_score", "mean", "median".
#' "log": using log transformation to normalize the data, which means take the logarithm of each value to base X. X should be set as parameter "logBase".
#' "z-score": using z-score to normalize the data.
#' "mean": each value minus the mean value of its feature or sample.
#' "median": each value minus the median value of its feature or sample.
#'
#' @param object Choose of "feature" and "sample" which means normalize for each feature or sample respectively.
#' @param logBase The base of log when using log transformation to normalize data.
#' @examples
#' data(COAD_Methy)
#' COAD_Methy <- normData(COAD_Methy, method = 'z_score', object = 'feature')
#' @export normData
normData <- function(data, method = 'z_score', object = 'feature', logBase = 2) {

  if(object == "feature") data <- t(data)

  if (method == "log") {
    if (is.null(logBase)) dt <- log2(data)
    else dt <- log(data, base = logBase)
  }
  if(method == "z_score") {
    dt <- apply(data, 2, scale)
    row.names(dt) <- row.names(data)

  }else if(method == "mean") {
    dd <- apply(data, 2, mean, na.rm = T)
    dt <- data
    for (i in c(1:ncol(dt))) {
      dt[,i] <- dt[,i] - dd[i]
    }

  }else if(method == "median") {
    dd <- apply(data, 2, median, na.rm = T)
    dt <- data
    for (i in c(1:ncol(dt))) {
      dt[,i] <- dt[,i] - dd[i]
    }
  }
  if(object == "feature") dt <- t(dt)
  return(dt)
}
