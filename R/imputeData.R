#' Impute the missing values of Data.
#'
#' Using different methods to impute the missing values of the data. CEPIM provides three imputation methods.
#'
#' @param data The data matrix needs to be imputed.
#' @param method Choose from "knn", "mean", and "median". "knn": using K nearest neighbors methods to impute the missing value.
#' "mean": using the mean value of each feature to impute the missing value.
#' "median": using the median value of each feature to impute the missing value.

#' @param k An integer value which means the number of the nearest neighbors when using KNN method to impute missing values.
#' @param trans_zero_to_NA A logical value. If TRUE, all zeros in the data will be converted to NA.
#'
#' @examples
#' data(COAD_Methy)
#' COAD_Methy <- imputeData(COAD_Methy, method = 'mean')
#' @references
#' Troyanskaya,O. et al. (2001) Missing value estimation methods for DNA microarrays. Bioinformatics, 17, 520-525.
#'
#' @export imputeData
imputeData <- function(data, method = "knn", k = 10, trans_zero_to_NA = FALSE) {

  if(method == "mean") {
    imp <- function(data_imp) mean(data_imp, na.rm = T)
  }else if(method == "median") {
    imp <- function(data_imp) median(data_imp, na.rm = T)
  }

  dt <- data

  if (trans_zero_to_NA) {
    dt[dt == 0] <- NA
  }

  if (method == "knn") {
    dt <- knnImputation(dt, k)
  }else {
    for(i in c(1:nrow(dt))) {
      val_imp <- imp(dt[i,])
      for (j in c(1:ncol(dt))) {
        if (is.na(dt[i,j])) dt[i,j] <- val_imp
      }
    }
  }

  return(dt)
}
