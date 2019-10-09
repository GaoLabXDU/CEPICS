#' Select features of dataset
#'
#' Using different strategies to select features of data.
#'
#' @param data A data matrix.
#' @param method A string means  the method of feature selection, choosing from "MAD", "VAR" and "PCA".
#' "MAD": using median absolute deviation to select features.
#' "VAR": select features according to the variance of each feature.
#' "PCA": using principal component analysis to select features.
#' @param select_number A real value from 0 to 1 which means the ratio of features we want to select. For "PCA", this parameter indicates the number of features to be selected instead of the ratio.
#'
#' @return
#' The data matrix after feature selection.
#' @examples
#' data(COAD_Methy)
#' COAD_Methy <- selectFeature(COAD_Methy, method = 'VAR', select_number = 0.7)
#'
#' @export selectFeature
selectFeature <- function(data, method = "VAR", select_number = 1) {
  dt <- data

  if(method == "MAD"){
    d_mad <- apply(dt, 1, mad)
    select_num <- select_number * nrow(dt)

    dt <- data[order(d_mad)[1:select_num],]
  }else if(method == "VAR") {
    d_var <- apply(dt, 1, var)
    select_num <- select_number * nrow(dt)

    dt <- data[order(d_var)[1:select_num],]
  }else if(method == "PCA") {
    d_pca <- prcomp(dt, rank. = select_number)
    d_pca <- t(d_pca$rotation)


    return(d_pca)
  }

  return(dt)
}
