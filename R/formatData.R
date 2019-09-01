#' Format Data
#'
#' Read data from a file.The first row and the first column are the names of the data.
#'
#' @param data File path of the data that needs to be read.
#' @param trans Whether you need to transpose the results.
#' @param sep The separator of the data in the file.
#'
#' @return
#' Data in matrix format read from file.
#'
#' @examples
#' data(COAD_Methy)
#' write.csv(COAD_Methy, 'temp.csv')
#' mydata <- formatData('temp.csv', sep = ',')
#'
#' @import pheatmap
#' @importFrom flexclust comPart
#' @import cluster
#' @import ggplot2
#' @import R.matlab
#' @import iClusterPlus
#' @import DMwR
#' @import survival
#' @import PINSPlus
#' @import LRAcluster
#' @import SNFtool
#' @import rmarkdown
#' @import stats
#' @import grid
#' @import reshape2
#' @import grDevices
#' @import graphics
#' @import utils
#' @import parallel
#'
#' @export formatData
formatData <- function(data, trans=FALSE ,sep=',') {

  data1 <- read.table(file = data, sep = sep, header = T, row.names = 1)


  if (trans == TRUE) data1 <- t(data1)

  dmat <- as.matrix(data1)
  return(dmat)
}
