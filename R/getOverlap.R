#' Get the overlap samples of  datasets.
#'
#' Select samples from different omics datasets that belong to the same patient cohort.
#'
#' @param datalist A list object. Different omics datasets should be gathered as a datalist.
#'
#' @examples
#' data(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' datalist <- list(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' datalist <- getOverlap(datalist)
#' @export getOverlap
getOverlap <- function(datalist) {
  datalen <- length(datalist)
  if (datalen == 1) return(datalist)

  lapColnames <- colnames(datalist[[1]])
  for (i in 2:datalen) {
    tcol <- colnames(datalist[[i]])
    lapColnames <- intersect(lapColnames, tcol)
  }

  deleteCol <- function(data, colname) {
    tncol <- ncol(data)
    for (i in tncol:1) {
      if (!(colnames(data)[i] %in% colname)) {
        data <- data[,-i]
        if (!is.matrix(data)) {
          data <- as.matrix(data)
          data <- t(data)
        }
      }
    }
    return(data)
  }

  if(length(lapColnames) == 0) {
    print('There are no overlap data in this datalist.')
    return()
  }
  datalist <- lapply(datalist, deleteCol, colname = lapColnames)
  return(datalist)
}
