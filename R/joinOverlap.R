joinOverlap <- function(datalist) {
  datalen <- length(datalist)
  if (datalen == 1) return(datalist)
  
  lapRownames <- rownames(datalist[[1]])
  for (i in 2:datalen) {
    trow <- rownames(datalist[[i]])
    lapRownames <- intersect(lapRownames, trow)
  }
  
  deleteRow <- function(data, rowname) {
    tnrow <- nrow(data)
    print(rownames(data)[1] %in% rowname)
    for (i in tnrow:1) {
      if (!(rownames(data)[i] %in% rowname)) {
        print(i)
        data <- data[,-i]
      }
    }
    return(data)
  }
  
  datalist <- lapply(datalist, deleteRow, rowname = lapRownames)
  return(datalist)
}