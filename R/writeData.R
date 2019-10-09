writeData <- function(startTime, data, funcName, clusNum, evaName) {

  strPic <- 'f'
  last <- 'RData'

  if (!dir.exists(startTime)) dir.create(startTime)

  fileName <- paste(strPic, funcName, clusNum, evaName, sep = '_')
  finfile <- paste(startTime, fileName, sep='/')
  finfile <- paste(finfile, last, sep='.')
  print(finfile)

  save(data, file = finfile)

}

getPath <- function(front, startTime, funcName, clusNum, evaName, last) {


  if (!dir.exists(startTime)) dir.create(startTime)

  fileName <- paste(front, funcName, clusNum, evaName, sep = '_')
  finfile <- paste(startTime, fileName, sep='/')
  finfile <- paste(finfile, last, sep='.')

  return(finfile)
}

writeToFileFromList <- function(data, pre, dir) {

  if (is.list(data)) {
    for (i in 1:length(data)) {
      tname <- names(data)[i]
      tname <- paste(pre, tname, sep = '_')
      writeToFileFromList(data[[i]], tname, dir)
    }
  }
  else {
    if (!dir.exists(dir)) dir.create(dir)

    if (!is.null(data)) {
      finfile <- paste(dir, pre, sep = '/')
      finfile <- paste(finfile, 'csv', sep = '.')

      write.csv(data, finfile)
    }

  }
}
