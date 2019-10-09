evaluateRes <- function(res, funcName, startTime, evalList, aff = FALSE, KMdata=NULL, trueLabel = NULL) {

  cluster <- res$cluster
  rawdata <- res$data
  rawdatas <- res$datas

  SIFinalRes <- NULL
  RIFinalRes <- NULL
  ARIFinalRes <- NULL
  NMIFinalRes <- NULL

  nogoldRIRes <- NULL
  nogoldARIRes <- NULL
  nogoldNMIRes <- NULL

  if ('KM' %in% evalList) {
    print(cluster)

    for (i in 1:ncol(cluster)) {
      if (!is.null(rawdatas)) data <- rawdatas[[i]]
      else data <- rawdata


      if (!aff) {
        data <- affinity(data)
      }

      picFilePath <- getPath(front = 'F', startTime = startTime, funcName = funcName, clusNum = colnames(cluster)[i], evaName = 'KM', last = 'RData')
      saveKMCurveData(cluster[,i], data = KMdata, fileName = picFilePath)
    }

  }

  SIAns <- NULL
  if ('SI' %in% evalList) {
    #the value of type is 0,1 or 2, respectly means row data matrix, similar matrix , distance matrix
    for (i in 1:ncol(cluster)) {
      if (!is.null(rawdatas)) data <- rawdatas[[i]]
      else data <- rawdata


      if (aff)
        SIRes <- calSI(data, cluster[,i],type = 1)
      else
        SIRes <- calSI(data, cluster[,i],type = 0)

      SIAns <- c(SIAns, SIRes)
    }
    names(SIAns) <- colnames(cluster)
    SIFinalRes <- SIAns
  }

  if ('Heatmap' %in% evalList) {
    for (i in 1:ncol(cluster)) {
      if (!is.null(rawdatas)) data <- rawdatas[[i]]
      else data <- rawdata

      if (!aff) data <- affinity(data)

      picFilePath <- getPath(front = 'F', startTime = startTime, funcName = funcName, clusNum = colnames(cluster)[i], evaName = 'heatmap', last = 'RData')
      drawHeatmap(list(affmat = data, clusters = cluster[,i]), fileName = picFilePath)
    }

  }

  if ('RI' %in% evalList) {
    if (!is.null(trueLabel)) {
     for (i in 1:ncol(cluster)) {
       RIRes <- calRI(cluster[,i], trueLabel)
       RIFinalRes <- c(RIFinalRes, RIRes)
      }
      names(RIFinalRes) <- colnames(cluster)
    }

    nogoldRIRes <- matrix(nrow = ncol(cluster), ncol = ncol(cluster))
    for (i in 1:ncol(cluster)) {
      for (j in 1:ncol(cluster)) {
        nogoldRIRes[i,j] = calRI(cluster[,i], cluster[,j])
      }
    }
    rownames(nogoldRIRes) <- colnames(cluster)
    colnames(nogoldRIRes) <- colnames(cluster)
  }


  if ('ARI' %in% evalList) {
    if (!is.null(trueLabel)) {
       for (i in 1:ncol(cluster)) {
         ARIRes <- calARI(cluster[,i], trueLabel)
          ARIFinalRes <- c(ARIFinalRes, ARIRes)
        }
      names(ARIFinalRes) <- colnames(cluster)
    }

    nogoldARIRes <- matrix(nrow = ncol(cluster), ncol = ncol(cluster))
    for (i in 1:ncol(cluster)) {
      for (j in 1:ncol(cluster)) {
        nogoldARIRes[i,j] = calARI(cluster[,i], cluster[,j])
      }
    }
    rownames(nogoldARIRes) <- colnames(cluster)
    colnames(nogoldARIRes) <- colnames(cluster)
  }


  if ('NMI' %in% evalList) {
    if (!is.null(trueLabel)) {
        for (i in 1:ncol(cluster)) {
          NMIRes <- calNMI(cluster[,i], trueLabel)
          NMIFinalRes <- c(NMIFinalRes, NMIRes)
       }
      names(NMIFinalRes) <- colnames(cluster)
    }

    nogoldNMIRes <- matrix(nrow = ncol(cluster), ncol = ncol(cluster))
    for (i in 1:ncol(cluster)) {
      for (j in 1:ncol(cluster)) {
        nogoldNMIRes[i,j] = calNMI(cluster[,i], cluster[,j])
      }
    }
    rownames(nogoldNMIRes) <- colnames(cluster)
    colnames(nogoldNMIRes) <- colnames(cluster)
  }

  return(list(SIFinalRes = SIFinalRes, RIFinalRes = RIFinalRes, ARIFinalRes = ARIFinalRes, NMIFinalRes = NMIFinalRes,
              nogoldRIRes = nogoldRIRes, nogoldARIRes = nogoldARIRes, nogoldNMIRes = nogoldNMIRes))
}
