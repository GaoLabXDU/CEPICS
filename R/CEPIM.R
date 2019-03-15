#' CEPIM, a comparison and evaluation platform for integration methods in cancer subtyping
#'
#' The main function of CEPIM R package Please refer to examples for details. When running CEPIM, you don't need to initialize the MATLAB server since  CEPIM will initialize it automatically.
#'
#' @param datalist A list of different omics data. Each data should be formatted as a data matrix with rows representing features and columns representing  samples.
#' @param datatype A list of string means the corresponding datatype of datalist. Each type can choose from "binary", "gaussian" or "poisson".
#' @param functionList A list of functions that need to be ran, choosing from "iCluster" (means iClusterBayes), "SNF", "LRA" (means LRAcluster), "PFA", and "PINS". We set SNF, LRA, iClusterBayes and PINS as default.evalList	The evaluation metrices that needs to be evaluated and compared, choosing from "KM", "SI", "Heatmap", "RI", "ARI", and "NMI". We set all the metrices as default.
#' @param evalList The evaluation method that needs to be ran. Should be some of "KM", "SI", "Heatmap", "RI", "ARI", "NMI".
#' @param kMax An integer value means the maximize number of clusters we will try.Res	A list contains the result of user's method. Please see Res for more details.
#' @param Res User's result.Put your data into \code{Res$data}, and put name into \code{Res$name}, put the time you cost into \code{Res$time}. If this parameter is not null, this data will appear at the top of the report.
#' @param KMdata The survival data of the sample. If this parameter is null, CEPIM will try to use built-in data including BRCA, COAD, KIRC, LUAD and LUSC downloaded from TCGA.
#' @param trueLabel A list of the true labels of the samples. If this parameter is null, CEPIM will compare different metrices between each two methods. Otherwise the results of each methods will be compared with the true labels. Please see turelabel for more details.
#' @param SNFPara Parameters of SNF. Should be a list of following parameters, \code{SNFPara$k}(default is one tenth of the sample size), \code{SNFPara$t}(default is 20). Please see runSNF for more details. If the parameter is null, CEPIM will use default values.
#' @param LRAPara Parameters of LRA. Should be a list of following parameter, \code{LRAPara$maxdimension} (default is 10) is the max number of dimension CEPIM will try to reduce. If the parameter is null, CEPIM will use default values.
#' @param iClusterPara Parameters of iClusterBayes. Should be a list of following parameters, \code{iClusterPara$n.burnin}(default is 1000), \code{iClusterPara$n.draw}(default is 1200), \code{iClusterPara$prior.gamma}(default is rep(0.1,6)), \code{iClusterPara$sdev}(default is 0.5), \code{iClusterPara$beta.var.scale}(default is 1), \code{iClusterPara$thin}(default is 1), \code{iClusterPara$pp.cutoff}(default is 0.5). Please see iClusterPlus R package for more details. If the parameter is null, CEPIM will use default values.
#' @param PINSPara Parameters of PINS. Should be a list of following parameter, \code{PINSPara$agreementCutoff}(default is 0.5) means agreement threshold to be considered consistent. If the parameter is null, CEPIM will use default values.
#' @param filename The file name of the report.
#' @param title The title of the report.
#' @param author The author of thereport.
#' @return Output the evaluation and comparion report and some experimental data files.
#' @examples
#' data(COAD_Methy)
#' data(COAD_miRNA)
#' data(COAD_mRNA)
#'
#'
#' datalist <- list(COAD_Methy,COAD_miRNA,COAD_mRNA)
#'
#' \donttest{
#' CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
#'     functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'),
#'     evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'), kMax=5)
#'
#' #Use Survival data
#' data(KMSURDATA)
#' CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
#'     functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'),
#'     evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'), kMax=5, KMdata = KMSURDATA)
#'
#' #Use your result
#' data(Res)
#'
#' CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
#'     functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'),
#'     evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'), kMax=5, Res = Res)
#'
#' #Use gold standard
#' data(truelabel)
#' CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
#'     functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'),
#'     evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'), kMax=5, trueLabel = truelabel)
#'
#' #Pass parameter to function
#' CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
#'     functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'),
#'     evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'),
#'     kMax=5, SNFPara = list(k = 8, t = 20))
#'
#' CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
#'     functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'),
#'     evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'),
#'     kMax=5, LRAPara = list(maxdimension=10))
#' }
#'
#'
#' @export CEPIM
CEPIM <- function(datalist, datatype, functionList = list('PINS', 'LRA', 'SNF', 'iCluster'),
                evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'), kMax=5, Res = NULL, KMdata=NULL, trueLabel = NULL, SNFPara=NULL, LRAPara=NULL, iClusterPara=NULL,
                  PINSPara=NULL, filename = 'CEPIM', title='Evaluation and Comparison Report',author="") {

  if (!('iCluster' %in% functionList | 'SNF' %in% functionList | 'PINS' %in% functionList | 'LRA' %in% functionList | 'PFA' %in% functionList)) {
    print('There are no methods to be implemented.')
    return()
  }
  if (!('Heatmap' %in% evalList | 'KM' %in% evalList | 'SI' %in% evalList | 'ARI' %in% evalList | 'RI' %in% evalList | 'NMI' %in% evalList)) {
    print('There are no evaluation indicators to be implemented.')
    return()
  }

  timeList <- NULL

  SI <- NULL
  NMI <- NULL
  RI <- NULL
  ARI <- NULL

  nogoldRIRes <- NULL
  nogoldARIRes <- NULL
  nogoldNMIRes <- NULL

  sampleClu <- NULL

  errorFunctionList <- NULL
  errorMessage <- NULL

  startTime <- as.character.Date(Sys.time())
  startTime <- gsub("\\D", '', startTime)
  startTime <- paste('CEPIM', startTime, sep = '_')
  if (!dir.exists(startTime)) dir.create(startTime)


  if (!is.null(Res)) {
    resName <- Res$name

    if (is.null(resName)) resName <- 'user'

    functionList <- c(resName, functionList)

    userEva <- evaluateRes(Res$data, funcName = resName, startTime = startTime, aff = T, evalList = evalList, trueLabel = trueLabel)

    writeToFileFromList(userEva, resName, startTime)

    if (!is.null(userEva$SIFinalRes)) {
      SI <- c(SI, list(user = userEva$SIFinalRes))
      names(SI)[length((SI))] <- resName
    }
    if (!is.null(userEva$NMIFinalRes)) {
      NMI <- c(NMI, list(user = userEva$NMIFinalRes))
      names(NMI)[length((NMI))] <- resName
    }
    if (!is.null(userEva$RIFinalRes)) {
      RI <- c(RI, list(user = userEva$RIFinalRes))
      names(RI)[length((RI))] <- resName
    }
    if (!is.null(userEva$ARIFinalRes)) {
      ARI <- c(ARI, list(user = userEva$ARIFinalRes))
      names(ARI)[length((ARI))] <- resName
    }

    if (!is.null(userEva$nogoldRIRes)) {
      nogoldRIRes <- c(nogoldRIRes, list(user = userEva$nogoldRIRes))
      names(nogoldRIRes)[length((nogoldRIRes))] <- resName
    }
    if (!is.null(userEva$nogoldARIRes)) {
      nogoldARIRes <- c(nogoldARIRes, list(user = userEva$nogoldARIRes))
      names(nogoldARIRes)[length((nogoldARIRes))] <- resName
    }
    if (!is.null(userEva$nogoldNMIRes)) {
      nogoldNMIRes <- c(nogoldNMIRes, list(user = userEva$nogoldNMIRes))
      names(nogoldNMIRes)[length((nogoldNMIRes))] <- resName
    }

    timeList <- c(timeList, Res$time)
    names(timeList)[length(timeList)] <- resName

  }

  if ('iCluster' %in% functionList) {
    print('iCluster')

    flag <- FALSE
    errorMes <- try({
      t1 <- Sys.time()
      if (is.null(iClusterPara)) {
         iCRes <- runiClusterBayes(datalist, kMax=kMax)


      }
      else {
        if (is.null(iClusterPara$n.burnin))iClusterPara$n.burnin = 1000
        if (is.null(iClusterPara$n.draw))iClusterPara$n.draw = 1200
        if (is.null(iClusterPara$prior.gamma))iClusterPara$prior.gamma = rep(0.1,6)
        if (is.null(iClusterPara$sdev))iClusterPara$sdev = 0.5
        if (is.null(iClusterPara$beta.var.scale))iClusterPara$beta.var.scale = 1
        if (is.null(iClusterPara$thin))iClusterPara$thin = 1
        if (is.null(iClusterPara$pp.cutoff))iClusterPara$pp.cutoff = 0.5

        iCRes <- runiClusterBayes(datalist, kMax=kMax, type = datatype,n.burnin=iClusterPara$n.burnin,n.draw=iClusterPara$n.draw,
                                  prior.gamma=iClusterPara$prior.gamma,sdev=iClusterPara$sdev,beta.var.scale=iClusterPara$beta.var.scale,
                                  thin=iClusterPara$thin,pp.cutoff=iClusterPara$pp.cutoff)
      }

      writeToFileFromList(iCRes, 'iCluster', startTime)

      iCEva <- evaluateRes(iCRes, funcName = 'iCluster', startTime = startTime, aff = F, evalList = evalList, trueLabel = trueLabel)

      writeToFileFromList(iCEva, 'iCluster', startTime)

      if (!is.null(iCEva$SIFinalRes)) SI <- c(SI, list(iCluster = iCEva$SIFinalRes))
      if (!is.null(iCEva$NMIFinalRes)) NMI <- c(NMI, list(iCluster = iCEva$NMIFinalRes))
      if (!is.null(iCEva$RIFinalRes)) RI <- c(RI, list(iCluster = iCEva$RIFinalRes))
      if (!is.null(iCEva$ARIFinalRes)) ARI <- c(ARI, list(iCluster = iCEva$ARIFinalRes))

      if (!is.null(iCEva$nogoldRIRes)) nogoldRIRes <- c(nogoldRIRes, list(Cluster = iCEva$nogoldRIRes))
      if (!is.null(iCEva$nogoldARIRes)) nogoldARIRes <- c(nogoldARIRes, list(iCluster = iCEva$nogoldARIRes))
      if (!is.null(iCEva$nogoldNMIRes)) nogoldNMIRes <- c(nogoldNMIRes, list(iCluster = iCEva$nogoldNMIRes))
      #output
      t2 <- Sys.time()
      timeList <- c(timeList, as.numeric(difftime(t2,t1, units = 'secs')))
      names(timeList)[length(timeList)] <- 'iCluster'
      flag <- TRUE
      })

    if (!flag) {
      tempMessage <- paste("Error in iCluster: ", errorMes, sep = '')
      errorMessage <- paste(errorMessage, tempMessage, sep = '')
      errorFunctionList <- c(errorFunctionList, 'iCluster')
      functionList <- functionList[functionList != 'iCluster']
    }

  }

  if ('SNF' %in% functionList) {
    print('SNF')

    flag <- FALSE
    errorMes <- try({
      t1 <- Sys.time()
      if (is.null(SNFPara)) {
        SNFRes <- runSNF(datalist, kMax = kMax)
      }
      else {
        if (is.null(SNFPara$t)) SNFPara$t = 20
        SNFRes <- runSNF(datalist,SNFPara$k, SNFPara$t, kMax)
      }

      nSNFRes <- list(cluster = SNFRes$clusters, data = SNFRes$affmat)

      SNFEva <- evaluateRes(nSNFRes, funcName = 'SNF', startTime = startTime, aff = T, evalList = evalList, trueLabel = trueLabel)

      writeToFileFromList(nSNFRes, 'SNF', startTime)
      writeToFileFromList(SNFEva, 'SNF', startTime)

      if (!is.null(SNFEva$SIFinalRes)) SI <- c(SI, list(SNF = SNFEva$SIFinalRes))
      if (!is.null(SNFEva$NMIFinalRes)) NMI <- c(NMI, list(SNF = SNFEva$NMIFinalRes))
      if (!is.null(SNFEva$RIFinalRes)) RI <- c(RI, list(SNF = SNFEva$RIFinalRes))
      if (!is.null(SNFEva$ARIFinalRes)) ARI <- c(ARI, list(SNF = SNFEva$ARIFinalRes))

      if (!is.null(SNFEva$nogoldRIRes)) nogoldRIRes <- c(nogoldRIRes, list(SNF = SNFEva$nogoldRIRes))
      if (!is.null(SNFEva$nogoldARIRes)) nogoldARIRes <- c(nogoldARIRes, list(SNF = SNFEva$nogoldARIRes))
      if (!is.null(SNFEva$nogoldNMIRes)) nogoldNMIRes <- c(nogoldNMIRes, list(SNF = SNFEva$nogoldNMIRes))
      #output
      t2 <- Sys.time()
      timeList <- c(timeList, as.numeric(difftime(t2,t1, units = 'secs')))
      names(timeList)[length(timeList)] <- 'SNF'
      flag <- TRUE
    })

    if (!flag) {
      tempMessage <- paste("Error in SNF: ", errorMes, sep = '')
      errorMessage <- paste(errorMessage, tempMessage, sep = '')
      errorFunctionList <- c(errorFunctionList, 'SNF')
      functionList <- functionList[functionList != 'SNF']
    }
  }

  if ('LRA' %in% functionList) {
    print('LRA')


    errorMes <- try({


      flag <- FALSE


      t1 <- Sys.time()
      if (is.null(LRAPara)) LRAPara <- list(maxdimension = 10)

      LRARes <- runLRA(datalist, type = datatype, LRAPara$maxdimension, kMax)


      nLRARes <- list(cluster = LRARes$clst, data = LRARes$x)
      writeToFileFromList(nLRARes, 'LRA', startTime)
      LRAEva <- evaluateRes(nLRARes, funcName = 'LRA', startTime = startTime, aff = F, evalList = evalList, trueLabel = trueLabel)
      writeToFileFromList(LRAEva, 'LRA', startTime)

      if (!is.null(LRAEva$SIFinalRes)) SI <- c(SI, list(LRA = LRAEva$SIFinalRes))
      if (!is.null(LRAEva$NMIFinalRes)) NMI <- c(NMI, list(LRA = LRAEva$NMIFinalRes))
      if (!is.null(LRAEva$RIFinalRes)) RI <- c(RI, list(LRA = LRAEva$RIFinalRes))
      if (!is.null(LRAEva$ARIFinalRes)) ARI <- c(ARI, list(LRA = LRAEva$ARIFinalRes))

      if (!is.null(LRAEva$nogoldRIRes)) nogoldRIRes <- c(nogoldRIRes, list(LRA = LRAEva$nogoldRIRes))
      if (!is.null(LRAEva$nogoldARIRes)) nogoldARIRes <- c(nogoldARIRes, list(LRA = LRAEva$nogoldARIRes))
      if (!is.null(LRAEva$nogoldNMIRes)) nogoldNMIRes <- c(nogoldNMIRes, list(LRA = LRAEva$nogoldNMIRes))
      #output
      t2 <- Sys.time()
      timeList <- c(timeList, as.numeric(difftime(t2,t1, units = 'secs')))
      names(timeList)[length(timeList)] <- 'LRA'
      flag <- TRUE

    })

    if (!flag) {
      message(errorMes)
      tempMessage <- paste("Error in LRA: ", errorMes, sep = '')
      errorMessage <- paste(errorMessage, tempMessage, sep = '')
      errorFunctionList <- c(errorFunctionList, 'LRA')
      functionList <- functionList[functionList != 'LRA']
    }

  }

  if ('PFA' %in% functionList) {
    print('PFA')


    flag <- FALSE
    errorMes <- try({
      m <- PFAinit(getwd())
      t1 <- Sys.time()

      PFARes <- runPFA(datalist, maxk = kMax, matlab = m)


      nPFARes <- list(cluster = PFARes$clst, data = PFARes$ans)
      nPFARes$data <- Re(nPFARes$data)

      PFAEva <- evaluateRes(nPFARes, funcName = 'PFA', startTime = startTime, aff = F, evalList = evalList, trueLabel = trueLabel)

      writeToFileFromList(nPFARes, 'PFA', startTime)
      writeToFileFromList(PFAEva, 'PFA', startTime)

      if (!is.null(PFAEva$SIFinalRes)) SI <- c(SI, list(PFA = PFAEva$SIFinalRes))
      if (!is.null(PFAEva$NMIFinalRes)) NMI <- c(NMI, list(PFA = PFAEva$NMIFinalRes))
      if (!is.null(PFAEva$RIFinalRes)) RI <- c(RI, list(PFA = PFAEva$RIFinalRes))
      if (!is.null(PFAEva$ARIFinalRes)) ARI <- c(ARI, list(PFA = PFAEva$ARIFinalRes))

      if (!is.null(PFAEva$nogoldRIRes)) nogoldRIRes <- c(nogoldRIRes, list(PFA = PFAEva$nogoldRIRes))
      if (!is.null(PFAEva$nogoldARIRes)) nogoldARIRes <- c(nogoldARIRes, list(PFA = PFAEva$nogoldARIRes))
      if (!is.null(PFAEva$nogoldNMIRes)) nogoldNMIRes <- c(nogoldNMIRes, list(PFA = PFAEva$nogoldNMIRes))
      #output
      t2 <- Sys.time()
      timeList <- c(timeList, as.numeric(difftime(t2,t1, units = 'secs')))
      names(timeList)[length(timeList)] <- 'PFA'

      flag <- TRUE
    })

    close(m)
    delfile("Algorithm1.m")
    delfile("Algorithm2.m")
    delfile("Algorithm4.m")
    delfile("computeerr.m")
    delfile("MainPFA.m")
    delfile("FindKMinEigen.m")
    delfile("InputStreamByteWrapper.class")
    delfile("MatlabServer.m")

    if (!flag) {
      tempMessage <- paste("Error in PFA: ", errorMes, sep = '')
      errorMessage <- paste(errorMessage, tempMessage, sep = '')
      errorFunctionList <- c(errorFunctionList, 'PFA')
      functionList <- functionList[functionList != 'PFA']
    }



  }

  if ('PINS' %in% functionList) {
    print('PINS')
    flag <- FALSE

    errorMes <- try({
      t1 <- Sys.time()
      if (is.null(PINSPara)) {
        PINSLabel <- runPINSPlus(datalist)
      }
      else PINSLabel <- runPINSPlus(datalist, kMax = kMax, agreementCutoff = PINSPara$agreementCutoff)

      PINSEvalist <- evalList[evalList != 'SI']
      PINSEvalist <- PINSEvalist[PINSEvalist != 'Heatmap']

      PINSEva <- evaluateRes(list(cluster = PINSLabel), funcName = 'PINS', startTime = startTime, aff = T, evalList = PINSEvalist, trueLabel = trueLabel)

      writeToFileFromList(PINSLabel, 'PINS_cluster', startTime)
      writeToFileFromList(PINSEva, 'PINS', startTime)

      if (!is.null(PINSEva$SIFinalRes)) SI <- c(SI, list(PINS = PINSEva$SIFinalRes))
      if (!is.null(PINSEva$NMIFinalRes)) NMI <- c(NMI, list(PINS = PINSEva$NMIFinalRes))
      if (!is.null(PINSEva$RIFinalRes)) RI <- c(RI, list(PINS = PINSEva$RIFinalRes))
      if (!is.null(PINSEva$ARIFinalRes)) ARI <- c(ARI, list(PINS = PINSEva$ARIFinalRes))

      if (!is.null(PINSEva$nogoldRIRes)) nogoldRIRes <- c(nogoldRIRes, list(PINS = PINSEva$nogoldRIRes))
      if (!is.null(PINSEva$nogoldARIRes)) nogoldARIRes <- c(nogoldARIRes, list(PINS = PINSEva$nogoldARIRes))
      if (!is.null(PINSEva$nogoldNMIRes)) nogoldNMIRes <- c(nogoldNMIRes, list(PINS = PINSEva$nogoldNMIRes))

      t2 <- Sys.time()
      timeList <- c(timeList, as.numeric(difftime(t2,t1, units = 'secs')))
      names(timeList)[length(timeList)] <- 'PINS'

      flag <- TRUE
    })

    if (!flag) {
      tempMessage <- paste("Error in PINS: ", errorMes, sep = '')
      errorMessage <- paste(errorMessage, tempMessage, sep = '')
      errorFunctionList <- c(errorFunctionList, 'PINS')
      functionList <- functionList[functionList != 'PINS']
    }
  }






  funcNMI <- NULL
  funcRI <- NULL
  funcARI <- NULL
  for (i in 2:kMax) {
    labelMatrix <- NULL
    if (!is.null(Res)){
      labelMatrix <- cbind(labelMatrix, Res$data$cluster[,(i-1)])
      colnames(labelMatrix)[(ncol(labelMatrix))] = resName
    }
    if ('SNF' %in% functionList){
      labelMatrix <- cbind(labelMatrix, SNFRes$clusters[,(i-1)])
      colnames(labelMatrix)[(ncol(labelMatrix))] = 'SNF'
    }
    if ('LRA' %in% functionList){
      labelMatrix <- cbind(labelMatrix, LRARes$clst[,(i-1)])
      colnames(labelMatrix)[(ncol(labelMatrix))] = 'LRA'
    }
    if ('iCluster' %in% functionList){
      labelMatrix <- cbind(labelMatrix, iCRes$cluster[,(i-1)])
      colnames(labelMatrix)[(ncol(labelMatrix))] = 'iCluster'
    }
    if ('PFA' %in% functionList){
      labelMatrix <- cbind(labelMatrix, nPFARes$cluster[,(i-1)])
      colnames(labelMatrix)[(ncol(labelMatrix))] = 'PFA'
    }

    if ('PINS' %in% functionList){
      if (colnames(PINSLabel)[1] == i) {
        labelMatrix <- cbind(labelMatrix, PINSLabel[,1])
        colnames(labelMatrix)[(ncol(labelMatrix))] = 'PINS'
      }
      else if (ncol(PINSLabel) == 2) {
        if (colnames(PINSLabel)[2] == i) {
          labelMatrix <- cbind(labelMatrix, PINSLabel[,2])
          colnames(labelMatrix)[(ncol(labelMatrix))] = 'PINS'
        }
      }
    }
    if (is.null(labelMatrix)) next

    colnum <- ncol(labelMatrix)
    rownum <- nrow(labelMatrix)
    tfuncNMI <- matrix(nrow = colnum, ncol = colnum)
    tfuncRI <- matrix(nrow = colnum, ncol = colnum)
    tfuncARI <- matrix(nrow = colnum, ncol = colnum)

    for (indi in 1:colnum) {
      for (indj in 1:colnum) {
        tfuncNMI[indi, indj] <- calNMI(labelMatrix[,indi], labelMatrix[,indj])
        tfuncRI[indi, indj] <- calRI(labelMatrix[,indi], labelMatrix[,indj])
        tfuncARI[indi, indj] <- calARI(labelMatrix[,indi], labelMatrix[,indj])
      }
    }

    colnames(tfuncNMI) <- colnames(labelMatrix)
    rownames(tfuncNMI) <- colnames(labelMatrix)
    colnames(tfuncRI) <- colnames(labelMatrix)
    rownames(tfuncRI) <- colnames(labelMatrix)
    colnames(tfuncARI) <- colnames(labelMatrix)
    rownames(tfuncARI) <- colnames(labelMatrix)

    funcNMI <- c(funcNMI, list(tfuncNMI))
    names(funcNMI)[i-1] <- i
    funcRI <- c(funcRI, list(tfuncRI))
    names(funcRI)[i-1] <- i
    funcARI <- c(funcARI, list(tfuncARI))
    names(funcARI)[i-1] <- i

    tsampleClu <- matrix(0, nrow = rownum, ncol = rownum)
    for (indi in 1:rownum) {
      for (indj in 1:rownum) {
        for (indc in 1:colnum) {
          if (labelMatrix[indi,indc] == labelMatrix[indj, indc])
            tsampleClu[indi, indj] <- tsampleClu[indi, indj] + 1
        }
      }
    }
    colnames(tsampleClu) <- rownames(labelMatrix)
    rownames(tsampleClu) <- rownames(labelMatrix)
    sampleClu <- c(sampleClu, list(list(x = tsampleClu)))
    names(sampleClu)[length(sampleClu)] <- i
  }


  nogold_method <- list(RI = funcRI, ARI = funcARI, NMI = funcNMI)
  save(nogold_method, file = as.character(paste(startTime, 'nogold_method.RData', sep='/')))

  if (!is.null(SI)) save(SI, file = as.character(paste(startTime, 'SI.RData', sep='/')))
  if (!is.null(NMI)) save(NMI, file = as.character(paste(startTime, 'NMI.RData', sep='/')))
  if (!is.null(RI)) save(RI, file = as.character(paste(startTime, 'RI.RData', sep='/')))
  if (!is.null(ARI)) save(ARI, file = as.character(paste(startTime, 'ARI.RData', sep='/')))


  if (!is.null(nogoldRIRes)) save(nogoldRIRes, file = as.character(paste(startTime, 'nogoldRIRes.RData', sep='/')))
  if (!is.null(nogoldARIRes)) save(nogoldARIRes, file = as.character(paste(startTime, 'nogoldARIRes.RData', sep='/')))
  if (!is.null(nogoldNMIRes)) save(nogoldNMIRes, file = as.character(paste(startTime, 'nogoldNMIRes.RData', sep='/')))


  if (!is.null(timeList)) {
    timeList <- timeList[names(timeList) %in% functionList]
    save(timeList, file = as.character(paste(startTime, 'timeList.RData', sep='/')))
    writeToFileFromList(timeList, 'timeList', startTime)
    drawTimeFig(timeList = timeList, dir = startTime)
  }

  if (!is.null(sampleClu)) save(sampleClu, file = as.character(paste(startTime, 'sampleClu.RData', sep='/')))
  writeToFileFromList(sampleClu, 'sampleClu', startTime)

  print(date())

  pic_eval <- NULL
  if ('Heatmap' %in% evalList) pic_eval <- c(pic_eval, 'HeatMap')
  if ('KM' %in% evalList) pic_eval <- c(pic_eval, 'KM')

  numeric_eval <- NULL
  if ('NMI' %in% evalList) numeric_eval <- c(numeric_eval, 'NMI')
  if ('ARI' %in% evalList) numeric_eval <- c(numeric_eval, 'ARI')
  if ('SI' %in% evalList) numeric_eval <- c(numeric_eval, 'SI')
  if ('RI' %in% evalList) numeric_eval <- c(numeric_eval, 'RI')

  if (is.null(trueLabel)) is_gold <- FALSE
  else is_gold <- TRUE


  outputReport(method=functionList,pic_eval=pic_eval,numeric_eval=numeric_eval,
               need_gold_eval=intersect(numeric_eval, c('NMI','ARI', 'RI')),is_gold=is_gold,n_pic_row=2,kmax=kMax,path=startTime,
               filename=paste(filename, 'Rmd', sep = '.'), title = title, author = author)

  if (!is.null(errorMessage)) {
    message(errorMessage)
  }
  message('Finish!')
}

drawTimeFig <- function(timeList, dir) {
  # Create data for the graph.
  x <- timeList
  myLabel <- as.vector(names(timeList))

  myLabel = paste(myLabel, " (", round(x / sum(x) * 100, 2), "%)", sep = "")

  svg(filename = as.character(paste(dir, 'timePie.svg', sep = '/')), width = 5, height = 5)


  timeList <- as.data.frame(timeList)
  timename <- rownames(timeList)
  FunctionName <- paste(timename, "(", round(timeList$timeList / sum(timeList$timeList) * 100, 2), "%)  ", sep = "")
  gp = ggplot(data=timeList, mapping=aes(x="",y = timeList ,fill=FunctionName))+
    geom_bar(stat="identity")+coord_polar(theta = "y") + labs(x = "", y = "", title = "")+theme(axis.ticks = element_blank()) + theme(axis.text.x = element_blank())+theme(panel.grid=element_blank())+theme(panel.border=element_blank())+theme(axis.line = element_blank())

  print(gp)

  dev.off(which = dev.cur())

}

renameCluRes <- function(clu) {
  fac <- factor(clu)
  fac <- levels(fac)
  for (i in 1:length(fac)) {
    for (j in 1:length(clu)) {
      if (fac[i] == clu[j]) clu[j] = i
    }
  }
  return(clu)
}

delfile <- function(x) {
  if(file.exists(x))
    file.remove(x)
  return()
}
