
saveKMCurveData <- function(clusters, data = NULL, fileName = NULL) {

  if(is.null(data)) {
    data <- KMSURDATA
  }

  clusters <- as.data.frame(clusters)
  clusters <- cbind(rownames(clusters), clusters)
  colnames(clusters)[1] <- 'name'

  finData <- merge(clusters, data)
  colnames(finData)[2] <- 'clusters'
  colnames(finData)[3] <- 'OS_Days'
  colnames(finData)[4] <- 'OS_Status'


  OS_DAYS<-finData$OS_Days
  OS_DAYS<-as.numeric(OS_DAYS)

  OS_STATUS<-finData$OS_Status
  OS_STATUS <- as.character(OS_STATUS)

  GROUP<-finData$clusters
  GROUP <- as.numeric(GROUP)

  maxGroup <- max(GROUP)

  cox = coxph(Surv(OS_DAYS,OS_STATUS=='Dead') ~GROUP)
  kmsurviaval<-survfit(Surv(OS_DAYS,OS_STATUS=='Dead') ~GROUP)


  if (!is.null(fileName)) {
    kmData <- list(kmsurviaval = kmsurviaval, cox = cox, maxGroup = maxGroup)
    save(kmData, file = as.character(fileName))
  }
  else {
    strGroup <- 'Group1'
    if (kmData$maxGroup >= 2) for (i in 2:kmData$maxGroup) {
      strGroup <-c(strGroup,paste('Group', i, sep=''))
    }
    palette(rainbow(kmData$maxGroup))
    plot(kmData$kmsurviaval,col = palette(),lty=1,lwd=2)

    title(main=paste("cox p:", round(summary(kmData$cox)$sctest[3],digits = 5),"  ",sep=" "),xlab = "Days",ylab = "Survival Probability")
    legend(("topright"),strGroup,fill= palette(), inset = -0.08, xpd = TRUE)

    p2=round(summary(kmData$cox)$sctest[3],digits = 5)
  }

}
