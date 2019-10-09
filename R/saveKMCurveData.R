
saveKMCurveData <- function(clusters, data = NULL, fileName = NULL) {

  if(is.null(data)) {
    data <- KMSURDATA
  }

  clusters <- as.data.frame(clusters)
  clusters <- cbind(rownames(clusters), clusters)
  colnames(clusters)[1] <- 'name'

  finData <- merge(clusters, data)

  OS_DAYS<-finData$OS_Days
  OS_DAYS<-as.numeric(OS_DAYS)

  OS_STATUS<-finData$OS_Status
  OS_STATUS <- as.character(OS_STATUS)

  GROUP<-finData$clusters
  GROUP <- as.numeric(GROUP)

  misv <- 0
  Gender<-finData$Gender
  Gender <- as.numeric(Gender)
  if (length(Gender) == 0) {
    Gender <- rep(0, length(GROUP))
    misv <- misv + 1
  }

  Age<-finData$Age
  Age <- as.numeric(Age)
  if (length(Age) == 0) {
    Age <- rep(0, length(GROUP))
    misv <- misv + 1
  }


  maxGroup <- max(GROUP)

  if (misv == 2){
	GROUP <- as.factor(GROUP)
	cox = coxph(Surv(OS_DAYS,OS_STATUS=='Dead') ~GROUP)
  }else{
	cox = coxph(Surv(OS_DAYS,OS_STATUS=='Dead') ~GROUP + Gender + Age)
  }
  kmsurviaval<-survfit(Surv(OS_DAYS,OS_STATUS=='Dead') ~GROUP)
  kmData <- list(kmsurviaval = kmsurviaval, cox = cox, maxGroup = maxGroup, misv = misv)

  if (!is.null(fileName)) {
    save(kmData, file = as.character(fileName))
  }
  else {
    strGroup <- 'Group1'
    if (kmData$maxGroup >= 2) for (i in 2:kmData$maxGroup) {
      strGroup <-c(strGroup,paste('Group', i, sep=''))
    }
    palette(rainbow(kmData$maxGroup))
    plot(kmData$kmsurviaval,col = palette(),lty=1,lwd=2)

    if (misv == 2) pv <- summary(kmData$cox)$sctest[3]
    else  pv <- summary(kmData$cox)$coefficients[1,5]

    title(main=paste("P-value:", round(pv, digits = 5),"  ",sep=" "),xlab = "Days",ylab = "Survival Probability")
    legend(("topright"),strGroup,fill= palette(), inset = -0.08, xpd = TRUE)
  }

}
