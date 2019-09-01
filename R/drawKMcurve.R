drawKMcurve <- function(kmData,tt) {
  strGroup <- 'Group1'
  if (kmData$maxGroup >= 2) for (i in 2:kmData$maxGroup) {
    strGroup <-c(strGroup,paste('Group', i, sep=''))
  }
  palette(rainbow(kmData$maxGroup))
  plot(kmData$kmsurviaval,col = palette(),lty=1,lwd=2)
  title(main=tt,xlab = "Days",ylab = "Survival Probability")
  legend(("topright"),strGroup,fill= palette(), inset = -0.08, xpd = TRUE)

  if (kmData$misv == 2) pv <- summary(kmData$cox)$sctest[3]
  else  pv <- summary(kmData$cox)$coefficients[1,5]

  legend("bottomleft",legend = paste("cox p:", round(pv, digits = 5),"  ",sep=" "), xpd = TRUE)
  p2=round(summary(kmData$cox)$sctest[3],digits = 5)

}
