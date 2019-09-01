drawPvalues <- function(method,k,ncol=1) {

  tryCatch({
    ans=matrix(0,length(method),length(k))
    rownames(ans)=method
    colnames(ans)=k
    for (m in method) {
      for(c in k)
      {
        load(paste('F_',m,'_',c,'_KM.RData', sep=''))
        if (kmData$misv == 2) pv <- summary(kmData$cox)$sctest[3]
        else  pv <- summary(kmData$cox)$coefficients[1,5]
        ans[m,as.character(c)]=round(pv,digits = 10)
      }
    }
    drawHeatmapRmd(list(t(ans)),k=k,x_continuous=TRUE,color=FALSE,ncol=ncol,x="Number of Clusters")
  },error=function(e){

  })
}
