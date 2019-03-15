drawPvalues <- function(method,k,ncol=1) {

  tryCatch({
    ans=matrix(0,length(method),length(k))
    rownames(ans)=method
    colnames(ans)=k
    for (m in method) {
      for(c in k)
      {
        load(paste('F_',m,'_',c,'_KM.RData', sep=''))
        ans[m,as.character(c)]=round(summary(kmData$cox)$sctest[3],digits = 10)
      }
    }
    drawHeatmapRmd(list(t(ans)),k=k,x_continuous=TRUE,color=FALSE,ncol=ncol,x="Number of Clusters")
  },error=function(e){

  })
}
