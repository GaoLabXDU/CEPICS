
kmeansn <- function(x,mink=2,maxk=10)
{

  x=as.matrix(x)
  ans=matrix(0,nrow=maxk-mink+1,ncol=nrow(x))
  for( k in c(mink:maxk) )
  {
    ans[k-1,]=kmeans(x,k)[[1]]
  }
  rownames(ans)=c(mink:maxk)
  return(ans)
}
