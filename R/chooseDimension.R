
chooseDimension <- function(data,type,maxdimension=10, cores = 1)
{
  if(length(data)!=length(type))
  {
    print("the dimensions of data and type is different!!!")
    return()
  }
  r = array()
  r[1] = 0
  t = array()
  t[1] = 0
  jiang=list()
  # multi-core
  cl <- makeCluster(cores)
  names = as.character(1:length(data))
  clusterExport(cl, "names")
  res = parLapplyLB(cl,1:maxdimension,function(i) LRAcluster::LRAcluster(data, type, i, names), chunk.size=1)
  stopCluster(cl)
  for(i in 1:maxdimension){
    r[i+1] = res[[i]]$potential #ev
    jiang[[i]]=t(res[[i]]$coordinate)
  }

  n=selectev(r)
  x=jiang[[n]]

  rec = r

  return(list(x=x,n=n,rec=rec,jiang=jiang))
}
