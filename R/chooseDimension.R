
chooseDimension <- function(data,type,maxdimension=10)
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
  for(i in 1:maxdimension){
    ins = paste(i, "calculating...", sep=" ")
    print(ins)
    res = LRAcluster(data, type, dimension = i, names = as.character(1:length(data)))
    r[i+1] = res$potential #ev
    jiang[[i]]=t(res$coordinate)
  }


  n=selectev(r)
  x=jiang[[n]]

  rec = r

  return(list(x=x,n=n,rec=rec,jiang=jiang))
}
