
splitFilename <- function(names,n,e) {
  yy=list()
  for (name in names) {
    split_name=unlist(strsplit(name,split='[_.]'))
      n0=n
      if(n<0)
        n0=length(split_name)+1+n
      if(split_name[n0]==e)
        yy=append(yy,name)
  }
  return(unlist(yy))
}
