
writeCode <- function(name,Option,code='',filename='xuhan.Rmd',closechunk=TRUE) {
  s= paste('```{r',name)
  if(length(Option)>0)
  {

    for(e in c(1:length(Option)))
    {
      s= paste(s,',',names(Option[e]),'=',Option[[e]])
    }

  }
  s=paste(s,'}\n')
  s=paste(s,code)
  if(closechunk)
    s=paste(s,'\n```\n')
  else
    s=paste(s,'\n')
  writeRmd(s,filename)
}
