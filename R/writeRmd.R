
writeRmd <- function(s,filename='report.Rmd',ishead=FALSE) {
  if(ishead)
  {
    cat(s,file = filename,sep='\n',append = TRUE)
  }
  else
  {
    cat(paste(s,'\n',sep=''),file = filename,sep='\n',append = TRUE)
  }
}
