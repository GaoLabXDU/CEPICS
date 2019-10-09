
writeHeader <- function(header,filename='report.Rmd') {
  if(file.exists(filename))
    file.remove(filename)
  writeRmd('---',filename,ishead=TRUE)
  writeRmd(paste('title: \"',header$title,'\"',sep=''),filename,ishead=TRUE)
  writeRmd(paste('author: \"',header$author,'\"',sep=''),filename,ishead=TRUE)
  writeRmd(paste('date: \"',header$date,'\"',sep=''),filename,ishead=TRUE)
  writeRmd(paste('output: \nhtml_document:\ndf_print: paged',sep=''),filename,ishead=TRUE)
  writeRmd('---\n\n\n\n',filename,ishead=TRUE)
}
