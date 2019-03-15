
insertPicture <- function(s ,filename) {
  writeRmd(paste('![](',s,')\n' ,sep = '' ) ,filename=filename  )
}
