
h <- function(s,n,filename) {
  hn=list('#','##','###','####','#####','######')
  writeRmd(paste(hn[[n]],s),filename)
}
