sampleClusterFig <- function(sampleClu) {
  
  sc <- sampleClu[[1]]$x
  if (length(sampleClu) >= 2) for (i in 2:length(sampleClu)) {
    sc <- sc + sampleClu[[1]]$x
  }
  norm <- function(x) {
    nu <- x[1,1]
    xg <- x/(nu+1)
    xg <- 1 - xg
    xg <- as.dist(xg)
    return(xg)
  }
  scgy <- norm(sc)
  hc <- hclust(scgy)
  ord <- hc$order
  sco <- sc[ord,ord]
  
  nu <- sco[1,1]
  sco <- sco/nu
  drawHeatmapRmd(sco,F)
  
}