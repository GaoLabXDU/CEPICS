sampleClusterFig <- function(sampleClu) {

  sc <- sampleClu[[1]]$x
  if (length(sampleClu) >= 2) for (i in 2:length(sampleClu)) {
    sc <- sc + sampleClu[[1]]$x
  }
  write.csv(sc, file = 'sample_overall.csv');

  norm <- function(x) {
    nu <- x[1,1]
    xg <- x/(nu+1)
    xg <- 1 - xg
    xg <- as.dist(xg)
    return(xg)
  }
  scgy <- norm(sc)
  hc <- hclust(scgy)
  cls <- NULL
  kMax <- length(sampleClu)+1
  for (i in 2:kMax) {
    tc <- cutree(hc, i)
    cls <- cbind(cls, tc)
  }
  colnames(cls) <- 2:kMax
  write.csv(cls, file = 'sample_overall_cluster.csv')

  ord <- hc$order
  sco <- sc[ord,ord]

  nu <- sco[1,1]
  sco <- sco/nu
  drawHeatmapRmd(sco,F)

}
