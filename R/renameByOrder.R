renamelable <- function(clusters) {
  for (i in 1:ncol(clusters)) {
    clusters[,i] <- renameByOrder(clusters[,i])
  }
  return(clusters)
}

renameByOrder <- function(cluster) {
  cluster <- paste0('a',cluster)
  uni_clu <- unique(cluster)

  for (j in 1:length(cluster)) {
    for (i in 1:length(uni_clu)) {
      if (uni_clu[i] == cluster[j]) {
        cluster[j] = i
      }
    }
  }
  return(as.numeric(cluster))
}
