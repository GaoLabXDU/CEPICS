
affinity <- function(data, K = 20, alpha = 0.5) {
  data <- (dist2(as.matrix(data),as.matrix(data)))^(1/2)
  data <- affinityMatrix(data, K, alpha)

  return(data)
}
