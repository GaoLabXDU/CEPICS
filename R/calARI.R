#'Calculate Adjusted Rand Index
#'
#'Calculate the Adjusted Rand Index of the clustering results.
#'
#'
#'
#' @param a A vector that contains the cluster label of each sample.
#' @param b A vector that contains the true label of each sample.
#' @return A real number that means the adjusted rand index between the two clustering results.
#' @examples
#' a=c(1,3,2,2,4,6,5)
#' b=c(3,1,2,5,6,4,3)
#' calARI(a,b)
#' @export
calARI <- function(a,b)
{
  a=as.vector(a)
  b=as.vector(b)
  return(flexclust::comPart(a,b,type=c('ARI')))
}
