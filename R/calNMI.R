#'Calculate the Normalized Mutual Information
#'
#'Calculate the Normalized Mutual Information of the clustering results.
#'
#'
#'
#' @param a A vector that contains the cluster label of each sample.
#' @param b A vector that contains the true label of each sample.
#' @return A real number that means the normalized mutual information between the two clustering results.
#' @examples
#' a=c(1,3,2,2,4,6,5)
#' b=c(3,1,2,5,6,4,3)
#' calNMI(a,b)
#' @export
calNMI <- function(a,b)
{
  len_a=length(a)
  xx=matrix(1,len_a,1)
  pa=aggregate(xx,list(a),sum)[,2]/len_a
  pb=aggregate(xx,list(b),sum)[,2]/len_a
  pab=aggregate(xx,list(a,b),sum)[3]/len_a

  ha=-sum(pa*log(pa))
  hb=-sum(pb*log(pb))
  hab=-sum(pab*log(pab))

  nmi=2*(ha+hb-hab)/(ha+hb)
  return(nmi)
}
