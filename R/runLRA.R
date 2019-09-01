#' Run LRAcluster to get the subtyping results.
#'
#' We use LRAcluster R package to apply dimension reduction and data integration, then use k-means to cluster the samples and get subtyping results.
#'
#' ev is a metric that provided by the author of LRAcluster to determine the dimension number of best low dimension. We determine the best low dimensional data by the metric in LRAcluster. When using LRAcluster R package, the user must decide the best dimension number manually according to the explained variation plot in the dimension reduction step so that LRAcluster cannot automatically get the integrated results. Therefore, we propose a method to choose the number of dimensions automatically.
#'
#' @param data A list of different omics data. Each data in data list should be format as a data matrix with rows representing features and  columns representing samples.
#' @param type A list of the types of different omics data.  Each element of the list can be one of "binary","gaussian","poisson"
#' @param maxdimension An integer value means the maximize dimension LRAcluster will tried from one, default 10.
#' @param maxk the max number of clustering we will try.
#' @param cores An integer value means the number of cores for parallel computing.
#' @return Return a list with the following elements:
#' \itemize{
#'    \item{\code{x} is the best low dimension answer}
#'    \item{\code{n} is best low dimension data we advise by ev values, the metric provide by author of LRAcluster.}
#'    \item{\code{rec} is all the ev values}
#'    \item{\code{result} is all of the low dimension data we tried from 1 to \code{maxdimension}}
#' }
#'
#' @examples
#' data(COAD_Methy)
#' data(COAD_miRNA)
#' data(COAD_mRNA)
#'
#' \donttest{
#' res=runLRA(data=list(COAD_Methy, COAD_miRNA, COAD_mRNA),
#'      type=list("gaussian","gaussian","gaussian"))
#'
#' res=runLRA(data=list(COAD_Methy, COAD_miRNA, COAD_mRNA),
#'      type=list("gaussian","gaussian","gaussian"),maxdimension=10,maxk=10,iskmeans=TRUE)
#' }
#' @references
#' Wu,D. et al. (2015) Fast dimension reduction and integrative clustering of multi-omics data using low-rank approximation: application to cancer molecular classification. BMC Genomics, 16, 1022.
#' @export
runLRA<- function(data,type,maxdimension=10,maxk=10, cores = 1)
{

  ans=chooseDimension(data,type,maxdimension, cores)
  clst=kmeansn(ans$x,2,maxk)
  colnames(clst)<-rownames(ans$x)
  clst <- t(clst)
  return(list(clst=clst,x=ans$x,n=ans$n,rec=ans$rec,jiang=ans$jiang))
}
