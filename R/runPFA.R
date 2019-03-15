#'Run PFA method to get the subtyping results.
#'
#'Using PFA to apply dimension reduction and data integration, then using k-means to get the cluster results.
#'
#' The author just provides the MATLAB code of the method. When we tried to transfer it to R, we found that the results of left division operation for non-square matrices in R and MATLAB were totally different. In order to compare and evaluate the original version of PFA, CEPIM uses R.matlab R package to initialize a MATLAB server to run the original MATLAB code. In addition, the original version code can only integrate three omics datasets. We improve the code and make it possible to integrate omics data without the number of datasets limit. You should execute \code{PFAinit()} to initialize a MATLAB server before running PFA, and close server after running.
#'
#' @param data A list of different omics data. Each data in data list should be format as a data matrix with rows representing features and columns representing samples.
#' @param maxk An integer value means the maximize number of clusters we will try.
#' @param matlab MATLAB server variable when execute \code{PFAinit()} function.
#' @examples
#' data(COAD_Methy)
#' data(COAD_miRNA)
#' data(COAD_mRNA)
#'
#' \donttest{
#' matlab=PFAinit(path=getwd())
#' res=runPFA(data=list(COAD_Methy, COAD_miRNA, COAD_mRNA),maxk=10,matlab=matlab)
#' close(matlab)
#' }
#'
#' @seealso
#' \code{\link{PFAinit}}
#'
#' @references
#' Shi,Q. et al. (2017) Pattern fusion analysis by adaptive alignment of multiple heterogeneous omics data. Bioinformatics, 33, 2706-2714.
#' @export
runPFA <- function(data,maxk=10,matlab)
{
  name=colnames(data[[1]])
  s="Y=MainPFA({"
  for (i in c(1:length(data))) {
    st=paste("setVariable(matlab,d",i,"=data[[i]])",sep='')
    eval(parse(text = st))
    s=paste(s,"d",i,",",sep='')
  }
  s=substr(s,1,nchar(s)-1)
  s=paste(s,"});",sep='')

  evaluate(matlab, s)

  ans=getVariable(matlab,"Y")$Y
  ans=as.matrix(ans)
  rownames(ans)=name

  clst=kmeansn(ans,2,maxk)
  clst=t(clst)
  rownames(clst)=name
  return(list(clst=clst,ans=ans))
}
