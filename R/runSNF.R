#' Run SNF method to get the subtyping results.
#'
#' Using SNFTool R package to get the subtyping results.
#'
#' @param Wall A list of different omics data. Each data in data list should be format as a data matrix with rows representing features and columns representing samples.
#' @param k An integer number means the number of neighbors in K-nearest neighbors.Default \code{NULL}.
#' @param t An integer number means the number of iterations for the diffusion process. Default \code{20}.
#' @param kMax An integer value means the maximize number of clusters we will try. Default \code{5}.
#'
#' @return
#' The sample clustering results matrix including 2 to kMax clustering results, and sample similarity matrix.
#' @examples
#' data(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' datalist <- list(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' res <- runSNF(datalist, k = 10, t = 20, kMax = 5)
#'
#' @references
#' Wang,B. et al. (2014) Similarity network fusion for aggregating data types on a genomic scale. Nature Methods, 11, 333-337.
#'
#' Concise description can be found here: http://compbio.cs.toronto.edu/SNF/SNF/Software.html
#'
#' @export runSNF
runSNF <- function(Wall, k = NULL, t = 20, kMax = 5) {

  if(is.null(k)) k <- ncol(Wall[[1]])/10

  K = k
  alpha = 0.5
  T = t

  pre <- function(X) {
    X <- t(X)
    Dist = (dist2(as.matrix(X),as.matrix(X)))^(1/2)
    Wx = affinityMatrix(Dist, K, alpha)
    return(Wx)
  }


  w_list <- lapply(Wall, pre)

  W = SNFtool::SNF(w_list, K, T)

  name <- colnames(W)

  if (kMax <= 2) {
    labels = spectralClustering(W, 2)
    names(labels) <- name
  }
  else {
    labels <- matrix(nrow = nrow(W), ncol = (kMax-1))

    for (i in 2:kMax) {
      labels[,i-1] = spectralClustering(W, i)

    }
    colnames(labels) <- c(2:kMax)
    rownames(labels) <- name
  }
  finres <- list(affmat = W, clusters = labels)

  return(finres)
}
