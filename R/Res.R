#' An example result of our own method.
#'
#' @description
#' CEPICS allows users to upload their own subtyping results so that it is very convenient for researchers to evaluate and compare their own methods to the state-of-the-art ones. "Res" is an example result of our own method which is a list with following elements:
#' \itemize{
#'    \item{
#'    data: a list with the following elements:
#'        \itemize{
#'              \item{data: a sample-sample similarity matrix. We calculate the Euclid distance between every two samples Then we divide 1 by the distance matrix and scaled it by the maximum of the matrix. Finally, we set the diagonal of the matrix to 1.}
#'              \item{cluster: the cluster result matrix. The n-th row is the cluster label of each sample when cluster samples into  (n+1) clusters.}
#'         }
#'
#'    }
#'    \item{name: "Ours", the name of our own method.}
#'    \item{time: 66, the time consumption of our own method.}
#' }
#'
#' @details
#' Res dataset can be only used in funtion \code{\link{CEPICS}}.
#'
#' @seealso
#' \code{\link{CEPICS}}
#'
#' @examples
#' data(Res)
#'
#' #we get the similarity matrix by following code:
#' data("COAD_miRNA")
#' sim=dist(t(COAD_miRNA))
#' sim=as.matrix(1/sim)
#' sim=sim/max(sim)
#' sim=sim+diag(1,nrow(sim),ncol(sim))
"Res"
