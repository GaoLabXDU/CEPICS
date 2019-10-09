#' Run the iClusterBayes method
#'
#' @param dl The datalist needs to be processed.
#' @param type Data type corresponding to dl1-6, which can be gaussian, binomial, or poisson.
#' @param kMax The maximum number of clusters.
#' @param n.burnin Please see "iClusterPlus" package for more details.
#' @param n.draw Please see "iClusterPlus" package for more details.
#' @param prior.gamma Please see "iClusterPlus" package for more details.
#' @param sdev Please see "iClusterPlus" package for more details.
#' @param beta.var.scale Please see "iClusterPlus" package for more details.
#' @param thin Please see "iClusterPlus" package for more details.
#' @param pp.cutoff Please see "iClusterPlus" package for more details.
#' @param cores An integer value means the number of cores for parallel computing.
#'
#' @return
#' 2 to kMax clustering results and sample similarity matrix.
#' @examples
#' \donttest{
#' data(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' datalist <- list(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' res <- runiClusterBayes(datalist, type = c("gaussian","gaussian","gaussian"),kMax=5,n.burnin=1000,
#'      n.draw=1200,prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5)
#' }
#'
#' @references
#' Mo,Q. et al. (2017) A fully Bayesian latent variable model for integrative clustering analysis of multi-type omics data. Biostatistics, 19, 71-86.
#'
#' @export
runiClusterBayes <- function(dl, type = c("gaussian","gaussian","gaussian","gaussian","gaussian","gaussian"),kMax=4,n.burnin=1000,n.draw=1200,
                                                                                     prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5, cores = 1) {
  dt1 <- NULL;
  dt2 <- NULL;
  dt3 <- NULL;
  dt4 <- NULL;
  dt5 <- NULL;
  dt6 <- NULL;
  if (!is.null(dl[1][[1]])) {
    dt1 <- t(dl[[1]])
  }
  if (!is.null(dl[2][[1]])) {
    dt2 <- t(dl[[2]])
  }
  if (!is.null(dl[3][[1]])) {
    print(3)
    dt3 <- t(dl[[3]])
  }
  if (!is.null(dl[4][[1]])) {
    dt4 <- t(dl[[4]])
  }
  if (!is.null(dl[5][[1]])) {
    dt5 <- t(dl[[5]])
  }
  if (!is.null(dl[6][[1]])) {
    dt6 <- t(dl[[6]])
  }

  #a matrix with rows and columns representing samples and genomic features, respectively.
  #result <- tune.iClusterBayes(cpus = 1, dt1 = dt1, dt2 = dt2, dt3 = dt3, dt4 = dt4, dt5 = dt5, dt6 = dt6, type = type, K = 1:(kMax-1), n.burnin = n.burnin, n.draw = n.draw,
  #                                      prior.gamma=prior.gamma ,sdev=sdev,beta.var.scale=beta.var.scale,thin=thin,pp.cutoff=pp.cutoff)
  cl <- makeCluster(cores)
  #clusterExport(cl, "dt1")
  #clusterExport(cl, "dt2")
  #clusterExport(cl, "dt3")
  #clusterExport(cl, "dt4")
  #clusterExport(cl, "dt5")
  #clusterExport(cl, "dt6")
  result <- parLapplyLB(cl, 1:(kMax-1), fun = function(x) {
    library(iClusterPlus)
    iClusterPlus::iClusterBayes(dt1 = dt1, dt2 = dt2, dt3 = dt3, dt4 = dt4, dt5 = dt5, dt6 = dt6, type = type, K = x, n.burnin = n.burnin, n.draw = n.draw,
                                prior.gamma=prior.gamma ,sdev=sdev,beta.var.scale=beta.var.scale,thin=thin,pp.cutoff=pp.cutoff)},
    chunk.size=1)
  stopCluster(cl)
  cat("End parallel computation\n")

  result <- list(fit = result)

  name <- rownames(dt1)
  finRes <- NULL
  clu <- NULL
  for (i in 1:(kMax-1)) {
    resm <- result$fit[[i]]$meanZ
    rownames(resm) <- name
    finRes <- c(finRes, list(resm))
    names(finRes)[length(finRes)] <- (i+1)

    clu <- cbind(clu, result$fit[[i]]$clusters)
  }


  rownames(clu) <- name
  colnames(clu) <- c(2:kMax)

  clu <- renamelable(clu)
  return(list(datas = finRes, cluster = clu))
}
