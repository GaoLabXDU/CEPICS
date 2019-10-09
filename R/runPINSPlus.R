#' Run PINSPlus method to get the subtyping results.
#'
#' Using PINSPlus R package to get the subtyping results.
#'
#' @param datalist A list of different omics data. Each data in data list should be format as a data matrix with rows representing features and columns representing samples.
#' @param agreementCutoff agreement threshold to be considered consistent. Default value is 0.5. Please see PINSPlus package for more details.
#' @param kMax An integer value means the maximize number of clusters we will try.
#' @param cores An integer value means the number of cores for parallel computing.
#'
#' @return
#' The sample clustering results matrix.
#' @examples
#' data(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' datalist <- list(COAD_Methy, COAD_miRNA, COAD_mRNA)
#' res <- runPINSPlus(datalist, 5, 0.5)
#'
#' @references
#' Nguyen,T. et al. (2017) A novel approach for data integration and disease subtyping. Genome Res., 27, 2025-2039.
#' @export runPINSPlus
runPINSPlus <- function(datalist, kMax = 5, agreementCutoff = 0.5, cores = 1) {
  pre <- function(x) {
    tx <- t(x)
    return(tx)
  }
  d_list <- lapply(datalist, pre)

  PINSRes <- PINSPlus::SubtypingOmicsData(d_list, kMax = kMax, agreementCutoff = agreementCutoff, ncore = cores)

  cname <- names(PINSRes$cluster1)
  label <- matrix(nrow = length(PINSRes[[1]]), ncol = (length(PINSRes)-1))
  colnames(label) <- (1:(length(PINSRes)-1))
  for (i in 1:(length(PINSRes)-1) ) {
    label[,i] <- PINSRes[[i]]
    label[,i] <- renameCluRes(label[,i])
  }
  label <- apply(label, 2, as.numeric)
  for (i in 1:(length(PINSRes)-1) ) {
    colnames(label)[i] <- max(label[,i])
  }
  rownames(label) <- cname
  if (colnames(label)[1] == colnames(label)[2])
    label <- label[,-2]

  label <- cbind(NULL, label)
  for (i in 1:ncol(label)) {
    colnames(label)[i] <- max(label[,i])
  }
  label <- renamelable(label)
  return(label)
}

renameCluRes <- function(clu) {
  fac <- factor(clu)
  fac <- levels(fac)
  for (i in 1:length(fac)) {
    for (j in 1:length(clu)) {
      if (fac[i] == clu[j]) clu[j] = i
    }
  }
  return(clu)
}
