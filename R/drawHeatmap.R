drawHeatmap <- function(result, fileName = NULL, startTime = NULL) {

   cc <- result$clusters
   sor <-order(cc, decreasing = F)
   jz <- result$affmat

   normalize <- function(X) X/rowSums(X)
   ind <- sort(as.vector(cc), index.return = TRUE)
   ind <- ind$ix

   diag(jz) <- median(as.vector(jz))
   jz <- normalize(jz)
   jz <- jz + t(jz)
   jz <- jz[ind, ind]

   if (!is.null(fileName)) {
     hmData <- list(jz = jz)
     save(hmData, file = as.character(fileName))
   }
   else {
     pheatmap(jz,cluster_rows = FALSE, cluster_cols = FALSE,color = colorRampPalette(c("white", "black"))(50), show_rownames = FALSE, show_colnames = FALSE)
   }


}
