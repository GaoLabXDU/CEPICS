dependence=c("pheatmap",
             "flexclust",
             "cluster",
             "ggplot2",
             "R.matlab",
#             "iClusterPlus", install from bioconductor.org
             "DMwR",
             "survival",
             "PINSPlus",
             "SNFtool",
             "rmarkdown",
             "grid",
             "reshape2",
             "stats",
             "grDevices",
             "graphics",
             "utils",
             "devtools"
              )
installed=rownames(installed.packages())

need=c()
for (e in dependence) {
  if(!(e %in% installed))
    need=append(need,e)
}
if(!is.null(need))
  install.packages(need)

if(!('LRAcluster' %in% installed))
{
  download.file('http://lifeome.net/software/lracluster/LRAcluster_1.0.tgz')
  install.packages('LRAcluster_1.0.tgz', repos = NULL, type = "source")
  if(file.exists('LRAcluster_1.0.tgz'))
    file.remove('LRAcluster_1.0.tgz')
}

if(!('iClusterPlus' %in% installed))
  BiocManager::install("iClusterPlus")
