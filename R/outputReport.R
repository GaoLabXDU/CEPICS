#' @export
outputReport <- function(method=list('SNF','LRA','iCluster','PINS'),pic_eval=list('HeatMap','KM'),pinsk=c(2,3),numeric_eval=list('NMI','ARI','SI'),
                    need_gold_eval=list('NMI','ARI'),is_gold=TRUE,n_pic_row=2,kmax=7,path='',
                    filename='report.Rmd',title='Evaluation and Comparison Report',author="",is_report=TRUE,num_truelabel=0) {

  eval_full_name=list("NMI"="Normalized Mutual Information","ARI"="Adjusted Rand Index","SI"="Silhouette Coefficient","RI"="Rand Index","KM"="Kaplan-Meier Survival Curves","P-value"="-Log10(Cox P-value)","SC"="Silhouette Coefficient")
  method_full_name=list("PFA"="PFA","iCluster"="iClusterBayes","LRA"="LRAcluster","PINS"="PINS","SNF"="SNF")
  pf=paste(path,filename,sep='/')
  has_PINS=FALSE
  newmethod=list()
  for(e in method)
    if(e!="PINS")
      newmethod=append(newmethod,e)
    else
      has_PINS=TRUE
  method=newmethod


  writeHeader(header = list(title=title,author=author,
                            date=paste(strsplit(as.character(Sys.time()),split=' ')[[1]][1],
                                       strsplit(as.character(Sys.time()),split=' ')[[1]][2])),filename = pf)

if(!is_report)
    writeCode('',list(echo=FALSE,inlcude=FALSE,warning = FALSE,message = FALSE,error = FALSE),code='knitr::opts_chunk$set(dev="pdf",fig.path="fig/")' ,filename = pf)


  writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code='
options(warn =-1)
library(pheatmap)
library(survival)
library(ggplot2)
library(grid)
library(reshape2)
library(ggthemr)
ggthemr("fresh")
' ,filename = pf)

  h('Time Consumption',1,pf)
  writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code='load("timeList.RData")\ndrawtime(timeList)' ,filename = pf)
  writeRmd('\n\n',pf)



  concat <- function(xx,sep='',n=1) {
    s=""
    for(e in xx)
    {
      s=paste(s,e,sep=sep)
    }
    return(substr(s,n,nchar(s)))
  }

  if(is_gold)
  {
    h('Comparison based on True Labels',1,pf)
    h(paste('There are',num_truelabel,'clusters in true label'),3,pf)
    lst=list(c(1,1),c(1,2),c(2,1),c(2,2))
    xxx=1
    code="pushViewport(viewport(layout = grid.layout(2,2)))"
    for (e in list("P-value","NMI","ARI","SC")) {
      fnm=paste(path,"/goldMatrix_",e,".RData",sep="")
      if(file.exists(fnm))
      {
        fnm1=paste("goldMatrix_",e,sep="")
        code=paste(code,'\ndrawBarRmd("', fnm1 ,'",ylabel="',eval_full_name[[e]],'",vp=c(',lst[[xxx]][1] ,',' ,lst[[xxx]][2] ,'))'  ,sep="")
      }
      xxx=xxx+1
    }
    writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=code ,filename = pf)

    writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code="
drawzmat('Metrices','Methods',key='Score')
" ,filename = pf)
  }



  if(length(method)!=0 && "KM" %in% pic_eval)
  {
    h('Survival Analysis',1,pf)
    writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('drawPvalues(method=list(',concat(method,sep='","',n=3),'"),k=c(',2,':',kmax,'))',sep='') ,filename = pf)
  }

  if(length(method)!=0)
  {
    if(is_gold)
    {
      for(e in numeric_eval)
      {
        if(e %in% need_gold_eval )
        {
          h(eval_full_name[[e]],1,pf)
          writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('drawBarRmd("',e,'",',kmax,',"',eval_full_name[e],'")',sep='') ,filename = pf)
        }
      }
      for(e in numeric_eval)
      {
        if(!(e %in% need_gold_eval) )
        {
          h(eval_full_name[[e]],1,pf)
          writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('drawBarRmd("',e,'",',kmax,',"',eval_full_name[e],'")',sep='') ,filename = pf)
        }
      }
    }
    else
    {
      h(paste(eval_full_name[['NMI']],"and",eval_full_name[['ARI']]) ,1,pf)
      writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('
load("meanNMI.RData")
load("meanARI.RData")
drawHeatmapRmd(list(meanNMI,meanARI),title=list("NMI","ARI"),xlabel="Methods",ylabel="Number of Clusters")',sep='') ,filename = pf)
      for(i in (2:kmax))
      {
        h(paste('Comparison Based on',i,'Clusters') ,2,pf)
        writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('load("nogold_method.RData")
                                                              drawHeatmapRmd(list(nogold_method[["NMI"]][["',i,'"]],nogold_method[["ARI"]][["',i,'"]]),title=list("NMI","ARI"))',sep='') ,filename = pf)
      }

      if('SI' %in% numeric_eval)
      {
        h(eval_full_name[['SI']],1,pf)
        writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('drawBarRmd("','SI','",',kmax,',"',eval_full_name["SI"],'")',sep='') ,filename = pf)
      }
    }
    h('Samples Similarity HeatMap',1,pf)
    writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('
load("sampleClu.RData")
sampleClusterFig(sampleClu)',sep='') ,filename = pf)

  }

  writeRmd('========================================================================================',filename = pf)
  writeRmd('========================================================================================',filename = pf)

  h('Performance of Each Method',1,pf)
  for(e in method)
  {
      s=""
      if(is.null(method_full_name[[e]]))
        s=e
      else
        s=method_full_name[[e]]
      h(s,2,pf)
      if(is_gold)
      {
        writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('table=drawtable_all_eval("',e,'",list(',
                                                              concat(numeric_eval,'","',n=3),'"),list(',concat(pic_eval,'","',n=3),'"),c(',2,':',kmax,')',')' ,sep='') ,filename = pf)
        h('Results of  Different Metrics',3,filename = pf)

        ss=""
        i=0
        for(ee in numeric_eval)
        {
          if(ee=='SI')
            aa=paste("SC","-----",eval_full_name[[ee]])
          else
            aa=paste(ee,"-----",eval_full_name[[ee]])
          for(nbb in c(nchar(aa):45))
            aa=paste(aa,"&nbsp;",sep="")
          ss=paste(ss,aa)
          if(i%%2==1)
          {
            writeRmd(ss,filename=pf)
            writeRmd("\n",filename=pf)
            ss=""
          }
          i=i+1
        }
        if(length(numeric_eval)%%2==1)
        {
          writeRmd(ss,filename=pf)
          writeRmd("\n",filename=pf)
        }
        writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('
drawHeatmapRmd(list(t(table)),manyeval=TRUE,x_continuous=TRUE,k=as.integer(colnames(table)),x="Number of Clusters")',sep='') ,filename = pf)

      }

      if('KM' %in% pic_eval)
      {
        h('Kaplan-Meier Survival Curves
          ',3,pf)
        writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('\n') ,closechunk = FALSE,filename = pf)
        writeRmd(paste('par(mfrow=c(1,',n_pic_row,'))\n'),filename = pf)
        i=1
        for (c in c(2:kmax))
        {

          if(c!=2 && c%%n_pic_row==0)
          {
            writeRmd(paste('```\n'),filename = pf)
            writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('\n') ,closechunk = FALSE,filename = pf)
            i=i+1
            writeRmd(paste('par(mfrow=c(1,',n_pic_row,'))\n'),filename = pf)
          }
          writeRmd(paste('load("F_',e,'_',c,'_KM.RData")\n',sep='') ,filename = pf)
          writeRmd(paste('drawKMcurve(kmData,"Number of Clusters: ',  c ,'")\n') ,filename = pf)
        }
        writeRmd(paste('```\n'),filename = pf)
      }
      if('HeatMap' %in% pic_eval)
      {
        h('HeatMaps',3,pf)
        for (c in c(2:kmax))
        {
          h(paste('Number of Clusters: ',c,sep=""),4,pf)
          writeCode('',list(echo=FALSE,inlcude=TRUE,fig.width=5.8,warning = FALSE,message = FALSE,error = FALSE),code=paste('\n'),closechunk = FALSE,filename = pf)
          writeRmd(paste('load("F_',e,'_',c,'_heatmap.RData")\n',sep='') ,filename = pf)
          writeRmd('pheatmap(hmData$jz,cluster_rows = FALSE, cluster_cols = FALSE,color = colorRampPalette(c("#294998", "#0672ba", "#0ba9de", "#48c4ed", "#9ad8f0", "#ffffff", "#f7d6ad", "#f9ae82", "#ef837b", "#dc515a", "#ae3c4a"))(50), show_rownames = FALSE, show_colnames = FALSE)\n',filename = pf)
          writeRmd(paste('```\n'),filename = pf)
        }
      }

  }
  if(has_PINS)
  {
    if(!(!("KM" %in% pic_eval) && !is_gold))
      h('PINS',2,pf)
    if(is_gold)
    {
      e="PINS"
      if(length(pinsk)>1)
        writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('table=drawtable_all_eval("',e,'",list(',
                                                              concat(numeric_eval,'","',n=3),'"),list(',concat(pic_eval,'","',n=3),'"),c(',pinsk[1],',',pinsk[2],')',')' ,sep='') ,filename = pf)
      else
        writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('table=drawtable_all_eval("',e,'",list(',
                                                              concat(numeric_eval,'","',n=3),'"),list(',concat(pic_eval,'","',n=3),'"),c(',pinsk[1],')',')' ,sep='') ,filename = pf)

      h('Results of  Different Metrics',3,filename = pf)
      ss=""
      i=0
      for(ee in numeric_eval)
      {
        if(ee=='SI')
          aa=paste("SC","-----",eval_full_name[[ee]])
        else
          aa=paste(ee,"-----",eval_full_name[[ee]])
        for(nbb in c(nchar(aa):45))
          aa=paste(aa,"&nbsp;",sep="")
        ss=paste(ss,aa)
        if(i%%2==1)
        {
          writeRmd(ss,filename=pf)
          writeRmd("\n",filename=pf)
          ss=""
        }
        i=i+1
      }
      if(length(numeric_eval)%%2==1)
      {
        writeRmd(ss,filename=pf)
        writeRmd("\n",filename=pf)
      }
      writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('
                                                            drawHeatmapRmd(list(t(table)),manyeval=TRUE,x_continuous=TRUE,k=as.integer(colnames(table)),x="Number of Clusters")',sep='') ,filename = pf)

    }
    if("KM"%in% pic_eval)
    {
      h('Kaplan-Meier Survival Curves',3,pf)
      writeCode('',list(echo=FALSE,inlcude=TRUE,warning = FALSE,message = FALSE,error = FALSE),code=paste('\n'),closechunk = FALSE,filename = pf)
      writeRmd(paste('par(mfrow=c(1,',n_pic_row,'))\n'),filename = pf)
      s=paste("F_PINS_",pinsk[[1]],"_KM.RData",sep='')
      writeRmd(paste('load("',s,'")\n',sep='') ,filename = pf)
      writeRmd(paste('drawKMcurve(kmData,"Number of Clusters: ', pinsk[[1]]  ,'")\n') ,filename = pf)
      if(length(pinsk)>1)
      {
        s=paste("F_PINS_",pinsk[[2]],"_KM.RData",sep='')
        writeRmd(paste('load("',s,'")\n',sep='') ,filename = pf)
        writeRmd(paste('drawKMcurve(kmData,"Number of Clusters: ',  pinsk[[2]]  ,'")\n') ,filename = pf)
      }
      writeRmd(paste('```\n'),filename = pf)
    }
  }

  render(pf, html_document())
  if(is_report)
    delfile(pf)
  else
  {
    delfile(pf)
    pf=paste(path,'/',strsplit(filename,".",fixed =TRUE)[[1]][1],".html",sep='')
    delfile(pf)
  }

}
