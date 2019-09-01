drawtable_all_eval <- function(e,numeric_eval,pic_eval,k) {
  if(e=="PINS")
  {
    x=list()
    for(ss in numeric_eval)
      if(ss!="SI")
        x=append(x,ss)
    numeric_eval=x
  }

  nrow_table=length(numeric_eval)
  if("KM" %in% pic_eval)
    nrow_table=nrow_table+1
  table=matrix(0,nrow_table,length(k))

  if("KM" %in% pic_eval)
  {
    rownames(table)=append(numeric_eval,"P-value",0)
  }else
  {
    rownames(table)=numeric_eval
  }

  aa=rownames(table)
  bb=list("P-value","NMI","RI","ARI","SI")
  cc=list()
  for (ee in bb) {
    if(ee %in% aa)
      cc=append(cc,ee)
  }
  if(length(aa)==length(cc))
    rownames(table)=cc


  colnames(table)=k

  for(r in rownames(table))
  {

    if(r=="P-value")
    {
      for(x in k)
      {
        load(paste('F_',e,'_',x,'_KM.RData', sep=''))
         if (kmData$misv == 2) pv <- summary(kmData$cox)$sctest[3]
          else  pv <- summary(kmData$cox)$coefficients[1,5]
        table[r,as.character(x)]=round(pv, digits = 10)
      }
    }
    if(r=="NMI")
    {
      load('NMI.RData')
      ans=NMI[[e]]
      for(x in k)
      {

        table[r,as.character(x)]=ans[as.character(x)]
      }
    }
    if(r=="ARI")
    {
      load('ARI.RData')
      ans=ARI[[e]]
      for(x in k)
      {

        table[r,as.character(x)]=ans[as.character(x)]
      }
    }
    if(r=="RI")
    {
      load('RI.RData')
      ans=RI[[e]]
      for(x in k)
      {

        table[r,as.character(x)]=ans[as.character(x)]
      }
    }
    if(r=="SI")
    {
      load('SI.RData')
      ans=SI[[e]]
      for(x in k)
      {

        table[r,as.character(x)]=ans[as.character(x)]
      }
    }
  }
  aa=as.list(rownames(table))
  for (i in c(1:length(aa))) {
    if(aa[[i]]=="SI")
      aa[i]="SC"
  }
  rownames(table)=aa
  return(table)
}
