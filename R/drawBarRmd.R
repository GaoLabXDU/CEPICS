
drawBarRmd <- function(s,kmax=10,ylabel="Num",vp=NULL,gmcp="gold_comp") {
  load(paste(s,'.RData',sep=''))
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  if(is.null(vp))
  {
    data=get(s)
    if (!is.null(data$PINS)) {
      if (as.numeric(names(data$PINS)[1]) > kmax) {
        data$PINS <- NULL
      }
    }
    mat <- matrix(nrow = kmax-1, ncol = length(data))
    rn <- c(2:kmax)
    cn <- c("K")
    for (i in 1:length(data)) {
      cn <- c(cn,names(data)[[i]])
      #names(mat)[i] <- data[[i]]$name
      for (j in 1:length(data[[i]])) {
        if(as.numeric(names(data[[i]])[j])<=kmax)
          mat[as.numeric(names(data[[i]])[j])-1,i] <- data[[i]][j]
      }
    }
    dat <- data.frame(rn, mat)
    names(dat) <- cn

    mdat <- melt(dat,id.vars = "K",variable.name="Function",value.name = "Num")
    g=ggplot(mdat,aes(x=K,y=Num,fill=Function)) + theme(text = element_text(size = 18)) +
      labs(x = "Number of Clusters", y = ylabel, fill = "Methods")+geom_bar(position="dodge",stat="identity")+ scale_x_continuous(breaks = 2:kmax)
    print(g)
  }
  else
  {
    data=get(gmcp)
    data <- round(data, 3)
    funcname <- names(data)
    golddata <- data
    data <- data.frame(funcname, golddata)
    g=ggplot(data, aes(x = funcname, y = golddata))+ geom_bar(aes(fill = funcname), stat = 'identity', width = 0.5) + geom_text(aes(label = data$golddata))+theme(legend.position = 'none',text = element_text(size = 10)) + labs(x = '', y = ylabel)
    print(g,vp=vplayout(vp[1],vp[2]))
  }


}
