drawHeatmapRmd <- function(data,text=TRUE,k,ncol=1,x=NULL,x_continuous=FALSE,sig=3,title=list(),color=TRUE,manyeval=FALSE,xlabel="",ylabel="") {
  Var1=0
  Var2=0
  value=0
  size=4.4
  xticksize=18
  if(!is.null(x))
    xlabel=x
  if(text)
  {
    pushViewport(viewport(layout = grid.layout(ncol,length(data))))
    vplayout <- function(x,y){
      viewport(layout.pos.row = x, layout.pos.col = y)
    }
    for(i in c(1:length(data)))
    {
      if(length(data)>1)
      {

        w1=11
        w2=-0.08
        nonadata=data[[i]]
        nonadata[is.na(nonadata)]=0
        npre=max(nchar(as.character(nonadata)))


        a1=6
        a2=-0.5
        size=w1*par("pin")[1]/3/nrow(nonadata)+w2*npre+0.1

        xticksize=a1*par("pin")[1]/3/nrow(nonadata)+a2*max(nchar(rownames(nonadata)))+10
      }

      if (x_continuous) {
        tdata <- as.matrix(data[[i]])
        rownames(tdata) <- 1:nrow(tdata)
        m <- melt(tdata,na.rm = TRUE)
      }
      else m = melt(as.matrix(data[[i]]),na.rm = TRUE)




      g = ggplot(m, aes(x=Var1, y=Var2, fill=value)) + theme(text = element_text(size = 18))
      g <- g + theme(axis.text.x = element_text(size = xticksize))
      if(manyeval)
        g <- g +labs(x = "Number of Clusters", y = "", title = "") + geom_tile() + coord_fixed() + geom_text(aes(label=signif(value, sig)), colour = 'black', size = 4.4) + guides(fill = F)
      else
        g <- g +labs(x = "", y = "", title = "") + geom_tile() + coord_fixed() + geom_text(aes(label=signif(value, sig)), colour = 'white', size = size) + guides(fill = F)
      if(color && !manyeval)
        g<-g+ scale_fill_gradient(low = "#56b1f7", high = "#122a42")
      else if(!manyeval)
        g<-g+ scale_fill_gradient(low = "#122a42", high = "#56b1f7")
      else
        g<-g+ scale_fill_gradient(low = "#FFFFFF", high = "#FFFFFF")
      if(length(title)!=0)
      {
        g <- g +labs(x = "", y = "", title = title[[i]])+theme(plot.title = element_text(hjust = 0.5))
      }
      if(x_continuous==TRUE)
      {
        g <- g +scale_x_continuous(breaks = 1:(length(k)), labels = k)
      }
      if(!is.null(x))
        g <- g +labs(x = x)
      g=g+labs(x = xlabel, y = ylabel)
      print(g,vp = vplayout(1,i))
    }

  }
  else
  {
    m = melt(data)
    g = ggplot(m, aes(x=Var1, y=Var2, fill=value))
    g <- g + labs(x = "", y = "", title = "", fill = "")+ geom_tile()  + coord_fixed() + theme(axis.text = element_blank(), axis.line = element_blank()) + theme(axis.ticks = element_blank())
    if(color)
      g<-g+ scale_fill_gradient(low = "#fffdcf", high = "#eb7170")
    else
      g<-g+ scale_fill_gradient(low = "#122a42", high = "#56b1f7")
    print(g)
  }
}
