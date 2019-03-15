drawtime <- function(timeList) {
  Var1=0
  Var2=0
  value=0
  x<-timeList
  vplayout <- function(x,y){
    viewport(layout.pos.row = x, layout.pos.col = y)
  }
  pushViewport(viewport(layout = grid.layout(20,20)))
  myLabel <- as.vector(names(timeList))
  myLabel = paste(myLabel, " (", round(x / sum(x) * 100, 2), "%)", sep = "")
  timeList <- as.data.frame(timeList)
  timename <- rownames(timeList)
  Methods <- paste(timename, "(", round(timeList$timeList / sum(timeList$timeList) * 100, 2), "%)  ", sep = "")
  gp = ggplot(data=timeList, mapping=aes(x="",y = timeList ,fill=Methods))+
    geom_bar(stat="identity")+coord_polar(theta = "y") + labs(x = "", y = "", title = "")+theme(axis.ticks = element_blank()) + theme(axis.text.x = element_blank())+theme(panel.grid=element_blank())+theme(panel.border=element_blank())+theme(axis.line = element_blank())
  print(gp,vp = vplayout(1:13,6:19))

  timeList=as.matrix(timeList)
  colnames(timeList)="Time(s)"
  m = melt(as.matrix(timeList))
  g = ggplot(m, aes(x=Var1, y=Var2, fill=value))
  g <- g +labs(x = "", y = "", title = "") + geom_tile() + coord_fixed(ratio = 2/5) + geom_text(aes(label=round(value, 2)), colour = 'white') + guides(fill = F)
  g<-g+ scale_fill_gradient(low = "#122a42", high = "#56b1f7")
  print(g,vp = vplayout(14:20,1:20))
}
