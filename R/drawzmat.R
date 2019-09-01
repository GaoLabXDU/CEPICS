drawzmat <- function(x="",y="",t="",key="score") {
  load("zmat.RData")
  maxName <- names(sort(zmat[key,], decreasing = TRUE))[1]
  m <- melt(zmat,na.rm = TRUE)
  m[["sign"]] = ifelse(m[["Var2"]] == maxName, "Yes", "No")

  ggplot(m, aes(x=Var1, y=Var2, fill=sign)) + theme(text = element_text(size = 18)) + geom_tile()+
    coord_fixed() + geom_text(aes(label=signif(value, 4)), colour = 'black', size = 4.4) + guides(fill = F) +
    scale_fill_manual(values = c("Yes" = "#87D3EF", "No" = "white"))+labs(x = x, y = y, title = t)
}
