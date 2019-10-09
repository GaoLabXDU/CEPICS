selectev <- function(ev)
{
  len=length(ev)
  mid=(1+len)/2
  r=matrix(0,len+1,len+1)
  v=var(c(1:len))

  for (i in c(1:len))
    for (j in c(1:i-1))
      r[i,j]=(ev[i]-ev[j])/(i-j)**2.2

  a=matrix(0,1,len)
  for (i in c(2:len))
    a[i]=(sum(r[i,0:i])-sum(r[(i+1):(len+1),i]))/exp(((i-mid)/v/5)**2)

  rank=order(a,decreasing = TRUE)
  #print(paste("the best number of dimension after reducing may be ",str(rank[1]-1), sep = ""))
  #plot(c(0:(len-1)),ev,xlim=c(0,10),xlab="Dimension",ylab="Explained Variation",type='l')
  #lines(rank[1]-1,ev[rank[1]],'p')
  #lines(rank[2]-1,ev[rank[2]],'o')
  #lines(rank[3]-1,ev[rank[3]],'o')

  return(rank[1]-1)
}
#
