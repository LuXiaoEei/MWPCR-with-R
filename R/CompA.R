# compute A(d,hs,W,X)
CompA <- function(location,h,W,Dis,data,Init){
  X <- cbind(1,t(data[nrow(data),-c(1:3)]))
  x <- as.numeric(location[1])
  y <- as.numeric(location[2])
  z <- as.numeric(location[3])
  index <- (z-1)*400+(y-1)*20+x
  w <- W[Dis[,index]<=h,index]
  res <- Init[Dis[,index]<=h,3]
  return(c(matrix(t(X)%*%X*sum(w/res),nrow=1)))
}
