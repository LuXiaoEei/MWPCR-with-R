# computr espn(i,w,hs)
epsn <- function(location,W,Dis,h,Init,beta,data){
  X <- cbind(1,t(data[nrow(data),-c(1:3)]))
  x <- as.numeric(location[1])
  y <- as.numeric(location[2])
  z <- as.numeric(location[3])
  index <- (z-1)*400+(y-1)*20+x
  Beta <- t(beta[Dis[,index]<=h,])
  Y <- t(data[-nrow(data),][Dis[,index]<=h,-c(1:3)])
  w <- W[Dis[,index]<=h,index]
  res <- Init[Dis[,index]<=h,3]
  return(t(w/res)%*%t(Y-X%*%Beta))
}
