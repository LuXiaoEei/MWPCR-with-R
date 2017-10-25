# updta beta
updatebeta <- function(location,h,W,Dis,data,Init,A){
  X <- cbind(1,t(data[nrow(data),-c(1:3)]))
  x <- as.numeric(location[1])
  y <- as.numeric(location[2])
  z <- as.numeric(location[3])
  index <- (z-1)*400+(y-1)*20+x
  temp <- matrix(A[index,],nrow = 2)
  w <- W[Dis[,index]<=h,index]
  res <- Init[Dis[,index]<=h,3]
  Y <- t(data[-nrow(data),][Dis[,index]<=h,-c(1:3)])%*%as.matrix(w/res)
  # print(c(solve(temp)%*%t(X)%*%Y))
  return(c(solve(temp)%*%t(X)%*%Y))
}
