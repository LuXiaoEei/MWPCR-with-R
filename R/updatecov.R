#  update cov
updatecov <- function(location,espn,h,Dis,A,data){
  X <- cbind(1,t(data[nrow(data),-c(1:3)]))
  x <- as.numeric(location[1])
  y <- as.numeric(location[2])
  z <- as.numeric(location[3])
  index <- (z-1)*400+(y-1)*20+x
  temp <- solve(matrix(A[index,],nrow = 2))
  cov <- temp%*%t(X)%*%diag(espn[index,]^2)%*%X%*%temp
  return(c(matrix(cov,nrow = 1)))
}
