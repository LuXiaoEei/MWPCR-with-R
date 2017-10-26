# compute WE X belong data:rbind(warm+Dis)
WE <- function(X,h){
  temp <- rep(0,4000)
  warm <- X[1:4000]
  dis <- X[4001:8000]
  index <- dis<=h
  temp[index] <- warm[index]/sum(warm[index])
  return(temp)
}
