# importance score weights
WI <- function(data,method='wald.test'){
  if(method=='wald.test'){
    data <- data[,-c(1:3)]
    X <- t(data[nrow(data),])
    dat <- data[-nrow(data),]
    return(apply(dat,1,waldP,X))
    }
}