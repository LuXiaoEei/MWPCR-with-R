# WARM计算h0时刻的系数，残差，估计值的协方差阵
OLS_loc<- function(location,
                   data,
                   Dis,
                   h
                     # SearchPoints0,
                     # xmin=1,
                     # xmax=20,
                     # ymin=1,
                     # ymax=20,
                     # zmin=1,
                     # zmax=10
){
  # print(location[1])
  x <- as.numeric(location[1])
  y <- as.numeric(location[2])
  z <- as.numeric(location[3])
  index <- (z-1)*400+(y-1)*20+x
  # print(index)
  loc <- t(data[-nrow(data),][Dis[,index]<=h,-c(1:3)])
  # Area <- SearchPoints3D(x,y,z,SearchPoints0,xmin,xmax,ymin,ymax,zmin,zmax)
  # loc <- t(merge(Area,data,by=c('x','y','z'))[,-c(1:3)])
  X <- rep(t(data[nrow(data),-c(1:3)]),ncol(loc))
  Y <- matrix(loc,ncol = 1)
  df <- lm(Y~X)
  beta <- coef(df)
  res <- sum(resid(df)^2)/(ncol(loc)*nrow(loc)-2)
  cov <- matrix(vcov(df),1)
  return(c(beta,res,cov))
}
