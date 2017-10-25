# 计算系数相似矩阵
Smatrix <- function(X,Init){
    XM <- solve(matrix(X[4:7],nrow = 2))
    diff <- t(apply(Init[,1:2],1,`-`,X[1:2]))
    return(apply(diff,1,function(Z,XM){t(Z)%*%XM%*%matrix(Z)},XM))
}


