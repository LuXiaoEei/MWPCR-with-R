# 生成多尺度的权重矩阵Q....Q1 Q2 Q3
createQ <- function(We,Wi,Dis,criter){
  Q <- array(dim = c(dim(We)[1],dim(We)[1],nrow(criter)))
  for (i in 1:nrow(criter)){
    We0 <- We
    We0[We<criter[i,2]|Dis>criter[i,3]] <- 0
    Q[,,i] <- We0%*%diag(sign(Wi>criter[i,1]))
  }
  # names(Q) <- paste('Q',1:nrow(criter),sep = '')
  return(Q)
}
