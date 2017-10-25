# 计算停止准则的数值
stopcriter <- function(Init,Init3){
  temp1 <- Init[,1:2]-Init3[,1:2]
  temp2 <- data.frame(x1=Init3[,7],x2=-Init3[,5],x3=-Init3[,6],x4=Init3[,4])/(Init3[,4]*Init3[,7]-Init3[,5]*Init3[,6])
  j=temp1[,1]^2*temp2[,1]+temp1[,1]*temp1[,2]*(temp2[,2]+temp2[,3])+temp1[,2]^2*temp2[,4]
}

