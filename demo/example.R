require(EBImage)
require(RMWPCR)
# short range spatial correlation noise
EspShort <- function(location,
                      noise,
                      SearchPoints0,
                      xmin=1,
                      xmax=20,
                      ymin=1,
                      ymax=20,
                      zmin=1,
                      zmax=10
                      ){
  # print(location[1])
  x <- as.numeric(location[1])
  y <- as.numeric(location[2])
  z <- as.numeric(location[3])
  Area <- SearchPoints3D(x,y,z,SearchPoints0,xmin,xmax,ymin,ymax,zmin,zmax)
  loc <- merge(Area,noise,by=c('x','y','z'))
  return(sum(loc$esp_iid)/nrow(loc))
}

# long range spatial correlation noise
EspLong <- function(g,avg,std,noise){
  zeta1 <- rnorm(nrow(noise),avg,std)
  zeta2 <- rnorm(nrow(noise),avg,std)
  zeta3 <- rnorm(nrow(noise),avg,std)
  return(2*sin(pi*g[1])*zeta1+2*cos(pi*g[2])*zeta2+2*sin(pi*g[3]/0.5)*zeta3+noise$esp_iid)
}


# x_ig=b_0(g)+b_1(g)l_i+esp_i(g)
# 3D image 20*20*20

dat2D <- data.frame(x=rep(c(1:20),20),y=rep(c(1:20),each=20),z=1,value=0.8)
dat2D$value[dat2D$x%in%c(7:14)|dat2D$y%in%c(7:14)] <- 0.2
# display(matrix(data2D$value,20,20))
dat3D <- dat2D
for (i in 1:9){
  dat3D <- rbind(dat3D,dat2D)
}
dat3D$z <- rep(c(1:10),each=400)
# display(matrix(dat3D[dat3D$z==1,]$value,20))
dat3D_zero <- dat3D_one <- dat3D
dat3D_one$value[dat3D_one$x%in%c(8:13)&dat3D_one$y%in%c(8:13)&dat3D_one$z%in%c(4:7)] <-0.5

# display(matrix(dat3d_two[dat3D$z==4,]$value,20))
rm(dat3D,dat2D)
## noise iid N(0,0.2^2)

# generate 100 subject with 60 from zero and the rest from one 
avg <- 0
std <- 0.1
g <- rnorm(3,avg,std)
SearchPoints0 <- SearchPoints3D0(h=1)
noise <- dat3D_zero[,c('x','y','z')]

esp_0 <- dat3D_zero[,c('x','y','z')]#不含噪声
esp_iid <- dat3D_zero[,c('x','y','z')]
esp_short <- dat3D_zero[,c('x','y','z')]
esp_long <- dat3D_zero[,c('x','y','z')]

for(i in 1:100){
  a <- Sys.time()
  noise$esp_iid <- rnorm(4000,0,0.1)
  noise$esp_short <- apply(noise,1,EspShort,noise,SearchPoints0)
  noise$esp_long <- EspLong(g,avg,std,noise)
  if (i %in% 1:60){
    esp_0 <- cbind(esp_0,dat3D_zero$value)
    esp_iid <- cbind(esp_iid,dat3D_zero$value+noise$esp_iid)
    esp_short <- cbind(esp_short,dat3D_zero$value+noise$esp_short)
    esp_long <- cbind(esp_long,dat3D_zero$value+noise$esp_long)
  }else{
    esp_0 <- cbind(esp_0,dat3D_one$value)
    esp_iid <- cbind(esp_iid,dat3D_one$value+noise$esp_iid)
    esp_short <- cbind(esp_short,dat3D_one$value+noise$esp_short)
    esp_long <- cbind(esp_long,dat3D_one$value+noise$esp_long)
  }
  print(c(Sys.time()-a,i))
}
esp_0 <- rbind(esp_0,c(NA,NA,NA,rep(0,60),rep(1,40)))
esp_iid <- rbind(esp_iid,c(NA,NA,NA,rep(0,60),rep(1,40)))
esp_short <- rbind(esp_short,c(NA,NA,NA,rep(0,60),rep(1,40)))
esp_long <- rbind(esp_long,c(NA,NA,NA,rep(0,60),rep(1,40)))
colnames(esp_0) <- NULL
colnames(esp_iid) <- NULL
colnames(esp_short) <- NULL
colnames(esp_long) <- NULL
# write.csv(esp_0,file = './data/esp_0.csv')
# write.csv(esp_iid,file = './data/esp_iid.csv')
# write.csv(esp_short,file = './data/esp_short.csv')
# write.csv(esp_long,file = './data/esp_long.csv')
# save(esp_0,file = './data/esp_0.Rdata')
# save(esp_iid,file = './data/esp_iid.Rdata')
# save(esp_short,file = './data/esp_short.Rdata')
# save(esp_long,file = './data/esp_long.Rdata')




# esp_iids
# compute importace score weights matrix
Dis <- as.matrix(dist(esp_iid[-nrow(esp_iid),c(1:3)],diag = TRUE,upper = TRUE)) 
W <- WI(data = esp_iid)
# chao jie
    W[W==0] <- exp(-745)
Wi <- -4000*log(W)/(-sum(log(W)))
# compute spatial weights matrix
ch=1.2
S=5
warm <- WARM(data = esp_iid,ch = ch,S = S,Dis = Dis)
We <- apply(rbind(warm,Dis),2,WE,h=ch^S) #ch=1.2
# 创建多尺度阈值数据集
criter <- data.frame(si=quantile(Wi)[2:4],
                     se1=quantile(We[We!=0])[2:4],
                     se2=ch^c(1,3,5))
# 生成多尺度的权重矩阵Q....Q1 Q2 Q3
Qmatrix <- createQ(We,Wi,Dis,criter)
# 数据加权
X <- t(esp_iid[-nrow(esp_iid),-c(1:3)])
# 减去均值
X <- X-matrix(rep(1,nrow(X)),ncol = 1)%*%colMeans(X)
Xweig <- array(dim=c(nrow(X),ncol(X),nrow(criter)))
for (i in 1:nrow(criter)) {
  Xweig[,,i] <- X%*%Qmatrix[,,i]
}
# svd 分解
K=5 #特征值个数
for (i in 1:nrow(criter)){
  message(i,'...',Sys.time())
  assign(paste('svd',i,sep=''),svd(Xweig[,,i],K,K))
}
# svd1 <- svd(Qmatrix[,,1],5,5)
# svd2 <- svd(Qmatrix[,,2],5,5)
# svd3 <- svd(Qmatrix[,,3],5,5)
U <- c()
for (i in 1:nrow(criter)){
  U <- cbind(U,eval(parse(text=paste('svd',i,sep = '')))[['u']])
}

Y <- t(esp_iid[nrow(esp_iid),-c(1:3)])
data_train <- data.frame(Y=c(Y),U)
df <- lm(Y~.,data = data_train)
df1 <- glm(Y~.,data = data_train,family = binomial(link='logit'),control = list(maxit = 100) )


# 生成新的esp_iid测试数据
esp_iid_test <- esp_iid[,-sample(3:103,50)]
esp_iid_test <- esp_iid[,c(1:3,4,5,6)]
esp_iid_test <- dat3D_zero[,c('x','y','z')]
for (i in 1:100){
  if (i<45){
    esp_iid_test <- cbind(esp_iid_test,dat3D_zero$value+rnorm(4000,0,0.01))
  }else{
    esp_iid_test <- cbind(esp_iid_test,dat3D_one$value+rnorm(4000,0,0.01))
  }
}
esp_iid_test <- rbind(esp_iid_test,c(NA,NA,NA,rep(0,10),rep(1,10)))
colnames(esp_iid_test) <- NULL

# solve new u
X_test <- t(esp_iid_test[-nrow(esp_iid_test),-c(1:3)])
# 减去均值
X_test <- X_test-matrix(rep(1,nrow(X_test)),ncol = 1)%*%colMeans(X_test)
U_test <- c()
for (i in 1:nrow(criter)){
  svd <- eval(parse(text=paste('svd',i,sep = '')))
  u <- svd[['u']]
  d <- diag(svd[['d']][1:K])
  v <- svd[['v']]
  U_test <- cbind(U_test,X_test%*%Qmatrix[,,i]%*%v%*%solve(d))
}
Y_test <- t(esp_iid_test[nrow(esp_iid_test),-c(1:3)])
data_test <- data.frame(Y=c(Y_test),U_test)
# 预测
Y_pred <- predict(df,data_test)
Y_pred1 <- predict(df1,data_test,type = 'response')

