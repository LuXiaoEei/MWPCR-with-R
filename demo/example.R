# x_ig=b_0(g)+b_1(g)l_i+esp_i(g)
# 3D image 20*20*20
require(EBImage)
dat2D <- data.frame(x=rep(c(1:20),20),y=rep(c(1:20),each=20),z=1,value=1)
dat2D$value[dat2D$x%in%c(7:14)|dat2D$y%in%c(7:14)] <- 0
# display(matrix(data2D$value,20,20))
dat3D <- dat2D
for (i in 1:9){
  dat3D <- rbind(dat3D,dat2D)
}
dat3D$z <- rep(c(1:10),each=400)
# display(matrix(dat3D[dat3D$z==1,]$value,20))
dat3D_one <- dat3d_two <- dat3D
dat3d_two$value[dat3d_two$x%in%c(8:13)&dat3d_two$y%in%c(8:13)&dat3d_two$z%in%c(4:7)] <-0.5
# display(matrix(dat3d_two[dat3D$z==4,]$value,20))
rm(dat3D,dat2D)
## noise iid N(0,0.2^2)
noise <- dat3D_one[,c('x','y','z')]
noise$esp1 <- rnorm(4000,0,0.2)
