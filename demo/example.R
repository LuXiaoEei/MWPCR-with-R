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
  x <- location[1]
  y <- location[2]
  z <- location[3]
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

esp_iid <- dat3D_zero[,c('x','y','z')]
esp_short <- dat3D_zero[,c('x','y','z')]
esp_long <- dat3D_zero[,c('x','y','z')]

for(i in 1:100){
  a <- Sys.time()
  noise$esp_iid <- rnorm(4000,0,0.2)
  noise$esp_short <- apply(noise,1,EspShort,noise,SearchPoints0)
  noise$esp_long <- EspLong(g,avg,std,noise)
  if (i %in% 1:60){
    esp_iid=cbind(esp_iid,dat3D_zero$value+noise$esp_iid)
    esp_short=cbind(esp_short,dat3D_zero$value+noise$esp_short)
    esp_long=cbind(esp_long,dat3D_zero$value+noise$esp_long)
  }else{
    esp_iid=cbind(esp_iid,dat3D_one$value+noise$esp_iid)
    esp_short=cbind(esp_short,dat3D_one$value+noise$esp_short)
    esp_long=cbind(esp_long,dat3D_one$value+noise$esp_long)
  }
  print(c(Sys.time()-a,i))
}
esp_iid[,nrow(esp_iid)+1] <- c(rep(0,60),rep(1,40))
esp_short[,nrow(esp_short)+1] <- c(rep(0,60),rep(1,40))
esp_long[,nrow(esp_long)+1] <- c(rep(0,60),rep(1,40))
x=1
y=1
z=1
xmin=1
xmax=20
ymin=1
ymax=20
zmin=1
zmax=10