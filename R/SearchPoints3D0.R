### 根据窗宽寻找球形邻域内的点的坐标
#(0,0)中心点坐标，h窗宽,interval坐标间隔,无边界

SearchPoints3D0 <- function(h,interval=1){
  #require(data.table)
  #require(dplyr)
  n <- c(-floor(h/interval):floor(h/interval))
  PointsPostotion0 <- data.frame(
      x = rep(n, length(n) ^ 2),
      y = rep(rep(n, each = length(n)), length(n)),
      z = rep(n, each = length(n) ^ 2)
    )
  PointsPostotion0 <- PointsPostotion0[apply(PointsPostotion0,1,function(x){sum(x^2)})<=h^2,]
  return(PointsPostotion0)
}
