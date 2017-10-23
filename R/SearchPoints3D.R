# 根据初始的范围 在加上中心坐标，边界来寻找球形邻域
SearchPoints3D <- function(x=0,
                         y=0,
                         z=0,
                         SearchPoints0,
                         xmin=-1000,
                         xmax=1000,
                         ymin=-1000,
                         ymax=1000,
                         zmin=-1000,
                         zmax=1000){
  #require(data.table)
  PointsPosition <- SearchPoints0
  PointsPosition <- data.frame(x=PointsPosition$x+x,y=PointsPosition$y+y,z=PointsPosition$z+z)
  PointsPosition <- subset(PointsPosition,x%in%c(xmin:xmax)&y%in%c(ymin:ymax)&z%in%c(zmin:zmax)) # 边界处理
  rownames(PointsPosition) <- c(1:nrow(PointsPosition))
  return(PointsPosition)
}
