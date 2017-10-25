# WARM
# data <- esp_iid
# names(data)[4:ncol(data)] <- NA
WE <- function(data,
               ch = 1.2,
               Cn = log(100) * qchisq(df = 2, p = 0.95),
               Kloc = TRUE,
               Kst = TRUE,
               S0 = 3,
               S = 5){
    message('初始化...',Sys.time())
    #距离矩阵
    message('初始化距离矩阵...',Sys.time())
    Dis <- as.matrix(dist(data[-nrow(data),c(1:3)],diag = TRUE,upper = TRUE)) 
    # SearchPoints0 <- SearchPoints3D0(h = ch^0)
    message('初始化beta...',Sys.time())
    Init <- t(apply(data[-nrow(data),],1,OLS_loc,data,Dis,ch^0))
    Init0 <- Init
    colnames(Init) <- c('beta0','beta1','res','cov00','cov01','cov10','cov11')
    message('初始化相似矩阵...',Sys.time())
    Sim <- apply(Init[-nrow(data),],1,Smatrix,Init)#相似矩阵
    W_loc <- diag(nrow(data)-1)
    W_st <- diag(nrow(data)-1)
    message('迭代开始...',Sys.time())
    for (i in 1:S){
      message('开始第',i,'此迭代')
      message('updata W...',Sys.time())
      if (Kloc){
        W_loc <- (1-Dis/ch^i)*sign((1-Dis/ch^i)>=0)
      }
      if (Kst){
        W_st <- exp(-Sim/Cn)
      }
      W <- W_loc*W_st
      if (i<S){
        message('updata A...',Sys.time())
        A <- t(apply(data[-nrow(data),],1,CompA,ch^i,W,Dis,data,Init))
        # updata para
        message('updata beta...',Sys.time())
        dict1 <- rep(FALSE,4000)
        beta <- t(apply(data[-nrow(data),],1,updatebeta,h=ch^i,W,Dis,data,Init,A))
        beta[dict1,] <- Init0[dict1,1:2]
        if (i>S0){
          D <- stopcriter(Init,Init3)
          dict0 <- D>qchisq(df = 2, p = 0.8)
          dict1 <- dict1|dict0 #保留所有超过的
          dict <- !dict1&dict0 #保留上一步超过的
          beta[dict,] <- Init0[dict,1:2]
        }
        Init[,1:2] <- beta
        message('updata cov...',Sys.time())## update cov
        espn <- t(apply(data[-nrow(data),],1,epsn,W,Dis,h=ch^i,Init,beta,data))
        Cov <- t(apply(data[-nrow(data),],1,updatecov,espn,h=ch^i,Dis,A,data))
        Init[,4:7] <- Cov
        if ( i==3) {Init3 <- Init}
        Init0 <- Init
        message('updata similar matrix...',Sys.time()) # update similar matrix
        Sim <- apply(Init[-nrow(data),],1,Smatrix,Init)  #相似矩阵
      }
    }
    return(W)
    }

