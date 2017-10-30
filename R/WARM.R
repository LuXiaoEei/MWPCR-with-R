# WARM
# data <- esp_iid
# names(data)[4:ncol(data)] <- NAs
WARM <- function(data,
               ch = 1.2,
               Cn = log(100) * qchisq(df = 2, p = 0.95),
               Kloc = TRUE,
               Kst = TRUE,
               S0 = 3,
               S = 5,
               Dis){
    message('初始化...',Sys.time())
    #距离矩阵
    message('初始化距离矩阵...',Sys.time())
    # Dis <- as.matrix(dist(data[-nrow(data),c(1:3)],diag = TRUE,upper = TRUE)) 
    # SearchPoints0 <- SearchPoints3D0(h = ch^0)
    message('初始化beta...',Sys.time())
    Init <- t(apply(data[-nrow(data),],1,OLS_loc,data,Dis,ch^0))
    Init0 <- Init
    colnames(Init) <- c('beta0','beta1','res','cov00','cov01','cov10','cov11')
    W_loc <- diag(nrow(data)-1)
    W_st <- diag(nrow(data)-1)
    message('迭代开始...',Sys.time())
    dict1 <- rep(FALSE,4000)
    for (i in 1:S){
      message('开始第',i,'此迭代\n')
      message('updata WW...',Sys.time())
      if (Kloc){#距离
        W_loc <- (1-Dis/ch^i)*sign((1-Dis/ch^i)>=0)
      }
      if (Kst){#相似度
        message('update similar matrix...',Sys.time())
        Sim <- apply(Init,1,Smatrix,Init)#相似矩阵
        W_st <- exp(-Sim/Cn)
      }
      WW <- W_loc*W_st
      if (i<S){
        message('updata A...',Sys.time())
        A <- t(apply(data[-nrow(data),],1,CompA,ch^i,WW,Dis,data,Init))
        # updata para
        message('updata beta...',Sys.time())
 
        beta <- t(apply(data[-nrow(data),],1,updatebeta,h=ch^i,WW,Dis,data,Init,A))
        beta[dict1,] <- Init0[dict1,1:2]#停止的不进行迭代
        if (i>S0){
          D <- stopcriter(Init,Init3)
          dict0 <- D>qchisq(df = 2, p = 0.8)
          dict <- (!dict1)&dict0 #保留上一步超过的
          dict1 <- dict1|dict0 #保留所有超过的
          beta[dict,] <- Init0[dict,1:2]#停止
        }
        Init[,1:2] <- beta
        message('updata cov...',Sys.time())## update cov
        espn <- t(apply(data[-nrow(data),],1,epsn,WW,Dis,h=ch^i,Init,beta,data))
        Cov <- t(apply(data[-nrow(data),],1,updatecov,espn,A,data))
        Init[,4:7] <- Cov
        if ( i==S0) {Init3 <- Init}
        Init0 <- Init
      }
    }
    message('迭代结束...',Sys.time()) 
    return(WW)
    }

