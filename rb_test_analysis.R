library('nloptr')
library(rgl)

optimal_weight_rb<- function(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf){
  rb_c_weight2<- function(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf){
    eval_f <- function(x) {
      sigma<- sqrt(t(x)%*%cov%*%x)
      return(-1* (ret%*%x-rf)/sigma) #sigma^2)
    }
    # Gradient
    eval_grad_f <- function(x) {
      sigma<- sqrt(t(x)%*%cov%*%x)
      return(-1* (ret/sigma - as.numeric((ret%*%x-rf)/sigma^3) * cov%*%x )) #2*cov%*%x)
    }
    # equality constraint function
    eval_g0 <- function(x) {
      return(sum(x)-1)
    }
    # jacobian of equality constraint
    eval_jac_g0 <- function(x) {
      return(rep(1,length(x)))
    }
    # inequality constraint function
    eval_g1 <- function(x) {
      var<- as.numeric(t(x)%*%cov%*%x)
      return(c(rw_lb - x*(cov%*%x)/var,x*(cov%*%x)/var - rw_ub))
    }
    # jacobian of inequality constraint
    eval_jac_g1 <- function(x) {
      #return n*n matrix
      #x<-w0
      var<- as.numeric(t(x)%*%cov%*%x)
      #i= 1~length(x)
      mat=NULL
      for (i in 1:length(x)){
        vec1<- cov[i,]*x[i]
        vec1[i]<- vec1[i]+(cov%*%x)[i]
        vec1<- vec1/var
        vec1<- t(vec1 - 2*x[i]*(cov%*%x)[i]/var^2 * (cov%*%x))
        mat= rbind(mat,vec1)
      }
      return(rbind(-mat,mat))
    }
    # solve fortfolio function
    res <- nloptr( x0=w0,
                   eval_f=eval_f,
                   eval_grad_f=eval_grad_f,
                   eval_g_eq=eval_g0,
                   eval_jac_g_eq=eval_jac_g0,
                   eval_g_ineq = eval_g1,
                   eval_jac_g_ineq = eval_jac_g1,
                   lb = w_lb,
                   ub = w_ub,
                   opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1.0e-8))
    #solution
    return(res$solution)
  }
  ###------------------------------------optimizing module-#
  w0<- rep(1/nrow(cov),nrow(cov))
  #start constraint
  dy_rw_lb<- rep(-1,nrow(cov))
  dy_rw_ub<- rep(1,nrow(cov))
  num<- 50
  for(i in 1:num){
    dy_rw_lb<- dy_rw_lb + (rw_lb - rep(-1,nrow(cov)))/num
    dy_rw_ub<- dy_rw_ub + (rw_ub - rep(1,nrow(cov)))/num
    w0<- rb_c_weight2(cov,w_lb,w_ub,dy_rw_lb,dy_rw_ub,w0,ret,rf)
    print(w0)
  }
  ###------------------------------------------------------#
  return(w0)
}

cov<- matrix(c( 1.0, 0.3,
                0.3, 1.5),nrow=2)
#2-d x1 x2
plot_sharpe<- function(cov,min,max){
    
  # initial values
  w0<- rep(1/nrow(cov),nrow(cov))
  # weight bounds
  w_lb<- rep(min,nrow(cov))
  w_ub<- rep(max,nrow(cov))
  # risk weight bounds
  rw_lb<- rep(0.2,nrow(cov))
  rw_ub<- rep(0.8,nrow(cov))
  #expected return
  ret<- c(0.1, 0.2)
  rf<- 0.01
  
  x1<- seq(w_lb[1],w_ub[1],length.out = 100)
  x2<- seq(w_lb[1],w_ub[1],length.out = 100)
  
  f<-function(x1,x2){
    sharpe<- (ret%*%c(x1,x2)-rf)/sqrt(t(c(x1,x2))%*%cov%*%c(x1,x2))
    return(sharpe)
  }
  
  mat<- matrix(nrow = length(x1),ncol = length(x2))
  for (i in 1:length(x1)){
    for (j in 1:length(x2))
      {mat[i,j] <- f(x1[i],x2[j])}
  }
  mat[is.infinite(mat)]<- 2*min(mat[is.finite(mat)])
  
  getmax<- function(mat){
    for (i in 1:length(x1)){
      for (j in 1:length(x2)){
          if(mat[i,j]==max(mat)){
            return(c(i,j))
          }
        }
    }
  }
  
  p_max<- list(x1.max=x1[getmax(mat)[1]],x2.max=x2[getmax(mat)[2]])
  print(p_max)
  plot3d(p_max[[1]],p_max[[2]],max(mat),col="red",size=5)
  surface3d(x1,x2,mat,col="grey")
}

plot_sharpe(cov,0,8)

#2-d    x1+x2=1
y= NULL
for (i in 1:length(x1)){
   y <- c(y,f(x1[i],1-x1[i]))
}
plot(x1,y,type="l",ylab="sharp ratio")


