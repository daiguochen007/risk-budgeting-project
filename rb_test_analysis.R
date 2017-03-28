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
    #print(w0)
  }
  ###------------------------------------------------------#
  return(w0)
}

f<-function(w,type="sharpe"){
  sharpe<- (ret%*%w-rf)/sqrt(t(w)%*%cov%*%w)
  RC<-w*(cov%*%w)/as.numeric(t(w)%*%cov%*%w)
  var<- t(w)%*%cov%*%w
  if(type=="sharpe"){return(sharpe)}
  if(type=="RC1"){return(RC[1])}
  if(type=="RC2"){return(RC[2])}
  if(type=="RC3"){return(RC[3])}
  if(type=="var"){return(var)}
}

################################################### 2D
cov<- matrix(c( 0.02, 0.01,
                0.01, 0.05),nrow=2)
ret<- c(0.1, 0.2)
rf<- 0.01
#2-d x1 x2
plot_2d<- function(min,max,type="sharpe"){
  # initial values
  w0<- rep(1/nrow(cov),nrow(cov))
  # weight bounds
  w_lb<- rep(min,nrow(cov))
  w_ub<- rep(max,nrow(cov))
  # risk weight bounds
  rw_lb<- rep(0.2,nrow(cov))
  rw_ub<- rep(0.8,nrow(cov))
  
  x1<- seq(w_lb[1],w_ub[1],length.out = 100)
  x2<- seq(w_lb[1],w_ub[1],length.out = 100)
  
  mat<- matrix(nrow = length(x1),ncol = length(x2))
  for (i in 1:length(x1)){
    for (j in 1:length(x2))
      {mat[i,j] <- f(c(x1[i],x2[j]),type)}
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
  #w0<- optimal_weight_rb(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf)
  
  p_max<- list(max.x1=x1[getmax(mat)[1]],max.x2=x2[getmax(mat)[2]])
  #pts=list(x1=c(p_max[[1]],w0[1]),x2=c(p_max[[2]],w0[2]),sharpe=c(max(mat),f(w0[1],w0[2])))
  #points3d(pts$x1,pts$x2,pts$sharpe, size=15,color = c("green","red"),xlab="x1",ylab="x2",zlab="z") 
  if(type=="sharpe"){
    print(p_max)
    points3d(p_max[[1]],p_max[[2]],max(mat), size=15,color ="green")
  }
  surface3d(x1,x2,mat,col="grey")
  axes3d()
  title3d(xlab="x1",ylab="x2",zlab=type)
}

plot_2d(-2,5,"sharpe")

#2-d    x1+x2=1
y= NULL
rc=NULL
for (i in 1:length(x1)){
   y <- c(y,f(c(x1[i],1-x1[i])))
   rc<- c(rc,f(c(x1[i],1-x1[i]),type = "RC1"))
}
plot(x1,y,type="l",xlab="x1",ylab="sharp ratio")
plot(x1,rc,type="l",xlab="x1",ylab="RC1")


################################################### 3D
cov<- matrix(c( 0.02 , 0.01 , 0.005,
                0.01 , 0.05 , 0.015,
                0.005, 0.015, 0.06 ),nrow=3)
ret<- c(0.1, 0.15, 0.2)
rf<- 0.01
#3-d x1+x2+x3=1
plot_3d<- function(w_min=-2,w_max=2,rw_min=0.2,rw_max=0.8,type="sharpe"){
  # initial values
  w0<- rep(1/nrow(cov),nrow(cov))
  # weight bounds
  w_lb<- rep(w_min,nrow(cov))
  w_ub<- rep(w_max,nrow(cov))
  # risk weight bounds
  rw_lb<- rep(rw_min,nrow(cov))
  rw_ub<- rep(rw_max,nrow(cov))
  
  x1<- seq(w_lb[1],w_ub[1],length.out = 100)
  x2<- seq(w_lb[1],w_ub[1],length.out = 100)
  
  mat<- matrix(nrow = length(x1),ncol = length(x2))
  for (i in 1:length(x1)){
    for (j in 1:length(x2))
    {mat[i,j] <- f(c(x1[i],x2[j],1-x1[i]-x2[j]),type)}
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
  getmin<- function(mat){
    for (i in 1:length(x1)){
      for (j in 1:length(x2)){
        if(mat[i,j]==min(mat)){
          return(c(i,j))
        }
      }
    }
  }
  w0<- optimal_weight_rb(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf)
   
  if(type=="sharpe"){
    p_max<- list(max.x1=x1[getmax(mat)[1]],max.x2=x2[getmax(mat)[2]])
    print(p_max)
    pts=list(x1=c(p_max[[1]],w0[1]),x2=c(p_max[[2]],w0[2]),sharpe=c(max(mat),f(w0)))
    points3d(pts$x1,pts$x2,pts$sharpe, size=15,color = c("green","red")) 
  }else if(type=="var"){
    p_min<- list(min.x1=x1[getmin(mat)[1]],min.x2=x2[getmin(mat)[2]])
    print(p_min)
    points3d(p_min[[1]],p_min[[2]],min(mat), size=15,color = c("green"))
  }else{
    surface3d(x1,x2,matrix(rw_min,ncol = length(x1),nrow = length(x2)),col="blue",alpha=0.5)
    surface3d(x1,x2,matrix(rw_max,ncol = length(x1),nrow = length(x2)),col="purple",alpha=0.5)
  }
  surface3d(x1,x2,mat,col="grey")
  axes3d()
  title3d(xlab="x1",ylab="x2",zlab=type)
}

plot_3d(-2,2,0.2,0.8,"var")

