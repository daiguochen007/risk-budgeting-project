library('nloptr')

################################## risk weight inequality constraints
cov<- matrix(c( 1.0, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 4.0),nrow=3)
# initial values
w0<- rep(1/nrow(cov),nrow(cov))
# weight bounds
w_lb<- rep(-3,nrow(cov))
w_ub<- rep(3,nrow(cov))
# risk weight bounds
rw_lb<- rep(0.15,nrow(cov))
rw_ub<- rep(0.5,nrow(cov))
#expected return
ret<- c(0.1, 0.2, 0.5)
rf<- 0.01


rb_c_weight<- function(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf){
  eval_f <- function(x) {
    sigma<- sqrt(t(x)%*%cov%*%x)
    return(-1* (ret%*%x-rf)/sigma )
  }
  # Gradient
  eval_grad_f <- function(x) {
    sigma<- sqrt(t(x)%*%cov%*%x)
    return( -1* (ret/sigma - as.numeric((ret%*%x-rf)/sigma^3) * cov%*%x ))
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

w0<- rb_c_weight(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf)

#risk contribution of each asset
w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)
#volatility
sqrt(t(w0)%*%(cov%*%w0))
#sharpe
(ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0)

############################### start from no bounds and max sharpe ratio portfolio
#### not stable with fixed constraint moving from min var
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

w0<- c(-0.5,1.75,-0.25)
#start constraint
dy_rw_lb<- rep(-1,nrow(cov))
dy_rw_ub<- rep(1,nrow(cov))

#rb_c_weight2(cov,w_lb,w_ub,rw_lb,rw_ub,rep(1/nrow(cov),nrow(cov)),ret,rf)
num<- 50
for(i in 1:num){
  dy_rw_lb<- dy_rw_lb + (rw_lb - rep(-1,nrow(cov)))/num
  dy_rw_ub<- dy_rw_ub + (rw_ub - rep(1,nrow(cov)))/num
  w0<- rb_c_weight2(cov,w_lb,w_ub,dy_rw_lb,dy_rw_ub,w0,ret,rf)
  print(w0)
}

#risk contribution of each asset
w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)
#volatility
sqrt(t(w0)%*%(cov%*%w0))
#sharpe
(ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0)



