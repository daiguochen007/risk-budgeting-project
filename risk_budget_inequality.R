library('nloptr')

################################## risk weight inequality constraints
cov<- matrix(c(   0.01,  0.009, -0.005,
                 0.009, 0.0225,   0.01,
                -0.005,   0.01,   0.04),nrow=3)
# initial values
w0<- rep(1/nrow(cov),nrow(cov))
# weight bounds
w_lb<- rep(-3,nrow(cov))
w_ub<- rep(3,nrow(cov))
# risk weight bounds
rw_lb<- rep(0.15,nrow(cov))
rw_ub<- rep(0.6,nrow(cov))
#expected return
ret<- c(0.1, 0.15, 0.25)
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

w0<- c(1.75,-0.5,-0.25)
#start constraint
dy_rw_lb<- rep(-1,nrow(cov))
dy_rw_ub<- rep(1,nrow(cov))

#rb_c_weight2(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf)

num<- 50
for(i in 1:num){
  dy_rw_lb<- dy_rw_lb + (rw_lb - rep(-1,nrow(cov)))/num
  dy_rw_ub<- dy_rw_ub + (rw_ub - rep(1,nrow(cov)))/num
  w0<- rb_c_weight2(cov,w_lb,w_ub,dy_rw_lb,dy_rw_ub,w0,ret,rf)
  cat(i)
  print(w0)
  print((ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0))
}

#risk contribution of each asset
w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)
#volatility
sqrt(t(w0)%*%(cov%*%w0))
#sharpe
(ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0)


############################### start from no bounds and min variance portfolio
#### 
rb_c_weight3<- function(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf){
  eval_f <- function(x) {
    return(t(x)%*%cov%*%x)
  }
  # Gradient
  eval_grad_f <- function(x) {
    return(2*cov%*%x)
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

w0<- c(1.75,-0.5,-0.25)
#start constraint
dy_rw_lb<- rep(-1,nrow(cov))
dy_rw_ub<- rep(1,nrow(cov))

#rb_c_weight3(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf)

num<- 50
for(i in 1:num){
  dy_rw_lb<- dy_rw_lb + (rw_lb - rep(-1,nrow(cov)))/num
  dy_rw_ub<- dy_rw_ub + (rw_ub - rep(1,nrow(cov)))/num
  w0<- rb_c_weight3(cov,w_lb,w_ub,dy_rw_lb,dy_rw_ub,w0,ret,rf)
  print(w0)
}

#risk contribution of each asset
w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)
#volatility
sqrt(t(w0)%*%(cov%*%w0))
#sharpe
(ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0)

#################################### 
# use global assets
# normal cov
cov<- cov(ret_xts(all_data[which(index(all_data)<"2008-09-01"&index(all_data)>"2006-09-01"),]))*252
# crisis cov
cov<- cov(ret_xts(all_data[which(index(all_data)>"2008-09-01"&index(all_data)<"2009-09-01"),]))*252
diag(cov)^0.5 # the stock is too volatile but bonds are not

w_lb<- rep(0,nrow(cov))
w_ub<- rep(1,nrow(cov))

opt_weight_rp2<- function(cov,w_lb,w_ub,adjust=F){
  opt_weight_rp<- function(cov,w_lb,w_ub){
    rb_weight2<- function(cov,rho,w_lb,w_ub,w0){
      #risk parity/b can be adjusted to risk weights
      b<- rep(1/nrow(cov),nrow(cov))
      eval_f <- function(x) {
        return(sum((x*(cov%*%x) - as.numeric(t(x)%*%cov%*%x)*b)^2) + rho*t(x)%*%cov%*%x)
      }
      # Gradient
      eval_grad_f <- function(x) {
        e<- diag(length(x))
        res<-0
        temp<- x*(cov%*%x) - as.numeric(t(x)%*%cov%*%x)*b
        for (i in 1:length(x)){
          res<- res+2*temp[i]*(e[i,]%*%t(cov[i,])+t(e[i,]%*%t(cov[i,])))%*%x
        }
        return(res + rho*2*cov%*%x)
      }
      # equality constraint function
      eval_g0 <- function(x) {
        return(sum(x)-1)
      }
      # jacobian of equality constraint
      eval_jac_g0 <- function(x) {
        return(rep(1,length(x)))
      }
      # solve fortfolio function
      res <- nloptr( x0=w0,
                     eval_f=eval_f,
                     eval_grad_f=eval_grad_f,
                     eval_g_eq=eval_g0,
                     eval_jac_g_eq=eval_jac_g0,
                     lb = w_lb,
                     ub = w_ub,
                     opts=list("algorithm"="NLOPT_LD_SLSQP","xtol_rel"=1.0e-8))
      #solution
      return(res$solution)
    }
    w0<- rep(1/nrow(cov),nrow(cov))
    rho<-100
    beta<-0.6
    while (rho >= 10^-6){
      rho<- rho*beta
      w0<- rb_weight2(cov,rho,w_lb,w_ub,w0)
      #print(w0)
    }
    return (w0)
  }
  w0<- opt_weight_rp(cov,w_lb,w_ub)
  if (adjust==F){return(w0)}
  
  b<- rep(1/nrow(cov),nrow(cov))
  #evaluate deviation 2% allowed
  dev = sum(((w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)/b)-1)^2)
  benchmk = sum(rep(0.05,nrow(cov))^2)
  while(dev > benchmk && w_ub[1]<3){
    w_lb<- w_lb-1
    w_ub<- w_ub+1
    w0<- opt_weight_rp(cov,w_lb,w_ub)
    dev = sum(((w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)/b)-1)^2)
  }
  cat(w_lb[1])
  cat(w_ub[1])
  return(w0)
}

opt_weight_rp2(cov,w_lb,w_ub)




