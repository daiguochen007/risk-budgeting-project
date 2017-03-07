
library('nloptr')

#################### risk parity with x in [a,b], sum(x)=1 and least variance solution among all RP portfolios
cov<- matrix(c( 1.0,-0.9, 0.6,
               -0.9, 1.0,-0.2,
                0.6,-0.2, 4.0),nrow=3)

w_lb = rep(-1,nrow(cov))
w_ub = rep(2,nrow(cov))

rp_weight<- function(cov,rho,w_lb,w_ub,w0){
  eval_f <- function(x) {
    return(sum((x*(cov%*%x) - as.numeric(t(x)%*%cov%*%x/length(x)))^2) + rho*t(x)%*%cov%*%x)
  }
  # Gradient
  eval_grad_f <- function(x) {
    e<- diag(length(x))
    res<-0
    temp<- x*(cov%*%x) - as.numeric(t(x)%*%cov%*%x/length(x))
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

# initial values
w0<- rep(1/nrow(cov),nrow(cov))

rho<-100
beta<-0.6
while (rho >= 10^-6){
  rho<- rho*beta
  w0<- rp_weight(cov,rho,w_lb,w_ub,w0)
  print(w0)
}

#risk contribution of each asset
w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)

#volatility
sqrt(t(w0)%*%(cov%*%w0))




#################### risk budget with general bounds [a,b], sum(x)=1
cov<- matrix(c( 1.0,-0.9, 0.6,
               -0.9, 1.0,-0.2,
                0.6,-0.2, 4.0),nrow=3)
# initial values
w0<- rep(1/nrow(cov),nrow(cov))
b<- c(0.2,0.3,0.5)
w_lb = rep(-3,nrow(cov))
w_ub = rep(3,nrow(cov))

rb_weight<- function(cov,w0,b){
  eval_f <- function(x) {
    return(sum((x*(cov%*%x) - as.numeric(t(x)%*%cov%*%x)*b)^2)) #+ rho*t(x)%*%cov%*%x)
  }
  # Gradient
  eval_grad_f <- function(x) {
    e<- diag(length(x))
    res<-0
    temp<- x*(cov%*%x) - as.numeric(t(x)%*%cov%*%x)*b
    for (i in 1:length(x)){
      res<- res+2*temp[i]*(e[i,]%*%t(cov[i,])+t(e[i,]%*%t(cov[i,])))%*%x
    }
    return(res) # + rho*2*cov%*%x)
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
#1   0.422 0.475 0.104
w0<- c(1/3,1/3,1/3)
w0<- rb_weight(cov,w0,b)
#2   0.587 0.541 -0.128
w0<- c(0.7,0.4,-0.1)
w0<- rb_weight(cov,w0,b)
#3   1.36 -2.01 1.64
w0<- c(1.5,-2,1.5)
w0<- rb_weight(cov,w0,b)

#problem: too sensitive
w0<- c(0.567,0.561,-0.128)
w0<- c(0.587,0.541,-0.128)
w0<- c(0.587,0.521,-0.108)
#risk contribution of each asset
w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)

#volatility
sqrt(t(w0)%*%(cov%*%w0))

######################### least var rb start from min-variance
rb_weight2<- function(cov,rho,w_lb,w_ub,w0){
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
rho<-100
beta<-0.6
while (rho >= 10^-6){
  rho<- rho*beta
  w0<- rb_weight2(cov,rho,w_lb,w_ub,w0)
  print(w0)
}



################################## another test of sensitivity
cov<- matrix(c( 1.0, 0.3, 0.2,
                0.3, 1.5, 0.4,
                0.2, 0.4, 4.0),nrow=3)

w0<- c(0.344,0.35,0.306)
w0<- c(0.324,0.36,0.316)
#risk contribution of each asset
w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)



