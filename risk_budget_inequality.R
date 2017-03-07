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
rw_lb<- rep(0,nrow(cov))
rw_ub<- rep(0.5,nrow(cov))
#expected return
ret<- c(0.1, 0.2, 0.5)
rf<- 0.01


rb_c_weight<- function(cov,rho,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf){
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
    sigma<- as.numeric(sqrt(t(x)%*%cov%*%x))
    return(rw_lb-x*(cov%*%x)/sigma)
  }
  # jacobian of inequality constraint
  eval_jac_g1 <- function(x) {
    #
    var<- as.numeric(t(x)%*%cov%*%x)
    return((cov%*%x+diag(cov)*x)/var-2*x*(cov%*%x)*(cov%*%x)/var^2)
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
                 opts=list("algorithm"="NLOPT_LD_MMA","xtol_rel"=1.0e-8))
  #solution
  return(res$solution)
}

w0<- rb_c_weight(cov,rho,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf)





