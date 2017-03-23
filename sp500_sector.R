library(quantmod)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(nloptr)

sectors<- c("finan","discre","health","infotech","utility","industry","material","staple","telecom","energy","realestate")

all_data<-NULL
for (i in 1:11){
  data<- read.csv(paste("/Users/guochendai/Desktop/4th semester/7043 capstone project/data/sec_",
                        sectors[i],".csv",sep=""))
  data<- xts(data[,2:3],order.by=as.Date(data[,1], format="%m/%d/%y"))
  #plot(data[,2],main=sectors[i])
  all_data<- cbind(all_data, data[,2]/as.numeric(data[1,2]))
}
colnames(all_data)<- sectors
all_data<- na.omit(all_data)
for (i in 1:ncol(all_data)){
  all_data[,i]<- all_data[,i]/as.numeric(all_data[1,i])
}
#autoplot(all_data,facets = FALSE)

############################################# ranged risk budget tool

w0<- rep(1/nrow(cov),nrow(cov))
# weight bounds
w_lb<- rep(-1,nrow(cov))
w_ub<- rep(2,nrow(cov))
# risk weight bounds
rw_lb<- rep(0.05,nrow(cov))
rw_ub<- rep(0.3,nrow(cov))
#rf
rf<- 0.03
#tool func
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
# return func from price
ret_xts<- function(data){
  for (i in 1:ncol(data)){
    data[,i]<- data[,i]/lag(data[,i])-1
  } 
  return(na.omit(data))
}

############################################# parameters
start_date = "2001/10/9"
end_date = "2017/3/19"
insam_length = 90
outsam_length = 30

all_data<- all_data[which(index(all_data)>=start_date&index(all_data)<=end_date),]

if (nrow(all_data)<insam_length+outsam_length){
  cat("The length of all_data is shorter than insam + outsam!")
}else{
  outsam_index_list <- seq(0,(nrow(all_data)-insam_length)%/%outsam_length-1,1)*outsam_length+insam_length+1
  outsam_index_list <- cbind(outsam_index_list,outsam_index_list+outsam_length-1)
  insam_index_list <- outsam_index_list - outsam_length
  insam_index_list[,1] <- insam_index_list[,1] - insam_length + outsam_length
}

# main loop
opt_port=NULL
for (j in 1:nrow(insam_index_list)){
  ############################################# insample data
  #seperate time period
  insam_data<- all_data[insam_index_list[j,1]:insam_index_list[j,2],]
  for (i in 1:ncol(insam_data)){
    insam_data[,i]<- insam_data[,i]/as.numeric(insam_data[1,i])
  }
  #autoplot(insam_data,facets = FALSE)
  
  # #simple mean ret(very bad!!! should not use)
  # ret<- (apply(insam_data,2,mean)-1)*252
  #geo ret
  ret<- as.vector((insam_data[nrow(insam_data),]^(1/nrow(insam_data))-1)*252)
  cov<- cov(insam_data)*252     #result will become bad with cov of ret
  #pheatmap(cov,cluster_rows = F,cluster_cols = F)
  
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
  
  # #risk contribution of each asset
  # w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)
  # #volatility
  # sqrt(t(w0)%*%(cov%*%w0))
  # #sharpe
  # (ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0)
  # 
  # #plot insample optimal
  # port<- as.matrix(insam_data)%*%w0
  # port<- xts(cbind(port, apply(insam_data,1,mean)),order.by = index(insam_data))
  # colnames(port)<- c("insam_optimal","insam_average")
  # autoplot(port,facets = F)
  
  
  ############################################# outsample data
  #seperate time period
  outsam_data<- all_data[outsam_index_list[j,1]:outsam_index_list[j,2],]
  for (i in 1:ncol(outsam_data)){
    outsam_data[,i]<- outsam_data[,i]/as.numeric(outsam_data[1,i])
  }
  #autoplot(outsam_data,facets = FALSE)
  
  # ret<- as.vector((outsam_data[nrow(outsam_data),]^(1/nrow(outsam_data))-1)*252)
  # cov<- cov(outsam_data)*252
  # rf<- 0.03
  # #risk contribution of each asset
  # w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)
  # #volatility
  # sqrt(t(w0)%*%(cov%*%w0))
  # #sharpe
  # (ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0)
  
  ##plot insample optimal
  # port<- as.matrix(outsam_data)%*%w0
  # port<- xts(cbind(port, apply(outsam_data,1,mean)),order.by = index(outsam_data))
  # colnames(port)<- c("outsam_optimal","outsam_average")
  # print(autoplot(port,facets = F,main=j))
  port<- xts(as.matrix(outsam_data)%*%w0,order.by = index(outsam_data))
  if (is.null(opt_port)){
    opt_port = port
  }else{
    opt_port = rbind(opt_port,port*as.numeric(opt_port[length(opt_port)]))
  }
}

#plot result
result<- cbind(opt_port,xts(apply(all_data,1,mean),order.by = index(all_data)))
result<- na.omit(result)
result[,2]<- result[,2]/as.numeric(result[1,2])
colnames(result)<-c("opt_pf","average")
autoplot(result, facets = F)

################### perf analysis
get_perf<- function(vec,rf){
    id<- c("annual_ret","annual_vol","sharpe_ratio","maxdrawdown","maxdown_time(y)")
    mat<- matrix(NA,length(id),1)
    rownames(mat)<-id
    colnames(mat)<- colnames(vec)
    mat[1,1]<- as.numeric((vec[length(vec)]^(1/length(vec))-1)*252)
    mat[2,1]<- sd(vec)*sqrt(252/length(vec))
    mat[3,1]<- (mat[1,1]-rf)/mat[2,1]
    
    v_max<- 0
    mdd<- 1
    mdt<- 0
    for (i in 1:length(vec)){
      v_max<- max(vec[i],v_max)
      mdd<- min(mdd,vec[i]/v_max)
      if (vec[i]<v_max){
        dt<- dt+1
      }else{
        dt=0
      }
      mdt<- max(mdt,dt)
    }
    mat[4,1]<- 1-mdd
    mat[5,1]<- mdt/252
    return(mat)
}

cbind(get_perf(result[,1],rf),get_perf(result[,2],rf))

