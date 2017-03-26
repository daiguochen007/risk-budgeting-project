library(quantmod)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(nloptr)

############################################# tool functions

#get data func
get_alldata<- function(namelist,path){
  all_data<-NULL
  for (i in 1:length(namelist)){
    data<- read.csv(paste(path, namelist[i],".csv",sep=""))
    data<- xts(data[,2:3],order.by=as.Date(data[,1], format="%m/%d/%y"))
    all_data<- cbind(all_data, data[,2]/as.numeric(data[1,2]))
  }
  colnames(all_data)<- namelist
  all_data<- na.omit(all_data)
  # for (i in 1:ncol(all_data)){
  #   all_data[,i]<- all_data[,i]/as.numeric(all_data[1,i])
  # }
  return(all_data)
}

#risk budget optimization func
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

# get return from price
ret_xts<- function(data,log=FALSE){
  if (log==FALSE){
    for (i in 1:ncol(data)){
      data[,i]<- data[,i]/lag(data[,i])-1
    }
    data[is.na(data)]<- 0
    return(data)
  }else{
    for (i in 1:ncol(data)){
      data[,i]<- log(data[,i]/lag(data[,i]))
    }
    data[is.na(data)]<- 0
    return(data)
  }
}

# get average risk free rate within a period
get_rf<- function(start,end){
  rf_vec<- rf_ts[which(index(rf_ts)>=start&index(rf_ts)<=end),]
  return(mean(rf_vec)/100)
}

# get performance analysis
get_perf<- function(xts_df,rf_ts){
  #vec func
  get_perf_vec<- function(vec,rf_ts){
    id<- c("annual_ret","annual_vol","sharpe_ratio","maxdrawdown","maxdown_time(y)")
    mat<- matrix(NA,length(id),1)
    rownames(mat)<-id
    colnames(mat)<- colnames(vec)
    mat[1,1]<- sum(ret_xts(vec,log=T))*252/length(vec)
    mat[2,1]<- sd(ret_xts(vec,log=T))*sqrt(252)
    mat[3,1]<- (mat[1,1] - get_rf(index(vec)[1],index(vec)[length(vec)]))/mat[2,1]
    
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
  
  if (ncol(xts_df)==1){
    return(get_perf_vec(xts_df,rf_ts))
  }else{
    res<-NULL
    for (i in 1:ncol(xts_df)){
      res<- cbind(res,get_perf_vec(xts_df[,i],rf_ts))
    }
    return(res)
  }
}

# back-test program
rb_back_test<- function(all_data,start_date,end_date,insam_length,outsam_length,
                        w_lb,w_ub,rw_lb,rw_ub,rf_ts){
  all_data<- all_data[which(index(all_data)>=start_date&index(all_data)<=end_date),]
  avg_port<- xts(cumprod(apply(ret_xts(all_data),1,mean)+1),order.by = index(all_data))
  
  if (nrow(all_data)<insam_length+outsam_length){
    cat("The length of all_data is shorter than insam + outsam!")
    return("Function terminated")
  }else{
    outsam_index_list <- seq(0,(nrow(all_data)-insam_length)%/%outsam_length-1,1)*outsam_length+insam_length+1
    outsam_index_list <- cbind(outsam_index_list,outsam_index_list+outsam_length-1)
    insam_index_list <- outsam_index_list - outsam_length
    insam_index_list[,1] <- insam_index_list[,1] - insam_length + outsam_length
  }
  
  cat("start backtesting...\n")
  # core loop
  # serval things we are interested: weights, risk contribution(in/out), performance(in/out), sharpe(i/o), benchmark sharpe(i/o)
  outsam_perf=NULL
  insam_perf=NULL
  weight_ts=NULL
  insamRC_ts = NULL
  outsamRC_ts = NULL
  insam_sharpe_ts = NULL
  outsam_sharpe_ts = NULL
  avg_insam_sharpe_ts = NULL
  avg_outsam_sharpe_ts = NULL
  for (j in 1:nrow(insam_index_list)){
    ############################################# insample data
    #seperate time period
    insam_data<- all_data[insam_index_list[j,1]:insam_index_list[j,2],]
    insam_benchmrk<- avg_port[insam_index_list[j,1]:insam_index_list[j,2],]
    #get insam rf rate
    rf<- get_rf(index(all_data[insam_index_list[j,1],]),index(all_data[insam_index_list[j,2],]))
    #geo insam ret
    #ret<- ((as.vector(insam_data[nrow(insam_data),])/as.vector(insam_data[1,]))^(1/nrow(insam_data))-1)*252
    # log ret
    ret<- apply(ret_xts(insam_data,log = T),2,mean)*252
    
    cov<- cov(ret_xts(insam_data,log = T))*252
    #pheatmap(cov,cluster_rows = F,cluster_cols = F)
    
    ###------------------------------------optimizing module-#
    w0<- optimal_weight_rb(cov,w_lb,w_ub,rw_lb,rw_ub,w0,ret,rf)
    ###------------------------------------------------------#
    
    #weight time series
    weight_ts<- rbind(weight_ts,w0)
    #insam risk contribution ts
    insamRC_ts<- rbind(insamRC_ts, t(w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)))
    #insam sharp ratio ts
    insam_sharpe_ts<- c(insam_sharpe_ts,(ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0))
    #avg insam benchmrk sharp ratio ts
    avg_insam_sharpe_ts<- c(avg_insam_sharpe_ts,(sum(ret_xts(insam_benchmrk,log = T))-rf)/(sd(ret_xts(insam_benchmrk,log = T))*sqrt(252)))
    
    ############################################# outsample data
    #seperate time period
    outsam_data<- all_data[outsam_index_list[j,1]:outsam_index_list[j,2],]
    outsam_benchmrk<- avg_port[outsam_index_list[j,1]:outsam_index_list[j,2],]
    #outsam risk contribution ts
    cov<- cov(ret_xts(outsam_data,log = T))*252 
    outsamRC_ts<- rbind(outsamRC_ts, t(w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)))
    #outsam sharp ratio ts
    rf<- get_rf(index(all_data[outsam_index_list[j,1],]),index(all_data[outsam_index_list[j,2],]))
    ret<- ((as.vector(outsam_data[nrow(outsam_data),])/as.vector(outsam_data[1,]))^(1/nrow(outsam_data))-1)*252
    outsam_sharpe_ts<- c(outsam_sharpe_ts,(ret%*%w0-rf)/sqrt(t(w0)%*%cov%*%w0))
    #avg outsam benchmrk sharp ratio ts
    avg_outsam_sharpe_ts<- c(avg_outsam_sharpe_ts,(sum(ret_xts(outsam_benchmrk,log = T))-rf)/(sd(ret_xts(outsam_benchmrk,log = T))*sqrt(252)))
    
    #insam portfolio performance
    in_perf<- xts(cumprod(as.matrix(ret_xts(insam_data))%*%w0+1),order.by = index(insam_data))
    #outsam portfolio performance
    out_perf<- xts(cumprod(as.matrix(ret_xts(outsam_data))%*%w0+1),order.by = index(outsam_data))
    if (is.null(outsam_perf)){
      outsam_perf = out_perf
      insam_perf = in_perf
    }else{
      outsam_perf = rbind(outsam_perf,out_perf*as.numeric(outsam_perf[length(outsam_perf)]))
      insam_perf = rbind(insam_perf,in_perf*as.numeric(insam_perf[length(insam_perf)]))
    }
  }
  cat("Back test finished!\n")
  
  # time series of weights, RC, sharpe
  weight_ts<- xts(weight_ts,order.by = index(all_data)[outsam_index_list[,1]])
  colnames(weight_ts)<- colnames(all_data)

  insamRC_ts<- xts(insamRC_ts,order.by = index(all_data)[outsam_index_list[,1]])
  colnames(insamRC_ts)<- colnames(all_data)

  outsamRC_ts<- xts(outsamRC_ts,order.by = index(all_data)[outsam_index_list[,1]])
  colnames(outsamRC_ts)<- colnames(all_data)

  insam_sharpe_ts<- xts(insam_sharpe_ts,order.by = index(all_data)[outsam_index_list[,1]])
  outsam_sharpe_ts<- xts(outsam_sharpe_ts,order.by = index(all_data)[outsam_index_list[,1]])
  sharpe_ts<- cbind(insam_sharpe_ts,avg_insam_sharpe_ts,outsam_sharpe_ts,avg_outsam_sharpe_ts)
  colnames(sharpe_ts)<- c("insample","insample benchmark","outsample","outsample benchmark")
  #time series of insam result
  in_result<- cbind(insam_perf,avg_port)
  in_result<- na.omit(in_result)
  in_result[,2]<- in_result[,2]/as.numeric(in_result[1,2])
  colnames(in_result)<-c("optimal portfolio","average")

  #time series of outsam result
  out_result<- cbind(outsam_perf,avg_port)
  out_result<- na.omit(out_result)
  out_result[,2]<- out_result[,2]/as.numeric(out_result[1,2])
  colnames(out_result)<-c("optimal portfolio","average")

  return(list("weight_ts"=weight_ts,
              "insamRC_ts"=insamRC_ts,
              "outsamRC_ts"=outsamRC_ts,
              "in_result"=in_result,
              "out_result"=out_result,
              "sharpe"=sharpe_ts))
}

#show back test result: plot&table
show_result<- function(res, rf_ts){
  print(autoplot(res$weight_ts,facets = F,main="time series of optimized weights"))
  print(autoplot(res$insamRC_ts,facets = F,main="time series of insample risk contribution"))
  print(autoplot(res$outsamRC_ts,facets = F,main="time series of outsample risk contribution"))
  print(autoplot(res$in_result, facets = F,main="time series of insample overall performance"))
  print(autoplot(res$out_result, facets = F,main="time series of outsample overall performance"))
  print(autoplot(res$sharpe[,1:2],facets = F,main="time series of insample sharpe ratio"))
  print(autoplot(res$sharpe[,3:4],facets = F,main="time series of outsample sharpe ratio"))
  
  print("insample perf analysis:")
  print(get_perf(res$in_result,rf_ts))
  print("outsample perf analysis:")
  print(get_perf(res$out_result,rf_ts))
}

############################################# risk free rate
rf_ts<- read.csv("/Users/guochendai/Desktop/4th semester/7043 capstone project/data/tbill3m.csv")
rf_ts<- xts(rf_ts[,2],order.by=as.Date(rf_ts[,1], format="%m/%d/%Y"))

############################################# sp500 sector data
sectors<- c("finan","discre","health","infotech","utility","industry","material","staple","telecom","energy","realestate")
path<- "/Users/guochendai/Desktop/4th semester/7043 capstone project/data/sec_"
all_data<- get_alldata(sectors,path)

############################################# asset group data
assets<- c("spx","dow","crb","us_agg","reit","nasdaq","euro_50","ftse100","cac40","dax","nikkei","hsi","csi300","euro_agg","asianpacific_agg","global_highyield")
path<- "/Users/guochendai/Desktop/4th semester/7043 capstone project/data/"
all_data<- get_alldata(assets,path)
# autoplot(all_data,facets = FALSE)
# cov<- cov(all_data)/nrow(all_data)*252
# pheatmap(cov,cluster_rows = F,cluster_cols = F)

############################################# parameters
start_date = index(all_data)[1]
end_date = "2017/3/19"
insam_length = 90
outsam_length = 30
# weight bounds
w_lb<- rep(-2,ncol(all_data))
w_ub<- rep(2,ncol(all_data))
# risk weight bounds
rw_lb<- rep(0,ncol(all_data))
rw_ub<- rep(0.5,ncol(all_data))

############################################ run
res<- rb_back_test(all_data,start_date,end_date,insam_length,outsam_length,w_lb,w_ub,rw_lb,rw_ub,rf_ts)
show_result(res)


