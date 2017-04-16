######### for weekly and monthly tests
# to weekly sample
all_data<- all_data[endpoints(all_data,on = "weeks"),]
rf_ts<- rf_ts[endpoints(rf_ts,on = "weeks"),]

#back test program
rb_back_test<- function(all_data,start_date,end_date,insam_length,outsam_length,
                        w_lb,w_ub,rw_lb,rw_ub,rf_ts,mv_target=NULL){
  all_data<- all_data[which(index(all_data)>=start_date&index(all_data)<=end_date),]
  
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
  weight_ts=NULL
  insamRC_ts = NULL
  outsamRC_ts = NULL
  #outsam performance
  outsam_perf=NULL
  mv_outsam_perf = NULL
  avg_outsam_perf = NULL
  ms_outsam_perf= NULL
  rp_oursam_perf = NULL
  for (j in 1:nrow(insam_index_list)){
    ############################################# insample data
    #seperate time period
    insam_data<- all_data[insam_index_list[j,1]:insam_index_list[j,2],]
    #get insam rf rate
    rf<- get_rf(index(all_data[insam_index_list[j,1],]),index(all_data[insam_index_list[j,2],]))
    #geo insam ret
    #ret<- ((as.vector(insam_data[nrow(insam_data),])/as.vector(insam_data[1,]))^(1/nrow(insam_data))-1)*252
    # log ret / simple ret
    ret<- apply(ret_xts(insam_data),2,mean)*50 # log =T
    
    cov<- cov(ret_xts(insam_data))*50  # log = T
    
    ###------------------------------------core optimizing module-#
    w0<- optimal_weight_mv_rb(cov,w_lb,w_ub,rw_lb,rw_ub,ret,rf)           #optimal
    #w0<- opt_weight_rp2(cov,w_lb,w_ub,adjust = F)
    #w0<- rep(1/nrow(cov),nrow(cov))
    w1<- optimal_weight_minvar(cov,w_lb,w_ub,ret,target=mv_target)     #min-var
    w2<- optimal_weight_ms(cov,w_lb,w_ub,ret,rf)                       #max-sharpe
    wn<- rep(1/nrow(cov),nrow(cov))                                    #1/n
    wp<- opt_weight_rp2(cov,w_lb,w_ub,adjust = F)                      #risk parity
    ###-----------------------------------------------------------#
    #if problem with w0,replace
    #if(sum(w0)>1.1||sum(w0)<0.9){
    #  w0<- rep(1/nrow(cov),nrow(cov))
    #}
    
    #weight time series
    weight_ts<- rbind(weight_ts,w0)
    #insam risk contribution ts
    insamRC_ts<- rbind(insamRC_ts, t(w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)))

    ############################################# outsample data
    #seperate time period
    outsam_data<- all_data[outsam_index_list[j,1]:outsam_index_list[j,2],]
    #outsam risk contribution ts
    cov<- cov(ret_xts(outsam_data))*50  #log = T
    outsamRC_ts<- rbind(outsamRC_ts, t(w0*(cov%*%w0)/as.numeric(t(w0)%*%cov%*%w0)))
    #outsam sharp ratio ts
    rf<- get_rf(index(all_data[outsam_index_list[j,1],]),index(all_data[outsam_index_list[j,2],]))
    ret<- apply(ret_xts(outsam_data),2,mean)*50  # log = T
    
    #outsam portfolio performance
    out_perf<- xts(cumprod(as.matrix(ret_xts(outsam_data))%*%w0+1),order.by = index(outsam_data))
    mv_out_perf<- xts(cumprod(as.matrix(ret_xts(outsam_data))%*%w1+1),order.by = index(outsam_data))
    avg_out_perf<- xts(cumprod(as.matrix(ret_xts(outsam_data))%*%wn+1),order.by = index(outsam_data))
    ms_out_perf<- xts(cumprod(as.matrix(ret_xts(outsam_data))%*%w2+1),order.by = index(outsam_data))
    rp_out_perf<- xts(cumprod(as.matrix(ret_xts(outsam_data))%*%wp+1),order.by = index(outsam_data))
    
    
    if (is.null(outsam_perf)){
      outsam_perf = out_perf
      mv_outsam_perf = mv_out_perf
      avg_outsam_perf = avg_out_perf
      ms_outsam_perf = ms_out_perf
      rp_outsam_perf=rp_out_perf
    }else{
      outsam_perf = rbind(outsam_perf, out_perf*as.numeric(outsam_perf[length(outsam_perf)]))
      mv_outsam_perf = rbind(mv_outsam_perf, mv_out_perf*as.numeric(mv_outsam_perf[length(mv_outsam_perf)]))
      avg_outsam_perf = rbind(avg_outsam_perf, avg_out_perf*as.numeric(avg_outsam_perf[length(avg_outsam_perf)]))
      ms_outsam_perf = rbind(ms_outsam_perf, ms_out_perf*as.numeric(ms_outsam_perf[length(ms_outsam_perf)]))
      rp_outsam_perf = rbind(rp_outsam_perf, rp_out_perf*as.numeric(rp_outsam_perf[length(rp_outsam_perf)]))
    }
    cat("=")
  }
  cat("\nBack test finished!\n")
  
  # time series of weights, RC, sharpe
  weight_ts<- xts(weight_ts,order.by = index(all_data)[outsam_index_list[,1]])
  colnames(weight_ts)<- colnames(all_data)
  
  insamRC_ts<- xts(insamRC_ts,order.by = index(all_data)[outsam_index_list[,1]])
  colnames(insamRC_ts)<- colnames(all_data)
  
  outsamRC_ts<- xts(outsamRC_ts,order.by = index(all_data)[outsam_index_list[,1]])
  colnames(outsamRC_ts)<- colnames(all_data)
  
  #time series of outsam result
  out_result<- cbind(outsam_perf,avg_outsam_perf,mv_outsam_perf,ms_outsam_perf, rp_outsam_perf)
  colnames(out_result)<-c("optimal portfolio","1/n","min-var portfolio","max-sharpe portfolio","risk-parity")
  
  return(list("weight_ts"=weight_ts,
              "insamRC_ts"=insamRC_ts,
              "outsamRC_ts"=outsamRC_ts,
              "out_result"=out_result))
}

#show back test result: plot&table
show_result<- function(res, rf_ts){
  basic<- paste("(",name,"/insam=",insam_length,"w outsam=",outsam_length,"w/weight=[",w_lb[1],",",w_ub[1],"] risk weight=[",rw_lb[1],",",rw_ub[1],"])",sep="")
  #weights
  w_data<- cbind(index(res$weight_ts),as.data.frame(res$weight_ts))
  colnames(w_data)[1]<-"date"
  grid.arrange(ggplot(melt(w_data,id="date"),aes(date,value,colour=variable))+geom_line()+ggtitle(paste("time series of optimized weights",basic)),
               ggplot(melt(w_data,id="date"),aes(date,value))+geom_area(aes(fill=variable)))
  
  #risk contribution
  d1<- cbind(index(res$insamRC_ts),as.data.frame(res$insamRC_ts),rep("insample",nrow(res$insamRC_ts)))
  colnames(d1)[c(1,ncol(d1))]<-c("date","class")
  d2<- cbind(index(res$outsamRC_ts),as.data.frame(res$outsamRC_ts),rep("outsample",nrow(res$outsamRC_ts)))
  colnames(d2)[c(1,ncol(d2))]<-c("date","class")
  p<-ggplot(melt(rbind(d1,d2),id=c("date","class")), aes(date,value,colour=variable)) + geom_line() + 
    facet_wrap(~ class,scales = "free",ncol = 1)+ggtitle(paste("time series of risk contributions",basic))
  print(p)
  p<-ggplot(melt(rbind(d1,d2),id=c("date","class")), aes(date,value)) + geom_area(aes(fill=variable)) + 
    facet_wrap(~ class,scales = "free",ncol = 1)+ggtitle(paste("time series of risk contributions (filled)",basic))
  print(p)
  
  #out sample performance
  p<- autoplot(res$out_result,facets = F)+ggtitle(paste("time series of performance",basic))
  print(p)

  print("outsample perf analysis:")
  print(get_perf(res$out_result,rf_ts))
}


start_date = index(all_data)[1]
end_date = "2017/3/19"
insam_length = 36
outsam_length = 12
# weight bounds
w_lb<- rep(0,ncol(all_data))
w_ub<- rep(1,ncol(all_data))
# risk weight bounds
rw_lb<- rep(1/ncol(all_data),ncol(all_data))-0.05
rw_ub<- rep(1/ncol(all_data),ncol(all_data))+0.05

res<- rb_back_test(all_data,start_date,end_date,insam_length,outsam_length,w_lb,w_ub,rw_lb,rw_ub,rf_ts,NULL)
show_result(res)

