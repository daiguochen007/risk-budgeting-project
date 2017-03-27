##### test insample sharpe ratio time series with different weight range
##### start with asset_group_backtest
##### 1. load functions
##### 2. add data

############################################# parameters
start_date = index(all_data)[1]
end_date = "2017/3/19"
insam_length = 90
outsam_length = 30

# risk weight bounds
rw_lb<- rep(0,ncol(all_data))
rw_ub<- rep(0.5,ncol(all_data))

# range matrix
mat<- matrix(c(0,1,0,2,0,5,-2,2,-5,5),ncol=2,byrow = T)

# loop
ts <- NULL
for(i in 1:nrow(mat)){
  w_lb<- rep(mat[i,1],ncol(all_data))
  w_ub<- rep(mat[i,2],ncol(all_data))
  res<- rb_back_test(all_data,start_date,end_date,insam_length,outsam_length,w_lb,w_ub,rw_lb,rw_ub,rf_ts)
  ts<- cbind(ts,res$sharpe[,1])
  cat(paste(i,"th loop finished\n"))
}
colnames(ts)<- paste(mat[,1],"to",mat[,2])
autoplot(ts,facets = F,main = "insample sharpe ratio time series with different weight range")


