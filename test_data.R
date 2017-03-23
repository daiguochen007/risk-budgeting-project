library(quantmod)
library(ggplot2)
library(ggfortify)
#sp500 sectors
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
autoplot(all_data,facets = FALSE)

#start from 2001
all_data<- na.omit(all_data)
for (i in 1:ncol(all_data)){
  all_data[,i]<- all_data[,i]/as.numeric(all_data[1,i])
}
autoplot(all_data,facets = FALSE)

#all assets
assets<- c("spx","dow","crb","us_agg","reit","nasdaq","euro_50","ftse100","cac40","dax","nikkei","hsi","csi300","euro_agg","asianpacific_agg","global_highyield")
all_data<-NULL
for (i in 1:16){
  data<- read.csv(paste("/Users/guochendai/Desktop/4th semester/7043 capstone project/data/",
                        assets[i],".csv",sep=""))
  data<- xts(data[,2:3],order.by=as.Date(data[,1], format="%m/%d/%y"))
  all_data<- cbind(all_data, data[,2]/as.numeric(data[1,2]))
}
colnames(all_data)<- assets
autoplot(all_data[which(index(all_data)<"2990/01/01"),],facets = FALSE)


########### compare spx and equally weighted sector index
avg<- xts(apply(all_data,1,mean),order.by = index(all_data))
spx<- read.csv("/Users/guochendai/Desktop/4th semester/7043 capstone project/data/spx.csv")
spx<- xts(spx[,3],order.by=as.Date(spx[,1], format="%m/%d/%y"))
df<- na.omit(cbind(spx,avg))
df[,1]<- df[,1]/as.numeric(df[1,1])
colnames(df)<- c("spx","synthesis spx")
autoplot(df,facets = F)
  
  
  