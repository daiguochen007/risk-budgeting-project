library(quantmod)
library(ggplot2)
library(ggfortify)

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
for (i in 1:11){
  all_data[,i]<- all_data[,i]/as.numeric(all_data[1,i])
}
autoplot(all_data,facets = FALSE)
