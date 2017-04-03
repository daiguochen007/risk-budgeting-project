library(pheatmap)
library(rgl)

mlist<- seq(60,300,30)
nlist<- seq(10,180,20)
mat1<- matrix(NA,ncol=length(mlist),nrow=length(nlist),dimnames=list(nlist,mlist))
mat2<- matrix(NA,ncol=length(mlist),nrow=length(nlist),dimnames=list(nlist,mlist))
mat3<- matrix(NA,ncol=length(mlist),nrow=length(nlist),dimnames=list(nlist,mlist))
mat4<- matrix(NA,ncol=length(mlist),nrow=length(nlist),dimnames=list(nlist,mlist)) 
  
for (m in 1:length(mlist)){
  for (n in 1:length(nlist)){
    if(nlist[n] <= mlist[m]){
      insam_length = mlist[m]
      outsam_length = nlist[n]
      res<- rb_back_test(all_data,start_date,end_date,insam_length,outsam_length,w_lb,w_ub,rw_lb,rw_ub,rf_ts,NULL)
      
      mat1[n,m]<- get_perf(res$out_result[,1])[3] #optimal 
      mat2[n,m]<- get_perf(res$out_result[,2])[3] #1/n
      mat3[n,m]<- get_perf(res$out_result[,3])[3] #min-var 
      mat4[n,m]<- get_perf(res$out_result[,4])[3] #max sharpe
    }
  }
}

bks<- NULL
pheatmap(mat1,cluster_rows = F,cluster_cols = F, breaks = bks,main="optimal sharpe")
pheatmap(mat2,cluster_rows = F,cluster_cols = F,breaks = bks,main="1/n sharpe")
pheatmap(mat3,cluster_rows = F,cluster_cols = F,breaks = bks,main="min-var sharpe")
pheatmap(mat4,cluster_rows = F,cluster_cols = F,breaks = bks,main="max sharpe sharpe")


persp3d(nlist,mlist,mat1,col="lightblue",main="optimal sharpe")
#surface3d(nlist,mlist,mat4,col="red")
persp3d(nlist,mlist,mat2,col="lightblue",main="1/n sharpe")
persp3d(nlist,mlist,mat3,col="lightblue",main="min-var sharpe")
persp3d(nlist,mlist,mat4,col="lightblue",main="max sharpe sharpe")

# evaluate the sharpe matrix
eva<- matrix(NA,ncol = 4,nrow=2)
colnames(eva)<- c("optimal sharpe","1/n sharpe","min-var sharpe","max sharpe sharpe")
rownames(eva)<- c("mean","std")

eva[1,1]<- mean(na.omit(as.vector(mat1)))
eva[1,2]<- mean(na.omit(as.vector(mat2)))
eva[1,3]<- mean(na.omit(as.vector(mat3)))
eva[1,4]<- mean(na.omit(as.vector(mat4)))
eva[2,1]<- sd(na.omit(as.vector(mat1)))
eva[2,2]<- sd(na.omit(as.vector(mat2)))
eva[2,3]<- sd(na.omit(as.vector(mat3)))
eva[2,4]<- sd(na.omit(as.vector(mat4)))




