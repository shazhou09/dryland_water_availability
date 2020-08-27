library(ncdf4)
library(raster)
library(parallel)
library(abind)
library(lubridate)

##########################################
CAN<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/ECHAM6/',pattern='.nc',full.names=TRUE)
###################01#####################
ts<-function(x){
  dim(x)<-c(192,96)
  y<-apply(x,1,rev)
  z<-t(y)
  dim(z)<-18432
  return(z)
}

nc<-nc_open(CAN[12])
v1<-nc$var[[1]]
ssc<-ncvar_get(nc,v1)
ssc[ssc< 0]=NA
ssc<-apply(ssc,3,ts)
dim(ssc)<-c(192,96,12,151)#REF

nc<-nc_open(CAN[10])
v1<-nc$var[[1]]
ssa<-ncvar_get(nc,v1)
ssa[ssa< 0]=NA
dim(ssa)<-c(192,96,12,151)#expA

nc<-nc_open(CAN[11])
v1<-nc$var[[1]]
ssb<-ncvar_get(nc,v1)
ssb[ssb< 0]=NA
dim(ssb)<-c(192,96,12,151)#expB

nc<-nc_open(CAN[30])
v1<-nc$var[[1]]
tsc<-ncvar_get(nc,v1)
tsc[tsc< 0]=NA
tsc<-apply(tsc,3,ts)
dim(tsc)<-c(192,96,12,151)#REF

nc<-nc_open(CAN[28])
v1<-nc$var[[1]]
tsa<-ncvar_get(nc,v1)
tsa[tsa< 0]=NA
dim(tsa)<-c(192,96,12,151)#expA

nc<-nc_open(CAN[29])
v1<-nc$var[[1]]
tsb<-ncvar_get(nc,v1)
tsb[tsb< 0]=NA
dim(tsb)<-c(192,96,12,151)#expB

################################
projlat<-CRS("+proj=longlat +datum=WGS84 +no_defs")

Rmean<-raster(nrow=180,ncol=360)
extent(Rmean)<-c(-180,180,-90,90)
Rmean[is.na(Rmean)]<-1
projection(Rmean)<-projlat

reproject1d<-function(mat,Rmean){
  pr_mask<-raster(t(mat))
  extent(pr_mask)<-c(-180,180,-90,90)
  projection(pr_mask)<-projlat
    
  meanpr_proj<-projectRaster(from = pr_mask,to=Rmean)
  mat_pr<-t(as.matrix(meanpr_proj))[,180:1]
  return(mat_pr)
}

###############################
ssa60<-ssa[,,,c(22:51,122:151)]
ssa_mean<-apply(ssa60,c(1,2,4),mean,na.rm=T)
    
mat_ssa<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_ssa[,,m]<-reproject1d(ssa_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/ECHAM6_ssa.RData")
save(mat_ssa,file=outfile1)

###############################
ssc60<-ssc[,,,c(22:51,122:151)]
ssc_mean<-apply(ssc60,c(1,2,4),mean,na.rm=T)
    
mat_ssc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_ssc[,,m]<-reproject1d(ssc_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/ECHAM6_ssc.RData")
save(mat_ssc,file=outfile2)

###############################
ssb60<-ssb[,,,c(22:51,122:151)]
ssb_mean<-apply(ssb60,c(1,2,4),mean,na.rm=T)
    
mat_ssb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_ssb[,,m]<-reproject1d(ssb_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/ECHAM6_ssb.RData")
save(mat_ssb,file=outfile2)

###############################
tsa60<-tsa[,,,c(22:51,122:151)]
tsa_mean<-apply(tsa60,c(1,2,4),mean,na.rm=T)
    
mat_tsa<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsa[,,m]<-reproject1d(tsa_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/ECHAM6_tsa.RData")
save(mat_tsa,file=outfile1)

###############################
tsc60<-tsc[,,,c(22:51,122:151)]
tsc_mean<-apply(tsc60,c(1,2,4),mean,na.rm=T)
    
mat_tsc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsc[,,m]<-reproject1d(tsc_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/ECHAM6_tsc.RData")
save(mat_tsc,file=outfile2)

###############################
tsb60<-tsb[,,,c(22:51,122:151)]
tsb_mean<-apply(tsb60,c(1,2,4),mean,na.rm=T)
    
mat_tsb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsb[,,m]<-reproject1d(tsb_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/ECHAM6_tsb.RData")
save(mat_tsb,file=outfile2)
