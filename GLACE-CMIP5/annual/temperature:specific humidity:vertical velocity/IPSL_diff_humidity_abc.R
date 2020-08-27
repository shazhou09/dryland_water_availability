library(ncdf4)
library(raster)
library(parallel)
library(abind)
library(lubridate)

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

##########################################
CAN<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/',pattern='.nc',full.names=TRUE)
###################01#####################
nc<-nc_open(CAN[4])
v1<-nc$var[[3]]
tsc<-ncvar_get(nc,v1)
tsc[tsc< 0]=NA
dim(tsc)<-c(96,96,12,151)#REF

nc<-nc_open(CAN[5])
v1<-nc$var[[2]]
tsa<-ncvar_get(nc,v1)
tsa[tsa< 0]=NA
dim(tsa)<-c(96,96,12,151)#expA

nc<-nc_open(CAN[6])
v1<-nc$var[[3]]
tsb<-ncvar_get(nc,v1)
tsb[tsb< 0]=NA
dim(tsb)<-c(96,96,12,151)#expB

###############################
tsa60<-tsa[,,,c(22:51,122:151)]
tsa_mean<-apply(tsa60,c(1,2,4),mean,na.rm=T)
    
mat_tsa<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsa[,,m]<-reproject1d(tsa_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_tsa.RData")
save(mat_tsa,file=outfile1)

###############################
tsc60<-tsc[,,,c(22:51,122:151)]
tsc_mean<-apply(tsc60,c(1,2,4),mean,na.rm=T)
    
mat_tsc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsc[,,m]<-reproject1d(tsc_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_tsc.RData")
save(mat_tsc,file=outfile2)

###############################
tsb60<-tsb[,,,c(22:51,122:151)]
tsb_mean<-apply(tsb60,c(1,2,4),mean,na.rm=T)
    
mat_tsb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsb[,,m]<-reproject1d(tsb_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_tsb.RData")
save(mat_tsb,file=outfile2)


#########################
CAN1<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/expA/',pattern='.nc',full.names=TRUE)
CAN2<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/expB/',pattern='.nc',full.names=TRUE)
#########################
nc<-nc_open(CAN2[4])
v1<-nc$var[[3]]
ssb1<-ncvar_get(nc,v1)
ssb1[ssb1<0]<-NA

nc<-nc_open(CAN2[10])
v1<-nc$var[[3]]
ssb2<-ncvar_get(nc,v1)
ssb2[ssb2<0]<-NA

ssb<-abind(ssb1,ssb2,along=3)
dim(ssb)<-c(96,96,12,151)

nc<-nc_open(CAN1[4])
v1<-nc$var[[3]]
ssa1<-ncvar_get(nc,v1)
ssa1[ssa1<0]<-NA

nc<-nc_open(CAN1[10])
v1<-nc$var[[3]]
ssa2<-ncvar_get(nc,v1)
ssa2[ssa2<0]<-NA

ssa<-abind(ssa1,ssa2,along=3)
dim(ssa)<-c(96,96,12,151)

###############################
ssa60<-ssa[,,,c(22:51,122:151)]
ssa_mean<-apply(ssa60,c(1,2,4),mean,na.rm=T)
    
mat_ssa<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_ssa[,,m]<-reproject1d(ssa_mean[,,m],Rmean)
}
    
outfile3<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_ssa.RData")
save(mat_ssa,file=outfile3)

###############################
ssb60<-ssb[,,,c(22:51,122:151)]
ssb_mean<-apply(ssb60,c(1,2,4),mean,na.rm=T)
    
mat_ssb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_ssb[,,m]<-reproject1d(ssb_mean[,,m],Rmean)
}
    
outfile3<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_ssb.RData")
save(mat_ssb,file=outfile3)

