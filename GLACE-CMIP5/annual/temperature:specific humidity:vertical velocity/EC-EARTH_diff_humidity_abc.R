library(ncdf4)
library(raster)
library(parallel)
library(abind)
library(lubridate)

##########################################
CAN<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/EC-EARTH/',pattern='.nc',full.names=TRUE)
###################01#####################
nc<-nc_open(CAN[18])
v1<-nc$var[[2]]
tsc<-ncvar_get(nc,v1)
tsc[tsc< 0]=NA
dim(tsc)<-c(320,160,12,151)#REF

nc<-nc_open(CAN[16])
v1<-nc$var[[1]]
tsa<-ncvar_get(nc,v1)
tsa[tsa< 0]=NA
dim(tsa)<-c(320,160,12,151)#expA

nc<-nc_open(CAN[17])
v1<-nc$var[[1]]
tsb<-ncvar_get(nc,v1)
tsb[tsb< 0]=NA
dim(tsb)<-c(320,160,12,151)#expB

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
tsa60<-tsa[,,,c(22:51,122:151)]
tsa_mean<-apply(tsa60,c(1,2,4),mean,na.rm=T)
    
mat_tsa<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsa[,,m]<-reproject1d(tsa_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/EC-EARTH_tsa.RData")
save(mat_tsa,file=outfile1)

###############################
tsc60<-tsc[,,,c(22:51,122:151)]
tsc_mean<-apply(tsc60,c(1,2,4),mean,na.rm=T)
    
mat_tsc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsc[,,m]<-reproject1d(tsc_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/EC-EARTH_tsc.RData")
save(mat_tsc,file=outfile2)

###############################
tsb60<-tsb[,,,c(22:51,122:151)]
tsb_mean<-apply(tsb60,c(1,2,4),mean,na.rm=T)
    
mat_tsb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_tsb[,,m]<-reproject1d(tsb_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/EC-EARTH_tsb.RData")
save(mat_tsb,file=outfile2)
