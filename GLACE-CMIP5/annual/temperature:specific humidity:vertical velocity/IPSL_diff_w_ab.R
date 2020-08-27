#########################
#########################
#########################
#########################
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

#########################
CAN1<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/expA/',pattern='.nc',full.names=TRUE)
CAN2<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/expB/',pattern='.nc',full.names=TRUE)

#########################
nc<-nc_open(CAN2[6])
v1<-nc$var[[3]]
wb1<-ncvar_get(nc,v1)

nc<-nc_open(CAN2[13])
v1<-nc$var[[3]]
wb2<-ncvar_get(nc,v1)

wb<-abind(wb1,wb2,along=4)
dim(wb)<-c(96,96,39,12,151)

nc<-nc_open(CAN1[6])
v1<-nc$var[[3]]
wa1<-ncvar_get(nc,v1)

nc<-nc_open(CAN1[13])
v1<-nc$var[[3]]
wa2<-ncvar_get(nc,v1)

wa<-abind(wa1,wa2,along=4)
dim(wa)<-c(96,96,39,12,151)

###############################
###############################
mat_wa<-array(NA,dim=c(360,180,117))

wa30<-apply(wa[,,,,22:51],c(1,2,3),mean,na.rm=T)
for(i in 1:39){
    mat_wa[,,i]<-reproject1d(wa30[,,i],Rmean)
}

wa30<-apply(wa[,,,,122:151],c(1,2,3),mean,na.rm=T)
for(i in 1:39){
    mat_wa[,,i+39]<-reproject1d(wa30[,,i],Rmean)
}

mat_wa[,,79:117]<-mat_wa[,,40:78]-mat_wa[,,1:39]
outfile3<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_wa_pressure.RData")
save(mat_wa,file=outfile3)

###############################
###############################
mat_wb<-array(NA,dim=c(360,180,117))

wb30<-apply(wb[,,,,22:51],c(1,2,3),mean,na.rm=T)
for(i in 1:39){
    mat_wb[,,i]<-reproject1d(wb30[,,i],Rmean)
}

wb30<-apply(wb[,,,,122:151],c(1,2,3),mean,na.rm=T)
for(i in 1:39){
    mat_wb[,,i+39]<-reproject1d(wb30[,,i],Rmean)
}

mat_wb[,,79:117]<-mat_wb[,,40:78]-mat_wb[,,1:39]
outfile3<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_wb_pressure.RData")
save(mat_wb,file=outfile3)

