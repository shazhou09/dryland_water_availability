library(ncdf4)
library(raster)
library(parallel)
library(abind)
library(lubridate)

##########################################
CAN<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/GFDL/',pattern='.nc',full.names=TRUE)
###################01#####################
nc<-nc_open(CAN[7])
v1<-nc$var[[4]]
lhc<-ncvar_get(nc,v1)/2.499*0.0864
dim(lhc)<-c(144,90,12,150)

nc<-nc_open(CAN[8])
v1<-nc$var[[4]]
lha<-ncvar_get(nc,v1)/2.499*0.0864
dim(lha)<-c(144,90,12,150)

nc<-nc_open(CAN[9])
v1<-nc$var[[4]]
lhb<-ncvar_get(nc,v1)/2.499*0.0864
dim(lhb)<-c(144,90,12,150)

nc<-nc_open(CAN[22])
v1<-nc$var[[1]]
prc<-ncvar_get(nc,v1)*86400
dim(prc)<-c(144,90,12,150)

nc<-nc_open(CAN[23])
v1<-nc$var[[1]]
pra<-ncvar_get(nc,v1)*86400
dim(pra)<-c(144,90,12,150)

nc<-nc_open(CAN[24])
v1<-nc$var[[1]]
prb<-ncvar_get(nc,v1)*86400
dim(prb)<-c(144,90,12,150)

nc<-nc_open(CAN[19])
v1<-nc$var[[3]]
msc<-ncvar_get(nc,v1)
msc[msc<0]<-NA
dim(msc)<-c(144,90,12,150)

nc<-nc_open(CAN[21])
v1<-nc$var[[3]]
msb<-ncvar_get(nc,v1)
msb[msb<0]<-NA
dim(msb)<-c(144,90,12,150)

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
msb60<-msb[,,,c(21:50,121:150)]
msb_mean<-apply(msb60,c(1,2,4),mean,na.rm=T)
    
mat_msb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_msb[,,m]<-reproject1d(msb_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_msb.RData")
save(mat_msb,file=outfile1)

###############################
msc60<-msc[,,,c(21:50,121:150)]
msc_mean<-apply(msc60,c(1,2,4),mean,na.rm=T)
    
mat_msc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_msc[,,m]<-reproject1d(msc_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_msc.RData")
save(mat_msc,file=outfile2)

###############################
pra60<-pra[,,,c(21:50,121:150)]
pra_mean<-apply(pra60,c(1,2,4),mean,na.rm=T)
    
mat_pra<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_pra[,,m]<-reproject1d(pra_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_pra.RData")
save(mat_pra,file=outfile1)

###############################
prb0<-prb[,,,c(21:50,121:150)]
prb_mean<-apply(prb60,c(1,2,4),mean,na.rm=T)
    
mat_prb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_prb[,,m]<-reproject1d(prb_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_prb.RData")
save(mat_prb,file=outfile2)

###############################
prc60<-prc[,,,c(21:50,121:150)]
prc_mean<-apply(prc60,c(1,2,4),mean,na.rm=T)
    
mat_prc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_prc[,,m]<-reproject1d(prc_mean[,,m],Rmean)
}
    
outfile3<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_prc.RData")
save(mat_prc,file=outfile3)

###############################
lha60<-lha[,,,c(21:50,121:150)]
lha_mean<-apply(lha60,c(1,2,4),mean,na.rm=T)
    
mat_lha<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_lha[,,m]<-reproject1d(lha_mean[,,m],Rmean)
}
    
outfile4<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_lha.RData")
save(mat_lha,file=outfile4)

###############################
lhb0<-lhb[,,,c(21:50,121:150)]
lhb_mean<-apply(lhb60,c(1,2,4),mean,na.rm=T)
    
mat_lhb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_lhb[,,m]<-reproject1d(lhb_mean[,,m],Rmean)
}
    
outfile5<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_lhb.RData")
save(mat_lhb,file=outfile5)

###############################
lhc60<-lhc[,,,c(21:50,121:150)]
lhc_mean<-apply(lhc60,c(1,2,4),mean,na.rm=T)
    
mat_lhc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_lhc[,,m]<-reproject1d(lhc_mean[,,m],Rmean)
}
    
outfile6<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/GFDL_lhc.RData")
save(mat_lhc,file=outfile6)
