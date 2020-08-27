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

nc<-nc_open(CAN[22])
v1<-nc$var[[1]]
prc<-ncvar_get(nc,v1)*86400
dim(prc)<-c(144,90,12,150)

nc<-nc_open(CAN[23])
v1<-nc$var[[1]]
pra<-ncvar_get(nc,v1)*86400
dim(pra)<-c(144,90,12,150)

nc<-nc_open(CAN[19])
v1<-nc$var[[3]]
msc<-ncvar_get(nc,v1)
msc[msc<0]<-NA
dim(msc)<-c(144,90,12,150)

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
###############################
msc1<-msc[,,c(3:12,1:2),c(21:50)]
msc2<-msc[,,c(3:12,1:2),c(121:150)]

dim(msc1)<-c(144,90,3,4,30)
dim(msc2)<-c(144,90,3,4,30)

msc1_mean<-apply(msc1,c(1,2,4),mean,na.rm=T)
msc2_mean<-apply(msc2,c(1,2,4),mean,na.rm=T)
    
mat_msc<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_msc[,,m]<-reproject1d(msc1_mean[,,m],Rmean)
    mat_msc[,,m+4]<-reproject1d(msc2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/GFDL_msc_season.RData")
save(mat_msc,file=outfile2)

###############################
pra1<-pra[,,c(3:12,1:2),c(21:50)]
pra2<-pra[,,c(3:12,1:2),c(121:150)]

dim(pra1)<-c(144,90,3,4,30)
dim(pra2)<-c(144,90,3,4,30)

pra1_mean<-apply(pra1,c(1,2,4),mean,na.rm=T)
pra2_mean<-apply(pra2,c(1,2,4),mean,na.rm=T)
    
mat_pra<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_pra[,,m]<-reproject1d(pra1_mean[,,m],Rmean)
    mat_pra[,,m+4]<-reproject1d(pra2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/GFDL_pra_season.RData")
save(mat_pra,file=outfile2)

###############################
prc1<-prc[,,c(3:12,1:2),c(21:50)]
prc2<-prc[,,c(3:12,1:2),c(121:150)]

dim(prc1)<-c(144,90,3,4,30)
dim(prc2)<-c(144,90,3,4,30)

prc1_mean<-apply(prc1,c(1,2,4),mean,na.rm=T)
prc2_mean<-apply(prc2,c(1,2,4),mean,na.rm=T)
    
mat_prc<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_prc[,,m]<-reproject1d(prc1_mean[,,m],Rmean)
    mat_prc[,,m+4]<-reproject1d(prc2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/GFDL_prc_season.RData")
save(mat_prc,file=outfile2)


###############################
lha1<-lha[,,c(3:12,1:2),c(21:50)]
lha2<-lha[,,c(3:12,1:2),c(121:150)]

dim(lha1)<-c(144,90,3,4,30)
dim(lha2)<-c(144,90,3,4,30)

lha1_mean<-apply(lha1,c(1,2,4),mean,na.rm=T)
lha2_mean<-apply(lha2,c(1,2,4),mean,na.rm=T)
    
mat_lha<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_lha[,,m]<-reproject1d(lha1_mean[,,m],Rmean)
    mat_lha[,,m+4]<-reproject1d(lha2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/GFDL_lha_season.RData")
save(mat_lha,file=outfile2)

###############################
lhc1<-lhc[,,c(3:12,1:2),c(21:50)]
lhc2<-lhc[,,c(3:12,1:2),c(121:150)]

dim(lhc1)<-c(144,90,3,4,30)
dim(lhc2)<-c(144,90,3,4,30)

lhc1_mean<-apply(lhc1,c(1,2,4),mean,na.rm=T)
lhc2_mean<-apply(lhc2,c(1,2,4),mean,na.rm=T)
    
mat_lhc<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_lhc[,,m]<-reproject1d(lhc1_mean[,,m],Rmean)
    mat_lhc[,,m+4]<-reproject1d(lhc2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/GFDL_lhc_season.RData")
save(mat_lhc,file=outfile2)

