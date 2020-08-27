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
CAN1<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/expA/',pattern='.nc',full.names=TRUE)
CAN2<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/expB/',pattern='.nc',full.names=TRUE)
CAN3<-list.files('/rigel/glab/users/sz2766/GLACE_paper_four/data/IPSL/CTL/',pattern='.nc',full.names=TRUE)
###################01#####################
nc<-nc_open(CAN3[1])
v1<-nc$var[[2]]
lhc1<-ncvar_get(nc,v1)/2.499*0.0864
lhc1[lhc1>1e10]<-NA

nc<-nc_open(CAN3[3])
v1<-nc$var[[2]]
lhc2<-ncvar_get(nc,v1)/2.499*0.0864
lhc2[lhc2>1e10]<-NA

lhc<-abind(lhc1,lhc2,along=3)
dim(lhc)<-c(96,96,12,151)

nc<-nc_open(CAN3[2])
v1<-nc$var[[3]]
prc1<-ncvar_get(nc,v1)*86400

nc<-nc_open(CAN3[4])
v1<-nc$var[[3]]
prc2<-ncvar_get(nc,v1)*86400

prc<-abind(prc1,prc2,along=3)
dim(prc)<-c(96,96,12,151)

mean(lhc1,na.rm=T)
mean(lhc2,na.rm=T)
mean(prc1,na.rm=T)
mean(prc2,na.rm=T)

#########################
nc<-nc_open(CAN2[1])
v1<-nc$var[[2]]
lhb1<-ncvar_get(nc,v1)/2.499*0.0864
lhb1[lhb1>1e10]<-NA

nc<-nc_open(CAN2[8])
v1<-nc$var[[2]]
lhb2<-ncvar_get(nc,v1)/2.499*0.0864
lhb2[lhb2>1e10]<-NA

lhb<-abind(lhb1,lhb2,along=3)
dim(lhb)<-c(96,96,12,151)

nc<-nc_open(CAN2[3])
v1<-nc$var[[3]]
prb1<-ncvar_get(nc,v1)*86400

nc<-nc_open(CAN2[10])
v1<-nc$var[[3]]
prb2<-ncvar_get(nc,v1)*86400

prb<-abind(prb1,prb2,along=3)
dim(prb)<-c(96,96,12,151)

mean(lhb1,na.rm=T)
mean(lhb2,na.rm=T)
mean(prb1,na.rm=T)
mean(prb2,na.rm=T)

#########################
nc<-nc_open(CAN1[1])
v1<-nc$var[[2]]
lha1<-ncvar_get(nc,v1)/2.499*0.0864
lha1[lha1>1e10]<-NA

nc<-nc_open(CAN1[8])
v1<-nc$var[[2]]
lha2<-ncvar_get(nc,v1)/2.499*0.0864
lha2[lha2>1e10]<-NA

lha<-abind(lha1,lha2,along=3)
dim(lha)<-c(96,96,12,151)

nc<-nc_open(CAN1[3])
v1<-nc$var[[3]]
pra1<-ncvar_get(nc,v1)*86400

nc<-nc_open(CAN1[10])
v1<-nc$var[[3]]
pra2<-ncvar_get(nc,v1)*86400

pra<-abind(pra1,pra2,along=3)
dim(pra)<-c(96,96,12,151)

######################
nc<-nc_open(CAN[1])
v1<-nc$var[[2]]
msc<-ncvar_get(nc,v1)
dim(msc)<-c(96,96,12,151)

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
msc1<-msc[,,c(3:12,1:2),c(22:51)]
msc2<-msc[,,c(3:12,1:2),c(122:151)]

dim(msc1)<-c(96,96,3,4,30)
dim(msc2)<-c(96,96,3,4,30)

msc1_mean<-apply(msc1,c(1,2,4),mean,na.rm=T)
msc2_mean<-apply(msc2,c(1,2,4),mean,na.rm=T)
    
mat_msc<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_msc[,,m]<-reproject1d(msc1_mean[,,m],Rmean)
    mat_msc[,,m+4]<-reproject1d(msc2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/IPSL_msc_season.RData")
save(mat_msc,file=outfile2)

###############################
pra1<-pra[,,c(3:12,1:2),c(22:51)]
pra2<-pra[,,c(3:12,1:2),c(122:151)]

dim(pra1)<-c(96,96,3,4,30)
dim(pra2)<-c(96,96,3,4,30)

pra1_mean<-apply(pra1,c(1,2,4),mean,na.rm=T)
pra2_mean<-apply(pra2,c(1,2,4),mean,na.rm=T)
    
mat_pra<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_pra[,,m]<-reproject1d(pra1_mean[,,m],Rmean)
    mat_pra[,,m+4]<-reproject1d(pra2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/IPSL_pra_season.RData")
save(mat_pra,file=outfile2)

###############################
prb1<-prb[,,c(3:12,1:2),c(22:51)]
prb2<-prb[,,c(3:12,1:2),c(122:151)]

dim(prb1)<-c(96,96,3,4,30)
dim(prb2)<-c(96,96,3,4,30)

prb1_mean<-apply(prb1,c(1,2,4),mean,na.rm=T)
prb2_mean<-apply(prb2,c(1,2,4),mean,na.rm=T)
    
mat_prb<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_prb[,,m]<-reproject1d(prb1_mean[,,m],Rmean)
    mat_prb[,,m+4]<-reproject1d(prb2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/IPSL_prb_season.RData")
save(mat_prb,file=outfile2)

###############################
prc1<-prc[,,c(3:12,1:2),c(22:51)]
prc2<-prc[,,c(3:12,1:2),c(122:151)]

dim(prc1)<-c(96,96,3,4,30)
dim(prc2)<-c(96,96,3,4,30)

prc1_mean<-apply(prc1,c(1,2,4),mean,na.rm=T)
prc2_mean<-apply(prc2,c(1,2,4),mean,na.rm=T)
    
mat_prc<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_prc[,,m]<-reproject1d(prc1_mean[,,m],Rmean)
    mat_prc[,,m+4]<-reproject1d(prc2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/IPSL_prc_season.RData")
save(mat_prc,file=outfile2)


###############################
lha1<-lha[,,c(3:12,1:2),c(22:51)]
lha2<-lha[,,c(3:12,1:2),c(122:151)]

dim(lha1)<-c(96,96,3,4,30)
dim(lha2)<-c(96,96,3,4,30)

lha1_mean<-apply(lha1,c(1,2,4),mean,na.rm=T)
lha2_mean<-apply(lha2,c(1,2,4),mean,na.rm=T)
    
mat_lha<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_lha[,,m]<-reproject1d(lha1_mean[,,m],Rmean)
    mat_lha[,,m+4]<-reproject1d(lha2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/IPSL_lha_season.RData")
save(mat_lha,file=outfile2)

###############################
lhb1<-lhb[,,c(3:12,1:2),c(22:51)]
lhb2<-lhb[,,c(3:12,1:2),c(122:151)]

dim(lhb1)<-c(96,96,3,4,30)
dim(lhb2)<-c(96,96,3,4,30)

lhb1_mean<-apply(lhb1,c(1,2,4),mean,na.rm=T)
lhb2_mean<-apply(lhb2,c(1,2,4),mean,na.rm=T)
    
mat_lhb<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_lhb[,,m]<-reproject1d(lhb1_mean[,,m],Rmean)
    mat_lhb[,,m+4]<-reproject1d(lhb2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/IPSL_lhb_season.RData")
save(mat_lhb,file=outfile2)

###############################
lhc1<-lhc[,,c(3:12,1:2),c(22:51)]
lhc2<-lhc[,,c(3:12,1:2),c(122:151)]

dim(lhc1)<-c(96,96,3,4,30)
dim(lhc2)<-c(96,96,3,4,30)

lhc1_mean<-apply(lhc1,c(1,2,4),mean,na.rm=T)
lhc2_mean<-apply(lhc2,c(1,2,4),mean,na.rm=T)
    
mat_lhc<-array(NA,dim=c(360,180,8))
for (m in 1:4){
    mat_lhc[,,m]<-reproject1d(lhc1_mean[,,m],Rmean)
    mat_lhc[,,m+4]<-reproject1d(lhc2_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/season/IPSL_lhc_season.RData")
save(mat_lhc,file=outfile2)

