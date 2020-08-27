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

prc<-abind(prc1,prc2,along=3)#REF
dim(prc)<-c(96,96,12,151)

#########################
nc<-nc_open(CAN2[1])
v1<-nc$var[[2]]
lhb1<-ncvar_get(nc,v1)/2.499*0.0864
lhb1[lhb1>1e10]<-NA

nc<-nc_open(CAN2[3])
v1<-nc$var[[2]]
lhb2<-ncvar_get(nc,v1)/2.499*0.0864
lhb2[lhb2>1e10]<-NA

lhb<-abind(lhb1,lhb2,along=3)
dim(lhb)<-c(96,96,12,151)

nc<-nc_open(CAN2[2])
v1<-nc$var[[3]]
prb1<-ncvar_get(nc,v1)*86400

nc<-nc_open(CAN2[4])
v1<-nc$var[[3]]
prb2<-ncvar_get(nc,v1)*86400

prb<-abind(prb1,prb2,along=3)#expB
dim(prb)<-c(96,96,12,151)

#########################
nc<-nc_open(CAN1[1])
v1<-nc$var[[2]]
lha1<-ncvar_get(nc,v1)/2.499*0.0864
lha1[lha1>1e10]<-NA

nc<-nc_open(CAN1[3])
v1<-nc$var[[2]]
lha2<-ncvar_get(nc,v1)/2.499*0.0864
lha2[lha2>1e10]<-NA

lha<-abind(lha1,lha2,along=3)
dim(lha)<-c(96,96,12,151)

nc<-nc_open(CAN1[2])
v1<-nc$var[[3]]
pra1<-ncvar_get(nc,v1)*86400

nc<-nc_open(CAN1[4])
v1<-nc$var[[3]]
pra2<-ncvar_get(nc,v1)*86400

pra<-abind(pra1,pra2,along=3)#expA
dim(pra)<-c(96,96,12,151)

######################
nc<-nc_open(CAN[1])
v1<-nc$var[[2]]
msc<-ncvar_get(nc,v1)
dim(msc)<-c(96,96,12,151)

nc<-nc_open(CAN[2])
v1<-nc$var[[2]]
msb<-ncvar_get(nc,v1)
dim(msb)<-c(96,96,12,151)

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
msb60<-msb[,,,c(22:51,122:151)]
msb_mean<-apply(msb60,c(1,2,4),mean,na.rm=T)
    
mat_msb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_msb[,,m]<-reproject1d(msb_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_msb.RData")
save(mat_msb,file=outfile1)

###############################
msc60<-msc[,,,c(22:51,122:151)]
msc_mean<-apply(msc60,c(1,2,4),mean,na.rm=T)
    
mat_msc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_msc[,,m]<-reproject1d(msc_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_msc.RData")
save(mat_msc,file=outfile2)

###############################
pra60<-pra[,,,c(22:51,122:151)]
pra_mean<-apply(pra60,c(1,2,4),mean,na.rm=T)
    
mat_pra<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_pra[,,m]<-reproject1d(pra_mean[,,m],Rmean)
}
    
outfile1<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_pra.RData")
save(mat_pra,file=outfile1)

###############################
prb60<-prb[,,,c(22:51,122:151)]
prb_mean<-apply(prb60,c(1,2,4),mean,na.rm=T)
    
mat_prb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_prb[,,m]<-reproject1d(prb_mean[,,m],Rmean)
}
    
outfile2<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_prb.RData")
save(mat_prb,file=outfile2)

###############################
prc60<-prc[,,,c(22:51,122:151)]
prc_mean<-apply(prc60,c(1,2,4),mean,na.rm=T)
    
mat_prc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_prc[,,m]<-reproject1d(prc_mean[,,m],Rmean)
}
    
outfile3<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_prc.RData")
save(mat_prc,file=outfile3)

###############################
lha60<-lha[,,,c(22:51,122:151)]
lha_mean<-apply(lha60,c(1,2,4),mean,na.rm=T)
    
mat_lha<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_lha[,,m]<-reproject1d(lha_mean[,,m],Rmean)
}
    
outfile4<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_lha.RData")
save(mat_lha,file=outfile4)

###############################
lhb60<-lhb[,,,c(22:51,122:151)]
lhb_mean<-apply(lhb60,c(1,2,4),mean,na.rm=T)
    
mat_lhb<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_lhb[,,m]<-reproject1d(lhb_mean[,,m],Rmean)
}
    
outfile5<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_lhb.RData")
save(mat_lhb,file=outfile5)

###############################
lhc60<-lhc[,,,c(22:51,122:151)]
lhc_mean<-apply(lhc60,c(1,2,4),mean,na.rm=T)
    
mat_lhc<-array(NA,dim=c(360,180,60))
for (m in 1:60){
    mat_lhc[,,m]<-reproject1d(lhc_mean[,,m],Rmean)
}
    
outfile6<-paste0("/rigel/glab/users/sz2766/GLACE_paper_four/out/month/IPSL_lhc.RData")
save(mat_lhc,file=outfile6)

