library(ncdf4)
library(parallel)
library(plsRglm)
library(plsdof)
library(boot)

pls_reg<-function(x){
  if (sum(is.na(x))>400){
    return(rep(NA,each=4))
  }

  tt<-length(x)/2

  pr<-x[1:tt]
  sm<-x[(1+tt):(tt*2)]

  if (max(abs(sm),na.rm=T)<1e-5){
    return(rep(NA,each=4))
  }else{
  tryCatch({
    dim(pr)<-c(12,tt/12)
    pr_mm<-apply(pr,1,mean,na.rm=T)
    pr01<-rep(pr_mm,times=tt/12)
    pr1<-pr-pr01
    dim(pr1)<-tt

    dim(sm)<-c(12,tt/12)
    sm_mm<-apply(sm,1,mean,na.rm=T)
    sm01<-rep(sm_mm,times=tt/12)
    sm1<-sm-sm01
    dim(sm1)<-tt

    tm<-1:tt
    data<-cbind(pr1,tm)
    data<-na.omit(data)
    pr02<-pr1-tm*lm(data[,1]~data[,2])$coefficients[2]

    data<-cbind(sm1,tm)
    data<-na.omit(data)
    sm02<-sm1-tm*lm(data[,1]~data[,2])$coefficients[2]

    apr1<-pr02[1:(tt-1)]
    asm1<-sm02[1:(tt-1)]
    apr2<-pr02[2:tt]
    asm2<-sm02[2:tt]
  
    out<-array(NA,dim=4)
    data1<-cbind(asm2,apr1,asm1)
    data2<-cbind(apr2,asm1,apr1)
    data1<-na.omit(data1)
    data2<-na.omit(data2)

    reg1<-plsRglm(data1[,1],data1[,2:3],2)
    reg2<-plsRglm(data2[,1],data2[,2:3],2)

    reg1_boot<-bootpls(reg1, R=500, verbose=FALSE)
    reg2_boot<-bootpls(reg2, R=500, verbose=FALSE)

    conf1<-boot.ci(reg1_boot, conf = 0.95, type = "bca",index=2)$bca[4:5]
    conf2<-boot.ci(reg2_boot, conf = 0.95, type = "bca",index=2)$bca[4:5]

    out[1]<-reg1$Std.Coeffs[2]
    out[2]<-reg2$Std.Coeffs[2]
    if (conf1[2]*conf1[1]>0){
      out[3]<-1
    }
    if (conf2[2]*conf2[1]>0){
      out[4]<-1
    }

    return(out)
  },error=function(e){return(rep(NA,each=4))})
  }
}

args = (commandArgs(trailingOnly = F))
print(args)
myargs <-sub('-','',args[length(args)])
args = as.numeric(myargs)
print(args)

##########################################
CAN<-list.files('/rigel/glab/users/sz2766/paper_four/ERA5M/var/',pattern='.nc',full.names=TRUE)
###################01#####################
nc<-nc_open(CAN[4])
v1<-nc$var[[1]]
dvit<-ncvar_get(nc,v1)
dvit<-dvit[,args,]
dim(dvit)<-c(360,480)

v1<-nc$var[[2]]
dvid<-ncvar_get(nc,v1)
dvid<-dvid[,args,]
dim(dvid)<-c(360,480)

v1<-nc$var[[3]]
dvitd<-ncvar_get(nc,v1)
dvitd<-dvitd[,args,]
dim(dvitd)<-c(360,480)

nc<-nc_open(CAN[7])
v1<-nc$var[[1]]
sm1<-ncvar_get(nc,v1)
sm1<-sm1[,args,]
dim(sm1)<-c(360,480)

v1<-nc$var[[2]]
sm2<-ncvar_get(nc,v1)
sm2<-sm2[,args,]
dim(sm2)<-c(360,480)

v1<-nc$var[[3]]
sm3<-ncvar_get(nc,v1)
sm3<-sm3[,args,]
dim(sm3)<-c(360,480)

sm<-(sm1*7+sm2*21+sm3*72)/100

dvit<- -dvit
dvid<- -dvid
dvitd<- -dvitd

#################
input7<-cbind(dvit,sm)#SM-TH
input8<-cbind(dvid,sm)#SM-MCD
input9<-cbind(dvitd,sm)#SM-COV

cl<-makeCluster(getOption("cl.cores",5))
clusterEvalQ(cl,library(parallel))
clusterEvalQ(cl,library(ncdf4))
clusterEvalQ(cl,library(plsRglm))
clusterEvalQ(cl,library(plsdof))
clusterEvalQ(cl,library(boot))

out7<-parRapply(cl,input7,pls_reg)
out8<-parRapply(cl,input8,pls_reg)
out9<-parRapply(cl,input9,pls_reg)

dim(out7)<-c(4,360)
dim(out8)<-c(4,360)
dim(out9)<-c(4,360)

out01<-rbind(out7,out8,out9)

#lat<-ncvar_get(nc,varid = 'latitude')
long<-ncvar_get(nc,varid = 'longitude')
result1<-1:12

#dimlat<-ncdim_def('latitude','deg',lat)
dimlong<-ncdim_def('longitude','deg',long)
dimout1<-ncdim_def('result','',result1)

ncdem<-ncvar_def('sc','',list(dimout1,dimlong),-9999,longname="sensitivity coefficient",prec='double')
ncout<-nc_create(paste('/rigel/glab/users/sz2766/paper_four/ERA5M/out/2var2/ERA5M_plsr_boot_2var_month_',
                       formatC(args,width=4,flag="0"),'.nc',sep=""),ncdem)
ncvar_put(ncout,varid=ncdem,out01)
nc_close(ncout)

