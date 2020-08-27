library(ncdf4)
library(parallel)
library(plsRglm)
library(plsdof)
library(boot)

pls_reg<-function(x){
  if (sum(is.na(x))>800){
    return(rep(NA,each=4))
  }

  tt<-length(x)/2

  pr<-x[1:tt]
  sm<-x[(1+tt):(tt*2)]

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

    reg1_boot<-bootpls(reg1, R=2000, verbose=FALSE)
    reg2_boot<-bootpls(reg2, R=2000, verbose=FALSE)

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

args = (commandArgs(trailingOnly = F))
print(args)
myargs <-sub('-','',args[length(args)])
args = as.numeric(myargs)
print(args)

##########################################
CAN<-list.files('/rigel/glab/users/sz2766/paper_four/model/data/ECHAM6/',pattern='.nc',full.names=TRUE)
###################01#####################
nc<-nc_open(CAN[1])
v1<-nc$var[[1]]
lhc<-ncvar_get(nc,v1)*0.0864/2.499
lhc[lhc>1e19]=NA
lhc<-lhc[,args,]
dim(lhc)<-c(192,12*151)
lhc<- -lhc

nc<-nc_open(CAN[3])
v1<-nc$var[[1]]
prc<-ncvar_get(nc,v1)*86400
prc[prc>1e19]=NA
prc<-prc[,args,]
dim(prc)<-c(192,12*151)

nc<-nc_open(CAN[2])
v1<-nc$var[[1]]
msc<-ncvar_get(nc,v1)
msc[msc< 0]=NA
msc<-msc[,args,]
dim(msc)<-c(192,12*151)

mtc<-prc-lhc

pr<-prc[,253:1812]
et<-lhc[,253:1812]
mt<-mtc[,253:1812]
sm<-msc[,253:1812]

#################
input1<-cbind(mt,sm)#SM-(P-E)
input2<-cbind(pr,sm)#SM-P
input3<-cbind(et,sm)#SM-E

cl<-makeCluster(getOption("cl.cores",5))
clusterEvalQ(cl,library(parallel))
clusterEvalQ(cl,library(ncdf4))
clusterEvalQ(cl,library(plsRglm))
clusterEvalQ(cl,library(plsdof))
clusterEvalQ(cl,library(boot))

out1<-parRapply(cl,input1,pls_reg)
out2<-parRapply(cl,input2,pls_reg)
out3<-parRapply(cl,input3,pls_reg)

dim(out1)<-c(4,192)
dim(out2)<-c(4,192)
dim(out3)<-c(4,192)

out01<-rbind(out1,out2,out3)

#lat<-ncvar_get(nc,varid = 'lat')
long<-ncvar_get(nc,varid = 'lon')
result1<-1:12

#dimlat<-ncdim_def('latitude','deg',lat)
dimlong<-ncdim_def('longitude','deg',long)
dimout1<-ncdim_def('result','',result1)

ncdem<-ncvar_def('r','',list(dimout1,dimlong),-9999,longname="r",prec='double')
ncout<-nc_create(paste('/rigel/glab/users/sz2766/paper_four/model/out/ECHAM6/echam6_plsr_month_',
                       formatC(args,width=4,flag="0"),'.nc',sep=""),ncdem)
ncvar_put(ncout,varid=ncdem,out01)
nc_close(ncout)
