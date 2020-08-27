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

args = (commandArgs(trailingOnly = F))
print(args)
myargs <-sub('-','',args[length(args)])
args = as.numeric(myargs)
print(args)

##########################################
CAN<-list.files('/rigel/glab/users/sz2766/paper_four/model/data/GFDL/',pattern='.nc',full.names=TRUE)
###################01#####################
nc<-nc_open(CAN[1])
v1<-nc$var[[4]]
lhc<-ncvar_get(nc,v1)*0.0864/2.499
lhc[lhc>1e19]=NA
lhc<-lhc[,args,]
dim(lhc)<-c(144,12*150)

nc<-nc_open(CAN[3])
v1<-nc$var[[1]]
prc<-ncvar_get(nc,v1)*86400
prc[prc>1e19]=NA
prc<-prc[,args,]
dim(prc)<-c(144,12*150)

nc<-nc_open(CAN[2])
v1<-nc$var[[3]]
msc<-ncvar_get(nc,v1)
msc[msc< 0]=NA
msc<-msc[,args,]
dim(msc)<-c(144,12*150)

mtc<-prc-lhc

pr1<-prc[,337:816]
et1<-lhc[,337:816]
mt1<-mtc[,337:816]
sm1<-msc[,337:816]

pr2<-prc[,1321:1800]
et2<-lhc[,1321:1800]
mt2<-mtc[,1321:1800]
sm2<-msc[,1321:1800]

#################
input1<-cbind(mt1,sm1)#SM-(P-E):1979-2018
input2<-cbind(pr1,sm1)#SM-P:1979-2018
input3<-cbind(et1,sm1)#SM-E:1979-2018
input4<-cbind(mt2,sm2)#SM-(P-E):2061-2100
input5<-cbind(pr2,sm2)#SM-P:2061-2100
input6<-cbind(et2,sm2)#SM-E:2061-2100

cl<-makeCluster(getOption("cl.cores",5))
clusterEvalQ(cl,library(parallel))
clusterEvalQ(cl,library(ncdf4))
clusterEvalQ(cl,library(plsRglm))
clusterEvalQ(cl,library(plsdof))
clusterEvalQ(cl,library(boot))

out1<-parRapply(cl,input1,pls_reg)
out2<-parRapply(cl,input2,pls_reg)
out3<-parRapply(cl,input3,pls_reg)
out4<-parRapply(cl,input4,pls_reg)
out5<-parRapply(cl,input5,pls_reg)
out6<-parRapply(cl,input6,pls_reg)

dim(out1)<-c(4,144)
dim(out2)<-c(4,144)
dim(out3)<-c(4,144)
dim(out4)<-c(4,144)
dim(out5)<-c(4,144)
dim(out6)<-c(4,144)

out01<-rbind(out1,out2,out3,out4,out5,out6)

#lat<-ncvar_get(nc,varid = 'lat')
long<-ncvar_get(nc,varid = 'lon')
result1<-1:24

#dimlat<-ncdim_def('latitude','deg',lat)
dimlong<-ncdim_def('longitude','deg',long)
dimout1<-ncdim_def('result','',result1)

ncdem<-ncvar_def('r','',list(dimout1,dimlong),-9999,longname="r",prec='double')
ncout<-nc_create(paste('/rigel/glab/users/sz2766/paper_four/model/out/PLSR40/GFDL/gfdl_plsr_month_',
                       formatC(args,width=4,flag="0"),'.nc',sep=""),ncdem)
ncvar_put(ncout,varid=ncdem,out01)
nc_close(ncout)

