##this script tested packages: curvclust, funHDDC, and FunCluster on clustering time series curves
require(fda)
require(curvclust)
require(funHDDC)
require(FunCluster)
require(ggplot2)
require(R.matlab)
#example 1
## test on the data set Berkeley Growth Data
K=2
###curvclust {Wavelet}, not good design, not good performance
fdat=list()
data=growth[["hgtf"]]
size=dim(data)
dataarr=as.data.frame(matrix(NA,nrow=size[1]*size[2],ncol=4))
colnames(dataarr)=c("person","age","height","cluster")
for(j in 1:ncol(data)){
  vec=data[,j]
  names(vec)=NULL
  fdat[[j]]=vec
}
dataarr[,"height"]=unlist(fdat)
dataarr[,"person"]=rep(seq(size[2]),each=size[1])
dataarr[,"age"]=rep(growth[["age"]],times=size[2])
#data
CCD=new("CClustData",Y=fdat,filter.number=1)#smoothness
CCDred=getUnionCoef(CCD)
#options
CCO=new("CClustO")
CCO["nbclust"]=K
CCO["Gamma2.structure"]="none"
CCR=getFCMM(CCDred,CCO)#Functional Clustering Mixed Models {FCMM,FCM,FMM}
summary(CCR)
cluster=apply(CCR["Tau"],1,which.max)
dataarr[,"cluster"]=rep(cluster,each=size[1])
p<-ggplot(data=dataarr,aes(age,height,colour=cluster,group=person))+
      geom_line(alpha=0.5)+
      xlab("age")+
      ylab("height")+
      theme_bw()

###funHDDC
age=growth[["age"]]
heightmat=growth[["hgtf"]]
norder=6
nbasis=length(age)+norder-2
heightbasis=create.bspline.basis(c(1,18),nbasis,norder,age)
heightfdPar=fdPar(heightbasis,4,0.01)
heightfd=smooth.basis(age,heightmat,heightfdPar)$fd
res.uni<-funHDDC(heightfd,K=2,model="AkBkQkDk",init="random",threshold=0.2)
plot(heightfd,col=res.uni$class)
# slopeHeuristic(res.uni)
res.pca<-mfpca(heightfd)
plot.mfpca(res.pca)
#FunCluster: really domain specific(hard to code with)

#example 2
## test on ridge tracking data set
load("/Users/mikeaalv/Dropbox (Edison_Lab@UGA)/Projects/Bioinformatics_modeling/spectral.related/result/compound.quan.record.absmore.upd.RData")
experi=1
lenlist=sapply(list.quan,function(x){
  dim(x)[1]
})
dataarr=as.data.frame(cbind(Reduce(rbind,list.quan),rep(names(list.quan),times=lenlist)))
colnames(dataarr)=c("time","intensity","exp","compd")
dataarr=dataarr[dataarr[,"exp"]==experi,]
# timecount=table(dataarr[,"time"])
# allexitime=names(timecount)[timecount==max(timecount)]
dataarr=dataarr[dataarr[,"time"]<7.2 & dataarr[,"time"]>0.8,]
dataarr=dataarr[order(dataarr[,"compd"]),]
compds=unique(dataarr[,"compd"])
for(compd in compds){
  ind=dataarr[,"compd"]==compd
  dataarr[ind,"intensity"]=scale(dataarr[ind,"intensity"])
}
time=sort(unique(dataarr[,"time"]))
###curvclust
fdat=list()
for(compd_i in seq(length(compds))){
  compd=compds[compd_i]
  vec=dataarr[dataarr[,"compd"]==compd,"intensity"]
  names(vec)=NULL
  fdat[[compd_i]]=vec
}
#data
CCD=new("CClustData",Y=fdat,filter.number=1)#smoothness
CCDred=getUnionCoef(CCD)
#options
CCO=new("CClustO")
CCO["nbclust"]=3
CCO["Gamma2.structure"]="none"
CCR=getFCM(CCDred,CCO)#Functional Clustering Mixed Models {FCMM,FCM,FMM}
summary(CCR)
cluster=apply(CCR["Tau"],1,which.max)
dataarrtemp=dataarr
dataarrtemp[,"cluster"]=as.factor(rep(cluster,each=length(time)))
p<-ggplot(data=dataarrtemp,aes(time,intensity,colour=cluster,group=compd))+
      geom_line(alpha=0.5)+
      xlab("time")+
      ylab("intensity")+
      theme_bw()

###funHDDC
list_intensity=sapply(compds,simplify=FALSE,function(x){
  tempmat=dataarr[dataarr[,"compd"]==x,]
  tempmat[order(tempmat[,"time"]),"intensity"]
})
intenmat=Reduce(cbind,list_intensity)
norder=6
nbasis=length(time)+norder-2
intbasis=create.bspline.basis(c(min(time),max(time)),nbasis,norder,time)
intfdPar=fdPar(intbasis,4,0.01)
intfd=smooth.basis(time,intenmat,intfdPar)$fd
res.uni<-funHDDC(intfd,K=3,model="AkBkQkDk",init="random",threshold=0.2,nb.rep=50)
plot(intfd,col=res.uni$class)
# slopeHeuristic(res.uni)
res.pca<-mfpca(intfd)
plot.mfpca(res.pca)
