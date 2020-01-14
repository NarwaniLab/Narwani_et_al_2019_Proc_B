####################################################
#####  Phenotypic Trajectory Analysis	S.R. Dennis
####################################################

install.packages("car", dependencies=TRUE,type='binary')
library("car")
library("ellipse")
library("lubridate")
library("vegan")
library("dplyr")
require('tidyverse')

myColors<-c('grey50','#0570b0', '#7fbc41','#dd3497','#f46d43') #agreed upon colors

# CODE ASSUMES factor 1=lineage, factor 2=evol. levels
master<-read.csv("Narwani_et_al_2019_Proceedings_B_data.csv")

#factor
master[,"Treatment"]<-factor(master[,"Treatment"])
master[,"Dreisse"]<-factor(master[,"Dreisse"])
master[,"Macrophytes"]<-factor(master[,"Macrophytes"])
master[,"Nutrients"]<-factor(master[,"Nutrients"])
master[,"Block"]<-factor(master[,"Block"])
master$newdate<-mdy(as.character(master$Date))


exp1week = master[master$Date=="3/6/2017",]
exp1week$Week
#this tells us that for the first part of the experiment we only want to deal with data where Week<35, which is anything before March 2017

basedates<-c("2016-08-08","2016-08-22","2016-09-05","2016-09-19","2016-10-10")
nutdats<-as.Date(c("2016-08-12", "2016-08-26", "2016-09-09", "2016-09-22", "2016-10-10"))
#dates on which we added nutrients

datesPTAplot<-c("2016-08-08","2016-10-17","2017-02-13")
submaster = master[master$Week<35,] #should be week 35

submaster[,"Week"]<-factor(submaster[,"Week"])
master[,"Week"]<-factor(master[,"Week"])

completesp<-submaster[complete.cases(submaster$c_Cyanobacteria), ] 

which( colnames(completesp)=="v_Cyanobacteria" )

#columns with functional group volume
gpvol <- data.frame(completesp[,97:102])
gpcomp <- decostand(gpvol, method = "hellinger",na.rm=T)
colnames(gpcomp)<-paste("h",colnames(gpvol),sep = "_")
species<-colnames(gpcomp)

completesp<-cbind(completesp,gpcomp)


beforedates<- as.Date(c("2016-07-25","2016-08-02","2016-08-08"))
duringdates<- as.Date(c("2016-08-15","2016-08-22","2016-08-29","2016-09-05","2016-09-12","2016-09-19","2016-09-26","2016-10-03","2016-10-10"))
afterdates<- as.Date(c("2016-10-17","2016-10-24","2016-10-31","2016-11-14","2016-11-24","2016-11-28","2016-12-19","2017-01-16","2017-02-13"))


wrk<-data.frame(completesp[,c(1:10,97:102)])
wrk$Date<-mdy(as.character(wrk$Date))

wrk$period[ymd(wrk$Date) < "2016-07-25"]<-"A" #transient
wrk$period[ymd(wrk$Date) %in%  beforedates]<-"B" #before
wrk$period[ymd(wrk$Date) %in% duringdates]<-"C" #during
wrk$period[ymd(wrk$Date) %in% afterdates]<-"D" #after


taxa<-as.factor(wrk$Treatment) # a different trajectory will be calculated for each of these
# for NUDREM taxa = Treatment  
level<-as.factor(wrk$period)  # this is for the points along the trajectories from taxa
# for NUDREM level = date = period

head(wrk)
gpvol <- data.frame(wrk[,-c(1:10,17)])# remove descriptive cols
gpcomp <- decostand(gpvol, method = "hellinger",na.rm=T)                    

Y<-gpcomp	## Removes factor columns from rest of data - assumes only two factor columns
## for NUDREM this is the matrix of transformed community biovolume data
taxaBylevel<-as.factor(paste(taxa,level))  #for LS means at each level
Y<-prcomp(Y)$x	#TO USE IF > 2 variables: for plotting PC axes at the end #,center=T,scale=T
OO<-prcomp(gpcomp) #,center=T,scale=T
summary(OO)

row.names(Y)<-row.names(wrk)


n<-length(levels(taxa)) 
p<-length(levels(level))  #
k<-ncol(Y) #dimensions of phenotypic data


lm.full<-lm(Y~taxa*level,model=T,x=T,y=T,qr=T) 
yhat.full<-predict(lm.full)
summary(manova(lm.full))



#For resid randomization (used later)
lm.red<-lm(Y ~ taxa + level, model=T,x=T,y=T,qr=T)   #just the covariate
yhat.red<-predict(lm.red)
res.red<-resid(lm.red)



lsmeans.obs<-NULL
for (i in 1:k){
  mean.temp<-tapply(yhat.full[,i],taxaBylevel,mean)
  lsmeans.obs<-cbind(lsmeans.obs,mean.temp)
}



arrayspecs<-function(A){
  n.m<-NULL # n.m stands for "new" means
  for(i in 1:n){
    temp<-as.matrix(A[((1+(i-1)*p):(i*p)),1:k])
    n.m<-cbind(n.m,temp)}
  trajectories<-array(n.m,dim=c(p,k,n))
}



####Functions calculating Trajectory attributes (don't edit anything in here)

# Pathlength distance
pathdist<-function(M) {as.matrix(dist(M))}
trajsize<-function(M){
  traj.pathdist<-array(0,dim=c(n,1))   		#loop across trajectories
  for (i in 1:n){
    temp<-pathdist(M[,,i])
    for (j in 1:(p-1)){
      traj.pathdist[i]<-traj.pathdist[i]+temp[j,j+1]
    }
  }
  traj.size.dist<-as.matrix(dist(traj.pathdist))		#TrajSizeDiff
}

#trajectory direction 
orient<-function(M) {(svd(var(M))$v[1:k,1])} 	#find orientation
trajorient<-function(M){
  traj.orient<-array(NA,dim=c(n,k))   		#loop across trajectories
  check.1<-array(NA,dim=c(n))
  for (i in 1:n){
    temp<-orient(M[,,i])
    traj.orient[i,]<-temp
    check.1[i]<-M[1,,i]%*%traj.orient[i,]  #check startingpoint location
    check.1[i]<-check.1[i]/abs(check.1[i])
    if(check.1[i]==-1) traj.orient[i,]<--1*traj.orient[i,]
  }
  options(warn=-1)				#b/c acos of 1 (e.g, on diagonal) yields warning
  traj.ang.diff<-(180/pi)*acos(traj.orient%*%t(traj.orient))
  #diag(traj.ang.diff)<-0
}

#trajectory shape
### GPA: following J. Claude 2008: Morphometrics in R
trans<-function(A){scale(A,scale=F)} 	##TRANSLATION
csize<-function(A)				##CSIZE
{p<-dim(A)[1]
size<-sqrt(sum(apply(A,2,var))*(p-1))
list("centroid_size"=size,"scaled"=A/size)}
mshape<-function(A){apply(A,c(1,2),mean)}	#meanshape	

pPsup<-function(M1,M2){				## OPA rotation 1-->2
  k<-ncol(M1)
  Z1<-trans(csize(M1)[[2]])
  Z2<-trans(csize(M2)[[2]])
  sv<-svd(t(Z2)%*%Z1)
  U<-sv$v; V<-sv$u; Delt<-sv$d
  sig<-sign(det(t(Z1)%*%Z2))
  Delt[k]<-sig*abs(Delt[k]); V[,k]<-sig*V[,k]
  Gam<-U%*%t(V)
  beta<-sum(Delt)
  list(Mp1=beta*Z1%*%Gam,Mp2=Z2,rotation=Gam,scale=beta,
       df=sqrt(1-beta^2))}

pgpa<-function(A)
{p<-dim(A)[1]; k<-dim(A)[2]; n<-dim(A)[3]  
temp2<-temp1<-array(NA,dim=c(p,k,n)); Siz<-numeric(n)#; Qm2<-numeric(n)
for (i in 1:n)
{Acs<-csize(A[,,i])
Siz[i]<-Acs[[1]]
temp1[,,i]<-trans(Acs[[2]])}
Qm1<-dist(t(matrix(temp1,k*p,n)))
Q<-sum(Qm1); iter<-0
while (abs(Q)> 0.00001)
{for(i in 1:n){
  M<-mshape(temp1[,,-i])
  temp2[,,i]<-pPsup(temp1[,,i],M)[[1]]}
  Qm2<-dist(t(matrix(temp2,k*p,n)))
  Q<-sum(Qm1)-sum(Qm2)
  Qm1<-Qm2
  iter=iter+1
  temp1<-temp2}
list("rotated"=temp2,"it.number"=iter,"Q"=Q,"intereucl.dist"=Qm2,"mshape"=
       csize(mshape(temp2))[[2]],"cent.size"=Siz)
}

## loop for GPA and shape distances
trajshape<-function(M){
  x<-pgpa(M)
  traj.shape.dist<-as.matrix(x$intereucl.dist) 
}

## TrajectorySummaryStat
sumstat<-function(M){
  M<-as.dist(M)
  x<-if(length(M)>1)(x=var(M)) else 0
}



##################MAIN LOOP################
traj.specs.obs<-arrayspecs(lsmeans.obs)
trajsize.obs<-trajsize(traj.specs.obs)
trajdir.obs<-trajorient(traj.specs.obs)
diag(trajdir.obs)<-0 #b/c some NA/N values on diagonal)
trajshape.obs<-trajshape(traj.specs.obs)
sumstatsize.obs<-sumstat(trajsize.obs)
sumstatdir.obs<-sumstat(trajdir.obs)
sumstatshape.obs<-sumstat(trajshape.obs)

### PERMUTATION PROCEDURE
permute<-9999
line<-nrow(Y)
PSize<-POrient<-PShape<-array(1,dim=c(n,n))
PSumSize<-PSumOrient<-PSumShape<-1
for(i in 1:permute){
  line.rand<-sample(line,replace=FALSE)
  res.temp<-cbind(line.rand,res.red)
  z<-(order(line.rand))
  res.temp2<-as.matrix(res.temp[z,])
  res.p<-res.temp2[,-1]  # Rows of residuals are now randomized
  y.rand<-yhat.red+res.p
  yhat.rand<-predict(lm(y.rand~taxa*level,model=T,x=T,y=T,qr=T))
  
  lsmeans.rand<-NULL
  for (j in 1:k){
    mean.temp<-tapply(yhat.rand[,j],taxaBylevel,mean)
    lsmeans.rand<-cbind(lsmeans.rand,mean.temp)
  }
  traj.specs.rand<-arrayspecs(lsmeans.rand)
  trajsize.rand<-trajsize(traj.specs.rand)
  trajdir.rand<-trajorient(traj.specs.rand)
  diag(trajdir.rand)<-0 #b/c some NA/N values on diagonal)
  trajshape.rand<-trajshape(traj.specs.rand)
  sumstatsize.rand<-sumstat(trajsize.rand)
  sumstatdir.rand<-sumstat(trajdir.rand)
  sumstatshape.rand<-sumstat(trajshape.rand)
  for(j in 1:n){
    for(jj in 1:n){
      PSize[j,jj]<-if(trajsize.rand[j,jj]>=trajsize.obs[j,jj]) 
        (PSize[j,jj]+1) else PSize[j,jj]
      POrient[j,jj]<-if(trajdir.rand[j,jj]>=trajdir.obs[j,jj]) 
        (POrient[j,jj]+1) else POrient[j,jj]
      PShape[j,jj]<-if(trajshape.rand[j,jj]>=trajshape.obs[j,jj]) 
        (PShape[j,jj]+1) else PShape[j,jj]
    }
  }
  PSumSize<-if(sumstatsize.rand>=sumstatsize.obs) (PSumSize+1) else PSumSize
  PSumSOrient<-if(sumstatdir.rand>=sumstatdir.obs) (PSumOrient+1) else PSumOrient
  PSumShape<-if(sumstatshape.rand>=sumstatshape.obs) (PSumShape+1) else PSumShape
  print(paste(i,"/",permute))
}  #end permute


for(j in 1:n){
  for(jj in 1:n){
    PSize[j,jj]<-PSize[j,jj]/(permute+1)
    POrient[j,jj]<-POrient[j,jj]/(permute+1)
    PShape[j,jj]<-PShape[j,jj]/(permute+1)
  }
}
PSumSize<-PSumSize/(permute+1)
PSumOrient<-PSumOrient/(permute+1)
PSumShape<-PSumShape/(permute+1)

asize <- triu(trajsize.obs,1) ; psize <- tril(PSize,-1)
size_tab <- as.data.frame(as.matrix(asize + psize)) ; diag(size_tab)<-NA ; names(size_tab) <-row.names(size_tab)<- levels(taxa) ; size_tab
adir <- triu(trajdir.obs,1) ; porient <- tril(POrient,-1)
dir_tab <- as.data.frame(as.matrix(adir + porient)) ; diag(dir_tab)<-NA ; names(dir_tab) <-row.names(dir_tab)<- levels(taxa) ; dir_tab
ashape <- triu(trajshape.obs,1) ; pshape <- tril(PShape,-1)
shape_tab <- as.data.frame(as.matrix(ashape + pshape)) ; diag(shape_tab)<-NA ; names(shape_tab) <-row.names(shape_tab)<- levels(taxa) ; shape_tab

######NOTE: summary statistics only valid for 2+ trajectories!!
sumstatsize.obs
PSumSize
sumstatdir.obs
PSumOrient
sumstatshape.obs
PSumShape

summary(OO)

#########Plot Data and LS means 
#n number of treatments
#p number of time points

Rownames=rownames(Y)
rows = Rownames [which(wrk$period == "D")] 
gr=taxa[41:60]
algaegroups<-c("Cyanobacteria", "Golden algae", "Diatoms", "Dinoflagellates", "Cryptophytes", "Green algae")


tiff(paste("PTA_period.tiff"), width = 14, height = 7, pointsize = 1/2000, units = 'in', res = 450)

pcx=1
pcy=2
par(mfrow=c(1,2), cex.lab=1.6,cex.axis=1.6, pin=c(5,5))
plot(Y[,pcx:pcy],type="n",xlab=paste("PC",pcx," (34%)"), ylab=paste("PC",pcy," (33%)",sep=""), xlim=c(-0.2,0.2),ylim=c(-1,1), asp=1)

for (i in 1:n){		 	
  for (j in 1:(p-1)){		
    points(traj.specs.obs[(j:(j+1)), pcx,i], traj.specs.obs[(j:(j+1)), pcy,i],type="l", pch=21, col=myColors[i])    
  }
}

for (i in 1:n){		 	
  for (j in 2:(p-1)){		
    points(traj.specs.obs[j,pcx,i], traj.specs.obs[j,pcy,i], pch=21, bg=myColors[i], col=myColors[i], cex=0.5)
  }
}

#STARTING POINT
for (i in 1:n){
  points(traj.specs.obs[1,pcx,i], traj.specs.obs[1,pcy,i], pch=21,bg="white", cex=2,lwd=2, col=myColors[i])
}

#ENDING POINT
for (i in 1:n){
  points(traj.specs.obs[p,pcx,i], traj.specs.obs[p,pcy,i], pch=21, bg=myColors[i], cex=2, col=myColors[i])
}

A=as.numeric(Y[rows,pcx])
B=as.numeric(Y[rows,pcy])
groupmultiplier = length(A)/20
sp<-substr(species,start=5,stop=100)
grps=rep(gr,groupmultiplier)


dataEllipse(A,B,groups=grps,group.labels="",col=myColors,center.cex=1.5,add=T,levels=0.5,lwd=0,fill=T,fill.alpha=0.05,plot.points=F) 
arrows(rep(0,7),rep(0,7),OO$rotation[,pcx],OO$rotation[,pcy],col="black", length=0.1, angle=40,lwd=1)
text(x=OO$rotation[,pcx]+ c(0.23,-0.02,-0.15,0,0,0), y=OO$rotation[,pcy]+c(0.06,-0.04,-0.02,0.05,-0.05,0.06), label=algaegroups, font=3, cex=1.2)
text(x=0.9, y=0.9, label="a", cex=1.8, font=1.5)


pcx=2
pcy=3

plot(Y[,pcx:pcy],type="n",xlab=paste("PC",pcx," (33%)",sep=""), ylab=paste("PC",pcy," (26%)",sep=""), asp=1,  xlim=c(-1,1), ylim=c(-1,1), cex=1.2) 

for (i in 1:n){		 	
  for (j in 1:(p-1)){		
    points(traj.specs.obs[(j:(j+1)),pcx,i],traj.specs.obs[(j:(j+1)),pcy,i],type="l",pch=21, col=myColors[i])    
  }
}

for (i in 1:n){		 	
  for (j in 2:(p-1)){		
    points(traj.specs.obs[j,pcx,i],traj.specs.obs[j,pcy,i],pch=21,bg=myColors[i],col=myColors[i],cex=0.5)
  }
}

#STARTING POINT
for (i in 1:n){
  points(traj.specs.obs[1,pcx,i],traj.specs.obs[1,pcy,i],pch=21,bg="white", cex=2,lwd=2,col=myColors[i])
}

#ENDING POINT
for (i in 1:n){
  points(traj.specs.obs[p,pcx,i],traj.specs.obs[p,pcy,i],pch=21,bg=myColors[i],cex=2,col=myColors[i])
}

A=as.numeric(Y[rows,pcx])
B=as.numeric(Y[rows,pcy])
groupmultiplier = length(A)/20
sp<-substr(species,start=5,stop=100)
grps=rep(gr,groupmultiplier)

dataEllipse(A,B,groups=grps,group.labels="",col=myColors,center.cex=1.5,add=T,levels=0.5,lwd=0,fill=T,fill.alpha=0.05,plot.points = F) 
arrows(rep(0,7),rep(0,7),OO$rotation[,pcx],OO$rotation[,pcy],col="black",length=0.1, angle=40,lwd=1)
text(x=OO$rotation[,pcx]+  c(0,-0.25,-0.05,0,-0.25,0), y=OO$rotation[,pcy]+c(0.05,0,-0.05,-0.05,0,-0.05), label=algaegroups, font=3, cex=1.2)
text(x=0.9, y=0.9, label="b", cex=1.8, font=1.5)
dev.off()

