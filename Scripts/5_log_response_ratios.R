theme_set(theme_cowplot())
master<-read.csv("Narwani_et_al_2019_Proceedings_B_data.csv")

#factor
master[,"Treatment"]<-factor(master[,"Treatment"])
master[,"Dreisse"]<-factor(master[,"Dreisse"])
master[,"Macrophytes"]<-factor(master[,"Macrophytes"])
master[,"Nutrients"]<-factor(master[,"Nutrients"])
master[,"Block"]<-factor(master[,"Block"])
master$newdate<-mdy(as.character(master$Date))
master$TP<-as.numeric(master$TP)

exp1week = master[master$Date=="3/6/2017",]
exp1week$Week
#this tells us that for the first part of the experiment we only want to deal with data where Week<35, which is anything before March 2017

nutdats<-as.Date(c("2016-08-12", "2016-08-26", "2016-09-09", "2016-09-22", "2016-10-10"))
#dates on which we added nutrients

myColors_sub<-c('#0570b0', '#7fbc41','#dd3497','#000000') # colors
submaster = master[master$Week<35,] #should be week 35

submaster[,"Week"]<-factor(submaster[,"Week"])
master[,"Week"]<-factor(master[,"Week"])

subdates<-unique(submaster$newdate)

submaster$v_Lagerheimia_1 = submaster$v_Lagerheimia_1 + submaster$v_Lagerheimia_2
submaster = subset(submaster, select = -c(v_Lagerheimia_2) )

names(submaster)[names(submaster) == 'v_Lagerheimia_1'] <- 'v_Lagerheimia'
head(submaster)

which(colnames(submaster)=="v_Cyanobacteria")
which(colnames(submaster)=="v_Green_algae")

which(colnames(submaster)=="v_Achnthaceae")
which(colnames(submaster)=="v_Unknown_13")


#columns with functional group volume
gpvol <- data.frame(submaster[,97:102])
#columns with species volume
spvol <- data.frame(submaster[,103:167])

########
#hellinger transform of species biovolume
library(vegan)
spcomp <- decostand(spvol, method = "hellinger",na.rm=T)
colnames(spcomp)<-paste("h",colnames(spcomp),sep="_")
submaster<-cbind(submaster,spcomp)

#hellinger transform of functional group biovolume
gpcomp <- decostand(gpvol, method = "hellinger",na.rm=T)
colnames(gpcomp)<-paste("h",colnames(gpcomp),sep = "_")
submaster<-cbind(submaster,gpcomp)

#log transform of species biovolume
spcomp <- decostand(spvol, method = "log",na.rm=T)
colnames(spcomp)<-paste("l",colnames(spcomp),sep="_")
submaster<-cbind(submaster,spcomp)

#log transfrom of functional group biovolume
gpcomp <- decostand(gpvol, method = "log",na.rm=T)
colnames(gpcomp)<-paste("l",colnames(gpcomp),sep="_")
submaster <-cbind(submaster,gpcomp)

head(submaster)

pd <- position_dodge(width = 2)


####################
# Chlorophyll-a LRR
####################

completechl<-submaster[complete.cases(submaster$CHL), ] 
L<-length(unique(completechl$newdate))
dates<-unique(completechl$newdate)
LRR_Trmt_CHL<-matrix(NA,L*12,4)
LRR_Control_CHL<-matrix(NA,L*4,4)

count=1
for (i in 1:L){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12

  subdat<-subset(completechl,newdate==dates[i])
  controls<-subset(subdat,Treatment=="C")
  others<-filter(subdat,Treatment!="C"&Treatment!="O")
  mean_C_CHL<- mean(controls$CHL)

  LRR_Trmt_CHL[count,1]<-as.character(others$Pond[j])
  LRR_Trmt_CHL[count,2]<-as.character(dates[i])
  LRR_Trmt_CHL[count,3]<-as.character(others$Treatment[j])
  LRR_Trmt_CHL[count,4]<-log(others$CHL[j]/mean_C_CHL)
  count=count+1
  }
}
LRR_Trmt_CHL<-data.frame(LRR_Trmt_CHL)
colnames(LRR_Trmt_CHL)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_CHL$Treatment<-factor(LRR_Trmt_CHL$Treatment)
LRR_Trmt_CHL$newdate<-as.Date(LRR_Trmt_CHL$newdate)
LRR_Trmt_CHL$LRR<-as.numeric(as.character(LRR_Trmt_CHL$LRR))

measurevar= "LRR"
groupvars=(c("Treatment","newdate"))
conf.interval=.95

x <- summarySE(LRR_Trmt_CHL, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$LRR), ] 

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_CHL<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_CHL<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_CHL+mean_D_CHL
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0

tiff("Chlorophyll_LRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
chl_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+

  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR Chl-a")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9))+
        #,
        #legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("a")

plot(chl_LRR)
dev.off()


##################
# Cyan dominance
##################


rm(completesp)
#only taking the portion of the data where we have species composition estimates
completesp<-submaster[complete.cases(submaster$h_v_Cyanobacteria), ] 

#columns with functional group volume
gpvol <- data.frame(completesp[,101:106])

#columns with species volume
spvol <- data.frame(completesp[,107:171])

#hellinger transform of functional group biovolume
gpcomp <- decostand(gpvol, method = "hellinger",na.rm=T)
colnames(gpcomp)<-paste("h",colnames(gpcomp),sep="_")
completesp<-cbind(completesp,gpcomp)

#hellinger transform of species biovolume
spcomp <- decostand(spvol, method = "hellinger",na.rm=T)
colnames(spcomp)<-paste("h",colnames(spcomp),sep="_")
completesp<-cbind(completesp,spcomp)


L<-length(unique(completesp$newdate))
dates<-unique(completesp$newdate)
LRR_Trmt_cyan<-matrix(NA,L*12,4)
LRR_Control_cyan<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completesp,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_cyan<- mean(controls$h_v_Cyanobacteria)
    
    LRR_Trmt_cyan[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_cyan[count,2]<-as.character(dates[i])
    LRR_Trmt_cyan[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_cyan[count,4]<-log((others$h_v_Cyanobacteria[j]/mean_C_cyan)+1)
    count=count+1
  }
}

LRR_Trmt_cyan<-data.frame(LRR_Trmt_cyan)
colnames(LRR_Trmt_cyan)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_cyan$Treatment<-factor(LRR_Trmt_cyan$Treatment)
LRR_Trmt_cyan$newdate<-as.Date(LRR_Trmt_cyan$newdate)
LRR_Trmt_cyan$LRR<-as.numeric(as.character(LRR_Trmt_cyan$LRR))

x <- summarySE(LRR_Trmt_cyan, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$LRR), ] 

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_cyan<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_cyan<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_cyan+mean_D_cyan
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0

extramatrix = data.frame(matrix(NA,3,7))
extramatrix[,1] = as.character(c("D","M","MD"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)

extramatrix2 = data.frame(matrix(NA,3,7))
extramatrix2[,1] = as.character(c("D","M","MD"))
extramatrix2[,2] = ymd(c("2016-07-04","2016-07-04","2016-07-04"))
xnames2<-names(x)
names(extramatrix2)<-xnames2
x<-rbind(x,extramatrix2)
dev.off()

tiff("Cyano_dom_LRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
cyan_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0.6931472, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR cyano dominance")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("b")

plot(cyan_LRR)
dev.off()
rm(Additive)
##################
# Trait evenness
##################
rm(Additive,completesp,x,LRR_Trmt_TED,LRR_Control_TED)
L<-length(unique(completechl$newdate))
dates<-unique(completechl$newdate)
LRR_Trmt_TED<-matrix(NA,L*12,4)
LRR_Control_TED<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completechl,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_TED<- mean(controls$TED)
    
    LRR_Trmt_TED[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_TED[count,2]<-as.character(dates[i])
    LRR_Trmt_TED[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_TED[count,4]<-log((others$TED[j]/mean_C_TED)+1)
    count=count+1
  }
}


LRR_Trmt_TED<-data.frame(LRR_Trmt_TED)
colnames(LRR_Trmt_TED)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_TED$Treatment<-factor(LRR_Trmt_TED$Treatment)
LRR_Trmt_TED$newdate<-as.Date(LRR_Trmt_TED$newdate)
LRR_Trmt_TED$LRR<-as.numeric(as.character(LRR_Trmt_TED$LRR))

x <- summarySE(LRR_Trmt_TED, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$sd), ] 
#dim(x)

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_TED<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_TED<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_TED+mean_D_TED
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0

extramatrix = data.frame(matrix(NA,3,7))
extramatrix[,1] = as.character(c("D","M","MD"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)

tiff("TEDLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
TED_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept =  0.6931472, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR TED")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("c")

plot(TED_LRR)
dev.off()


##################
# CtoN
##################

completeCN<-submaster[complete.cases(submaster$CtoN), ] 

L<-length(unique(completeCN$newdate))
dates<-unique(completeCN$newdate)
LRR_Trmt_CN<-matrix(NA,L*12,4)
LRR_Control_CN<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completeCN,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_CN<- mean(controls$CtoN)
    
    LRR_Trmt_CN[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_CN[count,2]<-as.character(dates[i])
    LRR_Trmt_CN[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_CN[count,4]<-log(others$CtoN[j]/mean_C_CN)
    count=count+1
  }
}

LRR_Trmt_CN<-data.frame(LRR_Trmt_CN)
colnames(LRR_Trmt_CN)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_CN$Treatment<-factor(LRR_Trmt_CN$Treatment)
LRR_Trmt_CN$newdate<-as.Date(LRR_Trmt_CN$newdate)
LRR_Trmt_CN$LRR<-as.numeric(as.character(LRR_Trmt_CN$LRR))


x <- summarySE(LRR_Trmt_CN, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)

x <-x[complete.cases(x$sd), ] 

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_CN<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_CN<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_CN+mean_D_CN
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0

extramatrix = data.frame(matrix(NA,3,7))
extramatrix[,1] = as.character(c("D","M","MD"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)


tiff("CNLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
CN_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR C:N")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("d")

plot(CN_LRR)
dev.off()





##################
# CtoP
##################

completeCP<-submaster[complete.cases(submaster$CtoP), ] 

L<-length(unique(completeCP$newdate))
dates<-unique(completeCP$newdate)
LRR_Trmt_CP<-matrix(NA,L*12,4)
LRR_Control_CP<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completeCP,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_CP<- mean(controls$CtoP)
    
    LRR_Trmt_CP[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_CP[count,2]<-as.character(dates[i])
    LRR_Trmt_CP[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_CP[count,4]<-log(others$CtoP[j]/mean_C_CP)
    count=count+1
  }
}

LRR_Trmt_CP<-data.frame(LRR_Trmt_CP)
colnames(LRR_Trmt_CP)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_CP$Treatment<-factor(LRR_Trmt_CP$Treatment)
LRR_Trmt_CP$newdate<-as.Date(LRR_Trmt_CP$newdate)
LRR_Trmt_CP$LRR<-as.numeric(as.character(LRR_Trmt_CP$LRR))

x <- summarySE(LRR_Trmt_CP, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)

x <-x[complete.cases(x$sd), ] 

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_CP<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_CP<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_CP+mean_D_CP
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0

extramatrix = data.frame(matrix(NA,3,7))
extramatrix[,1] = as.character(c("D","M","MD"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)


tiff("CPLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
CP_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR C:P")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("e")

plot(CP_LRR)
dev.off()



tiff("Algal LRR response.tiff", width=30, height=12, pointsize=1/300, units='in', res = 300)
multiplot(chl_LRR, CN_LRR, cyan_LRR, CP_LRR, TED_LRR, cols=3)
dev.off()

































##################
# Nitrate
##################
completechl<-submaster[complete.cases(submaster$CHL), ] 
L<-length(unique(completechl$newdate))
dates<-unique(completechl$newdate)
LRR_Trmt_nit<-matrix(NA,L*12,4)
LRR_Control_nit<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completechl,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_nit<- mean(controls$Nitrate)
    
    LRR_Trmt_nit[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_nit[count,2]<-as.character(dates[i])
    LRR_Trmt_nit[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_nit[count,4]<-log(others$Nitrate[j]/mean_C_nit)
    count=count+1
  }
}

LRR_Trmt_nit<-data.frame(LRR_Trmt_nit)
colnames(LRR_Trmt_nit)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_nit$Treatment<-factor(LRR_Trmt_nit$Treatment)
LRR_Trmt_nit$newdate<-as.Date(LRR_Trmt_nit$newdate)
LRR_Trmt_nit$LRR<-as.numeric(as.character(LRR_Trmt_nit$LRR))

x <- summarySE(LRR_Trmt_nit, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$sd), ] 
dim(x)


Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_nit<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_nit<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_nit+mean_D_nit
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0



tiff("nitLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
nit_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR Nitrate-N")+
  
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+

  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("a")

plot(nit_LRR)
dev.off()








##################
# Phosphate
##################
completechl<-submaster[complete.cases(submaster$CHL), ] 
L<-length(unique(completechl$newdate))
dates<-unique(completechl$newdate)
LRR_Trmt_oP<-matrix(NA,L*12,4)
LRR_Control_oP<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completechl,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_oP<- mean(controls$orthoP)
    
    LRR_Trmt_oP[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_oP[count,2]<-as.character(dates[i])
    LRR_Trmt_oP[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_oP[count,4]<-log(others$orthoP[j]/mean_C_oP)
    count=count+1
  }
}

LRR_Trmt_oP<-data.frame(LRR_Trmt_oP)
colnames(LRR_Trmt_oP)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_oP$Treatment<-factor(LRR_Trmt_oP$Treatment)
LRR_Trmt_oP$newdate<-as.Date(LRR_Trmt_oP$newdate)
LRR_Trmt_oP$LRR<-as.numeric(as.character(LRR_Trmt_oP$LRR))


x <- summarySE(LRR_Trmt_oP, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$sd), ] 
dim(x)

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_oP<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_oP<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_oP+mean_D_oP
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0

tiff("oPLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
oP_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR Phosphate-P")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("b")

plot(oP_LRR)
dev.off()




##################
# Turbidity
##################
completeTurb<-submaster[complete.cases(submaster$Turbidity), ] 
L<-length(unique(completeTurb$newdate))
dates<-unique(completeTurb$newdate)
LRR_Trmt_turb<-matrix(NA,L*12,4)
LRR_Control_turb<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completeTurb,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_turb<- mean(controls$Turbidity)
    
    LRR_Trmt_turb[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_turb[count,2]<-as.character(dates[i])
    LRR_Trmt_turb[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_turb[count,4]<-log(others$Turbidity[j]/mean_C_turb)
    count=count+1
  }
}

LRR_Trmt_turb<-data.frame(LRR_Trmt_turb)
colnames(LRR_Trmt_turb)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_turb$Treatment<-factor(LRR_Trmt_turb$Treatment)
LRR_Trmt_turb$newdate<-as.Date(LRR_Trmt_turb$newdate)
LRR_Trmt_turb$LRR<-as.numeric(as.character(LRR_Trmt_turb$LRR))



x <- summarySE(LRR_Trmt_turb, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$sd), ] 
dim(x)

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_turb<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_turb<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_turb+mean_D_turb
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0


tiff("turbLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
turb_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR Turbidity")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("c")

plot(turb_LRR)
dev.off()




##################
# Dissolved oxygen
##################

completeDO<-submaster[complete.cases(submaster$O2), ] 
L<-length(unique(completeDO$newdate))
dates<-unique(completeDO$newdate)
LRR_Trmt_DO<-matrix(NA,L*12,4)
LRR_Control_DO<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completeDO,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_DO<- mean(controls$O2)
    
    LRR_Trmt_DO[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_DO[count,2]<-as.character(dates[i])
    LRR_Trmt_DO[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_DO[count,4]<-log(others$O2[j]/mean_C_DO)
    count=count+1
  }
}


LRR_Trmt_DO<-data.frame(LRR_Trmt_DO)
colnames(LRR_Trmt_DO)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_DO$Treatment<-factor(LRR_Trmt_DO$Treatment)
LRR_Trmt_DO$newdate<-as.Date(LRR_Trmt_DO$newdate)
LRR_Trmt_DO$LRR<-as.numeric(as.character(LRR_Trmt_DO$LRR))


x <- summarySE(LRR_Trmt_DO, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$sd), ] 
dim(x)

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_turb<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_turb<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_turb+mean_D_turb
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0

tiff("DOLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
DO_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR Dissolved oxygen")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("d")

plot(DO_LRR)
dev.off()


##################
# pH
##################

completepH<-submaster[complete.cases(submaster$pH), ] 
L<-length(unique(completepH$newdate))
dates<-unique(completepH$newdate)
LRR_Trmt_pH<-matrix(NA,L*12,4)
LRR_Control_pH<-matrix(NA,L*4,4)

count=1
for (i in 1:(L)){ #unique time points 
  for (j in 1:12){ #non-oligotrophic and non-control ponds=12
    
    subdat<-subset(completepH,newdate==dates[i])
    controls<-subset(subdat,Treatment=="C")
    others<-filter(subdat,Treatment!="C"&Treatment!="O")
    mean_C_pH<- mean(controls$pH)
    
    LRR_Trmt_pH[count,1]<-as.character(others$Pond[j])
    LRR_Trmt_pH[count,2]<-as.character(dates[i])
    LRR_Trmt_pH[count,3]<-as.character(others$Treatment[j])
    LRR_Trmt_pH[count,4]<-log(others$pH[j]/mean_C_pH)
    count=count+1
  }
}

LRR_Trmt_pH<-data.frame(LRR_Trmt_pH)
colnames(LRR_Trmt_pH)<-c("Pond","newdate","Treatment","LRR")
LRR_Trmt_pH$Treatment<-factor(LRR_Trmt_pH$Treatment)
LRR_Trmt_pH$newdate<-as.Date(LRR_Trmt_pH$newdate)
LRR_Trmt_pH$LRR<-as.numeric(as.character(LRR_Trmt_pH$LRR))


x <- summarySE(LRR_Trmt_pH, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
x <-x[complete.cases(x$sd), ] 
dim(x)

Additive = matrix(NA, L, 7)
for (j in 1:L){
  #j=5
  subdat<-subset(x,newdate==dates[j])
  mean_M_pH<- subdat$LRR[subdat$Treatment=="M"]
  mean_D_pH<- subdat$LRR[subdat$Treatment=="D"]
  Additive[j,1] = as.character("Add")
  Additive[j,2] = as.character(dates[j])
  Additive[j,3] = as.numeric(1)
  Additive[j,4] = mean_M_pH+mean_D_pH
}


Additive<-data.frame(Additive)
colnames(Additive) = colnames(x)
Additive$Treatment<-factor(Additive$Treatment)
Additive$newdate<-as.Date(Additive$newdate)
Additive$LRR<-as.numeric(as.character(Additive$LRR))
Additive$se=as.numeric(as.character(Additive$se))
Additive$ci=as.numeric(as.character(Additive$ci))
x=rbind(x,Additive)
x$ci[x$Treatment!="MD"]=0


tiff("pHLRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
pH_LRR<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-ci, ymax=LRR+ci), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey')+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='grey')+
  geom_hline (yintercept = 0, col='dark grey')+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors_sub)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("LRR pH")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("e")

plot(pH_LRR)
dev.off()



tiff("EF LRR effects.tiff", width=30, height=12, pointsize=1/300, units='in', res=300)
multiplot(nit_LRR, DO_LRR, oP_LRR, pH_LRR, turb_LRR,  NULL,cols=3)
dev.off()


