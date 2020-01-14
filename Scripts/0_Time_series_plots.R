rm(list=ls())

library('ggplot2')
library('Rmisc')
library('devtools')
library('easyGgplot2')
library('bindrcpp')
library('tidyverse')
library('codyn')
library('sp')
library('permute')
library('labeling')
library('reshape')
library('labdsv')
library('digest')
library('moments')
library('digest')
library('vegan')
library('backports')
library('checkmate')
library('lme4')
library('htmlTable')
library('lmerTest')
library("ggfortify")
library('cowplot')
library('stringi')
library('stringr')
library('lubridate')
library('corrplot')
library('RColorBrewer')
library('scales')
library('gridExtra')
library('zoo')
library('philentropy')
library('factoextra')

theme_set(theme_cowplot())
myColors<-c('#f46d43','#0570b0', '#7fbc41','#dd3497','grey50') #agreed upon colors

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

pd <- position_dodge(width = 1.5)


##########################
#Algal response
##########################
#temperature

measurevar= "Temperature"
groupvars=(c("Treatment","newdate"))
conf.interval=.95

completetemp<-submaster[complete.cases(submaster$Temperature), ] 

x <- summarySE(completetemp, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)


tiff("Temperature.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
Temp<-ggplot(x, aes(x=newdate, y=Temperature, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Temperature-se, ymax=Temperature+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("Temperature")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) 
  #ggtitle("a")

plot(Temp) 
dev.off()


tiff("Chlorophyll.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
chl<-ggplot(x, aes(x=newdate, y=CHL, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=CHL-se, ymax=CHL+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("Chlorophyll-a (ug/L)")+
  
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
  
plot(chl) 
dev.off()

################################# Investigating Rep 2C.
MDdat = subset(completechl,Treatment=="MD")

chlMD<-ggplot(MDdat, aes(x=newdate, y=CHL, group=Pond, colour=Pond)) + 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+

  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("Chlorophyll-a (ug/L)")+
  
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

plot(chlMD)


#####################################
### Plot the Cyanobacteria dominance
#####################################
rm(completesp)
#only taking the portion of the data where we have species composition estimates
completesp<-submaster[complete.cases(submaster$c_Cyanobacteria), ] 

which(colnames(completesp)=="v_Cyanobacteria")
which(colnames(completesp)=="v_Green_algae")

which(colnames(completesp)=="v_Achnthaceae")
which(colnames(completesp)=="v_Unknown_13")

#columns with functional group volume
gpvol <- data.frame(submaster[,97:102])
#columns with species volume
spvol <- data.frame(submaster[,103:167])

#hellinger transform of functional group biovolume
gpcomp <- decostand(gpvol, method = "hellinger",na.rm=T)
colnames(gpcomp)<-paste("h",colnames(gpcomp),sep="_")
completesp<-cbind(completesp,gpcomp)

#hellinger transform of species biovolume
spcomp <- decostand(spvol, method = "hellinger",na.rm=T)
colnames(spcomp)<-paste("h",colnames(spcomp),sep="_")
completesp<-cbind(completesp,spcomp)

#log transfrom of functional group biovolume
gpcomp <- decostand(gpvol, method = "log",na.rm=T)
colnames(gpcomp)<-paste("l",colnames(gpcomp),sep="_")
completesp<-cbind(completesp,gpcomp)

colnames(completesp)<- str_replace(colnames(completesp)," ","_")
head(completesp)

#log transform of species biovolume
spcomp <- decostand(spvol, method = "log",na.rm=T)
colnames(spcomp)<-paste("l",colnames(spcomp),sep="_")
completesp<-cbind(completesp,spcomp)

which( colnames(completesp)=="l_v_Cyanobacteria" )
which( colnames(completesp)=="l_v_Achnthaceae" )


###plot Cyano dominance

measurevar= "h_v_Cyanobacteria"

x <- summarySE(completesp, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)

extramatrix = data.frame(matrix(NA,5,7))
extramatrix[,1] = as.character(c("C","D","M","MD","O"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)

extramatrix2 = data.frame(matrix(NA,5,7))
extramatrix2[,1] = as.character(c("C","D","M","MD","O"))
extramatrix2[,2] = ymd(c("2016-07-04","2016-07-04","2016-07-04","2016-07-04","2016-07-04"))
xnames2<-names(x)
names(extramatrix2)<-xnames2
x<-rbind(x,extramatrix2)

tiff("Cyano_hel.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
cyano<-ggplot(x, aes(x=newdate, y=h_v_Cyanobacteria, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=h_v_Cyanobacteria-se, ymax=h_v_Cyanobacteria+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey', size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Cyanobacterial dominance")+
  
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
plot(cyano)
dev.off()


MDcyan = subset(completesp,Treatment=="MD")

MDcyano<-ggplot(MDcyan, aes(x=newdate, y=h_v_Cyanobacteria, group=Pond, colour=Pond)) + 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey', size = 1.2)+
  geom_point(size=3, position=pd)+
  
  ylab("Cyanobacterial dominance")+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9))+
#        legend.position="none")+
  
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("b")
plot(MDcyano)

tiff("MD cyan reps.tiff", width=10, height=6, pointsize=1/300, units='in', res = 300)
plot(MDcyano)
dev.off()

tiff("MD reps.tiff", width=20, height=6, pointsize=1/300, units='in', res = 300)
multiplot(chlMD, MDcyano, cols=2)
dev.off()

###########################
## Supporting Information
###########################

#h_v_Synechococcales

measurevar= "h_v_Synechococcales"
x <- summarySE(completesp, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)
tiff("Synechococcus.tiff", width = 16, height = 6, pointsize = 1/300, units = 'in', res = 300)
Syn<-ggplot(x, aes(x=newdate, y=h_v_Synechococcales, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=h_v_Synechococcales-se, ymax=h_v_Synechococcales+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Dominance of Synechococcus")+
  
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

plot(Syn)
dev.off()

#h_v_Lagerheimia

measurevar= "h_v_Lagerheimia"
x <- summarySE(completesp, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)
tiff("Lagerheimia.tiff", width = 16, height = 6, pointsize = 1/300, units = 'in', res = 300)
Lag<-ggplot(x, aes(x=newdate, y=h_v_Lagerheimia, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=h_v_Lagerheimia-se, ymax=h_v_Lagerheimia+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Dominance of Lagerheimia")+
  
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

plot(Lag)
dev.off()


#Trait evenness

measurevar= "TED"
x <- summarySE(completechl, measurevar, groupvars, na.rm=TRUE,
               conf.interval, .drop=TRUE)
extramatrix = data.frame(matrix(NA,5,7))
extramatrix[,1] = as.character(c("C","D","M","MD","O"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)

tiff("TED.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
ted<-ggplot(x, aes(x=newdate, y=TED, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=TED-se, ymax=TED+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("TED")+
  
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

  plot(ted)
dev.off()



MDdat = subset(completechl,Treatment == "MD")

tedMD<-ggplot(MDdat, aes(x=newdate, y=TED, group=Pond, colour=Pond)) + 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  ylab("TED")+
  
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

plot(tedMD)

##############
#Stoichiometry
##############

completeCN<-submaster[complete.cases(submaster$CtoN), ] 

measurevar= "CtoN"
x <- summarySE(completeCN, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)
extramatrix = data.frame(matrix(NA,5,7))
extramatrix[,1] = as.character(c("C","D","M","MD","O"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)


tiff("CtoN.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
CtoN<-ggplot(x, aes(x=c(newdate), y=c(CtoN), group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=CtoN-se, ymax=CtoN+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("C:N molar ratio")+
  
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
plot(CtoN)
dev.off()

#CtoP
completeCP<-submaster[complete.cases(submaster$CtoP), ] 
measurevar= "CtoP"
x <- summarySE(completeCP, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

extramatrix = data.frame(matrix(NA,5,7))
extramatrix[,1] = as.character(c("C","D","M","MD","O"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)


tiff("CtoP.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
CtoP<-ggplot(x, aes(x=newdate, y=CtoP, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=CtoP-se, ymax=CtoP+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("C:P molar ratio")+
  
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

plot(CtoP)
dev.off()


##########################
#Ecosystem and Food Web Effects
##########################

#Nitrate

measurevar= "Nitrate"
x <- summarySE(completechl, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)
extramatrix = data.frame(matrix(NA,5,7))
extramatrix[,1] = as.character(c("C","D","M","MD","O"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)

extramatrix2 = data.frame(matrix(NA,5,7))
extramatrix2[,1] = as.character(c("C","D","M","MD","O"))
extramatrix2[,2] = ymd(c("2016-07-04","2016-07-04","2016-07-04","2016-07-04","2016-07-04"))
xnames2<-names(x)
names(extramatrix2)<-xnames2
x<-rbind(x,extramatrix2)

tiff("Nitrate.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
Nitrate<-ggplot(x, aes(x=newdate, y=Nitrate, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Nitrate-se, ymax=Nitrate+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Nitrate-N (mg/L)")+
  
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

plot(Nitrate)
dev.off()

#orthoP

measurevar= "orthoP"
x <- summarySE(completechl, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

extramatrix = data.frame(matrix(NA,5,7))
extramatrix[,1] = as.character(c("C","D","M","MD","O"))
extramatrix[,2] = ymd(c("2017-02-27","2017-02-27","2017-02-27","2017-02-27","2017-02-27"))
xnames<-names(x)
names(extramatrix)<-xnames
x<-rbind(x,extramatrix)

extramatrix2 = data.frame(matrix(NA,5,7))
extramatrix2[,1] = as.character(c("C","D","M","MD","O"))
extramatrix2[,2] = ymd(c("2016-07-04","2016-07-04","2016-07-04","2016-07-04","2016-07-04"))
xnames2<-names(x)
names(extramatrix2)<-xnames2
x<-rbind(x,extramatrix2)

tiff("orthoP.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
phos<-ggplot(x, aes(x=newdate, y=orthoP, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=orthoP-se, ymax=orthoP+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Phosphate-P (ug/L)")+
  
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

plot(phos)
dev.off()

#########################
#Ecosystem-level effects
#########################

#Turbidity
head(submaster)
completeTurb<-submaster[complete.cases(submaster$Turbidity), ] 
measurevar= "Turbidity"
x <- summarySE(completeTurb, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

extramatrix2 = data.frame(matrix(NA,5,7))
extramatrix2[,1] = as.character(c("C","D","M","MD","O"))
extramatrix2[,2] = ymd(c("2016-07-04","2016-07-04","2016-07-04","2016-07-04","2016-07-04"))
xnames2<-names(x)
names(extramatrix2)<-xnames2
x<-rbind(x,extramatrix2)

tiff("Turbidity.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
Turbid<-ggplot(x, aes(x=newdate, y=Turbidity, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Turbidity-se, ymax=Turbidity+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Turbidity (FNU)")+
  
  
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

plot(Turbid)
dev.off()

#dO
head(submaster)
completeO2<-submaster[complete.cases(submaster$O2), ] 
measurevar= "O2"
x <- summarySE(completeO2, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

extramatrix2 = data.frame(matrix(NA,5,7))
extramatrix2[,1] = as.character(c("C","D","M","MD","O"))
extramatrix2[,2] = ymd(c("2016-07-04","2016-07-04","2016-07-04","2016-07-04","2016-07-04"))
xnames2<-names(x)
names(extramatrix2)<-xnames2
x<-rbind(x,extramatrix2)

tiff("O2.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
O2<-ggplot(x, aes(x=newdate, y=O2, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=O2-se, ymax=O2+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Dissolved oxygen (%)")+
  
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

plot(O2)
dev.off()


#pH
head(submaster)
completepH<-submaster[complete.cases(submaster$pH), ] 
measurevar= "pH"
x <- summarySE(completepH, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)


tiff("pH.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
pH<-ggplot(x, aes(x=newdate, y=pH, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=pH-se, ymax=pH+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("pH")+
  
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
plot(pH)
dev.off()


#Daphnia biomass
completeDaph<-submaster[complete.cases(submaster$Daphnia_adult), ] 
measurevar= "Daphnia_adult"
x <- summarySE(completeDaph, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

tiff("Daphnia.tiff", width = 16, height = 6, pointsize = 1/300, units = 'in', res = 300)
Daphnia<-ggplot(x, aes(x=newdate, y=Daphnia_adult, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Daphnia_adult-se, ymax=Daphnia_adult+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+
  ylab("Adult daphnia (per 10L)")+
  
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

plot(Daphnia)
dev.off()

#Copepods
head(submaster)
completeCope<-submaster[complete.cases(submaster$Copepod_adult), ] 
measurevar= "Copepod_adult"
x <- summarySE(completeCope, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

tiff("Copepods.tiff", width = 16, height = 6, pointsize = 1/300, units = 'in', res = 300)
Copepods<-ggplot(x, aes(x=newdate, y=Copepod_adult, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=Copepod_adult-se, ymax=Copepod_adult+se), width=1, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='dark grey',size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-10")), col='dark grey',size = 1.2)+
  geom_point(size=3, position=pd)+
  
  scale_color_manual(values=myColors)+

  ylab("Adult copepods (per 10L)")+
  
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

plot(Copepods)
dev.off()



###########
## Figures
###########


tiff("Algal response.tiff", width=30, height=12, pointsize=1/300, units='in', res = 300)
multiplot(chl, CtoN, cyano, CtoP, ted, cols=3)
dev.off()

tiff("EF effects.tiff", width=30, height=12, pointsize=1/300, units='in', res=300)
multiplot(Nitrate, O2, phos, pH, Turbid,  NULL,cols=3)
dev.off()

tiff("Grazers.tiff", width=20, height=6, pointsize=1/300, units='in', res=300)
multiplot(Daphnia, Copepods, cols=2)
dev.off()

tiff("Phyto.tiff", width=20, height=6, pointsize=1/300, units='in', res=300)
multiplot(Syn, Lag, cols=2)
dev.off()
