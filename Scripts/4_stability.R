library('lubridate')
library('RColorBrewer')
library('Rmisc')
library('scales')
library('multcomp')
library('car')
library('MASS')
library('vegan')
library('tidyverse')

myColors<-c('#f46d43','#0570b0', '#7fbc41','#dd3497','grey50') #colors
myColorsmO<-c('#f46d43','#0570b0', '#7fbc41','#dd3497') #minus oligos

master<-read.csv("Narwani_et_al_2019_Proceedings_B_data.csv")

#factor
master[,"Treatment"]<-factor(master[,"Treatment"])
master[,"Dreisse"]<-factor(master[,"Dreisse"])
master[,"Macrophytes"]<-factor(master[,"Macrophytes"])
master[,"Nutrients"]<-factor(master[,"Nutrients"])
master[,"Block"]<-factor(master[,"Block"])
master[,"Rep"]<-factor(master[,"Rep"])
master$newdate<-mdy(as.character(master$Date))


submaster = master[master$Week<35,] #should be week 35
submaster[,"Week"]<-factor(submaster[,"Week"])
submaster = submaster %>% filter(submaster$newdate > "2016-08-02")
totdates = length(unique(submaster$newdate))

basedates<-c("2016-08-08","2016-08-22","2016-09-05","2016-09-19","2016-10-10")
nutdats<-c("2016-08-12", "2016-08-26", "2016-09-09", "2016-09-22", "2016-10-10")
nutsadd<-c(10,20,30,40,50)


########################
#### Cumulative distance
########################
timexpond = 16*totdates
dates = unique(submaster$newdate)

oligos = subset(submaster,Treatment=="O")
oligochl <-matrix(NA, totdates, 3)

for(a in 1:totdates){
  olchl=subset(oligos, newdate==dates[a])
  oligochl[a,1] = as.character(dates[a])
  oligochl[a,2] = as.numeric(olchl$Week[1])
  oligochl[a,3] = mean(olchl$CHL,na.rm=T)
}

oligochl = data.frame(oligochl)
colnames(oligochl)=c("newdate_O","Week_O","CHL_O")
oligochl$newdate_O=as.Date(oligochl$newdate_O)

notoligos = subset(submaster, Treatment!="O")
ponds = unique(notoligos$Pond)

Diff <-matrix(NA, timexpond, 7)

count=1
for (p in 1:16) {
  #p=1
  for (i in 1:1){
    #i=1
    O_CHL = as.numeric(paste(subset(oligochl, newdate_O==dates[i])$CHL_O))
    pond = as.character(ponds[p])
    subdat = notoligos[which(notoligos$newdate == dates[i] & notoligos$Pond == pond), ]
    Diff[count,1] = as.character(dates[i])
    Diff[count,2] = pond
    Diff[count,3] = as.numeric(subdat$Week)
    Diff[count,4] = as.character(subdat$Treatment)
    Diff[count,5] = as.numeric(subdat$CHL- O_CHL)
    Diff[count,6] = as.numeric(subdat$Rep)
    Diff[count,7] = "NA"
    count = count + 1 
  }
  for (i in 2:totdates){
    #i=1
    O_CHL = as.numeric(paste(subset(oligochl, newdate_O==dates[i])$CHL_O))
    pond = as.character(ponds[p])
    subdat = notoligos[which(notoligos$newdate == dates[i] & notoligos$Pond == pond), ]
    Diff[count,1] = as.character(dates[i])
    Diff[count,2] = pond
    Diff[count,3] = as.numeric(subdat$Week)
    Diff[count,4] = as.character(subdat$Treatment)
    Diff[count,5] = as.numeric(subdat$CHL- O_CHL)
    Diff[count,6] = as.numeric(subdat$Rep)
    Diff[count,7] = as.numeric(Diff[count-1,5])
    count = count + 1 
    }
  }

Diff = data.frame(Diff)
colnames(Diff) =c("newdate","Pond","Week","Treatment","CHL_diff", "Rep","CHL_diff_1")
Diff$newdate = as.Date(Diff$newdate)
Diff$Treatment = factor(Diff$Treatment)
Diff$CHL_diff = as.numeric(paste(Diff$CHL_diff))
Diff$CHL_diff_1 = as.numeric(paste(Diff$CHL_diff_1))
Diff$CHL_diff = (Diff$CHL_diff + 5)
Diff$CHL_diff_1 = (Diff$CHL_diff_1 + 5)

#############################################
######### Autocorrelation on raw CHL
#############################################

rm(subshort)
subshort = submaster %>% filter(submaster$newdate > "2016-08-02", submaster$newdate < "2016-11-24")
ps = unique(subshort$Pond)
lp = length(unique(subshort$Pond))
CHL_0 = rep("NA",length(subshort$CHL))
CHL_1 = rep("NA",length(subshort$CHL))
newCHL = cbind(CHL_0,CHL_1)
subdate = length(unique(subshort$newdate))

count = 1
for (p in 1:lp){
  #p=1  
  pond = ps[p]
  dat = subshort %>% filter(subshort$Pond == pond)
  for (j in 1:1){
  newCHL[count,1] = dat$CHL[j]
  count=count+1
  }
  for (j in 2:subdate){
  newCHL[count,1] = dat$CHL[j]
  newCHL[count,2] = newCHL[count-1,1]
  count=count+1
  }
}

newCHL=data.frame(newCHL)
colnames(newCHL) = c("CHL_0","CHL_1")
newCHL$CHL_0 = as.numeric(paste(newCHL$CHL_0 ))
newCHL$CHL_1 = as.numeric(paste(newCHL$CHL_1 ))
subshort = subshort%>% arrange(subshort$Pond,subshort$newdate)
subshort = cbind(newCHL,subshort)

subshort = subshort %>% filter(Treatment != "O")
ponds = unique(subshort$Pond)
head(subshort)
tail(subshort)

ps = unique(subshort$Pond)

AC<-matrix(NA, 16, 4)
for (p in 1:16) {
  #p=1
  pond = ps[p]
  dat = subshort %>% filter(Pond == pond)
  AC[p,1] = as.character(paste(pond))  #Pond 
  AC[p,2] = as.character(paste(dat$Treatment[1])) #Treatment
  AC[p,3] = as.numeric(paste(dat$Rep[1])) #Rep
  z = cor.test(dat$CHL_0,dat$CHL_1,method="spearman",na.rm=T)
  AC[p,4] = as.numeric(paste(z$estimate))
  }

head(AC)
AC = data.frame(AC)
colnames(AC) =c("Pond","Treatment","Rep", "ac")
AC$Treatment = factor(AC$Treatment)
AC$Pond = factor(AC$Pond)
AC$ac = as.numeric(paste(AC$ac))

pd <- position_dodge(width = 0.3)
groupvars=c("Treatment")
conf.interval=.95

measurevar="ac"
x <- summarySE(AC, measurevar, groupvars, na.rm=TRUE, conf.interval, .drop=TRUE)



tiff("Temporal autocorrelation.tiff", width=16, height=8, pointsize=1/300, units='in', res = 300)
e<- ggplot(x, aes(x=Treatment, y=ac, colour=Treatment)) +
  scale_color_manual(values=myColorsmO)+
  geom_errorbar(aes(ymin=ac-se, ymax=ac+se), width=0.3, position=pd)+
  geom_point(size=5, position=pd)+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  ylab("Autocorrelation (lag=1)")+
  xlab("Treatment")+
  theme_classic()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.text.x=element_text(size=22),
        axis.title.x = element_text(size=26, margin = margin(t = 15)),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  ggtitle("e")
plot(e)
dev.off()

########################

head(Diff)
ps = unique(Diff$Pond)

AC<-matrix(NA, 16, 4)
for (p in 1:16) {
  #p=1
  pond = ps[p]
  dat = Diff %>% filter(Pond == pond)
  AC[p,1] = as.character(paste(pond))  #Pond 
  AC[p,2] = as.character(paste(dat$Treatment[1])) #Treatment
  AC[p,3] = as.numeric(paste(dat$Rep[1])) #Rep
  z = cor.test(dat$CHL_diff,dat$CHL_diff_1,method="spearman",na.rm=T)
  AC[p,4] = as.numeric(paste(z$estimate))
}

head(AC)
AC = data.frame(AC)
colnames(AC) =c("Pond","Treatment","Rep", "ac")
AC$Treatment = factor(AC$Treatment)
AC$Pond = factor(AC$Pond)
AC$ac = as.numeric(paste(AC$ac))


pd <- position_dodge(width = 0.3)
groupvars=c("Treatment")
conf.interval=.95

measurevar="ac"
x <- summarySE(AC, measurevar, groupvars, na.rm=TRUE, conf.interval, .drop=TRUE)


tiff("Temporal autocorrelation diff.tiff", width=16, height=8, pointsize=1/300, units='in', res = 300)
e<- ggplot(x, aes(x=Treatment, y=ac, colour=Treatment)) +
  scale_color_manual(values=myColorsmO)+
  geom_errorbar(aes(ymin=ac-se, ymax=ac+se), width=0.3, position=pd)+
  geom_point(size=5, position=pd)+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  ylab("Autocorrelation (lag=1)")+
  xlab("Treatment")+
  theme_classic()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.text.x=element_text(size=22),
        axis.title.x = element_text(size=26, margin = margin(t = 15)),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  ggtitle("e")
plot(e)
dev.off()

###########################
###########################

rm(Cumulative)
Cumulative <-matrix(NA,timexpond, 17)

count=1
for (p in 1:16) {
  #p=1
  pond = as.character(ponds[p])
  subdat = notoligos %>% filter(Pond == pond)
  base = subdat %>% filter(newdate == dates[1])
  
  for(i in 1:1){
    #i=1
    subdat1 = filter(subdat, newdate == dates[i])
    Cumulative[count,1] = as.character(paste(dates[i]))
    Cumulative[count,2] = pond
    Cumulative[count,3] = paste(subdat1$Week)
    Cumulative[count,4] = as.character(subdat1$Treatment)
    Cumulative[count,5] = base$CHL #what was the baseline
    Cumulative[count,6] = 0 #weeks from last point
    Cumulative[count,7] = subdat1$CHL - base$CHL #change of the time point relative to the time before
    Cumulative[count,8] = 0
    Cumulative[count,9] = Cumulative[count,8]
    Cumulative[count,10] = subdat1$CHL #difference from oligos
    Cumulative[count,11] = 0
    Cumulative[count,12] = 0
    Cumulative[count,13] = 0
    Cumulative[count,14] = subdat1$Rep
    Cumulative[count,15] = "NA"
    Cumulative[count,16] = "NA"
    Cumulative[count,17] = "NA"
    count = count + 1 
  }
  7
  for (i in 2:totdates){
    #i=15
    subdat1 = subdat %>% filter(newdate == dates[i])
    
    if(is.na(subdat1$CHL)==FALSE){
      subdat2 = subdat %>% filter(newdate == dates[i-1])
      
      if(is.na(subdat2$CHL)==TRUE){
        subdat2.na = subdat %>% filter(newdate == dates[i-2])
        Cumulative[count,1] = as.character(paste(dates[i]))
        Cumulative[count,2] = pond
        Cumulative[count,3] = paste(subdat1$Week)
        Cumulative[count,4] = as.character(subdat1$Treatment)
        Cumulative[count,5] = base$CHL #what was the baseline
        weeks = as.numeric(paste(subdat1$Week)) - as.numeric(paste(subdat2.na$Week)) #weeks from last point
        Cumulative[count,6] = weeks
        Cumulative[count,7] = subdat1$CHL - subdat2.na$CHL #change of the time point relative to the time before
        Cumulative[count,8] = (log(subdat1$CHL)-log(subdat2.na$CHL))/weeks #exponential
        Cumulative[count,9] = ((subdat1$CHL - subdat2.na$CHL)/2)*weeks + subdat2.na$CHL*weeks #integral of CHL over the interval
        Cumulative[count,10] = sum(as.numeric(paste(Cumulative[count,9])),as.numeric(paste(Cumulative[(count-2),10])),na.rm=T) #cumulative CHL
        Cumulative[count,11] =(subdat1$CHL - subdat2.na$CHL)/weeks #dN/dt
        Cumulative[count,12] = (subdat1$CHL-subdat2.na$CHL)/subdat2.na$CHL/weeks #dN/dt/N
        Cumulative[count,13] = log(subdat1$CHL/subdat2.na$CHL) #LRR
        Cumulative[count,14] = as.numeric(paste(subdat1$Rep))
        Cumulative[count,15] = Cumulative[count-1,8] #exponential from the time before
        Cumulative[count,16] = Cumulative[count-1,11] #dNdt from the time before
        Cumulative[count,17] = Cumulative[count-1,13] #LRR from the time before
      }
      
      else{
        Cumulative[count,1] = as.character(paste(dates[i]))
        Cumulative[count,2] = pond
        Cumulative[count,3] = paste(subdat1$Week)
        Cumulative[count,4] = as.character(subdat1$Treatment)
        Cumulative[count,5] = base$CHL #what was the baseline
        weeks = as.numeric(paste(subdat1$Week)) - as.numeric(paste(subdat2$Week)) #weeks from last point
        Cumulative[count,6] = weeks
        Cumulative[count,7] = subdat1$CHL - subdat2$CHL #change of the time point relative to the time before
        Cumulative[count,8] = (log(subdat1$CHL)-log(subdat2$CHL))/weeks #exponential
        Cumulative[count,9] = ((subdat1$CHL - subdat2$CHL)/2)*weeks + subdat2$CHL*weeks #integral of difference from oligos over the interval
        Cumulative[count,10] = sum(as.numeric(paste(Cumulative[count,9])),as.numeric(paste(Cumulative[(count-1),10])),na.rm=T) #cumulative difference from oligos
        Cumulative[count,11] =(subdat1$CHL - subdat2$CHL)/weeks #dN/dt
        Cumulative[count,12] = (subdat1$CHL-subdat2$CHL)/subdat2$CHL/weeks #dN/dt/N
        Cumulative[count,13] = log(subdat1$CHL/subdat2$CHL) #LRR
        Cumulative[count,14] = subdat1$Rep
        Cumulative[count,15] = Cumulative[count-1,8] #exponential from the time before
        Cumulative[count,16] = Cumulative[count-1,11] #dNdt from the time before
        Cumulative[count,17] = Cumulative[count-1,13] #LRR from the time before
      }}
    
    else{
      Cumulative[count,1] = as.character(paste(dates[i]))
      Cumulative[count,2] = pond
      Cumulative[count,3] = paste(subdat1$Week)
      Cumulative[count,4] = as.character(subdat1$Treatment)
      Cumulative[count,5] = base$CHL #what was the baseline
      weeks = as.numeric(paste(subdat1$Week)) - as.numeric(paste(subdat2$Week)) #weeks from last point
      Cumulative[count,6] = weeks
    }
    count = count + 1
  }
}


Cumulative = data.frame(Cumulative)
colnames(Cumulative) =c("newdate","Pond","Week","Treatment","base","weeks",
                        "temp_change","exponential","integ_interv","cumulative",
                        "dNdt","dNdtN","LRR","Rep","exp_1","dNdt_1","LRR_1")
Cumulative$newdate = as.Date(Cumulative$newdate)
Cumulative$Treatment = factor(Cumulative$Treatment)
Cumulative$weeks = as.numeric(paste(Cumulative$weeks))
Cumulative$temp_change = as.numeric(paste(Cumulative$temp_change))
Cumulative$integ_interv = as.numeric(paste(Cumulative$integ_interv))
Cumulative$cumulative = as.numeric(paste(Cumulative$cumulative))
Cumulative$exponential = as.numeric(paste(Cumulative$exponential))
Cumulative$dNdt = as.numeric(paste(Cumulative$dNdt))
Cumulative$dNdtN = as.numeric(paste(Cumulative$dNdtN))
Cumulative$LRR = as.numeric(paste(Cumulative$LRR))
Cumulative$exp_1 = as.numeric(paste(Cumulative$exp_1))
Cumulative$dNdt_1 = as.numeric(paste(Cumulative$dNdt_1))
Cumulative$LRR_1 = as.numeric(paste(Cumulative$LRR_1))


head(Cumulative)
tail(Cumulative)
Cumulative %>% filter (Pond == "1A")

write.csv(Cumulative,"Cumulative.csv")

#LRR and the exponential are the same (except exponential is divided by the number of weeks).
#dNdtN and LRR are exactly but non-linearly related.

plot(Cumulative$dNdtN,Cumulative$exponential)
plot(Cumulative$dNdtN,Cumulative$LRR)
plot(Cumulative$exponential,Cumulative$LRR)
plot(Cumulative$LRR,log(Cumulative$dNdt))




##########################################################################################
#plot histogram of various stability metric by rep and treatment
##########################################################################################

qplot(LRR, binwidth = 0.1,data=Cumulative, colour=Treatment) +
  facet_wrap(~Treatment, nrow = 4) + 
  scale_color_manual(values = myColorsmO) +
  theme_classic() +
  theme(legend.position = "none")  


################################################
### CHL_Diff
################################################
diff = Diff[complete.cases(Diff$CHL_diff),]
#diff = diff %>% filter(diff$newdate > "2016-08-02",diff$newdate < "2016-11-24")

conf.interval=.95
pd = position_dodge(width = 1.5)
measurevar= "CHL_diff"
groupvars=(c("Treatment","newdate"))

x <- summarySE(diff, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

tiff("Chlorophyll_diff.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
a<-ggplot(x, aes(x=newdate, y=CHL_diff, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=CHL_diff-se, ymax=CHL_diff+se), width=2, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey', size = 1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-11")), col='grey', size = 1.2)+
  geom_point(size=4, position=pd)+
  
  scale_color_manual(values=myColors)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("O-Std. Chlorophyll-a (ug/L)")+
  theme_classic()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+

  geom_segment(aes(x=as.Date("2016-08-08"),xend=as.Date("2016-11-14"),y=62,yend=62),linetype="dashed",color="grey", size = 1.2)+
  geom_hline(yintercept=0, color = "grey")+
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("a")

plot(a)
dev.off()



################################################
### Integral
################################################
Cumul= Cumulative[complete.cases(Cumulative$cumulative),]
#cum = Cumul %>% filter(Cumul$newdate > "2016-08-02")
#cum = cum %>% filter(Cumul$newdate < "2016-11-24")

conf.interval=.95
pd = position_dodge(width = 1.5)
measurevar= "cumulative"
groupvars=(c("Treatment","newdate"))

x <- summarySE(Cumul, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

tiff("Cumulative_diff.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
b<-ggplot(x, aes(x=newdate, y=cumulative, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=cumulative-se, ymax=cumulative+se), width=2, position=pd)+ 
  geom_line(position=pd) +
  
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-11")), col='grey', size=1.2)+
  geom_point(size=4, position=pd)+
  
  scale_color_manual(values=myColors)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("Cumulative Chlorophyll-a (ug/L)")+
  theme_classic()+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+

  
  geom_segment(aes(x=as.Date("2016-08-08"),xend=as.Date("2016-11-14"),y=475,yend=475),linetype="dashed",color="grey", size=1.2)+
  geom_hline(yintercept=0, color = "grey")+
  scale_x_date(breaks = date_breaks("months"), 
               labels = date_format("%m/%y")) +
  ggtitle("b")

plot(b)
dev.off()


ggplot(dat=cum, aes(x=newdate, y=log(cumulative), group=Treatment, colour=Treatment)) +
  geom_point(shape=1) +
  scale_color_manual(values=myColorsmO) +
  theme_classic()



c = subset(submaster, submaster$Treatment=="C") 
m = subset(submaster, submaster$Treatment=="M")
d = subset(submaster, submaster$Treatment=="D")
md = subset(submaster, submaster$Treatment=="MD")
plot(ecdf(c$CHL),xlim = c(0,80), col = myColorsmO[1])
plot(ecdf(m$CHL),xlim = c(0,80), col = myColorsmO[2])
plot(ecdf(d$CHL),xlim = c(0,80), col = myColorsmO[3])
plot(ecdf(md$CHL),xlim = c(0,80), col = myColorsmO[4])

ks.test(m$CHL,md$CHL)
ks.test(d$CHL,md$CHL)
ks.test(md$CHL,c$CHL)
ks.test(c$CHL,d$CHL)
ks.test(c$CHL,m$CHL)


################################################
### LRR over time
################################################
Cumul = Cumulative %>% filter(Cumulative$newdate > "2016-08-02", Cumulative$newdate< "2016-11-24")
Cumul = Cumul %>% arrange(Pond,newdate)
Cumul$Padd = rep(factor(c(0,10,10,20,20,30,30,40,40,40,50,50,50,50)),16)
lrr= Cumul[complete.cases(Cumul$LRR),]

head(lrr)
tail(lrr)

conf.interval=.95
pd = position_dodge(width = 1.5)
measurevar= "LRR"
groupvars=(c("Treatment","newdate"))

x <- summarySE(lrr, measurevar, groupvars, na.rm=FALSE,
               conf.interval, .drop=TRUE)

tiff("LRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
c<-ggplot(x, aes(x=newdate, y=LRR, group=Treatment, colour=Treatment)) + 
  geom_errorbar(aes(ymin=LRR-se, ymax=LRR+se), width=3, position=pd)+ 
  geom_line(position=pd) +
  geom_hline(yintercept=0, color = "grey", size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-12")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-08-26")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-09")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-09-22")), col='grey', size=1.2)+
  geom_vline (xintercept = as.numeric(as.Date("2016-10-11")), col='grey', size=1.2)+
  geom_point(size=5, position=pd)+
  
  scale_color_manual(values=myColors)+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  
  ylab("Log Response Ratio")+
  theme_classic()+
  
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=22,angle=45,hjust=1),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+

  scale_x_date(breaks = date_breaks("weeks"), 
               labels = date_format("%d/%m")) +
  ggtitle("c")

plot(c)
dev.off()



nuts = unique(lrr$Padd)
nutlevs = length(nuts)
ponds = unique(lrr$Pond)

Change<-matrix(NA,80, 10)

count=1
for (p in 1:16) {
  #p=1
  pond = ponds[p]
  subdat = lrr %>% filter(Pond==pond)
  for (i in 2:nutlevs){
    #i=2
    nut = nuts[i]
    subdat1 = subdat %>% filter (Padd==nut)
    Change[count,1] = as.character(paste(ponds[p]))  #Pond 
    Change[count,2] = as.character(paste(subdat1$Treatment[1])) #Treatment
    Change[count,3] = as.numeric(paste(subdat1$Rep[1])) #Rep
    Change[count,4] = as.numeric(paste(subdat1$Padd[1])) #Paddition
    Change[count,5] = as.numeric(paste(max(subdat1$dNdt,na.rm=T)))   #maxdNdt
    Change[count,6] = as.numeric(paste(min(subdat1$dNdt,na.rm=T))) #mindNdt
    Change[count,7] = as.numeric(paste(max(subdat1$exponential,na.rm=T))) #maxExp
    Change[count,8] = as.numeric(paste(min(subdat1$exponential,na.rm=T))) #minExp
    Change[count,9] = as.numeric(paste(max(subdat1$LRR,na.rm=T))) #maxLRR
    Change[count,10] = as.numeric(paste(min(subdat1$LRR,na.rm=T))) #minLRR
    count = count+1
  }
}

head(Change)
Change = data.frame(Change)
colnames(Change) =c("Pond","Treatment","Rep","Padd", "maxdNdt", "mindNdt", "maxExp", "minExp", "maxLRR", "minLRR")
Change$Treatment = factor(Change$Treatment)
Change$Pond = factor(Change$Pond)
Change$maxdNdt = as.numeric(paste(Change$maxdNdt))
Change$mindNdt = as.numeric(paste(Change$mindNdt))
Change$maxExp = as.numeric(paste(Change$maxExp))
Change$minExp = as.numeric(paste(Change$minExp))
Change$maxLRR = as.numeric(paste(Change$maxLRR))
Change$minLRR= as.numeric(paste(Change$minLRR))


Change$Padd = as.numeric(paste(Change$Padd))
modmax = lm(Change$maxLRR~Change$Treatment*Change$Padd)
anova(modmax)
modmax.1 = lm(Change$maxLRR~Change$Treatment)

Dunmodmax1<-multcomp::glht(modmax.1,Treatment=mcp(Group="Dunnett"))
summary(Dunmodmax1)

modmin = lm(Change$minLRR~Change$Treatment*Change$Padd)
anova(modmin)

acmod=lm (AC$ac~AC$Treatment)
anova(acmod)
n = length(AC$ac)


library(pwr)
acmod = aov(ac ~ Treatment, data = AC)
summary(acmod)
pwr.anova.test(k=4,n=n/4,f=2.7667,sig.level=0.05)
pwr.anova.test(k=4,f=2.7667,sig.level=0.05,power=0.9999)

pd <- position_dodge(width = 0.3)
groupvars=c("Treatment","Padd")
conf.interval=.95

measurevar="maxLRR"
x <- summarySE(Change, measurevar, groupvars, na.rm=TRUE, conf.interval, .drop=TRUE)

S.a<- ggplot(x, aes(x=Padd, y=maxLRR, colour=Treatment)) +
  scale_color_manual(values=myColorsmO)+
  geom_errorbar(aes(ymin=maxLRR-se, ymax=maxLRR+se), width=0.3, position=pd)+
  geom_point(size=5, position=pd)+
  ylim(-0.6,2.6) +
  ylab("Maximum LRR")+
  xlab("Nutrient Addition (ug/L of P)")+
  
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  theme_classic()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=16),
        axis.title.y = element_text(size=20, margin = margin(r = 15)),
        axis.text.x=element_text(size=16),
        axis.title.x = element_text(size=20, margin = margin(t = 15)),
        plot.title = element_text(size=24,hjust=0.9),
        legend.position="none")+
  ggtitle("a")
plot(S.a)


measurevar="minLRR"
x<- summarySE(Change, measurevar, groupvars, na.rm=TRUE, conf.interval, .drop=TRUE)

S.b<- ggplot(x, aes(x=Padd, y=minLRR, colour=Treatment)) +
  scale_color_manual(values=myColorsmO)+
  geom_errorbar(aes(ymin=minLRR-se, ymax=minLRR+se), width=0.3, position=pd)+
  geom_point(size=5, position=pd)+
  ylim(0.6,-2.6)+
  ylab("Minimum LRR")+
  xlab("Nutrient Addition (ug/L of P)")+
  
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  theme_classic()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=16),
        axis.title.y = element_text(size=20, margin = margin(r = 15)),
        axis.text.x=element_text(size=16),
        axis.title.x = element_text(size=20, margin = margin(t = 15)),
        plot.title = element_text(size=24,hjust=0.9),
        legend.position="none")+
  ggtitle("b")
plot(S.b)


tiff("Max_Min_LRR.tiff", width = 12, height = 5, pointsize = 1/300, units = 'in', res = 300)
multiplot(S.a,S.b,cols=2)
dev.off()










tiff("max_v_min_LRR.tiff", width = 16, height = 8, pointsize = 1/300, units = 'in', res = 300)
d<-ggplot(dat=Change, aes(x=maxLRR, y=minLRR, group=Treatment, colour=Treatment)) +#,shape=Padd
  geom_smooth(method=lm, alpha=0.2, aes(fill=Treatment))+
  scale_color_manual(values=myColorsmO) +
  scale_fill_manual(values=myColorsmO) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+
  geom_abline(intercept = 0, slope = -1, color="grey", linetype="dashed")+
  geom_point(shape=1,size=4) +
  xlab("Maximum LRR")+
  ylab("Minimum LRR")+
  theme_classic()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"),
        axis.text.y=element_text(size=24),
        axis.title.y = element_text(size=22, margin = margin(r = 15)),
        axis.text.x=element_text(size=22),
        axis.title.x = element_text(size=26, margin = margin(t = 15)),
        plot.title = element_text(size=30,hjust=0.9),
        legend.position="none")+
  ggtitle("d")
  
plot(d)
dev.off()

groupvars=(c("Treatment","Padd"))
measurevar= "maxLRR"
x <- summarySE(Change, measurevar, groupvars, na.rm=FALSE, conf.interval, .drop=TRUE)
measurevar= "minLRR"
y <- summarySE(Change, measurevar, groupvars, na.rm=FALSE, conf.interval, .drop=TRUE)
colnames(y)[colnames(y)=="se"] <- "y.se"
colnames(y)[colnames(y)=="sd"] <- "y.sd"
colnames(y)[colnames(y)=="ci"] <- "y.ci"
colnames(y)[colnames(y)=="minLRR"] <- "minLRR"
colnames(x)[colnames(x)=="se"] <- "x.se"
colnames(x)[colnames(y)=="sd"] <- "x.sd"
colnames(x)[colnames(x)=="ci"] <- "x.ci"
colnames(x)[colnames(x)=="maxLRR"] <- "maxLRR"
dat = cbind(x,y[,4:7])


ggplot(dat=dat, aes(x=maxLRR, y=minLRR, group=Treatment, colour=Treatment,shape=Padd)) +
  geom_point(size=3) +
  #geom_smooth(method=lm, alpha=0.2, aes(fill=Treatment))+
  geom_errorbar(aes(ymin = minLRR-y.se,ymax = minLRR+y.se)) + 
  geom_errorbarh(aes(xmin =  maxLRR-x.se,xmax = maxLRR+x.se)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+
  geom_abline(intercept = 0, slope = -1, color="grey", linetype="dashed")+
  scale_color_manual(values=myColorsmO) +
  scale_shape_manual(values =  c(16,15,17,6,8))+
  #scale_fill_manual(values=myColorsmO) +
  theme_classic()




tiff("Figure 4.tiff", width=30, height=12, pointsize=1/300, units='in', res = 300)
multiplot(a,d,b,e,c,cols = 3)
dev.off()



