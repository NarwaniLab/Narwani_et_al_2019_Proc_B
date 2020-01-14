rm(list=ls())
#set the working directory to wherever you have saved this data file.
getwd()
master<-read.csv("Narwani_et_al_2019_Proceedings_B_data.csv")

#factor
#install.packages("lubridate")
master[,"Treatment"]<-factor(master[,"Treatment"])
master[,"Dreisse"]<-factor(master[,"Dreisse"])
master[,"Macrophytes"]<-factor(master[,"Macrophytes"])
master[,"Nutrients"]<-factor(master[,"Nutrients"])
master[,"Block"]<-factor(master[,"Block"])
master$newdate<-lubridate::mdy(as.character(master$Date))
master$TP<-as.numeric(master$TP)
master$Week<-as.numeric(paste(master$Week))

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

#columns with species volume
spvol <- data.frame(submaster[,103:167])

#columns with functional group volume
gpvol <- data.frame(submaster[,97:102])

########
#hellinger transform of species biovolume

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


completechl<-submaster[complete.cases(submaster$CHL), ] 

beforedates<- as.Date(c("2016-07-25","2016-08-02","2016-08-08"))
duringdates<- as.Date(c("2016-08-15","2016-08-22","2016-08-29","2016-09-05","2016-09-12","2016-09-19","2016-09-26","2016-10-03","2016-10-10"))
afterdates<- as.Date(c("2016-10-17","2016-10-24","2016-10-31","2016-11-14","2016-11-24","2016-11-28","2016-12-19","2017-01-16","2017-02-13"))
afterdates2<- as.Date(c("2016-10-17","2016-10-24","2016-10-31","2016-11-14","2016-11-24","2016-11-28"))
afterdates3<- as.Date(c("2016-10-17","2016-10-24","2016-10-31"))

before = subset(submaster, newdate %in% beforedates)
during = subset(submaster, newdate %in% duringdates)
after = subset(submaster, newdate %in% afterdates)
after2 = subset(submaster, newdate %in% afterdates2)
after3 = subset(submaster, newdate %in% afterdates3)


## Chlorophyll a
#Before

completeCHLbef<-before[complete.cases(before$CHL), ] 
unique(completeCHLbef$newdate)

modBefCHL<-lm(CHL~Treatment,data=before)
summary(modBefCHL)

mixmodBefCHL<-nlme::lme(CHL~ Treatment,random=~1|Week/Pond,data=before,na.action=na.omit)# correlation=nlme::corAR1(form=~1|Pond), 
DunBefCHL<-multcomp::glht(mixmodBefCHL,Treatment=mcp(Group="Dunnett"))
summary(DunBefCHL)

Trmt = c("D","M","MD","O")
Var = rep("CHL",4) #replace "CHL" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefCHL)$test$coefficients[2:5]
Zval = summary(DunBefCHL)$test$tstat[2:5]
Pval = summary(DunBefCHL)$test$pvalues[2:5]

DunnettsCHL1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#During
completeCHLdur<-during[complete.cases(during$CHL), ] 
unique(completeCHLdur$newdate)

modDurCHL<-lm(CHL~Treatment,data=during)
DunDurCHL<-multcomp::glht(modDurCHL,Treatment=mcp(Group="Dunnett"))
summary(DunDurCHL)

mixmodDurCHL<-nlme::lme(CHL~ Treatment,random=~1|Week/Pond, data=during,na.action=na.omit)# correlation=nlme::corAR1(form=~Week|Pond)
DunDurCHL<-multcomp::glht(mixmodDurCHL,Treatment=mcp(Group="Dunnett"))
summary(DunDurCHL)

Trmt = c("D","M","MD","O")
Var = rep("CHL",4) #replace "CHL" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurCHL)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurCHL)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurCHL)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCHL2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#After
completeCHLaft<-after[complete.cases(after$CHL), ] 
unique(completeCHLaft$newdate)

modAftCHL<-lm(CHL~Treatment,data=after)
DunAftCHL<-multcomp::glht(modAftCHL,Treatment=mcp(Group="Dunnett"))
summary(DunAftCHL)
rm(DunAftCHL)

mixmodAftCHL<-nlme::lme(CHL~ Treatment, random=~1|Week/Pond, data=after, na.action=na.omit)
#mixmodAftCHL<-nlme::lme(CHL~ Treatment, random=~1|Pond, correlation=nlme::corAR1(form=~1|Pond), data=after, na.action=na.omit)
anova(mixmodAftCHL)
acf(resid(mixmodAftCHL))

DunAftCHL<-multcomp::glht(mixmodAftCHL,Treatment=mcp(Group="Dunnett"))
summary(DunAftCHL)


Trmt = c("D","M","MD","O")
Var = rep("CHL",4) #replace "CHL" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftCHL)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftCHL)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftCHL)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCHL3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsCHL = rbind(DunnettsCHL1,DunnettsCHL2,DunnettsCHL3)
write.csv(DunnettsCHL,"Dunnetts_mixedmod_CHL.csv")




## Cyanobacteria dominance.
#Before.
completeCyanbef<-before[complete.cases(before$h_v_Cyanobacteria), ] 
unique(completeCyanbef$newdate)

modBefCyan<-lm(h_v_Cyanobacteria ~ Treatment,data=before)
DunBefCyan<-multcomp::glht(modBefCyan,Treatment=mcp(Group="Dunnett"))
summary(DunBefCyan)
rm(DunBefCyan)

mixmodBefCyan<-nlme::lme(h_v_Cyanobacteria ~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefCyan<-multcomp::glht(mixmodBefCyan,Treatment=mcp(Group="Dunnett"))
summary(DunBefCyan)

Trmt = c("D","M","MD","O")
Var = rep("Cyan",4) #replace "Cyan" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefCyan)$test$coefficients[2:5]
Zval = summary(DunBefCyan)$test$tstat[2:5]
Pval = summary(DunBefCyan)$test$pvalues[2:5]

DunnettsCyan1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#During.
completeCyandur<-during[complete.cases(during$h_v_Cyanobacteria), ] 
unique(completeCyandur$newdate)

modDurCyan<-lm(h_v_Cyanobacteria~Treatment,data=during)
DunDurCyan<-multcomp::glht(modDurCyan,Treatment=mcp(Group="Dunnett"))
summary(DunDurCyan)
rm(DunDurCyan)

mixmodDurCyan<-nlme::lme(h_v_Cyanobacteria~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurCyan<-multcomp::glht(mixmodDurCyan,Treatment=mcp(Group="Dunnett"))
summary(DunDurCyan)

Trmt = c("D","M","MD","O")
Var = rep("Cyan",4) #replace "Cyan" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurCyan)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurCyan)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurCyan)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCyan2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#After.
completeCyanaft<-after[complete.cases(after$h_v_Cyanobacteria), ] 
unique(completeCyanaft$newdate)

modAftCyan<-lm(h_v_Cyanobacteria~Treatment,data=after)
DunAftCyan<-multcomp::glht(modAftCyan,Treatment=mcp(Group="Dunnett"))
summary(DunAftCyan)
rm(DunAftCyan)

mixmodAftCyan<-nlme::lme(h_v_Cyanobacteria~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftCyan<-multcomp::glht(mixmodAftCyan,Treatment=mcp(Group="Dunnett"))
summary(DunAftCyan)

Trmt = c("D","M","MD","O")
Var = rep("Cyan",4) #replace "Cyan" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftCyan)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftCyan)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftCyan)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCyan3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsCyan = rbind(DunnettsCyan1,DunnettsCyan2,DunnettsCyan3)
write.csv(Dunnetts,"Dunnetts_mixedmod_Cyan.csv")




##Trait evenness
#Before
completeTEDbef<-before[complete.cases(before$TED), ] 
unique(completeTEDbef$newdate)

modBefTED<-lm(TED ~Treatment,data=before)
DunBefTED<-multcomp::glht(modBefTED,Treatment=mcp(Group="Dunnett"))
summary(DunBefTED)
rm(DunBefTED)

mixmodBefTED<-nlme::lme(TED ~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefTED<-multcomp::glht(mixmodBefTED,Treatment=mcp(Group="Dunnett"))
summary(DunBefTED)

Trmt = c("D","M","MD","O")
Var = rep("TED",4) #replace "TED" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefTED)$test$coefficients[2:5]
Zval = summary(DunBefTED)$test$tstat[2:5]
Pval = summary(DunBefTED)$test$pvalues[2:5]

DunnettsTED1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#During
completeTEDdur<-during[complete.cases(during$TED), ] 
unique(completeTEDdur$newdate)

modDurTED<-lm(TED ~Treatment,data=during)
DunDurTED<-multcomp::glht(modDurTED,Treatment=mcp(Group="Dunnett"))
summary(DunDurTED)
rm(DunDurTED)

mixmodDurTED<-nlme::lme(TED ~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurTED<-multcomp::glht(mixmodDurTED,Treatment=mcp(Group="Dunnett"))
summary(DunDurTED)

Trmt = c("D","M","MD","O")
Var = rep("TED",4) #replace "TED" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurTED)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurTED)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurTED)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsTED2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#After
completeTEDaft<-after[complete.cases(after$TED), ] 
unique(completeTEDaft$newdate)

modAftTED<-lm(TED ~Treatment,data=after)
DunAftTED<-multcomp::glht(modAftTED,Treatment=mcp(Group="Dunnett"))
summary(DunAftTED)
rm(DunAftTED)

mixmodAftTED<-nlme::lme(TED ~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftTED<-multcomp::glht(mixmodAftTED,Treatment=mcp(Group="Dunnett"))
summary(DunAftTED)

Trmt = c("D","M","MD","O")
Var = rep("TED",4) #replace "TED" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftTED)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftTED)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftTED)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsTED3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsTED = rbind(DunnettsTED1,DunnettsTED2,DunnettsTED3)
write.csv(DunnettsTED,"Dunnetts_mixedmod_TED.csv")




##CtoN
#Before
completeCtoNbef<-before[complete.cases(before$CtoN), ] 
unique(completeCtoNbef$newdate)

modBefCN<-lm(CtoN ~Treatment,data=completeCtoNbef)
DunBefCN<-multcomp::glht(modBefCN,Treatment=mcp(Group="Dunnett"))
summary(DunBefCN)
rm(DunBefCN)

mixmodBefCN<-nlme::lme(CtoN ~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefCN<-multcomp::glht(mixmodBefCN,Treatment=mcp(Group="Dunnett"))
summary(DunBefCN)

Trmt = c("D","M","MD","O")
Var = rep("CN",4) #replace "CN" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefCN)$test$coefficients[2:5]
Zval = summary(DunBefCN)$test$tstat[2:5]
Pval = summary(DunBefCN)$test$pvalues[2:5]

DunnettsCN1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)



#During
completeCtoNdur<-during[complete.cases(during$CtoN), ] 
unique(completeCtoNdur$newdate)

modDurCN<-lm(CtoN~Treatment,data=during)
DunDurCN<-multcomp::glht(modDurCN,Treatment=mcp(Group="Dunnett"))
summary(DunDurCN)
rm(DunDurCN)

mixmodDurCN<-nlme::lme(CtoN~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurCN<-multcomp::glht(mixmodDurCN,Treatment=mcp(Group="Dunnett"))
summary(DunDurCN)

Trmt = c("D","M","MD","O")
Var = rep("CN",4) #replace "CN" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurCN)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurCN)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurCN)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCN2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


## After
completeCtoNaft<-after[complete.cases(after$CtoN), ] 
unique(completeCtoNaft$newdate)

modAftCN<-lm(CtoN~Treatment,data=after)
DunAftCN<-multcomp::glht(modAftCN,Treatment=mcp(Group="Dunnett"))
summary(DunAftCN)
rm(DunAftCN)

mixmodAftCN<-nlme::lme(CtoN~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftCN<-multcomp::glht(mixmodAftCN,Treatment=mcp(Group="Dunnett"))
summary(DunAftCN)

Trmt = c("D","M","MD","O")
Var = rep("CN",4) #replace "CN" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftCN)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftCN)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftCN)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCN3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsCN = rbind(DunnettsCN1,DunnettsCN2,DunnettsCN3)
write.csv(DunnettsCN3,"Dunnetts_mixedmod_CN.csv")


##CtoP
#Before
completeCtoPbef<-before[complete.cases(before$CtoP), ] 
unique(completeCtoNbef$newdate)

modBefCP<-lm(CtoP~Treatment,data=completeCtoPbef)
DunBefCP<-multcomp::glht(modBefCP,Treatment=mcp(Group="Dunnett"))
summary(DunBefCP)
rm(DunBefCP)

mixmodBefCP<-nlme::lme(CtoP~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefCP<-multcomp::glht(mixmodBefCP,Treatment=mcp(Group="Dunnett"))
summary(DunBefCP)

Trmt = c("D","M","MD","O")
Var = rep("CP",4) #replace "CP" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefCP)$test$coefficients[2:5]
Zval = summary(DunBefCP)$test$tstat[2:5]
Pval = summary(DunBefCP)$test$pvalues[2:5]

DunnettsCP1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#During
completeCtoPdur<-during[complete.cases(during$CtoP), ] 
unique(completeCtoPdur$newdate)

modDurCP<-lm(CtoP~Treatment,data=during)
DunDurCP<-multcomp::glht(modDurCP,Treatment=mcp(Group="Dunnett"))
summary(DunDurCP)
rm(DunDurCP)

mixmodDurCP<-nlme::lme(CtoP~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurCP<-multcomp::glht(mixmodDurCP,Treatment=mcp(Group="Dunnett"))
summary(DunDurCP)

Trmt = c("D","M","MD","O")
Var = rep("CP",4) #replace "CP" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurCP)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurCP)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurCP)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCP2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#After
completeCtoPaft<-after[complete.cases(after$CtoP), ] 
unique(completeCtoPaft$newdate)

modAftCP<-lm(CtoP~Treatment,data=after)
DunAftCP<-multcomp::glht(modAftCP,Treatment=mcp(Group="Dunnett"))
summary(DunAftCP)
rm(DunAftCP)

mixmodAftCP<-nlme::lme(CtoP~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftCP<-multcomp::glht(mixmodAftCP,Treatment=mcp(Group="Dunnett"))
summary(DunAftCP)

Trmt = c("D","M","MD","O")
Var = rep("CP",4) #replace "CP" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftCP)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftCP)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftCP)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsCP3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsCP = rbind(DunnettsCP1,DunnettsCP2,DunnettsCP3)
write.csv(DunnettsCP,"Dunnetts_mixedmod_CP.csv")



###############
### Ecosystem
###############

##Nitrate
#Before
completeNitratebef<-before[complete.cases(before$Nitrate), ] 
unique(completeNitratebef$newdate)

modBefnitrate<-lm(Nitrate~Treatment,data=before)
DunBefnitrate<-multcomp::glht(modBefnitrate,Treatment=mcp(Group="Dunnett"))
summary(DunBefnitrate)
rm(DunBefnitrate)

mixmodBefnitrate<-nlme::lme(Nitrate~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefnitrate<-multcomp::glht(mixmodBefnitrate,Treatment=mcp(Group="Dunnett"))
summary(DunBefnitrate)

Trmt = c("D","M","MD","O")
Var = rep("nitrate",4) #replace "nitrate" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefnitrate)$test$coefficients[2:5]
Zval = summary(DunBefnitrate)$test$tstat[2:5]
Pval = summary(DunBefnitrate)$test$pvalues[2:5]

Dunnettsnitrate1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#During
completeNitratedur<-during[complete.cases(during$Nitrate), ] 
unique(completeNitratedur$newdate)

modDurnitrate<-lm(Nitrate~Treatment,data=during)
DunDurnitrate<-multcomp::glht(modDurnitrate,Treatment=mcp(Group="Dunnett"))
summary(DunDurnitrate)
rm(DunDurnitrate)

mixmodDurnitrate<-nlme::lme(Nitrate~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurnitrate<-multcomp::glht(mixmodDurnitrate,Treatment=mcp(Group="Dunnett"))
summary(DunDurnitrate)

Trmt = c("D","M","MD","O")
Var = rep("nitrate",4) #replace "nitrate" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurnitrate)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurnitrate)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurnitrate)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnettsnitrate2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#After
completeNitrateaft<-after[complete.cases(after$Nitrate), ] 
unique(completeNitrateaft$newdate)

modAftnitrate<-lm(Nitrate~Treatment,data=after)
DunAftnitrate<-multcomp::glht(modAftnitrate,Treatment=mcp(Group="Dunnett"))
summary(DunAftnitrate)
rm(DunAftnitrate)

mixmodAftnitrate<-nlme::lme(Nitrate~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftnitrate<-multcomp::glht(mixmodAftnitrate,Treatment=mcp(Group="Dunnett"))
summary(DunAftnitrate)

Trmt = c("D","M","MD","O")
Var = rep("nitrate",4) #replace "nitrate" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftnitrate)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftnitrate)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftnitrate)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnettsnitrate3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

Dunnettsnitrate = rbind(Dunnettsnitrate1,Dunnettsnitrate2,Dunnettsnitrate3)
write.csv(Dunnettsnitrate,"Dunnetts_mixedmod_CHL.csv")


##OrthoP
#Before
head(before)
completeorthobef<-before[complete.cases(before$orthoP), ] 
unique(completeorthobef$newdate)

modBefortho<-lm(orthoP~Treatment,data=before)
DunBefortho<-multcomp::glht(modBefortho,Treatment=mcp(Group="Dunnett"))
summary(DunBefortho)
rm(DunBefortho)

mixmodBefortho<-nlme::lme(orthoP~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefortho<-multcomp::glht(mixmodBefortho,Treatment=mcp(Group="Dunnett"))
summary(DunBefortho)

Trmt = c("D","M","MD","O")
Var = rep("ortho",4) #replace "ortho" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefortho)$test$coefficients[2:5]
Zval = summary(DunBefortho)$test$tstat[2:5]
Pval = summary(DunBefortho)$test$pvalues[2:5]

Dunnettsortho1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#During
completeorthodur<-during[complete.cases(during$orthoP), ] 
unique(completeorthodur$newdate)

modDurortho<-lm(orthoP~Treatment,data=during)
DunDurortho<-multcomp::glht(modDurortho,Treatment=mcp(Group="Dunnett"))
summary(DunDurortho)
rm(DunDurortho)

mixmodDurortho<-nlme::lme(orthoP~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurortho<-multcomp::glht(mixmodDurortho,Treatment=mcp(Group="Dunnett"))
summary(DunDurortho)

Trmt = c("D","M","MD","O")
Var = rep("ortho",4) #replace "ortho" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurortho)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurortho)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurortho)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnettsortho2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#After
completeorthoaft<-after[complete.cases(after$orthoP), ] 
unique(completeorthoaft$newdate)

modAftortho<-lm(orthoP~Treatment,data=after)
DunAftortho<-multcomp::glht(modAftortho,Treatment=mcp(Group="Dunnett"))
summary(DunAftortho)
rm(DunDurortho)

mixmodAftortho<-nlme::lme(orthoP~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftortho<-multcomp::glht(mixmodAftortho,Treatment=mcp(Group="Dunnett"))
summary(DunAftortho)

Trmt = c("D","M","MD","O")
Var = rep("ortho",4) #replace "ortho" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftortho)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftortho)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftortho)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnettsortho3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

Dunnettsortho = rbind(Dunnettsortho1,Dunnettsortho2,Dunnettsortho3)
write.csv(Dunnettsortho,"Dunnetts_mixedmod_ortho.csv")


##Turbidity
#Before
completeTurbbef<-before[complete.cases(before$Turb), ] 
unique(completeTurbbef$newdate)

modBefTurb<-lm(Turbidity~Treatment,data=before)
DunBefTurb<-multcomp::glht(modBefTurb,Treatment=mcp(Group="Dunnett"))
summary(DunBefTurb)
rm(DunBefTurb)

mixmodBefTurb<-nlme::lme(Turbidity~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefTurb<-multcomp::glht(mixmodBefTurb,Treatment=mcp(Group="Dunnett"))
summary(DunBefTurb)

Trmt = c("D","M","MD","O")
Var = rep("Turb",4) #replace "Turb" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefTurb)$test$coefficients[2:5]
Zval = summary(DunBefTurb)$test$tstat[2:5]
Pval = summary(DunBefTurb)$test$pvalues[2:5]

DunnettsTurb1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#During
completeTurbdur<-during[complete.cases(during$Turb), ] 
unique(completeTurbdur$newdate)

modDurTurb<-lm(Turbidity~Treatment,data=during)
DunDurTurb<-multcomp::glht(modDurTurb,Treatment=mcp(Group="Dunnett"))
summary(DunDurTurb)
rm(DunDurTurb)

mixmodDurTurb<-nlme::lme(Turbidity~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurTurb<-multcomp::glht(mixmodDurTurb,Treatment=mcp(Group="Dunnett"))
summary(DunDurTurb)

Trmt = c("D","M","MD","O")
Var = rep("Turb",4) #replace "Turb" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurTurb)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurTurb)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurTurb)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsTurb2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#After
completeTurbaft<-after[complete.cases(after$Turb), ] 
unique(completeTurbaft$newdate)

modAftTurb<-lm(Turbidity~Treatment,data=after)
DunAftTurb<-multcomp::glht(modAftTurb,Treatment=mcp(Group="Dunnett"))
summary(DunAftTurb)
rm(DunAftTurb)

mixmodAftTurb<-nlme::lme(Turbidity~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
#mixmodAftTurb<-nlme::lme(Turbidity~ Treatment, random=~1|Pond, correlation=nlme::corAR1(form=~1|Pond), data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
acf(resid(mixmodAftTurb))
DunAftTurb<-multcomp::glht(mixmodAftTurb,Treatment=mcp(Group="Dunnett"))
summary(DunAftTurb)

Trmt = c("D","M","MD","O")
Var = rep("Turb",4) #replace "Turb" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftTurb)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftTurb)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftTurb)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsTurb3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


DunnettsTurb = rbind(DunnettsTurb1,DunnettsTurb2,DunnettsTurb3)
write.csv(DunnettsTurb,"Dunnetts_mixedmod_Turb.csv")

##Dissolved Oxygen
#Before
head(before)
completeDObef<-before[complete.cases(before$O2_1), ] 
unique(completeDObef$newdate)

modBefdO<-lm(O2_1~Treatment,data=before)
DunBefdO<-multcomp::glht(modBefdO,Treatment=mcp(Group="Dunnett"))
summary(DunBefdO)
rm(DunBefdO)

mixmodBefdO<-nlme::lme(O2_1~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefdO<-multcomp::glht(mixmodBefdO,Treatment=mcp(Group="Dunnett"))
summary(DunBefdO)

Trmt = c("D","M","MD","O")
Var = rep("dO",4) #replace "dO" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefdO)$test$coefficients[2:5]
Zval = summary(DunBefdO)$test$tstat[2:5]
Pval = summary(DunBefdO)$test$pvalues[2:5]

DunnettsdO1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#During
completeDOdur<-during[complete.cases(during$O2_1), ] 
unique(completeDOdur$newdate)

modDurdO<-lm(O2_1~Treatment,data=during)
DunDurdO<-multcomp::glht(modDurdO,Treatment=mcp(Group="Dunnett"))
summary(DunDurdO)
rm(DunDurdO)

mixmodDurdO<-nlme::lme(O2_1~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurdO<-multcomp::glht(mixmodDurdO,Treatment=mcp(Group="Dunnett"))
summary(DunDurdO)

Trmt = c("D","M","MD","O")
Var = rep("dO",4) #replace "dO" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurdO)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurdO)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurdO)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsdO2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#After
completeDOaft<-after3[complete.cases(after3$O2_1), ] 
unique(completeDOaft$newdate)

modAftdO<-lm(O2_1~Treatment,data=after3)
DunAftdO<-multcomp::glht(modAftdO,Treatment=mcp(Group="Dunnett"))
summary(DunAftdO)
rm(DunAftdO)

mixmodAftdO<-nlme::lme(O2_1~ Treatment, random=~1|Week/Pond, data=after3,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftdO<-multcomp::glht(mixmodAftdO,Treatment=mcp(Group="Dunnett"))
summary(DunAftdO)

Trmt = c("D","M","MD","O")
Var = rep("dO",4) #replace "dO" with name of variable of interest.
Period = rep("After3",4) #replace time period you're considering.
Coef = summary(DunAftdO)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftdO)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftdO)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettsdO3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsdO = rbind(DunnettsdO1,DunnettsdO2,DunnettsdO3)
write.csv(DunnettsdO,"Dunnetts_mixedmod_dO.csv")



##pH
#Before
completepHbef<-before[complete.cases(before$pH), ] 
unique(completepHbef$newdate)

modBefpH<-lm(pH~Treatment,data=completepHbef)
DunBefpH<-multcomp::glht(modBefpH,Treatment=mcp(Group="Dunnett"))
summary(DunBefpH)
rm(DunBefpH)

mixmodBefpH<-nlme::lme(pH~ Treatment, random=~1|Week/Pond, data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefpH<-multcomp::glht(mixmodBefpH,Treatment=mcp(Group="Dunnett"))
summary(DunBefpH)

Trmt = c("D","M","MD","O")
Var = rep("pH",4) #replace "pH" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefpH)$test$coefficients[2:5]
Zval = summary(DunBefpH)$test$tstat[2:5]
Pval = summary(DunBefpH)$test$pvalues[2:5]

DunnettspH1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#During
completepHdur<- during[complete.cases(during$pH), ] 
unique(completepHdur$newdate)

modDurpH<-lm(pH~Treatment,data=during)
DunDurpH<-multcomp::glht(modDurpH,Treatment=mcp(Group="Dunnett"))
summary(DunDurpH)
rm(DunDurpH)

mixmodDurpH<-nlme::lme(pH~ Treatment, random=~1|Week/Pond, data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurpH<-multcomp::glht(mixmodDurpH,Treatment=mcp(Group="Dunnett"))
summary(DunDurpH)

Trmt = c("D","M","MD","O")
Var = rep("pH",4) #replace "pH" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurpH)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurpH)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurpH)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettspH2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#After
completepHaft<- after[complete.cases(after$pH), ] 
unique(completepHaft$newdate)

modAftpH<-lm(pH~Treatment,data=after)
DunAftpH<-multcomp::glht(modAftpH,Treatment=mcp(Group="Dunnett"))
summary(DunAftpH)
rm(DunAftpH)

mixmodAftpH<-nlme::lme(pH~ Treatment, random=~1|Week/Pond, data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftpH<-multcomp::glht(mixmodAftpH,Treatment=mcp(Group="Dunnett"))
summary(DunAftpH)

Trmt = c("D","M","MD","O")
Var = rep("pH",4) #replace "pH" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftpH)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftpH)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftpH)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

DunnettspH3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettspH = rbind(DunnettspH1,DunnettspH2,DunnettspH3)
write.csv(DunnettspH,"Dunnetts_mixedmod_pH.csv")



#### others

#Daphnia Before
completeDaphbef<-before[complete.cases(before$Daphnia_adult), ] 
unique(completeDaphbef$newdate)

modBefdaphnia<-lm(Daphnia_adult~Treatment,data=before)
DunBefdaphnia<-multcomp::glht(modBefdaphnia,Treatment=mcp(Group="Dunnett"))
summary(DunBefdaphnia)
rm(DunBefdaphnia)

mixmodBefdaphnia<-nlme::lme(Daphnia_adult~ Treatment,random=~1|Week/Pond,data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefdaphnia<-multcomp::glht(mixmodBefdaphnia,Treatment=mcp(Group="Dunnett"))
summary(DunBefdaphnia)

Trmt = c("D","M","MD","O")
Var = rep("daphnia",4) #replace "daphnia" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefdaphnia)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunBefdaphnia)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunBefdaphnia)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnetts1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#Daphnia During
completeDaphdur<-during[complete.cases(during$Daphnia_adult), ] 
unique(completeDaphdur$newdate)

modDurdaphnia<-lm(Daphnia_adult~Treatment,data=during)
DunDurdaphnia<-multcomp::glht(modDurdaphnia,Treatment=mcp(Group="Dunnett"))
summary(DunDurdaphnia)
rm(DunDurdaphnia)

mixmodDurdaphnia<-nlme::lme(Daphnia_adult~ Treatment,random=~1|Week/Pond,data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurdaphnia<-multcomp::glht(mixmodDurdaphnia,Treatment=mcp(Group="Dunnett"))
summary(DunDurdaphnia)

Trmt = c("D","M","MD","O")
Var = rep("daphnia",4) #replace "daphnia" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurdaphnia)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurdaphnia)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurdaphnia)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnetts2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)


#Daphnia After
completeDaphaft<-after[complete.cases(after$Daphnia_adult), ] 
unique(completeDaphaft$newdate)

modAftdaphnia<-lm(Daphnia_adult~Treatment,data=after)
DunAftdaphnia<-multcomp::glht(modAftdaphnia,Treatment=mcp(Group="Dunnett"))
summary(DunAftdaphnia)
rm(DunAftdaphnia)

mixmodAftdaphnia<-nlme::lme(Daphnia_adult~ Treatment,random=~1|Week/Pond,data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftdaphnia<-multcomp::glht(mixmodAftdaphnia,Treatment=mcp(Group="Dunnett"))
summary(DunAftdaphnia)

Trmt = c("D","M","MD","O")
Var = rep("daphnia",4) #replace "daphnia" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftdaphnia)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftdaphnia)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftdaphnia)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnetts3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsDaph = rbind(Dunnetts1,Dunnetts2,Dunnetts3)
write.csv(DunnettsDaph,"Dunnetts_mixedmod_CHL.csv")


#Copepod Before
completeCopebef<-before[complete.cases(before$Copepod_adult), ] 
unique(completeCopebef$newdate)

modBefcope<-lm(Copepod_adult~Treatment,data=before)
DunBefcope<-multcomp::glht(modBefcope,Treatment=mcp(Group="Dunnett"))
summary(DunBefcope)
rm(DunBefcope)

mixmodBefcope<-nlme::lme(Copepod_adult~ Treatment,random=~1|Week/Pond,data=before,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunBefcope<-multcomp::glht(mixmodBefcope,Treatment=mcp(Group="Dunnett"))
summary(DunBefcope)

Trmt = c("D","M","MD","O")
Var = rep("cope",4) #replace "cope" with name of variable of interest.
Period = rep("Before",4) #replace time period you're considering.
Coef = summary(DunBefcope)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunBefcope)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunBefcope)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnetts1 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#Copepod During
completeCopedur<-during[complete.cases(during$Copepod_adult), ] 
unique(completeCopedur$newdate)

modDurcope<-lm(Copepod_adult~Treatment,data=during)
DunDurcope<-multcomp::glht(modDurcope,Treatment=mcp(Group="Dunnett"))
summary(DunDurcope)

mixmodDurcope<-nlme::lme(Copepod_adult~ Treatment,random=~1|Week/Pond,data=during,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunDurcope<-multcomp::glht(mixmodDurcope,Treatment=mcp(Group="Dunnett"))
summary(DunDurcope)

Trmt = c("D","M","MD","O")
Var = rep("cope",4) #replace "cope" with name of variable of interest.
Period = rep("During",4) #replace time period you're considering.
Coef = summary(DunDurcope)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunDurcope)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunDurcope)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnetts2 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

#Copepod After
completeCopeaft<-after[complete.cases(after$Copepod_adult), ] 
unique(completeCopeaft$newdate)

modAftcope<-lm(Copepod_adult~Treatment,data=after)
DunAftcope<-multcomp::glht(modAftcope,Treatment=mcp(Group="Dunnett"))
summary(DunAftcope)
rm(DunAftcope)

mixmodAftcope<-nlme::lme(Copepod_adult~ Treatment,random=~1|Week/Pond,data=after,na.action=na.omit)#,correlation=nlme::corAR1(form=~1|Pond)
DunAftcope<-multcomp::glht(mixmodAftcope,Treatment=mcp(Group="Dunnett"))
summary(DunAftcope)

Trmt = c("D","M","MD","O")
Var = rep("cope",4) #replace "cope" with name of variable of interest.
Period = rep("After",4) #replace time period you're considering.
Coef = summary(DunAftcope)$test$coefficients[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Zval = summary(DunAftcope)$test$tstat[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur
Pval = summary(DunAftcope)$test$pvalues[2:5] #replace 3 letters indicating model time period e.g. Bef becomes Dur

Dunnetts3 = cbind(Trmt,Var,Period,Coef,Zval,Pval)

DunnettsCope = rbind(Dunnetts1,Dunnetts2,Dunnetts3)
write.csv(DunnettsCope,"Dunnetts_mixedmod_cope.csv")


#####Data binding
DunnettsAll = rbind(DunnettsCHL,DunnettsCyan,DunnettsTED,DunnettsCN,DunnettsCP,Dunnettsnitrate,Dunnettsortho,DunnettsTurb,DunnettsdO,DunnettspH,DunnettsDaph,DunnettsCope)
write.csv(DunnettsAll,"Dunnetts_mixedmod_final_All.csv")
