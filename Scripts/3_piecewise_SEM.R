rm(list=ls())

library("vegan")
library("dplyr")
library('ape')
library('caper')
library('nlme')
library('lavaan')
library('curl')
library('mnormt')
library('piecewiseSEM')
library('diptest')
library('corrgram')
library('stringr')
library('lubridate')
library("ggplot2")
library("digest")
library('RColorBrewer')
library('lubridate')
library("tidyverse")

####### PIECEWISE SEM #######
master<-read.csv("Narwani_et_al_2019_Proceedings_B_data.csv")

#factor
master[,"Treatment"]<-factor(master[,"Treatment"])
master[,"Dreisse"]<-factor(master[,"Dreisse"])

submaster = master[master$Week<22,] #subsetting up to 7 weeks after the last nutrient addition.
submaster[,"Week"]<-factor(submaster[,"Week"])

master$Dreisse
submaster$Dreissena <- ifelse(submaster$Dreisse==1,c(1), c(-1))
submaster$Macrophytes
submaster$Myriophyllum <- ifelse(submaster$Macrophytes==1,c(1), c(-1))
is.numeric(submaster$Myriophyllum)

which(colnames(submaster)=="v_Cyanobacteria")
which(colnames(submaster)=="v_Green_algae")
#columns with functional group volume
gpvol <- data.frame(submaster[,97:102])

########


#log transfrom of functional group biovolume
gpcomp <- decostand(gpvol, method = "log",na.rm=T)
colnames(gpcomp)<-paste("l",colnames(gpcomp),sep="_")
submaster <-cbind(submaster,gpcomp)
head(submaster)

SEMvars<-c("Pond","Week","Block", "Rep", "Dreissena", "Myriophyllum","Nutrients", "TreatNR", "Treatment",
           "Turbidity", "pH", "orthoP", "Nitrate", "Temperature", "O2", "CHL", "TED",
           "l_v_Cyanobacteria","l_v_Golden_algae","l_v_Diatoms","l_v_Dinoflagellates","l_v_Cryptophytes","l_v_Green_algae",
           "Daphnia_adult","Copepod_adult","CtoN","CtoP")

SEMdata = submaster[,SEMvars]

#### Fixed categorical variables for mod1. Need to do this for the rest as well.
SEMdata_1 <- SEMdata %>% na.omit(CHL,l_v_Cyanobacteria,l_v_Green_algae,orthoP,Nitrate)
SEMdata_2<- subset(SEMdata_1, Treatment!="O")
head(SEMdata_2)
SEMdata_2$Pond<- as.factor(SEMdata_2$Pond)

dim(SEMdata_2)
 
# lme and temporal autocorrelation.

# Started with all phyto groups in explanation of chlorophyll and TED.
# Looked for effects of phytos on the nutrients.
# Looking for predictors of the grazer community (suggestion incorporated from E-L reviews).
# Added correlation structures for variables where a causal link shows up in d-sep tests, but I can't imagine the causal relation being valid.

psem0.dummy <- psem(
  lme(CHL~ l_v_Cyanobacteria + l_v_Golden_algae + l_v_Diatoms + l_v_Dinoflagellates + l_v_Cryptophytes + l_v_Green_algae,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria + l_v_Golden_algae + l_v_Diatoms + l_v_Dinoflagellates + l_v_Cryptophytes + l_v_Green_algae,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena*Myriophyllum  + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_Golden_algae ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Diatoms ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Dinoflagellates ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cryptophytes ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Green_algae ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ CHL + TED + l_v_Cyanobacteria, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~ CHL + TED + l_v_Cyanobacteria, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(Nitrate ~  CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(orthoP ~  CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~  CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ CHL+ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit)
  )

summary(psem0.dummy, .progressBar = F)

# first add variables from d-sep tests


psem1.dummy <- psem(
  lme(CHL~ TED + Temperature + l_v_Cyanobacteria + l_v_Golden_algae + l_v_Diatoms + l_v_Dinoflagellates + l_v_Cryptophytes + l_v_Green_algae,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria + l_v_Golden_algae + l_v_Diatoms + l_v_Dinoflagellates + l_v_Cryptophytes + l_v_Green_algae,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena*Myriophyllum  + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_Golden_algae ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Diatoms ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Dinoflagellates ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cryptophytes ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Green_algae ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ CHL + TED + l_v_Cyanobacteria, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~Temperature + CHL + TED + l_v_Cyanobacteria, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(Nitrate ~  l_v_Green_algae + l_v_Dinoflagellates + CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(orthoP ~  l_v_Green_algae + CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~  l_v_Dinoflagellates + l_v_Diatoms + Temperature + CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ Temperature + Myriophyllum + CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~  l_v_Green_algae + l_v_Cryptophytes + l_v_Dinoflagellates + Myriophyllum + CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  orthoP %~~% Nitrate,#
  orthoP %~~% Temperature,#
  orthoP %~~% Daphnia_adult,#
  O2 %~~% Nitrate,#
  pH  %~~%  Turbidity,#
  pH  %~~%  Nitrate, #
  pH  %~~%  orthoP,#
  pH  %~~%  O2,#
  Turbidity  %~~% Nitrate,#
  l_v_Dinoflagellates %~~% l_v_Cyanobacteria, #
  l_v_Green_algae %~~% l_v_Cyanobacteria, #
  l_v_Green_algae %~~% l_v_Dinoflagellates #
)
  
summary(psem1.dummy, .progressBar = F)

#now remove things that are not significant
psem2.dummy <- psem(
  lme(CHL~ TED + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena*Myriophyllum  + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_Golden_algae ~  Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Diatoms ~  Dreissena + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Dinoflagellates ~  Dreissena*Myriophyllum, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Green_algae ~ Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~Temperature + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(Nitrate ~  l_v_Dinoflagellates, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(orthoP ~  l_v_Green_algae + CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~  l_v_Diatoms + Temperature + CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ Temperature + Myriophyllum + CHL, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~  l_v_Green_algae + Myriophyllum, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  orthoP %~~% Nitrate,
  orthoP %~~% Daphnia_adult,
  pH  %~~%  Nitrate, 
  pH  %~~%  orthoP,
  pH  %~~%  O2,
  Turbidity %~~% Nitrate,
  l_v_Dinoflagellates %~~% l_v_Cyanobacteria, 
  l_v_Green_algae %~~% l_v_Cyanobacteria, 
  l_v_Green_algae %~~% l_v_Dinoflagellates 
)

summary(psem2.dummy, .progressBar = F)


#again return things with significant d-sep tests
#deal with strange directionality in orthoP
psem3.dummy <- psem(
  lme(CHL~   TED + l_v_Golden_algae + l_v_Green_algae + l_v_Dinoflagellates + orthoP + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria , random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena*Myriophyllum + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_Golden_algae ~  Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Diatoms ~  Dreissena + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Dinoflagellates ~  Dreissena*Myriophyllum, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Green_algae ~ Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~ TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(Nitrate ~  l_v_Dinoflagellates, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(orthoP ~ l_v_Green_algae + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  #due to implausible +ve effect of CHL orthoP, we reversed direction, i.e. effect of orthoP on CHL (ns).
  lme(Turbidity ~  CHL + TED + l_v_Diatoms + l_v_Dinoflagellates + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ CHL + Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~  l_v_Green_algae + Myriophyllum, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  orthoP %~~% Nitrate,
  orthoP %~~% Daphnia_adult,
  orthoP %~~% Temperature,
  pH  %~~%  Nitrate, 
  pH  %~~%  orthoP,
  pH  %~~%  O2,
  Turbidity %~~% Nitrate,
  Turbidity %~~% pH,
  l_v_Dinoflagellates %~~% l_v_Cyanobacteria, 
  l_v_Green_algae %~~% l_v_Cyanobacteria, 
  l_v_Green_algae %~~% l_v_Dinoflagellates 
)


summary(psem3.dummy, .progressBar = F)


















#### Proc B ESM model alternative

# lme and temporal autocorrelation.

# only cyanobacteria and eukaryotes in explanation of chlorophyll and TED.
# Looked for effects of phytos on the nutrients.
# Looking for predictors of the grazer community.
# Added correlation structures for variables where a causal link shows up in d-sep tests, but I can't imagine the causal relation being valid.


SEMdata_2$l_v_eukaryotes = SEMdata_2$l_v_Golden_algae + SEMdata_2$l_v_Diatoms + SEMdata_2$l_v_Dinoflagellates + SEMdata_2$l_v_Cryptophytes + SEMdata_2$l_v_Green_algae
head(SEMdata_2)

psem0.P.dummy <- psem(
  lme(CHL~ l_v_Cyanobacteria + l_v_eukaryotes,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria + l_v_eukaryotes,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena*Myriophyllum  + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_eukaryotes ~  Dreissena*Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(Nitrate ~  CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(orthoP ~  CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~  CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2~ CHL+ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit)
)

summary(psem0.P.dummy, .progressBar = F)


##### Adding

psem1.P.dummy <- psem(
  lme(CHL ~ TED + l_v_Cyanobacteria + l_v_eukaryotes + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria + l_v_eukaryotes,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena * Myriophyllum + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_eukaryotes ~ Dreissena * Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~ CHL + TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(Nitrate ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(orthoP ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~ CHL + TED + Nitrate + orthoP + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ CHL + TED + Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~ CHL + TED + Myriophyllum + Dreissena, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  orthoP %~~% Temperature,
  orthoP %~~% Nitrate,
  orthoP %~~% Daphnia_adult,
  O2 %~~% Nitrate,
  pH %~~%  Nitrate,
  pH %~~%  orthoP,
  pH %~~% Turbidity,
  pH %~~%  O2
)

summary(psem1.P.dummy, .progressBar = F)

##### Simplifying

psem2.P.dummy <- psem(
  lme(CHL ~ Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena * Myriophyllum + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_eukaryotes ~ Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~ TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  #lme(Nitrate ~ CHL + TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(orthoP ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~ CHL + TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ CHL+ Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~ CHL + TED + Myriophyllum + Dreissena, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  #orthoP %~~% Temperature,
  #orthoP %~~% Nitrate,
  orthoP %~~% Daphnia_adult,
  #O2 %~~% Nitrate,
  #pH %~~%  Nitrate,
  pH %~~%  orthoP,
  #pH %~~% Turbidity,
  pH %~~%  O2
)

summary(psem2.P.dummy, .progressBar = F)

##### Adding again

psem3.P.dummy <- psem(
  lme(CHL ~ l_v_Cyanobacteria + l_v_eukaryotes + TED + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena * Myriophyllum + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_eukaryotes ~ Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~ TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(orthoP ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~ CHL + TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ CHL+ Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~ CHL + TED + Turbidity + Myriophyllum + Dreissena, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  orthoP %~~% Temperature,
  orthoP %~~% Turbidity,
  orthoP %~~% Daphnia_adult,
  pH %~~%  orthoP,
  pH %~~% Turbidity,
  pH %~~%  O2
)

summary(psem3.P.dummy, .progressBar = F)

##### Simplifying again

psem4.P.dummy <- psem(
  lme(CHL ~ l_v_Cyanobacteria + l_v_eukaryotes + TED + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(TED ~ l_v_Cyanobacteria,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(l_v_Cyanobacteria ~ Dreissena * Myriophyllum + Temperature,random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit), 
  lme(l_v_eukaryotes ~ Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Daphnia_adult ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Copepod_adult ~ TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  lme(orthoP ~ TED, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(Turbidity ~ CHL + TED + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(O2 ~ CHL+ Myriophyllum + Temperature, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  lme(pH ~ l_v_eukaryotes + Turbidity + Myriophyllum + Dreissena, random=~1|Pond,correlation=corAR1(form=~1|Pond),data= SEMdata_2,na.action=na.omit),
  
  orthoP %~~% Temperature,
  orthoP %~~% Turbidity,
  orthoP %~~% Daphnia_adult,
  pH %~~%  orthoP,
  pH %~~% Turbidity,
  pH %~~%  O2
)

summary(psem4.P.dummy, .progressBar = F)

