#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN SEED CANOLA FIELDS (2016)
# Libraries and ggplot theme ----------------------------------------------
library(ggplot2)
library(MCMCglmm)
library(nlme)
library(mgcv)
library(MASS)
library(dplyr)
library(tidyr)
library(reshape2)

#Big-text theme, no grid lines (used for Bayer 2016 presentation)
prestheme=theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15), 
                axis.title=element_text(size=20),
                title=element_text(size=20),
                axis.line.x=element_line(colour='black'),
                axis.line.y=element_line(colour='black'),
                #panel.grid.major=element_line(size=0.5,colour='black',linetype='dotted'),
                #panel.grid.minor=element_blank,
                panel.border=element_blank())                
theme_set(theme_classic()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

# Load and organize data  ---------------------------------------------------------

# Load data
# load("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\seedfieldData2016.RData")

fields=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\field info 2016.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))

survey=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\survey data 2016.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))

plants=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\plant info 2016 seed.csv",stringsAsFactors=T,strip.white=T,na.strings = c("",'NA'))

#Nectar data

brix2mgul=function(B){ #B is corrected %brix read from refractometer
  1.177*((0.9988*B) +(0.389*(B^2)/100) +(0.113*(B^3)/10000) +(0.091*(B^4)/1000000) -(0.039*(B^5)/100000000))/100 } #Pyke's formula (from Ralph)

temp=fields[grep('obtained',fields$Other.notes),] %>% #Extra rows of pre-survey nectar levels
  mutate(FlwNetEnd=FlwNetStart)
temp[,c(26:35)]=rbind(c(0,0,0,0,0.5,1,0,0.5,0,0.5), #Nectar values taken before survey
                     c(0,0,0,0,NA,1,1,1,0.5,NA),
                     c(0,0,0,0,1,0,0,0,0.5,1))
fields$Surveyed[fields$Field=='Bateman 7'|fields$Field=='Bateman 8']=fields$Surveyed[fields$Field=='Bateman 9'] #Date on Bateman 7 and Bateman 8 fields (nectar surveyed same day as Bateman 9, but no insect survey)

nectar=rbind(fields,temp) %>%
  select(Field,Variety,Surveyed,
         StartTime=FlwNetStart,EndTime=FlwNetEnd,
         Temp=Temperature,WindSp=Wind.Speed.km.hr.,WindDir=Wind.Direction,RH=Humidity,Cloud=Cloud.Cover,
         NectarF1:MicrocapVol) %>%
  filter(!is.na(EndTime)) %>% #Removes unobserved fields
  mutate(Field=as.factor(Field),Variety=as.factor(Variety)) %>%
  mutate(StartTime=as.POSIXct(paste(Surveyed,StartTime),format='%b %d, %Y %H:%M'), #Recasts times/dates
    EndTime=as.POSIXct(paste(Surveyed,EndTime),format='%b %d, %Y %H:%M'),
    Surveyed=as.POSIXct(Surveyed,format='%b %d, %Y'),
    TotalTime=as.numeric(EndTime-StartTime,units='hours')) %>% #Total time (hrs)
  gather("Flower","Len",NectarF1:NectarM5) %>% #Gathers all nectar samples into 1 column
  mutate(Flower=substr(Flower,7,8)) %>% #Strips text from Flower column
  separate(Flower,c('Bay','Sample'),1) %>% #Creates Bay and Sample column
  mutate(Bay=as.factor(Bay),Sample=as.numeric(Sample)) %>% #Recasts Bay and Sample
  mutate(mguLF=brix2mgul(NectarPercF),mguLM=brix2mgul(NectarPercM)) %>% #Converts %sugar to mg/uL
  mutate(Vol=MicrocapVol*Len/32,Sugar=ifelse(Bay=='F',Vol*mguLF,Vol*mguLM)) %>% #Nectar volume (uL), sugar (mg)
  filter(!is.na(Vol)) %>% #Strips empty rows
  select(-NectarPercF,-NectarPercM,-MicrocapVol,-Len,-mguLF,-mguLM)  #Cleanup

rm(brix2mgul,temp)

#Organize field-level data
fields=rename(fields,Temp=Temperature,WindSp=Wind.Speed.km.hr.,WindDir=Wind.Direction,RH=Humidity,Cloud=Cloud.Cover) %>%
  select(Field:Cloud,-Other.notes,-Elevation,-FlwNetStart,-FlwNetEnd) %>%
  filter(!is.na(Harvested)) %>% #Removes unused fields
  mutate(Area=Acres*0.404686) %>% #Converts to hectares
  select(-Acres) %>%
  mutate(Start.Time=paste(Surveyed,Start.Time),End.Time=paste(Surveyed,End.Time)) %>%
  mutate(BowlStart=ifelse(is.na(BowlStart),NA,paste(Surveyed,BowlStart)),BowlEnd=ifelse(is.na(BowlEnd),NA,paste(Surveyed,BowlEnd))) %>%
  mutate(Start.Time=as.POSIXct(Start.Time,format='%b %d, %Y %H:%M'),End.Time=as.POSIXct(End.Time,format='%b %d, %Y %H:%M')) %>%
  mutate(BowlStart=as.POSIXct(BowlStart,format='%b %d, %Y %H:%M'),BowlEnd=as.POSIXct(BowlEnd,format='%b %d, %Y %H:%M'))

#Join data from field data
survey=left_join(select(fields,Field,Treatment,Variety,Surveyed,Temp:Area),
          survey,by='Field')

#Organize survey data
survey=rename(survey,FlDens=Fls.50cm.2) %>%
  rowwise() %>%
  #Minimum distance to leafcutter shelter
  mutate(Surveyed=as.Date(Surveyed,format='%b %d, %Y'),minDist=min(LeafShelter1,LeafShelter2,LeafShelter3,LeafShelter4,na.rm=T)) %>%
  select(-LeafShelter1:-LeafShelter4) %>%
  mutate(StartTime=as.POSIXct(paste(Surveyed,StartTime),format='%b %d, %Y %H:%M'),
         EndTime=as.POSIXct(paste(Surveyed,EndTime),format='%b %d, %Y %H:%M')) %>%
  mutate(Surveyed=as.Date(Surveyed,format='%b %d, %Y %H:%M'))

#Pollen data
pollen=filter(survey,Bay=='F') %>%
  select(Field:EdgeDir,minDist,HoneybeePollen:MountingOK,-Bay) %>%
  gather(Flower,Count,Pollen1:Pollen5) %>%
  mutate(FieldPlot=paste(Field,Distance,sep='_')) %>%
  select(-Flower)

#Visitation data
visitors=select(survey,Field:EdgeDir,StartTime:OtherFly,minDist)

#Plant data
seeds=transmute(plants,Field=FieldName,BagNum=Bag.,VegMass=TotalMass-SeedMass,SeedMass,Branch,
                Pods=Pods+Bag+BagTscar,Missing=Missing-Bag,
                Pod1,Pod2,Pod3,Pod4,Pod5,
                Weigh1,Weigh2,Weigh3,Weigh4,Weigh5) %>% #Selects correct columns
  gather('Pod','Value',Pod1:Weigh5) %>% #Melts pod count and pod weight
  separate(Pod,into=c('Param','PodNum'),sep=-1) %>%
  mutate(Param=as.factor(Param),PodNum=as.numeric(PodNum)) %>%
  spread(Param,Value) %>% #Casts pod count and pod weight back
  rename(PodCount=Pod,PodMass=Weigh,Pod=PodNum) %>%
  group_by(Field,Pod)

se=function(x) sd(x)/sqrt(length(x)) #Convenience SE function

#Mean and SE of metrics for each plant (all mean values except for pod measurements will equal plant-level metrics)
plants=summarise_all(group_by(seeds,Field,BagNum),funs(mean,se)) %>%
  select(Field,BagNum,VegMass=VegMass_mean,SeedMass=SeedMass_mean,
         Branch=Branch_mean,Pods=Pods_mean,Missing=Missing_mean,
         AvPodCount=PodCount_mean,SEPodCount=PodCount_se,AvPodMass=PodMass_mean,SEPodMass=PodMass_se)

#Data from plots to merge in
bagdat=filter(survey,Bay=='F') %>%
  mutate(AvgPol=mean(Pollen1:Pollen5,na.rm=T)) %>%
  select(-Pollen1:-Pollen5) %>%
  #Data to merge from survey dataframe
  transmute(Field,Distance,Variety,EdgeCent=Edge.Cent,EdgeDir,HoneybeePollen,HoneybeeNectar,Leafcutterbee,Bumblebee,Otherbee,Hoverfly,OtherFly,Other,Other,minShelter=minDist,AvgPol,PlDens=Plants.m.2,Hectares=Area,Treatment,NumHives=NA,BagNum=BagNumbers) %>%
  separate(BagNum,c('Bag1','Bag2','Bag3')) %>%
  gather('Bag','BagNum',Bag1:Bag3) %>%
  mutate(BagNum=as.numeric(BagNum)) %>%
  select(-Bag)

#Merges Distances and other data into seed and plant dataframes
seeds=left_join(seeds,select(bagdat,-Field),by='BagNum') %>%
  rename(Plant=BagNum) %>%
  arrange(Field,Distance,Plant) %>%
  ungroup() %>%
  unite(FieldPlot,Field,Distance,EdgeCent,remove=F)%>%
  mutate(FieldPlot=as.factor(FieldPlot)) #%>%
  #filter(!is.na(VegMass))

plants=left_join(plants,select(bagdat,-Field),by='BagNum') %>%
  rename(Plant=BagNum) %>%
  arrange(Field,Distance,Plant) %>%
  ungroup() %>%
  unite(FieldPlot,Field,Distance,EdgeCent,remove=F) %>%
  mutate(FieldPlot=as.factor(FieldPlot)) #%>%
  #filter(!is.na(VegMass))

# pollen=filter(flowers,Bay=='F') %>%
#   select(Field,Distance,EdgeCent,minShelter,Leafbee:NumLeafbee,Flower,Pollen,Treatment,Honeybee:NumHoneybee) %>%
#   mutate(FieldPlot=paste(Field,Distance,sep='_')) #Unique plot names (for random effect)

rm(bagdat) #Remove extra dataframe

# Save workspace
save.image("~/Projects/UofC/canola_yield_project/Seed field analysis/seedfieldData2016.RData")

# Functions ---------------------------------------------------------------



AICc=function(model,...){ #Second-order bias corrected AIC
  #require(plyr)
  corAIC=function(x){
    k=length(coef(x))
    n=length(residuals(x))
    return(AIC(x)+(2*k*k-1)/(n-k-1))
  }
  if(missing(...)){
    return(corAIC(model))
  }
  else{
    mlist=list(model,...)
    len=plyr::laply(mlist,function(x) length(residuals(x)))
    if(length(len)!=sum(len %in% len[1]))
      warning('Models use different # data points. Comparing AICc values is invalid')
    forms=sapply(mlist, function(x) paste(formula(x)))
    data.frame(row.names=paste(forms[2,],forms[1,],forms[3,]),AICc=plyr::laply(mlist,corAIC))
  }
}

deltaAIC=function(m1,m2,...,correct=TRUE){
  if(correct) a=AICc(m1,m2,...) else a=AIC(m1,m2,...)
  a$delta_AICc=round(a$AICc-min(a$AICc),4)
  a$w_i=round(exp(-0.5*a$delta_AICc)/sum(exp(-0.5*a$delta_AICc)),4) #Model Probabilities (Akaike weights)
  a
}

DIC=function(m1,...){ #Convenience function to extract DIC values from MCMCglmm models
  modlist=list(m1,...)
  dicvals=unlist(lapply(modlist,'[[',c('DIC')))
  fixeffs=as.character(lapply((lapply(modlist,'[[','Fixed')),'[[','formula'))
  raneffs=as.character(lapply((lapply(modlist,'[[','Random')),'[[','formula'))
  data.frame(Fixed=fixeffs,Random=raneffs,DIC=dicvals,deltaDIC=dicvals-min(dicvals))
}

plotFixedTerms=function(m1){ #Convenience function to plot 95% distributions for MCMCglmm fixed effects
  dat=as.data.frame(summary(m1)$solutions)
  dat$params=as.factor(rownames(dat))
  colnames(dat)[c(2,3)]=c('lwr','upr')
  ggplot(dat,aes(x=params,y=post.mean,ymax=upr,ymin=lwr))+geom_pointrange()+labs(x='Fixed Effects',y='Posterior Mean')+geom_hline(yintercept=0,colour='red')
}

zeroCheck=function(mod,dat) { #Posterior predictive checks for zero-inflation (MCMCglmm course notes pg. 104)
  require(MCMCglmm)
  require(ggplot2)
  try(if(is.null(mod$X)) stop("Needs fixed effect design matrix (use saveX=T in MCMCglmm arguments)"))
  nz=1:nrow(mod$Sol) #Number of samples to hold
  oz=sum(dat[as.character(mod$Fixed$formula[2])]==0) #Actual number of 0s in the dataset
  for (i in 1:nrow(mod$Sol)) {
    pred.l <- rnorm(nrow(dat), (mod$X %*% mod$Sol[i,])@x, sqrt(mod$VCV[i])) #Predicted mean
    nz[i] <- sum(rpois(nrow(dat), exp(pred.l)) == 0)} #Expected number under poisson
  qplot(nz,xlab='Expected number of zeros')+geom_vline(aes(xintercept=oz))
  #Expected zeros, given the overdispersed poisson model
  }

varComp=function(mod,units=T) { #Tabulates variance components of random effects for MCMCglmm models
  VCV=mod$VCV
  if(units==F) VCV=VCV[,c(1:ncol(VCV)-1)] #Strips out "units" column (for binomial models)
  round(data.frame(mode=posterior.mode(VCV/rowSums(VCV)),HPDinterval(VCV/rowSums(VCV))),3)
}

# Nectar production -------------------------------------------------------
ylab=expression(paste('Nectar(',mu,'L)',sep='')) #Nectar volume
ylabhr=expression(paste('Nectar(',mu,'L)/hr',sep='')) #Nectar volume/hr

#Nectar production over hours
ggplot(nectar,aes(TotalTime,Vol,col=Bay))+geom_point(position=position_jitter(width=0.2,height=0))+
  labs(x='Time (hrs)',y=ylab)+scale_colour_manual(values=c('red','blue'))+geom_smooth(method='lm',se=F)+
  ylim(0,1.2)

a=lm(Vol~TotalTime+TotalTime:Bay,data=nectar) #Nectar production over time
b=unname(c(coef(a)[2],coef(a)[2]+coef(a)[3])) #F = 0.04uL/hr, M = 0.07uL/hr

c=lm(Vol~TotalTime+Variety,data=nectar,subset=Bay=='F')

#Nectar production (uL/hr) over season
ggplot(nectar[with(nectar,is.finite(Vol/TotalTime)),], #Strips 0 hr measurements
       aes(Surveyed,Vol/TotalTime,col=Bay))+
  geom_point(position=position_jitter(width=40000,height=.01))+labs(x='Day',y=ylabhr)+
  scale_colour_manual(values=c('red','blue'))

#Nectar production with temperature
ggplot(nectar[with(nectar,is.finite(Vol/TotalTime)),], #Strips 0 hr measurements
       aes(Temp,Vol/TotalTime,col=Bay))+
  geom_point(position=position_jitter(width=0.2,height=0))+labs(x='Temperature',y=ylabhr)+
  scale_colour_manual(values=c('red','blue'))+geom_smooth(method='loess',span=1)

