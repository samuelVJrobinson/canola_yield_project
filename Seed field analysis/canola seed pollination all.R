#CODE USED IN ANALYSIS OF INSECT VISITATION, POLLEN DEPOSITION, AND YIELD IN SEED CANOLA FIELDS (2015+2016)

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
                panel.border=element_blank(),
                strip.text=element_text(size=15))                
theme_set(theme_classic()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

setwd("~/Projects/UofC/canola_yield_project/Seed field analysis")
# load("seedfieldDataAll.RData")

# Load 2015 data ----------------------------------------------------------

load("./seedfieldData2015.RData")

rm(lbees,hbees,flowers,figs)

fields2015=mutate(fields,Year=2015,Variety='3007C',Treatment=factor(Treatment,levels=c('control','double tent','double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees')))
plants2015=mutate(plants,Year=2015,Variety='3007C',Treatment=factor(Treatment,levels=c('control','double tent','double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees')))
pollen2015=mutate(pollen,Year=2015,Variety='3007C',Treatment=factor(Treatment,levels=c('control','double tent','double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees')))
seeds2015=mutate(seeds,Year=2015,Variety='3007C',Treatment=factor(Treatment,levels=c('control','double tent','double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees')))
behav2015=filter(hbeesWide,EdgeCent!='Center') %>% #Temporary dataframe with side/top, nectar/pollen working info
  select(Field:Bay,Surveyed,Treatment,SideNecVis,TopNecVis,TopPolVis)
survey2015=mutate(survey,Year=2015,Variety='3007C',Treatment=factor(Treatment,levels=c('control','double tent','double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees'))) %>%
  bind_cols(transmute(hbeesWide,hbeePol=TopPolVis,hbeeNec=SideNecVis+TopNecVis))

rm(fields,plants,hbeesWide,pollen,seeds,survey)

# Load 2016 data ----------------------------------------------------------

load("./seedfieldData2016.RData")

fields2016=mutate(fields,Year=2016,Treatment=factor(Treatment,levels=c('Control','Double tent','Double tent and bee'),labels=c('Control','Double tent','Double tent\n+bees')))
plants2016=mutate(plants,Year=2016,
                  Treatment=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees')))
seeds2016=mutate(seeds,Year=2016,
                 Treatment=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees')))
nectar2016=nectar
pollen2016=mutate(pollen,Year=2016,Treatment=factor(Treatment,levels=c('Control','Double tent','Double tent and bee'),labels=c('Control','Double tent','Double tent\n+bees')))
survey2016=mutate(survey,Year=2016,Treatment=factor(Treatment,levels=c('Control','Double tent','Double tent and bee'),labels=c('Control','Double tent','Double tent\n+bees')))
visitors2016=visitors %>%
  mutate(Treatment=factor(Treatment,levels=c('Control','Double tent','Double tent and bee'),labels=c('Control','Double tent','Double tent\n+bees')))

rm(fields,nectar,pollen,survey,visitors,plants,seeds)
rm(AICc,deltaAIC,DIC,plotFixedTerms,se,varComp,zeroCheck)

# Put 2015 and 2016 together ----------------------------------------------

#Field data
allFields=bind_rows(transmute(fields2015,Field,Area=Hectares,Variety,Treatment,Surveyed,Harvested,WindSp=Wind,AirTemp,RH,Year),
                    transmute(fields2016,Field,Area,Variety,Treatment,Surveyed,Harvested,WindSp,AirTemp=Temp,RH,Year)) %>%
  mutate(Field=factor(Field),Variety=factor(Variety),Surveyed=as.Date(Surveyed,format='%b %d, %Y'),
         Harvested=as.Date(Harvested,format='%b %d, %Y'))
rm(fields2015,fields2016)

#Survey data
allSurvey=bind_rows(transmute(survey2015,Field,Treatment,Variety,Distance,minDist=minShelter,Bay,EdgeCent,EdgeDir,StartTime,
                              TotalTime=as.double(TotalTime,units='mins'),FlDens,
                              lbee=Leafbee,hbee=Honeybee,hbeePol,hbeeNec,otherBee=Otherbee,hFly=Hoverfly,
                              AirTemp,WindSp=Wind,RH,Year,PlDens),
                    transmute(survey2016,Field,Treatment,Variety,Distance,minDist,Bay,EdgeCent=Edge.Cent,EdgeDir,StartTime,
                              TotalTime=10,FlDens,
                              lbee=Leafcutterbee,hbee=HoneybeePollen+HoneybeeNectar,hbeePol=HoneybeePollen,
                              hbeeNec=HoneybeeNectar,
                              otherBee=Otherbee,hFly=Hoverfly,
                              AirTemp=Temp,WindSp,RH,Year,PlDens=Plants.m.2)) %>%
  mutate(EdgeCent=sub('Center','Cent',EdgeCent),EdgeCent=sub('Cent','Center',EdgeCent)) %>%
  mutate(Field=factor(Field),Variety=factor(Variety),EdgeCent=factor(EdgeCent))
rm(survey2015,survey2016)

#Pollen data
allPollen=bind_rows(transmute(pollen2015,Field,Distance,Treatment,Variety,minDist=minShelter,EdgeCent,Pollen,Year),
                    transmute(pollen2016,Field,Distance,Treatment,Variety,minDist,EdgeCent=Edge.Cent,Pollen=Count,Year)) %>%
  mutate(FieldPlot=paste(Field,Distance,EdgeCent,sep='_')) %>%
  mutate(EdgeCent=sub('Center','Cent',EdgeCent),EdgeCent=sub('Cent','Center',EdgeCent)) %>%
  mutate(Field=factor(Field),Variety=factor(Variety),EdgeCent=factor(EdgeCent))
rm(pollen2015,pollen2016)

#Plant data
allPlants=bind_rows(transmute(plants2015,Year,Field,Distance,EdgeCent,EdgeDir,Variety,VegMass,SeedMass,Branch,Pods,Missing,AvPodCount,AvPodMass,lbee=Leafbee,hbee=Honeybee,TotalTime=5),
                    transmute(plants2016,Year,Field,Distance,EdgeCent,EdgeDir,Variety,VegMass,SeedMass,Branch,Pods,Missing,AvPodCount,AvPodMass,lbee=Leafcutterbee,hbee=HoneybeeNectar+HoneybeePollen,TotalTime=10)) %>%
  mutate(Field=factor(Field),EdgeCent=factor(ifelse(EdgeCent=='Cent','Center',EdgeCent))) %>%
  mutate(EdgeDir=factor(ifelse(nchar(EdgeDir)>2,substr(EdgeDir,1,1),EdgeDir))) %>%
  mutate(Variety=factor(Variety))
rm(plants2015,plants2016)

#Seed data
allSeeds=bind_rows(transmute(seeds2015,Year,Field,Distance,EdgeCent,EdgeDir,Variety,VegMass,SeedMass,Branch,Pods,Missing,Pod,PodCount,PodMass,lbee=Leafbee,hbee=Honeybee,TotalTime=5),
                   transmute(seeds2016,Year,Field,Distance,EdgeCent,EdgeDir,Variety,VegMass,SeedMass,Branch,Pods,Missing,Pod,PodCount,PodMass,lbee=Leafcutterbee,hbee=HoneybeeNectar+HoneybeePollen,TotalTime=10))
rm(seeds2015,seeds2016)

folder="C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\Figures"

# Save workspace
save.image("~/Projects/UofC/canola_yield_project/Seed field analysis/seedfieldDataAll.RData")

# Nectar production -------------------------------------------------------
#PREDICTIONS: 
#1) Males produce more nectar than females (TRUE)
#2) Nectar production differs by variety (FALSE)
#3) Nectar production changes across season (INCONCLUSIVE)
#4) Nectar production changes with temperature (TRUE)

ylab=expression(paste('Nectar(',mu,'L)',sep='')) #Nectar volume
ylabhr=expression(paste('Nectar(',mu,'L)/hr',sep='')) #Nectar volume/hr

#1) Nectar production over hours
p=
ggplot(nectar2016,aes(TotalTime,Vol,col=Bay))+geom_point(position=position_jitter(width=0.2,height=0))+
  labs(x='Time (hrs)',y=ylab)+scale_colour_manual(values=c('red','blue'))+geom_smooth(method='lm',se=T)+
  ylim(0,1.2)

ggplot(nectar2016,aes(TotalTime,Sugar,col=Bay))+geom_point(position=position_jitter(width=0.2,height=0))+
  labs(x='Time (hrs)',y='Sugar')+scale_colour_manual(values=c('red','blue'))+geom_smooth(method='lm',se=T)

ggplot(nectar2016,aes(TotalTime,Sugar/Vol,col=Bay))+geom_point(position=position_jitter(width=0.2,height=0))+
  labs(x='Time (hrs)',y='Sugar/Vol')+scale_colour_manual(values=c('red','blue'))+geom_smooth(method='lm',se=T)
  

a=lm(Vol~TotalTime+TotalTime:Bay,data=nectar2016) #Male flowers have significantly higher rate of nectar production (p<0.001)
b=unname(c(coef(a)[2],coef(a)[2]+coef(a)[3])) #F = 0.04uL/hr, M = 0.07uL/hr

ggsave(paste(folder,'nectarProdBay.png',sep='\\'),p,width=8,height=6)


#2) Nectar production with variety
ggplot(nectar2016,aes(TotalTime,Vol,col=Variety))+geom_point(position=position_jitter(width=0.2,height=0))+
  labs(x='Time (hrs)',y=ylab)+scale_colour_manual(values=c('darkorange','purple'))+geom_smooth(method='lm',se=T)+
  ylim(0,1.2)

c=lm(Vol~TotalTime+Variety,data=nectar2016) #No difference in nectar production b/w varieties (p=0.51)

#3) Nectar production (uL/hr) over season
ggplot(nectar2016[with(nectar2016,is.finite(Vol/TotalTime)),], #Strips 0 hr measurements
       aes(Surveyed,Vol/TotalTime,col=Bay))+
  geom_point(position=position_jitter(width=40000,height=.01))+labs(x='Day',y=ylabhr)+
  scale_colour_manual(values=c('red','blue'))

#4) Nectar production with temperature
p=ggplot(nectar2016[with(nectar2016,is.finite(Vol/TotalTime)),], #Strips 0 hr measurements
       aes(Temp,Vol/TotalTime,col=Bay))+
  geom_point(position=position_jitter(width=0.2,height=0))+labs(x=expression(paste('Temperature('~degree~C,')',sep='')),y=ylabhr)+
  scale_colour_manual(values=c('red','blue'))+geom_smooth(method='loess',span=1)

ggsave(paste(folder,'nectarProdTemp.png',sep='\\'),p,width=8,height=6)

# Honeybee visitation -----------------------------------------------------
#PREDICTIONS: 
#1) Visits decrease with distance from the colony (TRUE for 2015)
#2) Visits are lower in male bays (Visit~Distance intercept changes for male/female bays) (FALSE)
#3) If distance-based competition operates differently for pollen/nectar, slope of visit:distance changes with bay (FALSE)
#4) Visits decrease with distance to leafcutter shelter (INCONCLUSIVE - DIFFICULT TO MODEL)
#5) Visit~leafcutter shelter relationship changes with Treatment (competition is different across treatments)
#6) For female bays, visits decrease in the center of the bay (FALSE - opposite is true), and more so in the center of the field (FALSE)
#7) The majority of visitors to the female bay were nectar visitors, and vice-versa
#8) Edge visitors will be a mixture of pollen/nectar foragers

#1) Visits~Distance
#2) Visits~Bay
#3) Visits~ Distance + Distance:Bay
ggplot(allSurvey,aes(Distance,(hbee*60)/TotalTime))+geom_point()+facet_wrap(~Year,ncol=1)+
  labs(x='Distance to honeybees(m)',y='Visits/hr')+geom_smooth(method='lm',se=F) #Doesn't look like much difference b/w 

#zeroCheck indicates that for both 2015 and 2016 the count data is zero-inflated
#Prior
pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        #G = list(G1 = list(V = diag(2), nu = 0.002))) #Old prior (unused)
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

#For 2015 data: Dist:Bay(p=0.58,0.94) and Bay(p=0.36,0.97) not significant; Distance significant (p=0.83,0.01)
temp2015=filter(allSurvey,Year==2015,EdgeCent=='Edge')

hbeeVis2015=MCMCglmm(hbee~trait*(Distance),
                     rcov=~idh(trait):units,
                     random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                     family='zapoisson',prior=pr,
                     nitt=50000+5000,thin=50,burnin=5000,verbose=F,
                     data=temp2015)
pred2015=with(temp2015,data.frame(predict(hbeeVis2015,interval='confidence')*12,Distance,Year=2015)) %>%
  distinct() #Visits/hr predictions for 2015


#For 2016 data: Dist:Bay(p=0.23,0.42), Bay(p=0.91,0.50), and Dist(p=0.28,0.35) not significant
temp2016=filter(allSurvey,Year==2016,EdgeCent=='Edge')
hbeeVis2016=MCMCglmm(hbee~trait*(Distance),
                     rcov=~idh(trait):units,
                     random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                     family='zapoisson',prior=pr,
                     nitt=100000+5000,thin=100,burnin=5000,verbose=F,
                     data=temp2016)
pred2016=with(temp2016,data.frame(predict(hbeeVis2016,interval='confidence')*6,Distance,Year=2016)) %>%
  distinct() #Visits/hr predictions for 2016

pred=bind_rows(pred2015,pred2016) #Predictions for both years

p=ggplot(allSurvey)+geom_ribbon(data=pred,aes(x=Distance,ymax=upr,ymin=lwr),alpha=0.3)+
  geom_point(aes(Distance,(hbee*60)/TotalTime),position=position_jitter(width=5))+
  geom_line(data=pred,aes(x=Distance,y=fit),size=1)+
  facet_wrap(~Year,ncol=1)+labs(x='Distance to honeybees(m)',y='Visits/hr')+
  ylim(0,760)

ggsave(paste(folder,'hbeeVisDist.png',sep='\\'),p,width=8,height=6)

#4) Visits~minDist
#5) Visit~minDist*Treatment

#Not an obvious trend...
# ggplot(allSurvey,aes(minDist,hbee))+geom_point()+facet_wrap(~Treatment+Year,ncol=2)+
#   geom_smooth(method='gam',formula=y~s(x,k=5),se=F)
# 
# temp2015=filter(allSurvey,Year==2015,EdgeCent=='Edge')
# hbeeShelter2015=gam(hbee~s(minDist,by=Treatment),data=temp2015) #Spline by Treatment doesn't really look that significant. Not enough df to do Distance by Treatment interaction.
# anova(hbeeShelter2015)
# ggplot(temp2015,aes(minDist,hbee,col=Treatment))+geom_point()+geom_smooth(method='lm',se=F)
# 
# hbeeShelter2015=MCMCglmm(hbee~minDist*Treatment,family='poisson',data=temp2015,verbose=F,nitt=100000+5000,burnin=5000,thin=100)
# 
# pred2015=unique(with(temp2015,data.frame(minDist,Treatment,predict(hbeeShelter2015,interval='confidence'))))
# 
# ggplot(temp2015,aes(minDist,hbee))+geom_point()

#6) Visits~Distance*EdgeCent for female bays

temp=filter(allSurvey,Distance!=20,Distance!=100,Distance!=250,Bay=='F') %>%
  mutate(Distance=factor(ifelse(Distance==5,'Near','Far')))%>%
  mutate(Distance=factor(Distance,levels=c('Near','Far')))
ggplot(temp,aes(Distance,hbee*(60/TotalTime),fill=EdgeCent))+geom_boxplot()+facet_wrap(~Year,ncol=1)+
  labs(y='Visits/hr',fill='Bay Position') #Looks like low visits to the center of the field during 2015

temp2015=filter(temp,Year==2015)
#zeroCheck indicates that 2015 data isn't zero-inflated, but it's getting very large CIs and predicted means
pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        #G = list(G1 = list(V = diag(2), nu = 0.002))) #Old prior (unused)
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

#2015 data - Distance:EdgeCent (p=0.59,0.79) is not significant; Distance (p=0.05,<0.001) and EdgeCent (0.01,0.18) are significant, ~139 visits less in center, ~34 less at far plot
hbeeVisFbay2015=MCMCglmm(hbee~trait*(Distance*EdgeCent),
         rcov=~idh(trait):units,
         random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
         family='zapoisson',prior=pr,
         nitt=100000+5000,thin=100,burnin=5000,verbose=F,
         data=temp2015)
pred2015= with(temp2015,data.frame(predict(hbeeVisFbay2015b,interval='confidence')*12,Distance,EdgeCent,Year=2015)) %>%
  distinct() #Visits/hr predictions for 2015


temp2016=filter(temp,Year==2016)
#zeroCheck indicates that 2016 data is zero-inflated
pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

#2016 data - Distance:EdgeCent (p=0.09,0.15), Distance (p=0.26,0.30) and EdgeCent (p=0.14,0.67) are not significant
hbeeVisFbay2016=MCMCglmm(hbee~trait*(Distance*EdgeCent),
                         rcov=~idh(trait):units,
                         random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                         family='zapoisson',prior=pr,
                         nitt=200000+5000,thin=200,burnin=5000,verbose=F,
                         data=temp2016)

pred2016=with(temp2016,data.frame(predict(hbeeVisFbay2016,interval='confidence')*6,Distance,EdgeCent,Year=2016)) %>%
  distinct() #Visits/hr predictions for 2016

pred=bind_rows(pred2015,pred2016)

p=ggplot(temp,aes(Distance,hbee*(60/TotalTime),colour=EdgeCent))+
  geom_point(position=position_jitterdodge())+
  facet_wrap(~Year,ncol=2)+labs(y='Visits/hr',colour='Bay\nPosition')+
  scale_colour_manual(values=c('purple','forestgreen'))+
  geom_pointrange(data=pred,aes(x=Distance,y=fit,ymax=upr,ymin=lwr,group=EdgeCent),col='black',position=position_dodge(width=0.75))
ggsave(paste(folder,'hbeeVisFbay.png',sep='\\'),p,width=8,height=6)

#7) The majority of visitors to the female bay will be nectar visitors, and vice-versa (TRUE, but only for females)
#8) Edge visitors will be a mixture of pollen/nectar foragers (TRUE)

temp=ungroup(allSurvey) %>%
  select(Field,Distance,Bay,EdgeCent,TotalTime,hbeePol,hbeeNec,Year) %>%
  gather('Behav','count',hbeePol:hbeeNec) %>%
  mutate(Behav=factor(ifelse(Behav=='hbeePol','Pollen','Nectar')))



pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

#Temporary mixed effects model
a=lme(log(count+1)~Bay*Behav,random=~1|Field,data=filter(temp,EdgeCent!='Center',Year==2015),method='REML')
b=lme(log(count+1)~Bay*Behav,random=~1|Field,data=filter(temp,EdgeCent!='Center',Year==2016),method='REML')
#Fitted values (no CIs)
pred=data.frame(fit=c(exp(c(sum(fixef(a)[1]),sum(fixef(a)[1:2]),sum(fixef(a)[1:3]),sum(fixef(a)[1:4])))*12,exp(c(sum(fixef(b)[1]),sum(fixef(b)[1:2]),sum(fixef(b)[1:3]),sum(fixef(b)[1:4])))*6),
                Bay=factor(rep(c('F','M','F','M'),2)),
                Behav=factor(rep(rep(c('Nectar','Pollen'),each=2),2)),
                Year=rep(c(2015,2016),each=4))

p=ggplot(filter(temp,EdgeCent!='Center'),aes(Bay,count*(60/TotalTime),col=Behav))+
  geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.75))+
  facet_wrap(~Year)+
  scale_colour_manual(values=c('cyan','orange'))+labs(y='Visits/hr',col='Forager')+
  geom_point(data=pred,aes(x=Bay,y=fit,group=Behav),col='black',position=position_dodge(width=0.75),size=2)
ggsave(paste(folder,'hbeeBehav.png',sep='\\'),p,width=8,height=6)

#Temporary mixed effects model
a=lme(log(count+1)~EdgeCent*Behav,random=~1|Field,data=filter(temp,Bay!='M',Year==2015),method='REML')
b=lme(log(count+1)~EdgeCent*Behav,random=~1|Field,data=filter(temp,Bay!='M',Year==2016),method='REML')
#Fitted values (no CIs)
pred=data.frame(fit=c(exp(c(sum(fixef(a)[1]),sum(fixef(a)[1:2]),sum(fixef(a)[1:3]),sum(fixef(a)[1:4])))*12,exp(c(sum(fixef(b)[1]),sum(fixef(b)[1:2]),sum(fixef(b)[1:3]),sum(fixef(b)[1:4])))*6),
                EdgeCent=factor(rep(c('Center','Edge','Center','Edge'),2)),
                Behav=factor(rep(rep(c('Nectar','Pollen'),each=2),2)),
                Year=rep(c(2015,2016),each=4))

p=ggplot(filter(temp,Bay!='M'),aes(EdgeCent,count*(60/TotalTime),col=Behav))+
  geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.75))+
  facet_wrap(~Year)+
  scale_colour_manual(values=c('cyan','orange'))+labs(x='Bay Position',y='Visits/hr',col='Forager')+
  geom_point(data=pred,aes(x=EdgeCent,y=fit,group=Behav),col='black',position=position_dodge(width=0.75),size=2)
ggsave(paste(folder,'hbeeBehavFbay.png',sep='\\'),p,width=8,height=6)

#Temporary graph of worker types
temp=mutate(behav2015,SidePolVis=0) %>% #Adds side working pollen forager column
  gather('Clade','Visits',SideNecVis:SidePolVis) %>%
  mutate(Behav=ifelse(grepl('Side',Clade),'Side','Top'),Type=ifelse(grepl('Nec',Clade),'Nectar','Pollen')) %>%
  group_by(Bay,Type,Behav) %>% #,Distance
  summarize(sumVisits=sum(Visits),avgVisits=mean(Visits))

p=ggplot(temp,aes(x=Behav,y=sumVisits,fill=Type))+geom_bar(stat='identity',position=position_dodge())+
  facet_wrap(~Bay)+
  labs(y='Total Visits',fill='Forager\nType',x='Behaviour')+
  scale_fill_manual(values=c('forestgreen','orange'))
ggsave(paste(folder,'hbeeVisType.png',sep='\\'),p,width=8,height=6)



# Leafcutter visitation ---------------------------------------------------


#PREDICTIONS: 
#1) Visitation is higher closer to shelters (TRUE, but only for 2015 data)
#2) Slope will be lower at Double Tent treatment (less bees per tent) (appears TRUE, find way of testing this)
#3) Visitation will be higher in male bay (FALSE)
#4) Visitation will be higher at edge of field (FALSE)
#5) Visitation will be higher at edge of bay, not at center

#Incorporate Riley's extra seed field data (mainly for leafcutter distance)
#Data from 2015
temp2015=readRDS("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\rileyExtra2015.rds")
#Data from 2016
temp2016=readRDS("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\rileyExtra2016.rds")

temp2015=temp2015 %>%
  transmute(year=2015,site,date=as.Date(as.character(date),format='%m/%d/%Y'),treatment,ldist,hdist,
            visit=Visit,bay=Bay,pollinator) #%>%
# filter(treatment!='Double tent')

temp2016=temp2016 %>%
  transmute(year=2016,site=Site,date=as.character(Julian),ldist=Distance.from.leafcutter.tent..m.,
            hdist=Distance.from.honey.bee.hive..m.,visit=Visit,bay=Bay,pollinator=Pollinator) %>%
  mutate(date=as.Date(date,format='%d-%b-%y')) #%>%
# filter(treatment!='Double tent')

temp=bind_rows(temp2015,temp2016) %>%
  mutate(Riley=T,time=10)
rm(temp2015,temp2016)

temp=bind_rows(temp,
               gather(allSurvey,Pollinator,Visits,hbee:lbee) %>%
                 # filter(Treatment!='Double tent') %>%
                 transmute(year=Year,site=Field,date=as.Date(StartTime),treatment=Treatment,
                           ldist=minDist,hdist=Distance,visit=Visits,bay=Bay,pollinator=Pollinator,Riley=F,time=TotalTime)) %>%
  mutate(pollinator=ifelse(pollinator=='ALCB'|pollinator=='lbee','Leafcutter bee','Honey bee'))

visGam=gam(visit~s(ldist)+
             s(ldist,by=as.numeric(pollinator=='Leafcutter bee'))+
             offset(log(time)),
           family='nb',data=temp,subset=ldist<100)
summary(visGam) #Both leafcutter and honeybee splines are significant.

newdata=expand.grid(time=60,ldist=1:100,pollinator=unique(temp$pollinator))
pred=data.frame(newdata,predict(visGam,newdata=newdata,type='link',se.fit=T)) %>%
  mutate(upr=exp(fit+se.fit*1.96),lwr=exp(fit-se.fit*1.96),fit=exp(fit)) %>%
  select(-se.fit)

p=pred %>%
  ggplot(aes(ldist,fit,group=pollinator))+
  geom_ribbon(aes(y=NULL,ymin=lwr,ymax=upr,fill=pollinator),alpha=0.3)+
  geom_line(aes(col=pollinator),size=1)+
  geom_point(data=temp,aes(ldist,visit*(60/time),col=pollinator),size=1)+
  scale_colour_manual(values=c('forestgreen','orange'))+
  scale_fill_manual(values=c('forestgreen','orange'))+
  labs(x='Distance to leafcutter shelter (m)',y='Visits/hr',col='Visitor',fill='Visitor')+
  xlim(0,100)+ylim(0,700)

ggsave(paste(folder,'lbeeShelterDist.png',sep='\\'),p,width=8,height=6)

#1) Visits ~ minDist
#2) Slope will be lower at Double Tent treatment (less bees per tent)
require(mgcv)
#2015 data - Dist:Bay (0.11), Bay (0.822), and Distance (0.74) are not significant under glmm framework - s(minDist) is not significant (p=0.43)
temp2015=filter(allSurvey,Year==2015)

lbeeVis2015=gam(lbee~s(minDist,by=Treatment),family='nb',data=temp2015)
summary(lbeeVis2015)
anova(lbeeVis2015,lbeeVis2015b,test='Chisq') #Seems like using "by" argument is better. Not sure how well this test performs with GAMs
gam.check(lbeeVis2015)

#2016 data - Dist:Bay (0.53), Bay (0.44), and Distance (0.87) are not significant under glmm framework - s(minDist) is significant (p=0.002)
temp2016=filter(allSurvey,Year==2016)

#Overdispersion: 15.6 (high), so using quasipoisson
# sum(resid(lbeeVis2016,type='pearson')^2)/lbeeVis2016$df.res
lbeeVis2016=gam(lbee~s(minDist,bs='cr')+
                  s(minDist,by=as.numeric(Treatment=='Double tent'),bs='cr')+
                  s(minDist,by=as.numeric(Treatment=='Double tent\n+bees'),bs='cr')+Treatment,
                family=quasipoisson,data=temp2016)
summary(lbeeVis2016)

summary(lbeeVis2016) #All terms are "significant".
gam.check(lbeeVis2016)

#No pattern in 2015 because of distance, pattern exists in 2016

p=ggplot(allSurvey,aes(minDist,(lbee*60)/TotalTime,col=Treatment))+geom_point()+facet_wrap(~Year,ncol=1)+
  labs(x='Distance to shelter (m)',y='Visits/hr')+geom_smooth(method='gam',method.args=list(family=quasipoisson),formula=y~s(x,bs='cr'),se=F,size=1)+
  ylim(NA,600)+scale_colour_manual(values=c('blue','darkorange','red'))
ggsave(paste(folder,'lbeeShelterTreat.png',sep='\\'),p,width=8,height=6)

p=ggplot(allSurvey,aes(minDist,(lbee*60)/TotalTime))+geom_point()+facet_wrap(~Year,ncol=1)+
  labs(x='Distance to shelter (m)',y='Visits/hr')+geom_smooth(method='gam',formula=y~s(x),se=F,size=1)+
  ylim(NA,600)
ggsave(paste(folder,'lbeeShelter.png',sep='\\'),p,width=8,height=6)

#3) lbee~Bay
#4) lbee~Distance*Bay 
pr=list(R = list(V = 1, nu = 0.02),
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for poisson process
               
#2015 data - Bay:Distance is not significant, nor is Bay (0.87), nor Distance (0.74). When run as a factor, Distance100 is marginally lower than other distances (p=0.07)
temp2015=filter(allSurvey,Year==2015,EdgeCent=='Edge')
lbeeVisBay2015=MCMCglmm(lbee~Distance*Bay,random=~Field,family='poisson',data=temp2015,verbose=F,thin=25,prior=pr,nitt=50000+5000,burnin=5000)
summary(lbeeVisBay2015)

ggplot(allSurvey,aes(Distance,lbee*(60/TotalTime),col=Bay))+geom_point()+
  facet_wrap(~Year,ncol=1)+scale_colour_manual(values=c('red','blue'))+ylim(0,600)

#2016 data - Bay:Distance is not significant, nor is Bay (0.48), or Distance (0.84)
temp2016=filter(allSurvey,Year==2016,EdgeCent=='Edge')
lbeeVisBay2016=MCMCglmm(lbee~Distance*Bay,random=~Field,family='poisson',data=temp2016,verbose=F,thin=25,prior=pr,nitt=50000+5000,burnin=5000)
summary(lbeeVisBay2016)

#5) lbee~EdgeCent*Distance
temp=filter(allSurvey,Distance!=20,Distance!=100,Distance!=250,Bay=='F') %>%
  mutate(Distance=factor(ifelse(Distance==5,'Near','Far'))) %>%
  mutate(Distance=factor(Distance,levels=c('Near','Far')))

pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

#2015 data - Treatment, Distance:Edgecent (p=0.18,0.85), Distance(0.40,0.78), and EdgeCent (0.97,0.39) are not significant,
temp2015=filter(temp,Year==2015)

lbeeVisFbay2015=MCMCglmm(lbee~trait*(Treatment),
         rcov=~idh(trait):units,
         random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
         family='zapoisson',prior=pr,
         nitt=200000+20000,thin=200,burnin=5000,verbose=F,
         data=temp2015)

summary(lbeeVisFbay2015)
pred2015= with(temp2015,data.frame(predict(lbeeVisFbay2015,interval='confidence')*12,Distance,EdgeCent,Year=2015)) %>%
  distinct() #Visits/hr predictions for 2015

#2016 data - Treatment, Distance:EdgeCent (p=0.67,0.79), Distance (0.16,0.32), and EdgeCent (p=0.57,0.99) are not significant
temp2016=filter(temp,Year==2016)

lbeeVisFbay2016=MCMCglmm(lbee~trait*(Treatment),
                         rcov=~idh(trait):units,
                         random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                         family='zapoisson',prior=pr,
                         nitt=200000+20000,thin=200,burnin=5000,verbose=F,
                         data=temp2016)
summary(lbeeVisFbay2016)
pred2016= with(temp2016,data.frame(predict(lbeeVisFbay2016,interval='confidence')*6,Distance,EdgeCent,Year=2016)) %>%
  distinct() #Visits/hr predictions for 2016

pred=bind_rows(pred2015,pred2016)

p=ggplot(temp,aes(Distance,lbee*(60/TotalTime),col=EdgeCent))+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.width=0.4))+facet_wrap(~Year,ncol=2)+
  labs(y='Visits/hr',col='Bay\nPosition')+ylim(0,225)+
  scale_colour_manual(values=c('purple','forestgreen'))+
  geom_pointrange(data=pred,aes(x=Distance,y=fit,ymax=upr,ymin=lwr,group=EdgeCent),col='black',position=position_dodge(width=0.75))
ggsave(paste(folder,'lbeeVisFbay.png',sep='\\'),p,width=8,height=6)

#Overall comparison - lumping both years

#Multiplies the 2015 survey counts by 2 (to account for smaller observation time)
lbeetemp=with(allSurvey,ifelse(TotalTime==5,lbee*2,lbee))

#ZAPoisson version
pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

#Treatment:Year (none significant, DoubleTent+bee:Year2016 = 0.13), Year (p=0.04,0.87), Treatment=(none significant)

allSurvey$Treatment=factor(allSurvey$Treatment,levels=c('Control','Double tent','Double tent+bees'))

lbeeVisFbay=MCMCglmm(lbeetemp~trait*(Treatment*as.factor(Year)),
                     rcov=~idh(trait):units,
                     random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                     family='zapoisson',prior=pr,
                     nitt=200000+20000,thin=200,burnin=5000,verbose=F,
                     data=allSurvey)
#Poisson version
pr=list(R = list(V = diag(1), nu = 0.02),
        G=list(G1=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000)))

lbeeVisFbay=MCMCglmm(lbeetemp~Treatment*as.factor(Year),
                     random=~Field,
                     family='poisson',prior=pr,
                     nitt=50000+20000,thin=50,burnin=20000,verbose=F,
                     data=allSurvey)
summary(lbeeVisFbay)

pred=with(allSurvey,data.frame(Year,Treatment,predict(lbeeVisFbay,interval='confidence')*6)) %>%
  distinct()

p=ggplot(allSurvey,aes(Treatment,lbee*(60/TotalTime),col=as.factor(Year)))+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.width=0.4))+ylim(0,225)+
  labs(y='Visits/hr',col='Year')+
  scale_colour_manual(values=c('darkorange','red'))+
  geom_pointrange(data=pred,aes(x=Treatment,y=fit,ymax=upr,ymin=lwr,group=as.factor(Year)),col='black',position=position_dodge(width=0.75))
ggsave(paste(folder,'lbeeVisTreat.png',sep='\\'),p,width=8,height=6)

#Lumped 
ggplot(allSurvey,aes(minDist,lbee))+geom_point()


# Pollen deposition -------------------------------------------------------

#PREDICTIONS: 
temp=filter(allPollen,EdgeCent!='Center',!is.na(Pollen)) %>%
  mutate(Treatment=factor(Treatment,levels=c('Control','Double tent','Double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees')))

pr=list(R = list(V = diag(1), nu = 0.02), 
        G=list(G1=list(V=diag(2), nu=1, alpha.mu=c(0,0), alpha.V=diag(2)*1000)))
               
polDistAll=MCMCglmm(Pollen~Distance*factor(Year),
                    random=~idh(1+Distance):Field,
                    family='poisson',prior=pr,
                    #nitt=200000+30000,thin=200,burnin=30000,verbose=F,
                    data=temp)

pred=with(temp,data.frame(round(predict(polDistAll,interval='confidence'),2),Distance,Year)) %>%
  distinct()

p=ggplot(pred,aes(Distance,fit))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  facet_wrap(~Year,ncol=1)+
  geom_line(size=1)+geom_point(data=temp,aes(x=Distance,y=Pollen))+ylim(0,200)+
  labs(x='Distance from honeybees(m)',y='Pollen grains/stigma')
ggsave(paste(folder,'polDist.png',sep='\\'),p,width=8,height=6)


#1) Pollen deposition decreases with distance from bee colony
#4) Pollen changes with pollination treatment 
pr=list(R = list(V = diag(1), nu = 0.02),
        G=list(G1=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000)))


#2015 data - Distance (p<0.001) and Distance:Treatment Doubletent+bees (p=0.05) are significant 
temp2015=filter(temp,Year==2015)
polDist2015=MCMCglmm(Pollen~Distance*Treatment,random=~Field,family='poisson',data=temp2015,prior=pr,verbose=F)
summary(polDist2015)
pred2015=with(temp2015,data.frame(predict(polDist2015,interval='confidence'),Distance,Treatment,Year=2015)) %>%
  distinct()

#2016 data - Distance is significant (p<0.001)
temp2016=filter(temp,Year==2016)

#ZAPoisson version - Distance is significant (p=0.002), Distance:Doubletent, Distance:Doubletent+bees
pr=list(R = list(V = diag(2), nu = 0.02, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

polDist2016=MCMCglmm(Pollen~trait+(Distance*Treatment),
                     rcov=~idh(trait):units,
                     random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                     family='zapoisson',prior=pr,
                     nitt=200000+30000,thin=200,burnin=30000,verbose=F,
                     data=temp2016)
pred2016=with(temp2016,data.frame(predict(polDist2016,interval='confidence'),Distance,Treatment,Year=2016)) %>%
  distinct()
summary(polDist2016)

#Poisson version
# pr=list(R = list(V = diag(1), nu = 0.02),
#         G=list(G1=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000),
#                G2=list(V=1, nu=0.01, alpha.mu=0, alpha.V=1000)))
# polDist2016b=MCMCglmm(Pollen~Distance*Treatment,random=~Field+FieldPlot,family='poisson',data=temp2016,
#                      prior=pr,verbose=F,nitt=50000+5000,thin=50,burnin=5000)
# summary(polDist2016b)
# pred2016b=with(temp2016,data.frame(predict(polDist2016b,interval='confidence'),Distance,Treatment,Year=2016)) %>%
#   distinct()

pred=bind_rows(pred2015,pred2016) %>%
  arrange(Year,Treatment)

p=ggplot(pred,aes(x=Distance,y=fit,col=Treatment))+geom_ribbon(aes(ymax=upr,ymin=lwr,fill=Treatment,col=NULL),alpha=0.3)+
  geom_line(size=1)+facet_wrap(~Year,ncol=1)+
  scale_colour_manual(values=c('blue','darkorange','red'))+
  scale_fill_manual(values=c('blue','darkorange','red'))+
  geom_point(data=temp,aes(y=Pollen),position=position_jitterdodge(dodge.width=7,jitter.width=0.2))+
  labs(y='Pollen count')+ylim(0,250)
ggsave(paste(folder,'polDistTreat.png',sep='\\'),p,width=8,height=6)

#2) Pollen is lower at center of F bay
temp=filter(allPollen,Distance!=20,Distance!=100,Distance!=250,!is.na(Pollen)) %>%
  mutate(Distance=factor(ifelse(Distance==5,'Near','Far'))) %>%
  mutate(Distance=factor(Distance,levels=c('Near','Far')))

#2015 data
pr=list(R = list(V = diag(1), nu = 0.02),G=list(G1=list(V=1, nu=0.01)))
temp2015=filter(temp,Year==2015)
#2015 data - Dist:EdgeCent (p=0.53) and EdgeCent (0.13) are not significant; Treatment (and Treatment:Distance:EdgeCent combinations) are not significant; Distance is significant (p<0.001)
polDistFbay2015=MCMCglmm(Pollen~Distance*EdgeCent,random=~Field,data=temp2015,family='poisson',prior=pr,
                         nitt=50000+5000,burnin=5000,thin=50,verbose=F)

summary(polDistFbay2015)
pred2015=with(temp2015,data.frame(predict(polDistFbay2015,interval='confidence'),Distance,EdgeCent,Year=2015)) %>%distinct()

#2016 data - Distance:EdgeCent (p=0.09) is marginally significant; EdgeCent (p<0.001) and Distance (p<0.001) are significant
temp2016=filter(temp,Year==2016)
polDistFbay2016=MCMCglmm(Pollen~Distance*EdgeCent,random=~Field,data=temp2016,family='poisson',prior=pr,
                         nitt=50000+5000,burnin=5000,thin=50,verbose=F)
summary(polDistFbay2016)


pred2016=with(temp2016,data.frame(predict(polDistFbay2016,interval='confidence'),Distance,EdgeCent,Year=2016)) %>%distinct()

pred=bind_rows(pred2015,pred2016)

p=ggplot(temp,aes(Distance,Pollen,col=EdgeCent))+facet_wrap(~Year,ncol=1)+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.width=0.4))+facet_wrap(~Year,ncol=2)+
  labs(y='Pollen count',col='Bay\nPosition')+ylim(0,250)+
  scale_colour_manual(values=c('purple','forestgreen'))+
  geom_pointrange(data=pred,aes(x=Distance,y=fit,ymax=upr,ymin=lwr,group=EdgeCent),col='black',position=position_dodge(width=0.75))
ggsave(paste(folder,'polDistFbay.png',sep='\\'),p,width=8,height=6)

#3) Pollen is lower away from leafcutter shelters
temp=filter(allPollen,!is.na(Pollen))

#2015 data
temp2015=filter(temp,Year==2015)

polDistShelter2015=gam(Pollen~s(minDist,by=Treatment),family='nb',data=temp2015)
anova(polDistShelter2015)
gam.check(polDistShelter2015)

#2016 data
temp2016=filter(temp,Year==2016)

polDistShelter2016=gam(Pollen~s(minDist,by=Treatment),family='nb',data=temp2016)
anova(polDistShelter2016)
gam.check(polDistShelter2016)

p=ggplot(temp,aes(minDist,Pollen,col=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees'))))+geom_point()+facet_wrap(~Year,ncol=1)+labs(y='Pollen count',x='Distance to shelter',col='Treatment')+
  geom_smooth(method='gam',formula=y~s(x,k=5),se=F)+ylim(0,300)+
  scale_colour_manual(values=c('blue','darkorange','red'))
ggsave(paste(folder,'polDistShelter.png',sep='\\'),p,width=8,height=6)

#Both years together
temp=filter(allPollen,!is.na(Pollen),Treatment!='Double tent')

polDistShelterAll=gam(Pollen~s(minDist),family='nb',data=temp)
summary(polDistShelterAll)

newdata=expand.grid(minDist=min(temp$minDist):max(temp$minDist))
pred=data.frame(newdata,predict(polDistShelterAll,newdata=newdata,se.fit=T)) %>%
  mutate(upr=fit+se.fit,lwr=fit-se.fit) %>%
  select(-se.fit) %>%
  mutate(fit=exp(fit),upr=exp(upr),lwr=exp(lwr))

ggplot(pred,aes(minDist,fit))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+geom_line()+
  geom_point(data=temp,aes(minDist,Pollen))+
  labs(x='Distance to shelter (m)',y='Pollen grains per stigma')+
  ylim(NA,150)

#Trying nlm using Micalis-menten model
max(temp$Pollen)
Vm1=getInitial(Pollen~SSmicmen(minDist,578,1),data=temp)[1]
K1=getInitial(Pollen~SSmicmen(minDist,578,1),data=temp)[2]
mod1=nls(Pollen~SSmicmen(minDist,Vm1,K1),data=temp)

newdata=expand.grid(minDist=5:80)
newdata$fit=predict(mod1,newdata=newdata)
ggplot(newdata,aes(minDist,fit))+geom_line()+
  geom_point(data=temp,aes(minDist,Pollen))+
  labs(x='Distance to shelter (m)',y='Pollen grains per stigma')+
  ylim(NA,150)

#Looks alright, but there's still not very much data close to the shelter
  


p=ggplot(pred,aes(minDist,fit,group=Treatment))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=Treatment),alpha=0.3)+
  geom_line(aes(col=Treatment),size=1)+
  geom_point(data=temp,aes(y=Pollen,col=Treatment),size=1)+
  scale_colour_manual(values=c('blue','darkorange','red'))+
  scale_fill_manual(values=c('blue','darkorange','red'))+
  labs(x='Distance to shelter (m)',y='Pollen grains per stigma')+
  ylim(0,300) #Weird tails to the distribution
ggsave(paste(folder,'polDistShelterAll.png',sep='\\'),p,width=8,height=6)

#Both years together, with Honeybee distance instead of treatment
#Note: splines don't really work that well for this data because of the limited number of samples very close to the shelter. Trying nlm instead  with Michalis-Menton model 
Vm1=getInitial(Pollen~SSmicmen(minDist,max(Pollen),1),data=temp)[1]
K1=getInitial(Pollen~SSmicmen(minDist,max(Pollen),1),data=temp)[2]
mod1=nls(Pollen~SSmicmen(minDist,Vm1,K1),data=temp)

#Could also try: Asymptotic Regression Model (SSasymp)
# Vm2=getInitial(Pollen~SSasymp(minDist,20,200,-3),data=temp)[1]
# K2=getInitial(Pollen~SSmicmen(minDist,max(Pollen),1),data=temp)[2]
# mod2=nls(Pollen~SSmicmen(minDist,Vm2,K2),data=temp)

newdata=data.frame(minDist=5:88)
newdata$fit=predict(mod1,newdata=newdata,se.fit=T)

ggplot(newdata,aes(minDist,fit))+geom_line(size=1)+
  geom_point(data=temp,aes(minDist,Pollen),size=1)+
  ylim(NA,150)+
  labs(x='Distance from shelter (m)',y='Pollen per stigma',col='Distance from\nHoney bees(m)')





# Pod production ----------------------------------------------------------
#PREDICTIONS: 
#1) Floral success declines with distance from beehives/leafcutter tents

temp=mutate(plants2015,success=Pods/(Pods+Missing))

#Leafcutter distances
a=gam(cbind(Pods,Missing)~Distance+s(minShelter,by=as.factor(Distance),k=3),family='quasibinomial',data=temp)
summary(a)
newdata=expand.grid(Distance=c(5,20,100,400),minShelter=5:75)

#Looks weird
data.frame(newdata,fit=predict(a,newdata=newdata,type='response')) %>%
  ggplot(aes(minShelter,fit,col=as.factor(Distance)))+geom_line()+
  geom_point(data=plants2015,aes(minShelter,Pods/(Pods+Missing),col=as.factor(Distance)))

#Trying NLS with Michalis-Menton model - appears to work OK, but no way to do everything in one equation, or compare if different fits are significantly different
Vm1=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==5))[1]
K1=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==5))[2]
mod1=nls(success~SSmicmen(minShelter,Vm,K),data=filter(temp,Distance==5))

Vm2=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==20))[1]
K2=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==20))[2]
mod2=nls(success~SSmicmen(minShelter,Vm,K),data=filter(temp,Distance==20))

Vm3=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==100))[1]
K3=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==100))[2]
mod3=nls(success~SSmicmen(minShelter,Vm,K),data=filter(temp,Distance==100))

Vm4=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==400))[1]
K4=getInitial(success~SSmicmen(minShelter,max(success),1),data=filter(temp,Distance==400))[2]
mod4=nls(success~SSmicmen(minShelter,Vm,K),data=filter(temp,Distance==400))

newdata=data.frame(minShelter=5:75)
newdata$dist5=predict(mod1,newdata=newdata)
newdata$dist20=predict(mod2,newdata=newdata)
newdata$dist100=predict(mod3,newdata=newdata)
newdata$dist400=predict(mod4,newdata=newdata)

gather(newdata,'Distance','pred',2:5) %>%
  mutate(Distance=as.numeric(sub('dist','',Distance))) %>%
  ggplot(aes(minShelter,pred,col=as.factor(Distance)))+geom_line(size=1)+
  geom_point(data=temp,aes(minShelter,success),size=1)+
  scale_colour_manual(values=c('blue','steelblue','salmon','red'))+
  labs(x='Distance from shelter (m)',y='Pod success (%)',col='Distance from\nHoney bees(m)')+
  ylim(0.3,1)

#Pod success~Distance+VegMass -- no effect of Treatment (-1.21) or minShelter (-0.84, p=0.98)
temp=filter(plants2015,EdgeCent=='Edge',!is.na(VegMass))  #Removes center plots
  
pr=list(R=list(V=1,nu=0.2),
        G=list(G1=list(V = diag(1), nu = 1),G2=list(V = diag(1), nu = 1))) #Priors
#G=list(G1=list(V = diag(1), nu = 1))) #Priors
#flWinsDist1=MCMCglmm(cbind(Pods,Missing)~Distance,random=~Field+FieldPlot,

#VegMass (p<0.001) and Distance(<0.001) are significant
flWinsDist1=MCMCglmm(cbind(Pods,Missing)~VegMass+Distance,random=~Field+FieldPlot,
                     family='multinomial2',data=temp,verbose=F,prior=pr)
summary(flWinsDist1)

#Can't figure out how to do partial plots, so doing this instead:
flWinsDist1a=MCMCglmm(cbind(Pods,Missing)~log(VegMass),random=~Field+FieldPlot,
                      family='multinomial2',data=temp,verbose=F,prior=pr)
flWinsDist1b=MCMCglmm(cbind(Pods,Missing)~Distance,random=~Field+FieldPlot,
                      family='multinomial2',data=temp,verbose=F,prior=pr)

#Distance effects
pred=with(temp,data.frame(Distance,total=Pods+Missing,
                          predict(flWinsDist1b,interval='confidence',marginal=~Field+FieldPlot))) %>% 
  mutate(fitprob=fit/total,lwrprob=lwr/total,uprprob=upr/total)%>%
  select(Distance,fitprob,lwrprob,uprprob) %>%
  round(.,4) %>%
  distinct()

p=ggplot(temp,aes(Distance,Pods/(Pods+Missing)))+
  geom_point()+
  geom_ribbon(data=pred,aes(x=Distance,y=NULL,ymax=uprprob,ymin=lwrprob),alpha=0.3)+
  geom_line(data=pred,aes(x=Distance,y=fitprob))+
  labs(x='Distance to Honeybees (m)',y='Proportion flower success')
ggsave(paste(folder,'flWinsDist.png',sep='\\'),p,width=8,height=6)



#2) Floral success is lower at center of middle bay

temp=filter(plants2015,Distance!=20,Distance!=100) %>%
  mutate(Distance=factor(Distance,labels=c('Near','Far'))) #Convert distance to factor

#EdgeCent:Distance is not significant (p=0.36); EdgeCent (0.01) and Distance (<0.001) were significant
flWinsFbay=MCMCglmm(cbind(Pods,Missing)~EdgeCent+Distance,random=~Field+FieldPlot,
                     family='multinomial2',data=temp,verbose=F,prior=pr)
summary(flWinsFbay)

pred=with(temp,data.frame(Distance,EdgeCent,total=Pods+Missing,
                          predict(flWinsFbay,interval='confidence',marginal=~Field+FieldPlot))) %>% 
  mutate(fitprob=round(fit/total,4),lwrprob=round(lwr/total,4),uprprob=round(upr/total,4)) %>%
  select(-total:-upr) %>%
  distinct()

p=ggplot(temp,aes(Distance,Pods/(Pods+Missing),col=EdgeCent))+
  geom_point(position=position_jitterdodge(jitter.width=0.3,dodge.width=0.5),size=2)+
  labs(x='Distance',y='Proportion flower success',col='Bay\nposition')+
  scale_colour_manual(values=c('purple','forestgreen'))+ylim(0.4,0.85)+
  geom_pointrange(data=pred,aes(y=fitprob,ymin=lwrprob,ymax=uprprob,group=EdgeCent),position=position_dodge(width=0.5),col='black',size=1)
ggsave(paste(folder,'flWinsFbay.png',sep='\\'),p,width=8,height=6)

#3) Floral success is lower at center of middle bay

flWinsDistShelter=MCMCglmm(cbind(Pods,Missing)~minShelter*(Treatment),random=~Field+FieldPlot,
                    family='multinomial2',data=plants2015,verbose=F,prior=pr)

pred=with(temp,data.frame(minShelter,Treatment,total=Pods+Missing,
                          predict(flWinsDistShelter,interval='confidence',marginal=~Field+FieldPlot))) %>% 
  mutate(fitprob=round(fit/total,4),lwrprob=round(lwr/total,4),uprprob=round(upr/total,4)) %>%
  select(-total:-upr) %>%
  distinct()

p=ggplot(plants2015,aes(minShelter,Pods/(Pods+Missing),col=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees'))))+
  geom_ribbon(data=pred,aes(x=minShelter,y=NULL,ymax=uprprob,ymin=lwrprob,fill=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees')),col=NULL),alpha=0.3)+
  geom_point(position=position_jitterdodge(jitter.width=0.3,dodge.width=0.5),size=1)+
  geom_line(data=pred,aes(x=minShelter,y=fitprob),size=1)+
  labs(x='Distance to leafcutter shelter (m)',y='Proportion flower success',fill='Treatment',col='Treatment')+
  scale_colour_manual(values=c('blue','darkorange','red'))+
  scale_fill_manual(values=c('blue','darkorange','red'))
ggsave(paste(folder,'flWinsDistShelter.png',sep='\\'),p,width=8,height=6)


# Seed production ---------------------------------------------------------
#PREDICTIONS: 
#1) Seeds/pod declines with distance from beehives/leafcutter tents

#PodCount~Distance + Treatment -- No effect of Shelter distance (p=0.85), Shelter Dist:Honeybee Dist (p=0.96), or Treatment:Dist (p=0.15 & 0.56)
temp=filter(seeds2015,!is.na(PodCount),!is.na(PodMass),EdgeCent=='Edge')

pr=list(R=list(V=diag(1),nu=0.02),
        G=list(G1=list(V = diag(1),nu=0.02),
               G2=list(V = diag(1), nu =0.02, alpha.mu=0, alpha.V=1000),
               G3=list(V = diag(1), nu =0.02, alpha.mu=0, alpha.V=1000))) #Priors
seedsPod1=MCMCglmm(PodCount~Treatment+Distance,family='poisson',random=~Field+FieldPlot+Plant,data=temp,
                   prior=pr,nitt=50000+5000,burnin=5000,thin=50,verbose=F)
pred=unique(with(temp,data.frame(Distance,Treatment,predict(seedsPod1,interval='confidence'))))

p=ggplot(pred,aes(x=Distance,y=fit,colour=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees')),fill=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees')),ymax=upr,ymin=lwr))+
  geom_line(size=1)+geom_ribbon(aes(colour=NULL),alpha=0.3)+
  scale_colour_manual(values=c('blue','darkorange','red'))+
  scale_fill_manual(values=c('blue','darkorange','red'))+
  geom_point(data=temp,aes(x=Distance,y=PodCount,ymax=NULL,ymin=NULL),position=position_jitter(),size=1)+
  labs(x='Distance to honeybees (m)',y='Seeds per pod',col='Treatment',fill='Treatment')
ggsave(paste(folder,'seedsPodDistanceTreatment.png',sep='\\'),p,width=8,height=6)

# Yield changes -----------------------------------------------------------
#PREDICTIONS: 
#50 lbs/bushel * 453.492 g/lb = 22679.6 g/bushel = 4.409249e-05 bushels/g
#1 acre = 4046.86m2 = 0.0002471052 acres/m2

#Grams/m2 * (4.409249e-05/0.0002471052) = g/m2 * 0.1784361 = bushels/acre
conversion= (1/(50*453.492))/(1/4046.86)

#1) Seed yield/m2 higher in center of field
temp=filter(plants2015,EdgeCent!='Center') %>%
  mutate(yield=SeedMass*PlDens*conversion)%>%
  mutate(Treatment=factor(Treatment,levels=c('control','double tent','double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees')))

pr=list(R = list(V = diag(1), nu = 0.02),
        G=list(G1=list(V=diag(2), nu=0.01,alpha.mu=c(0,0), alpha.V=diag(2)*1000)))

#Distance is significant (0.03), but treatment is not (right direction of effects, but not significant)
yieldDist=MCMCglmm(yield~Distance,random=~idh(1+Distance):Field,data=temp,
                   verbose=F,prior=pr,nitt=50000+5000,burnin=5000,thin=50)
summary(yieldDist)

pred=with(temp,data.frame(predict(yieldDist,interval='confidence'),Distance)) %>%distinct()

p=ggplot(temp,aes(Distance,yield))+
  geom_ribbon(data=pred,aes(x=Distance,y=NULL,ymax=upr,ymin=lwr),alpha=.3)+
  geom_line(data=pred,aes(x=Distance,y=fit),size=1)+ylim(15,100)+
  geom_point()+labs(x='Distance from honeybees(m)',y='Bushels/acre')
ggsave(paste(folder,'yieldDist.png',sep='\\'),p,width=8,height=6)

#2) Seed yield/m2 higher at edge of bay
temp=filter(plants2015,Distance!=20,Distance!=100,Distance!=250,!is.na(SeedMass)) %>%
  mutate(Distance=factor(ifelse(Distance==5,'Near','Far'))) %>%
  mutate(Distance=factor(Distance,levels=c('Near','Far')),yield=SeedMass*PlDens*conversion)%>%
  mutate(Treatment=factor(Treatment,levels=c('control','double tent','double tent+bees'),labels=c('Control','Double tent','Double tent\n+bees')))

pr=list(R = list(V = diag(1), nu = 0.02),
        G=list(G1=list(V=1, nu=0.01,alpha.mu=0, alpha.V=1000),
               G2=list(V=1, nu=0.01,alpha.mu=0, alpha.V=1000)))

#Distance:EdgeCent is not significant (p=0.90), but Distance (0.02) and EdgeCent (<0.001) are. Treatment or any combination of Distance/EdgeCent are not significant
yieldDistFbay=MCMCglmm(yield~Distance*EdgeCent,random=~Field+FieldPlot,data=temp,
                   verbose=F,prior=pr,nitt=50000+5000,burnin=5000,thin=50)

summary(yieldDistFbay)

pred=with(temp,data.frame(predict(yieldDistFbay,interval='confidence'),Distance,EdgeCent)) %>%distinct()

p=ggplot(temp,aes(Distance,yield,col=EdgeCent))+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.width=0.4))+
  labs(y='Bushels/acre',col='Bay\nPosition')+ylim(0,90)+
  scale_colour_manual(values=c('purple','forestgreen'))+
  geom_pointrange(data=pred,aes(x=Distance,y=fit,ymax=upr,ymin=lwr,group=EdgeCent),col='black',position=position_dodge(width=0.75))

ggsave(paste(folder,'yieldDistFbay.png',sep='\\'),p,width=8,height=6)

#3) Seed yield/m2 higher towards leafcutter shelters, not really for control

pr=list(R = list(V = diag(1), nu = 0.02),
        G=list(G1=list(V=1, nu=0.01,alpha.mu=0, alpha.V=1000),
               G2=list(V=1, nu=0.01,alpha.mu=0, alpha.V=1000)))
temp=mutate(plants2015,yield=SeedMass*PlDens*conversion)
yieldDistShelter=MCMCglmm(yield~minShelter*Treatment,random=~Field+FieldPlot,data=temp,
                          verbose=F,prior=pr,nitt=50000+5000,burnin=5000,thin=100)
summary(yieldDistShelter)
pred=with(temp,data.frame(round(predict(yieldDistShelter,interval='confidence'),1),minShelter,Treatment)) %>%distinct()

p=ggplot(temp,aes(minShelter,yield,col=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees'))))+
  geom_ribbon(data=pred,aes(x=minShelter,y=NULL,ymax=upr,ymin=lwr,col=NULL,fill=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees'))),alpha=0.3)+
  geom_point()+#geom_smooth(method='gam',formula=y~s(x),se=F)+
  geom_line(data=pred,aes(x=minShelter,y=fit,col=factor(Treatment,labels=c('Control','Double tent','Double tent\n+bees'))),size=1)+
  labs(y='Bushels/acre',x='Distance to leafcutter shelter(m)',col='Treatment',fill='Treatment')+
  scale_colour_manual(values=c('blue','darkorange','red'))+ylim(0,125)+
  scale_fill_manual(values=c('blue','darkorange','red'))
ggsave(paste(folder,'yieldDistShelter.png',sep='\\'),p,width=8,height=6)


