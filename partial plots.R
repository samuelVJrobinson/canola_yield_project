# MAKES PARTIAL EFFECTS PLOTS FOR SUB-MODELS

# Load everything ---------------------------------------------------------
library(tidyverse)
library(ggpubr)
theme_set(theme_bw())
setwd('~/Documents/canola_yield_project/Figures') #Galpern machine path
# setwd('~/Projects/UofC/canola_yield_project/Models') #Multivac path

source('../helperFunctions.R') #Helper functions

logCent <- function(x) log(x)-mean(log(x)) #Log-centers data
summaryDat <- function(x){
  x <- na.omit(x)
  data.frame(min=min(x),mean=mean(x),med=median(x),max=max(x))
} 

load('../Models/datalist_commodity.Rdata') #Commodity data
commData <- datalist
load('../Models/datalist_seed.Rdata') #Seed data
seedData <- datalist; rm(datalist)

#Summary stats (min,median,max) to use for plots

#Commodity data
commData$logHbeeVis <- log((commData$hbeeVis/commData$totalTime)+1) #Visitation rate
commData$logHbeeDist <- log(commData$dist) #Distance
commData$clogHbeeDist <- commData$logHbeeDist-mean(commData$logHbeeDist) #Centered log distance
commData$logPodCount <- log(commData$podCount) #Log pod count
commData$logFlwCount <- log(commData$flwCount) #Log pod count
commData$propFlwSurv <- commData$podCount/commData$flwCount #Flower survival
commData$propFlwSurv[commData$propFlwSurv==1] <- 0.99 #Deal with 100% success rates
commData$logitFlwSurv <- logit(commData$propFlwSurv)
commData$clogitFlwSurv <- commData$logitFlwSurv-mean(commData$logitFlwSurv) #Centered logit
avgCommData <- lapply(commData,summaryDat) 

#Seed data
seedData$lbeeDistAll <- with(seedData,c(lbee_dist,lbee_dist_extra)) #Distances
seedData$hbeeDistAll <- with(seedData,c(hbee_dist,hbee_dist_extra))
seedData$clogLbeeDist <- with(seedData,logCent(lbeeDistAll)) #Log-centered distances
seedData$clogHbeeDist <- with(seedData,logCent(hbeeDistAll))
seedData$logLbeeVis <- log((seedData$lbeeVis/seedData$totalTime)+1) #Visitation rates
seedData$logHbeeVis <- log((seedData$hbeeVis/seedData$totalTime)+1) 
seedData$logFlwCount <- log(seedData$flwCount) #Log flw count

avgSeedData <- lapply(seedData,summaryDat) 

avgSeedData$logHbeeDistCent <- mean(log(with(seedData,c(hbee_dist,hbee_dist_extra)))) #Center
avgSeedData$logLbeeDistCent <- mean(log(with(seedData,c(lbee_dist,lbee_dist_extra))))

load('../Models/modSummaries_commodity.Rdata') #Commodity data models
load('../Models/modSummaries_seed.Rdata') #Seed data models

# Visitation plots --------------------------------------------------------

timeOffset <- 6 #Use 1 hr offset (1 = 10 mins, 6 = 1 hr)
lab <- c('HB Seed','LCB Seed','HB Commodity') 
labCols <- c('darkorange','darkgreen','black')

#Hbee stocking rate numbers
with(avgCommData,
           list('intHbeeVis'=1,
                'slopeHbeeDistHbeeVis'=log(1),
                'slopeNumHivesHbeeVis'= log(c(0,20,40)+1), 
                'slopeFlDensHbeeVis' = 0)) %>% 
  getPreds(modSummaries_commodity[[2]],parList = .,offset=timeOffset,
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(type='Commodity',numHives=c(0,20,40),mean,med,lwr,upr) 

#Hbee + Lbee visitation from edge of field
dists <- with(avgCommData,seq(min(log(dist)),max(log(400)),length=20)) #Distances
dists2 <- with(avgSeedData,seq(min(clogHbeeDist),max(clogHbeeDist),length=20)) #Centered log distances for seed fields

d1 <- with(avgCommData, #Hbee visitation in commodity
     list('intHbeeVis'=1,
     'slopeHbeeDistHbeeVis'=dists,
     'slopeNumHivesHbeeVis'= log(40+1), 
     'slopeFlDensHbeeVis' = 0)) %>% 
  getPreds(modSummaries_commodity[[2]],parList = .,offset=timeOffset,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(dists),mean,med,lwr,upr)

d2 <- with(avgSeedData, #Hbee visitation in seed fields
     list('intHbeeVis'=1,
          'slopeFlDensHbeeVis' = flDens_obs$mean,
          'slopeHbeeDistHbeeVis'=dists2,
          'slopeLbeeDistHbeeVis'= 0,
          'slopeCentHbeeVis'= 0
          )) %>% 
  getPreds(modSummaries_seed[[2]],parList = .,offset=timeOffset,
           ZIpar = 'thetaHbeeVis',trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(dists),mean,med,lwr,upr)

# bind_rows(d1,d2,.id='type') %>% mutate(type=factor(type,labels=c('Commodity','Seed'))) %>% 
#   ggplot(aes(x=dist,y=mean))+
#   geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
#   geom_line(aes(col=type))+
#   labs(x='Distance from apiary (m)',y='Visits per hour',
#        fill='Field\nType',col='Field\nType')

d3 <- with(avgSeedData,
                  list('intVisitLbeeVis'=1,
                       'slopeFlDensLbeeVis' = 0,
                       'slopeHbeeDistLbeeVis'=dists2,
                       'slopeLbeeDistLbeeVis'= 0,
                       'slopeCentLbeeVis'= 0
                  )) %>% 
         getPreds(modSummaries_seed[[3]],parList = .,offset=timeOffset,
                  ZIpar = 'thetaLbeeVis',trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(dists),mean,med,lwr,upr)

p1 <- bind_rows(d2,d3,d1,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  coord_cartesian(y=c(0,150))+
  labs(x='Distance from apiary (m)',y='Visits per hour',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols)


#Hbee and Lbee visitation from shelters

dists <- log(1:30)-avgSeedData$logLbeeDistCent #log Shelter Distances

d1 <- with(avgSeedData, #Hbee visits
           list('intHbeeVis'=1,
                'slopeFlDensHbeeVis' = 0,
                'slopeHbeeDistHbeeVis'=0,
                'slopeLbeeDistHbeeVis'=dists,
                'slopeCentHbeeVis'= 0
           )) %>% 
  getPreds(modSummaries_seed[[2]],parList = .,offset=timeOffset,ZIpar = 'thetaHbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=1:30,mean,med,lwr,upr)

d2 <- with(avgSeedData, #Lbee visits
           list('intVisitLbeeVis'=1,
                'slopeFlDensLbeeVis' = 0,
                'slopeHbeeDistLbeeVis'=0,
                'slopeLbeeDistLbeeVis'=dists,
                'slopeCentLbeeVis'= 0
           )) %>% 
  getPreds(modSummaries_seed[[3]],parList = .,offset=timeOffset,ZIpar = 'thetaLbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=1:30,mean,med,lwr,upr)

p2 <- bind_rows(d1,d2,.id='type') %>% 
  mutate(type=factor(type,labels=lab[1:2])) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Distance from shelter (m)',y='Visits per hour',
       fill=NULL,col=NULL)+
  coord_cartesian(y=c(0,300))+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols)

#Center and edge of bays
d1 <- with(avgSeedData, #Hbee visits
           list('intHbeeVis'=1,
                'slopeFlDensHbeeVis' = 0,
                'slopeHbeeDistHbeeVis'=0,
                'slopeLbeeDistHbeeVis'=0,
                'slopeCentHbeeVis'= c(0,1)
           )) %>% 
  getPreds(modSummaries_seed[[2]],parList = .,offset=timeOffset,ZIpar = 'thetaHbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(bayType=c('Edge','Centre'),mean,med,lwr,upr)

d2 <- with(avgSeedData, #Lbee visits
           list('intVisitLbeeVis'=1,
                'slopeFlDensLbeeVis' = 0,
                'slopeHbeeDistLbeeVis'=0,
                'slopeLbeeDistLbeeVis'=0,
                'slopeCentLbeeVis'= c(0,1)
           )) %>% 
  getPreds(modSummaries_seed[[3]],parList = .,offset=timeOffset,ZIpar = 'thetaLbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(bayType=c('Edge','Centre'),mean,med,lwr,upr)

p3 <- bind_rows(d1,d2,.id='type') %>% 
  mutate(type=factor(type,labels=lab[1:2])) %>% 
  mutate(bayType=factor(bayType,levels=c('Edge','Centre'))) %>% 
  ggplot(aes(x=bayType,col=type))+
  geom_pointrange(aes(y=mean,ymax=upr,ymin=lwr))+
  labs(x='Female bay position',y='Visits per hour',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols[1:2])

(p <- ggarrange(p1,p2,p3,ncol=3,nrow=1,common.legend = TRUE,legend='bottom',labels = 'auto'))

ggsave('allVisits.png',p,bg='white',width = 12,height=4)

# Pollen plots ------------------------------------------------------------

#Visitation effects

#median, mean, max:
#0,1.72,36 max for commodity, 4,22.6,215 for seed
hbeeVisComm <- log(seq(0,max(with(commData,hbeeVis/totalTime)),length=30)+1)
hbeeVisSeed <- log((seq(0,max(with(seedData,hbeeVis/totalTime)),length=100))+1)
lbeeVisSeed <- log((seq(0,max(with(seedData,lbeeVis/totalTime)),length=100))+1)
lab <- c('HB Seed','LCB Seed','HB Commodity') 
labCols <- c('darkorange','darkgreen','black')

d1 <- with(avgCommData, #Hbee visitation effect - commodity
     list('intPollen'=1,
          'slopeHbeeVisPollen'=hbeeVisComm,
          'slopeHbeeDistPollen'= 0
          )) %>% 
  getPreds(modSummaries_commodity[[3]],parList = .,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(vis=exp(hbeeVisComm)-1,mean,med,lwr,upr)

d2 <- with(avgSeedData, #Hbee vis effect - seed
           list('intPollen'=1,
                'slopeHbeeVisPollen'=hbeeVisSeed,
                'slopeLbeeVisPollen'=avgSeedData$logLbeeVis$mean,
                'slopeCentPollen'=0,
                'slopeHbeeDistPollen'=0,
                'slopeLbeeDistPollen'=0,
                'slopeFlDensPollen'=0
           )) %>% 
  getPreds(modSummaries_seed[[4]],parList = .,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(vis=exp(hbeeVisSeed)-1,mean,med,lwr,upr)

d3 <- with(avgSeedData, #Lbee vis effect - seed
           list('intPollen'=1,
                'slopeHbeeVisPollen'=avgSeedData$logHbeeVis$mean,
                'slopeLbeeVisPollen'=lbeeVisSeed,
                'slopeCentPollen'=0,
                'slopeHbeeDistPollen'=0,
                'slopeLbeeDistPollen'=0,
                'slopeFlDensPollen'=0
           )) %>% 
  getPreds(modSummaries_seed[[4]],parList = .,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(vis=exp(hbeeVisSeed)-1,mean,med,lwr,upr)

p1 <- bind_rows(d2,d3,d1,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=vis,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Visits per hour',y='Pollen grains per stigma',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols)+
  scale_y_log10()

#Distance effects
hbeeDistComm <- with(avgCommData,seq(clogHbeeDist$min,log(400)-avgCommData$logHbeeDist$mean,length=100))
hbeeDistSeed <- with(avgSeedData,seq(clogHbeeDist$min,clogHbeeDist$max,length=100))
lbeeDistSeed <- with(avgSeedData,seq(clogLbeeDist$min,log(100)-avgSeedData$logLbeeDist,length=100))

d1 <- with(avgCommData, #Hbee distance - commodity
           list('intPollen'=1,
                'slopeHbeeVisPollen'=avgCommData$logHbeeVis$mean,
                'slopeHbeeDistPollen'= hbeeDistComm
           )) %>%
  getPreds(modSummaries_commodity[[3]],parList = .,trans='exp',q=c(0.5,0.05,0.95)) %>%
  transmute(dist=exp(hbeeDistComm+avgCommData$logHbeeDist$mean),mean,med,lwr,upr)

d2 <- with(avgSeedData, #Hbee dist - seed
     list('intPollen'=1,
          'slopeHbeeVisPollen'=avgSeedData$logHbeeVis$mean,
          'slopeLbeeVisPollen'=avgSeedData$logLbeeVis$mean,
          'slopeCentPollen'=0,
          'slopeHbeeDistPollen'=hbeeDistSeed,
          'slopeLbeeDistPollen'=0,
          'slopeFlDensPollen'=0
     )) %>% 
  getPreds(modSummaries_seed[[4]],parList = .,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(hbeeDistSeed+avgSeedData$logHbeeDist),mean,med,lwr,upr)

d3 <- with(avgSeedData, #Lbee dist - seed
           list('intPollen'=1,
                'slopeHbeeVisPollen'=avgSeedData$logHbeeVis$mean,
                'slopeLbeeVisPollen'=avgSeedData$logLbeeVis$mean,
                'slopeCentPollen'=0,
                'slopeHbeeDistPollen'=0,
                'slopeLbeeDistPollen'=lbeeDistSeed,
                'slopeFlDensPollen'=0
           )) %>% 
  getPreds(modSummaries_seed[[4]],parList = .,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(lbeeDistSeed+avgSeedData$logLbeeDist),mean,med,lwr,upr)

p2 <- bind_rows(d2,d3,d1,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Distance from apiary or shelter (m)',y='Pollen grains per stigma',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols)+
  scale_y_log10()

p3 <- with(avgSeedData, #Bay position
     list('intPollen'=1,
          'slopeHbeeVisPollen'=avgSeedData$logHbeeVis$mean,
          'slopeLbeeVisPollen'=avgSeedData$logLbeeVis$mean,
          'slopeCentPollen'=c(0,1),
          'slopeHbeeDistPollen'=0,
          'slopeLbeeDistPollen'=0,
          'slopeFlDensPollen'=0
     )) %>% 
  getPreds(modSummaries_seed[[4]],parList = .,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(bayType=factor(c('Edge','Centre'),levels=c('Edge','Centre')),mean,med,lwr,upr) %>% 
  ggplot(aes(x=bayType))+
  geom_pointrange(aes(y=mean,ymax=upr,ymin=lwr))+
  labs(x='Female bay position',y='Pollen grains per stigma',fill=NULL,col=NULL)

(p <- ggarrange(p1,p2,p3,ncol=3,nrow=1,legend='bottom',common.legend = TRUE,labels='auto'))

ggsave('allPollen.png',p,bg='white',width = 12,height=4)

# Seed production - pollen ---------------------------------------------------

#pollenPlot means - Taken from model output
commPollenPlot <- seq(-1.2644033,0.9869057,length=20) 
modSummaries_commodity[[3]]$mean[1] #Intercept
seedPollenPlot <- seq(-3.573654,2.364599,length=20) 
modSummaries_seed[[4]]$mean[1] #Intercept

lab <- c('Commodity','Seed')
labCols <- c('blue','darkred')

#Flower survival
d1 <- with(avgCommData, 
     list('intFlwSurv'=1,
          'slopeHbeeVisFlwSurv'=avgCommData$logHbeeVis$mean,
          'slopePlSizeFlwSurv'= avgCommData$plantSize$mean,
          'slopePollenFlwSurv'=commPollenPlot
          )) %>% 
  getPreds(modSummaries_commodity[[5]],trans='invLogit',parList = .,q=c(0.5,0.05,0.95)) %>% 
  transmute(pol=exp(commPollenPlot+modSummaries_commodity[[3]]$mean[1]),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
     list('intFlwSurv'=1,
          'slopePollenFlwSurv'=seedPollenPlot,
          'slopePlSizeFlwSurv'=avgSeedData$plantSize$mean,
          'slopeCentFlwSurv'=0,
          'slopeHbeeDistFlwSurv'=0,
          'slopeFlDensFlwSurv'=0
     )) %>% 
  getPreds(modSummaries_seed[[6]],parList = .,trans='invLogit',q=c(0.5,0.05,0.95)) %>% 
  transmute(pol=exp(seedPollenPlot+modSummaries_seed[[4]]$mean[1]),mean,med,lwr,upr) 

p1 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  mutate(across(c(mean:upr),~.x*100)) %>% 
  ggplot(aes(x=pol,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Pollen grains per stigma',y='Flower survival (%)',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
  scale_x_log10()

#Seeds per pod
d1 <- with(avgCommData, 
           list('intSeedCount'=1,
                'slopeHbeeVisSeedCount'=avgCommData$logHbeeVis$mean,
                'slopePollenSeedCount'= commPollenPlot,
                'slopePlSizeSeedCount'=avgCommData$plantSize$mean,
                'slopeFlwSurvSeedCount'=0, #centered logit flw surv
                'slopeFlwCountSeedCount'=avgCommData$logFlwCount$mean #Log flower count
           )) %>% 
  getPreds(modSummaries_commodity[[6]],parList = .,q=c(0.5,0.05,0.95),otherPars = 'lambdaSeedCount') %>% 
  mutate(across(c(mean:upr),~.x+1/lambdaSeedCount)) %>% 
  transmute(pol=exp(commPollenPlot+modSummaries_commodity[[3]]$mean[1]),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
           list('intSeedCount'=1,
                'slopePollenSeedCount'=seedPollenPlot,
                'slopePlSizeSeedCount'=avgSeedData$plantSize$mean,
                'slopeCentSeedCount'=0,
                'slopeHbeeDistSeedCount'=0, #Centered log dist
                'slopeFlDensSeedCount'=0,
                'slopeFlwSurvSeedCount'=avgSeedData$logitFlwSurv$mean,
                'slopeFlwCountSeedCount'=avgSeedData$logFlwCount$mean
           )) %>% 
  getPreds(modSummaries_seed[[7]],parList = .,otherPars='lambdaSeedCount',q=c(0.5,0.05,0.95)) %>% 
  mutate(across(c(mean:upr),~.x+1/lambdaSeedCount)) %>% 
  transmute(pol=exp(seedPollenPlot+modSummaries_seed[[4]]$mean[1]),mean,med,lwr,upr) 

p2 <- bind_rows(d1,d2,.id = 'type') %>% 
    mutate(type=factor(type,labels=lab)) %>% 
    ggplot(aes(x=pol,y=mean))+
    geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
    geom_line(aes(col=type))+
    labs(x='Pollen grains per stigma',y='Seeds per pod',fill=NULL,col=NULL)+
    scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
    scale_x_log10()

#Seed size

d1 <- with(avgCommData, 
           list('intSeedWeight'=1,
                'slopeHbeeVisSeedWeight'=avgCommData$logHbeeVis$mean,
                'slopePollenSeedWeight'= commPollenPlot,
                'slopeSeedCountSeedWeight'=avgCommData$seedCount$mean,
                'slopePlSizeSeedWeight'=avgCommData$plantSize$mean,
                'slopePlDensSeedWeight'=avgCommData$plDens_obs$mean,
                'slopeHbeeDistSeedWeight'=0,
                'slopeFlwSurvSeedWeight'=0, #centered logit flw surv
                'slopeFlwCountSeedWeight'=avgCommData$logFlwCount$mean #Log flower count
           )) %>% 
  getPreds(modSummaries_commodity[[7]],parList = .,q=c(0.5,0.05,0.95),otherPars = 'lambdaSeedWeight') %>% 
  mutate(across(c(mean:upr),~.x+(1/lambdaSeedWeight))) %>% 
  transmute(pol=exp(commPollenPlot+modSummaries_commodity[[3]]$mean[1]),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
           list('intSeedWeight'=1,
                'slopePollenSeedWeight'=seedPollenPlot,
                'slopeSeedCountSeedWeight'=avgSeedData$seedCount_obs$mean,
                'slopePlSizeSeedWeight'=avgSeedData$plantSize$mean,
                'slopePlDensSeedWeight'=avgSeedData$plDens_obs$mean,
                'slopeLbeeDistSeedWeight'=0
           )) %>% 
  getPreds(modSummaries_seed[[8]],parList = .,otherPars='lambdaSeedWeight',q=c(0.5,0.05,0.95)) %>% 
  mutate(across(c(mean:upr),~.x+(1/lambdaSeedWeight))) %>% 
  transmute(pol=exp(seedPollenPlot+modSummaries_seed[[4]]$mean[1]),mean,med,lwr,upr) 

p3 <- bind_rows(d1,d2,.id = 'type') %>% 
    mutate(type=factor(type,labels=lab)) %>% 
    ggplot(aes(x=pol,y=mean))+
    geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
    geom_line(aes(col=type))+
    labs(x='Pollen grains per stigma',y='Seed size (mg)',fill=NULL,col=NULL)+
    scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
    scale_x_log10()

(p <- ggarrange(p1,p2,p3,ncol=3,common.legend = TRUE,legend='bottom'))
ggsave('allSeeds_pollen.png',p,bg='white',width = 12,height=4)

# Seed production - plant size ---------------------------------------------------

#plSize - both log-transformed
commPlSize <- with(avgCommData$plantSize,seq(min,max,length=20)) 
seedPlSize <- with(avgSeedData$plantSize,seq(min,max,length=20))

lab <- c('Commodity','Seed')
labCols <- c('blue','darkred')

#Flower survival
d1 <- with(avgCommData, 
           list('intFlwSurv'=1,
                'slopeHbeeVisFlwSurv'=avgCommData$logHbeeVis$mean,
                'slopePlSizeFlwSurv'= commPlSize,
                'slopePollenFlwSurv'=0
           )) %>% 
  getPreds(modSummaries_commodity[[5]],trans='invLogit',parList = .,q=c(0.5,0.05,0.95)) %>% 
  transmute(plSize=exp(commPlSize),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
           list('intFlwSurv'=1,
                'slopePollenFlwSurv'=0,
                'slopePlSizeFlwSurv'=seedPlSize,
                'slopeCentFlwSurv'=0,
                'slopeHbeeDistFlwSurv'=0,
                'slopeFlDensFlwSurv'=0
           )) %>% 
  getPreds(modSummaries_seed[[6]],parList = .,trans='invLogit',q=c(0.5,0.05,0.95)) %>% 
  transmute(plSize=exp(seedPlSize),mean,med,lwr,upr) 

p1 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  mutate(across(c(mean:upr),~.x*100)) %>% 
  ggplot(aes(x=plSize,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Plant size (g)',y='Flower survival (%)',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols) +
  scale_x_log10()

#Seeds per pod
d1 <- with(avgCommData, 
           list('intSeedCount'=1,
                'slopeHbeeVisSeedCount'=avgCommData$logHbeeVis$mean,
                'slopePollenSeedCount'= 0,
                'slopePlSizeSeedCount'=commPlSize,
                'slopeFlwSurvSeedCount'=0, #centered logit flw surv
                'slopeFlwCountSeedCount'=avgCommData$logFlwCount$mean #Log flower count
           )) %>% 
  getPreds(modSummaries_commodity[[6]],parList = .,q=c(0.5,0.05,0.95),otherPars = 'lambdaSeedCount') %>% 
  mutate(across(c(mean:upr),~.x+1/lambdaSeedCount)) %>% 
  transmute(plSize=exp(commPlSize),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
           list('intSeedCount'=1,
                'slopePollenSeedCount'=0,
                'slopePlSizeSeedCount'=seedPlSize,
                'slopeCentSeedCount'=0,
                'slopeHbeeDistSeedCount'=0, #Centered log dist
                'slopeFlDensSeedCount'=0,
                'slopeFlwSurvSeedCount'=avgSeedData$logitFlwSurv$mean,
                'slopeFlwCountSeedCount'=avgSeedData$logFlwCount$mean
           )) %>% 
  getPreds(modSummaries_seed[[7]],parList = .,otherPars='lambdaSeedCount',q=c(0.5,0.05,0.95)) %>% 
  mutate(across(c(mean:upr),~.x+1/lambdaSeedCount)) %>% 
  transmute(plSize=exp(seedPlSize),mean,med,lwr,upr) 

p2 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=plSize,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Plant size (g)',y='Seeds per pod',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
  scale_x_log10()

#Seed size

d1 <- with(avgCommData, 
           list('intSeedWeight'=1,
                'slopeHbeeVisSeedWeight'=avgCommData$logHbeeVis$mean,
                'slopePollenSeedWeight'= 0,
                'slopeSeedCountSeedWeight'=avgCommData$seedCount$mean,
                'slopePlSizeSeedWeight'=commPlSize,
                'slopePlDensSeedWeight'=avgCommData$plDens_obs$mean,
                'slopeHbeeDistSeedWeight'=0, #centered log hbee dist
                'slopeFlwSurvSeedWeight'=0, #centered logit flw surv
                'slopeFlwCountSeedWeight'=avgCommData$logFlwCount$mean #Log flower count
           )) %>% 
  getPreds(modSummaries_commodity[[7]],parList = .,q=c(0.5,0.05,0.95),otherPars = 'lambdaSeedWeight') %>% 
  mutate(across(c(mean:upr),~.x+1/lambdaSeedWeight)) %>% 
  transmute(plSize=exp(commPlSize),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
           list('intSeedWeight'=1,
                'slopePollenSeedWeight'=0,
                'slopeSeedCountSeedWeight'=avgSeedData$seedCount_obs$mean,
                'slopePlSizeSeedWeight'=seedPlSize,
                'slopePlDensSeedWeight'=avgSeedData$plDens_obs$mean,
                'slopeLbeeDistSeedWeight'=0 #centered log lbee dist
           )) %>% 
  getPreds(modSummaries_seed[[8]],parList = .,otherPars='lambdaSeedWeight',q=c(0.5,0.05,0.95)) %>% 
  mutate(across(c(mean:upr),~.x+(1/lambdaSeedWeight))) %>% 
  transmute(plSize=exp(seedPlSize),mean,med,lwr,upr) 

p3 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=plSize,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Plant size (g)',y='Seed size (mg)',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
  scale_x_log10()

(p <- ggarrange(p1,p2,p3,ncol=3,common.legend = TRUE,legend='bottom'))
ggsave('allSeeds_plSize.png',p,bg='white',width = 12,height=4)


# Total yield and seed size for commodity fields -----------------------------------------------------

#Fixed terms: plSize, plDens, flDens, hbeeDist, numHives

dat <- list(
  avgPlSize = log(mean(exp(commData$plantSize))), #log avg plant size
  avgPlDens = log(mean(exp(commData$plDens_obs))), #log avg plant dens
  avgFlDens = 0, #fl Dens (centered)
  hbeeDists = 1:100, #1-100 m
  numHives = log(40+1) #40 hives
)

dat$clogHbeeDists <- log(dat$hbeeDists)-avgCommData$logHbeeDist$mean

# datdf <- expand.grid(
#   avgPlSize = log(mean(exp(commData$plantSize))), #log avg plant size
#   avgPlDens = log(mean(exp(commData$plDens_obs))), #log avg plant dens
#   avgFlDens = 0, #fl Dens (centered)
#   hbeeDists = 0:100, #0-100 m
#   numHives = log(40+1) #40 hives
# )

#Visitation
dat$hbeeVis <- list('intHbeeVis'=1,
                    'slopeHbeeDistHbeeVis'=dat$clogHbeeDists,
                    'slopeNumHivesHbeeVis'= dat$numHives, 'slopeFlDensHbeeVis' = 0) %>% 
  getPreds(modSummaries_commodity[[2]],parList = .,offset=1,trans='exp',q=c(0.5,0.05,0.95)) %>% pull(mean)

#Pollen
dat$pollen <- list('intPollen'=1,'slopeHbeeVisPollen'=log(dat$hbeeVis+1),
     'slopeHbeeDistPollen'= dat$clogHbeeDists) %>% 
  getPreds(modSummaries_commodity[[3]],parList = .,trans='exp',q=c(0.5,0.05,0.95),expandGrid = FALSE) %>% pull(mean)

dat$clogPollen <- log(dat$pollen)-modSummaries_commodity[[3]]$mean[1]

#Flower count


#Seed count

#Seed size

#Yield
