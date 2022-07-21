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

# avgCommData$flDensMean$mean #Avg sqrt flw dens
# avgCommData$logHbeeDist$mean #avg log hbee dist
# avgCommData$logitFlwSurv$mean #avg logit flw surv

debugonce(genCommYield)
d <- list()
d_mean <- genCommYield(dat=d,simPar=TRUE)
d <- replicate(100,genCommYield(dat=d,simPar = TRUE),simplify = FALSE) %>% bind_rows()
hist(d$yield_tha); abline(v=d_mean$yield_tha)

#Distance and stocking effect
datList <- list(
  hbeeDist = rep(1:100,2), #1-100 m
  plDens=50, #Same plant density
  numHives = rep(c(0,40),each=100) #0,40 hives
) 

d_mean <- genCommYield(datList) %>% as.data.frame() %>% 
  mutate(numHives=factor(numHives))

#Parallel version
# detectCores()
{ #Takes about 1 min
  cl <- makeCluster(14) 
  clusterExport(cl=cl,c('genCommYield'))
  d <- parLapply(cl,1:1000,function(i,d,m){
    source('../helperFunctions.R') #Helper functions
    genCommYield(dat = d, mods = m, simPar=TRUE)},d=datList,m=modSummaries_commodity) %>% 
    bind_rows(.id='rep')
  stopCluster(cl)
}

d %>% mutate(numHives=factor(numHives)) %>% group_by(hbeeDist,numHives) %>% 
  summarize(med=median(yield_tha),lwr=quantile(yield_tha,0.25),upr=quantile(yield_tha,0.75),medVis=median(hbeeVis)) %>% 
  ggplot(aes(x=hbeeDist,y=med))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=numHives),alpha=0.3)+
  geom_line(data=d_mean,aes(x=hbeeDist,y=yield_tha,col=numHives))
  
#Visitation effect
datList <- list(
  hbeeVis=0:36,
  plDens=50 #Same plant density
) 

debugonce(genCommYield)
debugonce(getPreds)
d_mean <- genCommYield(datList) %>% as.data.frame() 

{
  cl <- makeCluster(14)
  clusterExport(cl=cl,c('genCommYield'))
  d <- parLapply(cl,1:1000,function(i,d,m){
    source('../helperFunctions.R') #Helper functions
    genCommYield(dat = d, mods = m, simPar=TRUE)},d=datList,m=modSummaries_commodity) %>% 
    bind_rows(.id='rep')
  stopCluster(cl)
}

d %>% group_by(hbeeVis) %>% 
  summarize(med=median(yield_tha),lwr=quantile(yield_tha,0.05),upr=quantile(yield_tha,0.95)) %>% 
  ggplot(aes(x=hbeeVis*6))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=med))+
  labs(x='Honey bee visits per hour',y='Simulated yield (t/ha)')


# simulation for seed fields - not done yet
genSeedYield <- function(dat,mods=modSummaries_seed,
                         avgLogHbeeDist=4.319893,
                         avgLogLbeeDist=3.181154,
                         
                         avgFlDens=21.03358,
                         avgLogitFlwSurv=0.7344765,
                         sim=FALSE){
  
  dat <- list()
  
  library(dplyr)
  #Exogenous variables
  if(!'hbeeDist' %in% names(dat)) dat$hbeeDist <- exp(avgLogHbeeDist) #Hbee dist
  dat$clogHbeeDist <- log(dat$hbeeDist) - avgLogHbeeDist #centered log dist
  
  if(!'numHives' %in% names(dat)) dat$numHives <- 40
  dat$logNumHives <- log(dat$numHives+1) 
  
  #Endogenous variables
  if((!'plDens' %in% names(dat))|(('plDens' %in% names(dat))&any(is.na(dat$plDens)))){
    m <- list('intPlDens'=1,'slopeHbeeDistPlDens'=dat$clogHbeeDist) %>%
      getPreds(mods[['flDens']],parList = .,trans='exp',expandGrid=FALSE) %>% 
      pull(mean)  
    if(sim){ #Simulate
      s <- with(mods[['flDens']], mean[param=='sigmaPlDens'])
      m <- exp(rnorm(length(m),log(m),s)) 
    }
    na <- is.na(dat$plDens) #NAs present?
    if(any(na)){
      warning('NAs filled in plDens')
      dat$plDens[na] <- m[na]
    } else {
      dat$plDens <- m
    }
  }
  dat$logPlDens <- log(dat$plDens)
  
  if(!'plSize' %in% names(dat)){
    m <- list('intPlSize'=1,'slopePlDensPlSize'=dat$logPlDens,
              'slopeHbeeDistPlSize'= dat$clogHbeeDist) %>%
      getPreds(mods[['flDens']],parList = .,trans='exp',expandGrid=FALSE) %>% 
      pull(mean)
    if(sim){ #Simulate
      s <- with(mods[['flDens']], mean[param=='sigmaPlSize'])
      m <- exp(rnorm(length(m),log(m),s)) 
    }
    dat$plSize <- m
  }
  dat$logPlSize <- log(dat$plSize)
  
  if(!'flDens' %in% names(dat)){
    m <- list('intFlDens'=1,'slopePlSizeFlDens'=dat$logPlSize,
              'slopeHbeeDistFlDens'= dat$clogHbeeDist,
              'slopePlDensFlDens'=dat$logPlDens) %>%
      getPreds(mods[['flDens']],parList = .,expandGrid=FALSE) %>% 
      pull(mean)
    if(sim){ #Simulate
      s <- with(mods[['flDens']], mean[param=='sigmaFlDens'])
      m <- rnorm(length(m),m,s)
    }
    dat$flDens <- (m + avgFlDens)^2
  }
  dat$csqrtFlDens <- sqrt(dat$flDens)-avgFlDens
  
  #Visitation
  if(!'hbeeVis' %in% names(dat)){
    m <- list('intHbeeVis'=1,
              'slopeHbeeDistHbeeVis'=dat$clogHbeeDist,
              'slopeNumHivesHbeeVis'= dat$logNumHives, 
              'slopeFlDensHbeeVis' = dat$csqrtFlDens) %>% 
      getPreds(mods[['hbeeVis']],parList = .,offset=1,trans='exp',expandGrid=FALSE) %>% 
      pull(mean)
    if(sim){ #Simulate
      s <- with(mods[['hbeeVis']], mean[param=='phiHbeeVis'])
      m <- rnbinom(length(m),mu=m,size=s)
    }
    dat$hbeeVis <- m
  }
  dat$logHbeeVis <- log(dat$hbeeVis+1) 
  
  #Pollen
  if(!'pollen' %in% names(dat)){
    m <- list('intPollen'=1,'slopeHbeeVisPollen'=dat$logHbeeVis,
              'slopeHbeeDistPollen'= dat$clogHbeeDist) %>% 
      getPreds(mods[['pollen']],parList = .,trans='exp',expandGrid=FALSE) %>% 
      pull(mean)  
    if(sim){ #Simulate
      s <- with(mods[['pollen']], mean[param=='pollenPhi'])
      m <- rnbinom(length(m),mu=m,size=s)
    }
    dat$pollen <- m
  }
  dat$clogPollen <- log(dat$pollen+1)-with(mods[['pollen']],mean[param=='intPollen'])
  
  #Flower survival to pod (proportion)
  if(!'flwSurv' %in% names(dat)){
    m <- list('intFlwSurv'=1,'slopeHbeeVisFlwSurv'=dat$logHbeeVis,
              'slopePlSizeFlwSurv'= dat$logPlSize,
              'slopePollenFlwSurv'=dat$clogPollen) %>% 
      getPreds(mods[['flwSurv']],parList = .,trans='invLogit',expandGrid=FALSE) %>% 
      pull(mean)    
    if(sim){ #Simulate  
      s <- exp(with(mods[['flwSurv']], mean[param=='intPhiFlwSurv']))
      a <- invLogit(m)*s; b <- (1-invLogit(m))*s #beta binomial parameters
      m <- rbeta(length(m),a,b)
    }
    dat$flwSurv <- m
  }
  dat$clogitFlwSurv <- logit(dat$flwSurv)-avgLogitFlwSurv #centered logit flw surv
  
  #Flower count
  if(!'flwCount' %in% names(dat)){
    m <- list('intFlwCount'=1,'slopePlSizeFlwCount'=dat$logPlSize,
              'slopeFlwSurvFlwCount'= dat$clogitFlwSurv) %>% 
      getPreds(mods[['flwCount']],parList = .,trans='exp',expandGrid=FALSE) %>% 
      pull(mean)  
    if(sim){ #Simulate
      s <- list('intPhiFlwCount'=1,'slopePlSizePhiFlwCount'=dat$logPlSize) %>% 
        getPreds(mods[['flwCount']],parList = .,trans='exp',expandGrid=FALSE) %>% 
        pull(mean)
      m <- rnbinom(length(m),mu=m,size=s)
    }
    dat$flwCount <- m
  }
  dat$logFlwCount <- log(dat$flwCount)
  
  #Seeds per pod
  if(!'seedCount' %in% names(dat)){
    m <- list('intSeedCount'=1,
              'slopeHbeeVisSeedCount'=dat$logHbeeVis,
              'slopePollenSeedCount'= dat$clogPollen,
              'slopePlSizeSeedCount'=dat$logPlSize,
              'slopeFlwSurvSeedCount'=dat$clogitFlwSurv,
              'slopeFlwCountSeedCount'=dat$logFlwCount) %>% 
      getPreds(mods[['avgCount']],parList = .,otherPars = 'lambdaSeedCount',expandGrid=FALSE) %>% 
      mutate(across(c(mean:upr),~.x+1/lambdaSeedCount)) %>% pull(mean)
    if(sim){ #Simulate  
      s <- with(mods[['avgCount']], mean[param=='sigmaSeedCount'])
      l <- with(mods[['avgCount']], mean[param=='lambdaSeedCount'])
      m <- rnorm(length(m),m,s) + rexp(length(m),l)
    }
    dat$seedCount <- m
  }
  
  #Seed size
  if(!'seedWeight' %in% names(dat)){
    m <- list('intSeedWeight'=1,
              'slopeHbeeVisSeedWeight'=dat$logHbeeVis,
              'slopePollenSeedWeight'= dat$clogPollen,
              'slopeSeedCountSeedWeight'=dat$seedCount,
              'slopePlSizeSeedWeight'=dat$logPlSize,
              'slopePlDensSeedWeight'=dat$logPlDens,
              'slopeHbeeDistSeedWeight'=dat$clogHbeeDist,
              'slopeFlwSurvSeedWeight'=dat$clogitFlwSurv, #centered logit flw surv
              'slopeFlwCountSeedWeight'=dat$logFlwCount #Log flower count
    ) %>% 
      getPreds(mods[['avgWeight']],parList = .,otherPars = 'lambdaSeedWeight',expandGrid=FALSE) %>% 
      mutate(across(c(mean:upr),~.x+(1/lambdaSeedWeight))) %>% pull(mean)
    if(sim){ #Simulate
      s <- with(mods[['avgWeight']], mean[param=='sigmaSeedWeight'])
      l <- with(mods[['avgWeight']], mean[param=='lambdaSeedWeight'])
      m <- rnorm(length(m),m,s) + rexp(length(m),l)
    }
    dat$seedWeight <- m
  }
  
  #Yield g/plant
  dat$calcYield <- with(dat,flwSurv*flwCount*seedCount*seedWeight/1000)
  dat$logCalcYield <- log(dat$calcYield)
  
  m <- list('intYield'=1,'slopeYield'=dat$logCalcYield) %>% 
    getPreds(mods[[8]],parList=.,trans='exp',expandGrid=FALSE) %>% pull(mean)
  
  if(sim){ #Simulate
    s <- with(mods[['yield']], mean[param=='sigmaYield'])
    m <- exp(rnorm(length(m),log(m),s)) 
  }
  dat$yield <- m
  
  #Yield tonnes/ha
  dat$yield_tha <- (dat$plDens*dat$yield)/100
  
  chooseThese <- names(dat) %in% c('hbeeDist','numHives','plDens','plSize','flDens','hbeeVis',
                                   'pollen','flwSurv','flwCount','seedCount','seedWeight','yield','yield_tha')
  dat <- dat[chooseThese]
  return(dat)
}

