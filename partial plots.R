# MAKES PARTIAL EFFECTS PLOTS FOR SUB-MODELS
# ALSO CREATES SUMMARY TABLES OF VARIABLES (DATA) AND PARAMETERS FROM MODELS

# Load everything ---------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(metR)
library(parallel)
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

#Hbee visitation at edge given different stocking rates
with(avgCommData,
           list('intHbeeVis'=1,
                'slopeHbeeDistHbeeVis'=log(1),
                'slopeNumHivesHbeeVis'= log(c(0,20,40)+1), 
                'slopeFlDensHbeeVis' = 0)) %>% 
  getPreds(modSummaries_commodity[['hbeeVis']],parList = .,offset=timeOffset,
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
  getPreds(modSummaries_commodity[['hbeeVis']],parList = .,offset=timeOffset,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(dists),mean,med,lwr,upr)

d2 <- with(avgSeedData, #Hbee visitation in seed fields
     list('intHbeeVis'=1,
          'slopeFlDensHbeeVis' = flDens_obs$mean,
          'slopeHbeeDistHbeeVis'=dists2,
          'slopeLbeeDistHbeeVis'= 0,
          'slopeCentHbeeVis'= 0
          )) %>% 
  getPreds(modSummaries_seed[['hbeeVis']],parList = .,offset=timeOffset,
           ZIpar = 'thetaHbeeVis',trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(dists),mean,med,lwr,upr)

# bind_rows(d1,d2,.id='type') %>% mutate(type=factor(type,labels=c('Commodity','Seed'))) %>% 
#   ggplot(aes(x=dist,y=mean))+
#   geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
#   geom_line(aes(col=type))+
#   labs(x='Distance from apiary (m)',y='Visits per hour',
#        fill='Field\nType',col='Field\nType')

d3 <- with(avgSeedData, #Leafcutters in seed fields
                  list('intVisitLbeeVis'=1,
                       'slopeFlDensLbeeVis' = 0,
                       'slopeHbeeDistLbeeVis'=dists2,
                       'slopeLbeeDistLbeeVis'= 0,
                       'slopeCentLbeeVis'= 0
                  )) %>% 
         getPreds(modSummaries_seed[['lbeeVis']],parList = .,offset=timeOffset,
                  ZIpar = 'thetaLbeeVis',trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=exp(dists),mean,med,lwr,upr)

(p1 <- bind_rows(d2,d3,d1,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  coord_cartesian(y=c(0,150))+
  labs(x='Distance from apiary (m)',y='Visits per hour',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols))


#Hbee and Lbee visitation from shelters
dists <- log(1:30)-avgSeedData$logLbeeDistCent #log Shelter Distances

d1 <- with(avgSeedData, #Hbee visits
           list('intHbeeVis'=1,
                'slopeFlDensHbeeVis' = 0,
                'slopeHbeeDistHbeeVis'=0,
                'slopeLbeeDistHbeeVis'=dists,
                'slopeCentHbeeVis'= 0
           )) %>% 
  getPreds(modSummaries_seed[['hbeeVis']],parList = .,offset=timeOffset,ZIpar = 'thetaHbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=1:30,mean,med,lwr,upr)

# debugonce(getPreds)

d2 <- with(avgSeedData, #Lbee visits
           list('intVisitLbeeVis'=1,
                'slopeFlDensLbeeVis' = 0,
                'slopeHbeeDistLbeeVis'=0,
                'slopeLbeeDistLbeeVis'=dists,
                'slopeCentLbeeVis'= 0
           )) %>% 
  getPreds(modSummaries_seed[['lbeeVis']],parList = .,offset=timeOffset,ZIpar = 'thetaLbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(dist=1:30,mean,med,lwr,upr)

(p2 <- bind_rows(d1,d2,.id='type') %>% 
  mutate(type=factor(type,labels=lab[1:2])) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Distance from shelter (m)',y='Visits per hour',
       fill=NULL,col=NULL)+
  coord_cartesian(y=c(0,300))+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols))

#Center and edge of bays
d1 <- with(avgSeedData, #Hbee visits
           list('intHbeeVis'=1,
                'slopeFlDensHbeeVis' = 0,
                'slopeHbeeDistHbeeVis'=0,
                'slopeLbeeDistHbeeVis'=0,
                'slopeCentHbeeVis'= c(0,1)
           )) %>% 
  getPreds(modSummaries_seed[['hbeeVis']],parList = .,offset=timeOffset,ZIpar = 'thetaHbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(bayType=c('Edge','Centre'),mean,med,lwr,upr)

d2 <- with(avgSeedData, #Lbee visits
           list('intVisitLbeeVis'=1,
                'slopeFlDensLbeeVis' = 0,
                'slopeHbeeDistLbeeVis'=0,
                'slopeLbeeDistLbeeVis'=0,
                'slopeCentLbeeVis'= c(0,1)
           )) %>% 
  getPreds(modSummaries_seed[['lbeeVis']],parList = .,offset=timeOffset,ZIpar = 'thetaLbeeVis',
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(bayType=c('Edge','Centre'),mean,med,lwr,upr)

(p3 <- bind_rows(d1,d2,.id='type') %>% 
  mutate(type=factor(type,labels=lab[1:2])) %>% 
  mutate(bayType=factor(bayType,levels=c('Edge','Centre'))) %>% 
  ggplot(aes(x=bayType,col=type))+
  geom_pointrange(aes(y=mean,ymax=upr,ymin=lwr))+
  labs(x='Female bay position',y='Visits per hour',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols[1:2]))

(p <- ggarrange(p1,p2,p3,ncol=3,nrow=1,common.legend = TRUE,legend='bottom',labels = 'auto'))

ggsave('allVisits.png',p,bg='white',width = 12,height=4)

# Pollen plots ------------------------------------------------------------

#Visitation effects

#median, mean, max:
#0,1.72,36 max for commodity, 4,22.6,215 for seed
hbeeVisComm <- log(seq(0,max(with(commData,hbeeVis/totalTime)),length=30)+1)
hbeeVisSeed <- log((seq(0,100,length=100))+1)
lbeeVisSeed <- log((seq(0,100,length=100))+1)
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

(p1 <- bind_rows(d2,d3,d1,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=vis,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Visits per hour',y='Pollen grains per stigma',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols)+
  scale_y_log10())

#Distance effects
hbeeDistComm <- with(avgCommData,seq(clogHbeeDist$min,log(100)-avgCommData$logHbeeDist$mean,length=100))
hbeeDistSeed <- with(avgSeedData,seq(clogHbeeDist$min,log(100)-avgSeedData$logHbeeDist,length=100))
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

(p2 <- bind_rows(d2,d3,d1,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Distance from apiary or shelter (m)',y='Pollen grains per stigma',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+
  scale_fill_manual(values=labCols)+
  scale_y_log10())

(p3 <- with(avgSeedData, #Bay position
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
    geom_pointrange(aes(y=mean,ymax=upr,ymin=lwr),col='darkolivegreen3')+
    labs(x='Female bay position',y='Pollen grains per stigma',fill=NULL,col=NULL))

(p <- ggarrange(p1,p2,p3,ncol=3,nrow=1,legend='bottom',common.legend = TRUE,labels='auto'))

ggsave('allPollen.png',p,bg='white',width = 12,height=4)

# data.frame(x=1,y=1:10) %>% ggplot(aes(x=x,y=y,col=y))+geom_point()+scale_color_gradient(low='darkorange',high='darkgreen')

# Seed production - pollen ---------------------------------------------------

#pollenPlot means - Taken from model output
commPollenPlot <- seq(-1.2644033,0.9869057,length=20) 
modSummaries_commodity[[3]]$summary$mean[1] #Intercept
seedPollenPlot <- seq(-3.573654,2.364599,length=20) 
modSummaries_seed[[4]]$summary$mean[1] #Intercept

lab <- c('Commodity','Seed')
labCols <- c('black','darkolivegreen3')

#Flower survival
d1 <- with(avgCommData, 
     list('intFlwSurv'=1,
          'slopeHbeeVisFlwSurv'=avgCommData$logHbeeVis$mean,
          'slopePlSizeFlwSurv'= avgCommData$plantSize$mean,
          'slopePollenFlwSurv'=commPollenPlot
          )) %>% 
  getPreds(modSummaries_commodity[[5]],trans='invLogit',parList = .,q=c(0.5,0.05,0.95)) %>% 
  transmute(pol=exp(commPollenPlot+modSummaries_commodity[[3]]$summary$mean[1]),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
     list('intFlwSurv'=1,
          'slopePollenFlwSurv'=seedPollenPlot,
          'slopePlSizeFlwSurv'=avgSeedData$plantSize$mean,
          'slopeCentFlwSurv'=0,
          'slopeHbeeDistFlwSurv'=0,
          'slopeFlDensFlwSurv'=0
     )) %>% 
  getPreds(modSummaries_seed[[6]],parList = .,trans='invLogit',q=c(0.5,0.05,0.95)) %>% 
  transmute(pol=exp(seedPollenPlot+modSummaries_seed[[4]]$summary$mean[1]),mean,med,lwr,upr) 

(p1 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  mutate(across(c(mean:upr),~.x*100)) %>% 
  ggplot(aes(x=pol,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Pollen grains per stigma',y='Flower survival (%)',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
  scale_x_log10()+theme(legend.position = c(0.75,0.2),legend.background = element_rect(fill='white',colour='black',size=0.1)))

#Seeds per pod
d1 <- with(avgCommData, 
           list('intSeedCount'=1,
                'slopeHbeeVisSeedCount'=avgCommData$logHbeeVis$mean,
                'slopePollenSeedCount'= commPollenPlot,
                'slopePlSizeSeedCount'=avgCommData$plantSize$mean,
                'slopeFlwSurvSeedCount'=0, #centered logit flw surv
                'slopeFlwCountSeedCount'=avgCommData$logFlwCount$mean #Log flower count
           )) %>% 
  getPreds(modSummaries_commodity[[6]],parList = .,q=c(0.5,0.05,0.95),ENpar = 'lambdaSeedCount') %>% 
  transmute(pol=exp(commPollenPlot+modSummaries_commodity[[3]]$summary$mean[1]),mean,med,lwr,upr) 

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
  getPreds(modSummaries_seed[[7]],parList = .,ENpar='lambdaSeedCount',q=c(0.5,0.05,0.95)) %>%
  transmute(pol=exp(seedPollenPlot+modSummaries_seed[[4]]$summary$mean[1]),mean,med,lwr,upr) 

(p2 <- bind_rows(d1,d2,.id = 'type') %>% 
    mutate(type=factor(type,labels=lab)) %>% 
    ggplot(aes(x=pol,y=mean))+
    geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
    geom_line(aes(col=type))+
    labs(x='Pollen grains per stigma',y='Seeds per pod',fill=NULL,col=NULL)+
    scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
    scale_x_log10()+theme(legend.position = 'none'))

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
  getPreds(modSummaries_commodity[[7]],parList = .,q=c(0.5,0.05,0.95),ENpar = 'lambdaSeedWeight') %>% 
  transmute(pol=exp(commPollenPlot+modSummaries_commodity[[3]]$summary$mean[1]),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
           list('intSeedWeight'=1,
                'slopePollenSeedWeight'=seedPollenPlot,
                'slopeSeedCountSeedWeight'=avgSeedData$seedCount_obs$mean,
                'slopePlSizeSeedWeight'=avgSeedData$plantSize$mean,
                'slopePlDensSeedWeight'=avgSeedData$plDens_obs$mean,
                'slopeLbeeDistSeedWeight'=0
           )) %>% 
  getPreds(modSummaries_seed[[8]],parList = .,ENpar='lambdaSeedWeight',q=c(0.5,0.05,0.95)) %>% 
  transmute(pol=exp(seedPollenPlot+modSummaries_seed[[4]]$summary$mean[1]),mean,med,lwr,upr) 

(p3 <- bind_rows(d1,d2,.id = 'type') %>% 
    mutate(type=factor(type,labels=lab)) %>% 
    ggplot(aes(x=pol,y=mean))+
    geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
    geom_line(aes(col=type))+
    labs(x='Pollen grains per stigma',y='Seed size (mg)',fill=NULL,col=NULL)+
    scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
    scale_x_log10()+theme(legend.position = 'none'))

# (p <- ggarrange(p1,p2,p3,ncol=3,common.legend = TRUE,legend='bottom'))
# ggsave('allSeeds_pollen.png',p,bg='white',width = 12,height=4)

# Seed production - plant size ---------------------------------------------------

#plSize - both log-transformed
commPlSize <- with(avgCommData$plantSize,seq(min,max,length=20)) 
seedPlSize <- with(avgSeedData$plantSize,seq(min,max,length=20))

lab <- c('Commodity','Seed')
labCols <- c('black','darkolivegreen3')

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

(p4 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  mutate(across(c(mean:upr),~.x*100)) %>% 
  ggplot(aes(x=plSize,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Plant size (g)',y='Flower survival (%)',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols) +
  scale_x_log10()+theme(legend.position = 'none'))

#Seeds per pod
d1 <- with(avgCommData, 
           list('intSeedCount'=1,
                'slopeHbeeVisSeedCount'=avgCommData$logHbeeVis$mean,
                'slopePollenSeedCount'= 0,
                'slopePlSizeSeedCount'=commPlSize,
                'slopeFlwSurvSeedCount'=0, #centered logit flw surv
                'slopeFlwCountSeedCount'=avgCommData$logFlwCount$mean #Log flower count
           )) %>% 
  getPreds(modSummaries_commodity[[6]],parList = .,q=c(0.5,0.05,0.95),ENpar = 'lambdaSeedCount') %>% 
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
  getPreds(modSummaries_seed[[7]],parList = .,ENpar='lambdaSeedCount',q=c(0.5,0.05,0.95)) %>% 
  transmute(plSize=exp(seedPlSize),mean,med,lwr,upr) 

(p5 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=plSize,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Plant size (g)',y='Seeds per pod',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
  scale_x_log10()+theme(legend.position = 'none'))

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
  getPreds(modSummaries_commodity[[7]],parList = .,q=c(0.5,0.05,0.95),ENpar = 'lambdaSeedWeight') %>% 
  transmute(plSize=exp(commPlSize),mean,med,lwr,upr) 

d2 <- with(avgSeedData, 
           list('intSeedWeight'=1,
                'slopePollenSeedWeight'=0,
                'slopeSeedCountSeedWeight'=avgSeedData$seedCount_obs$mean,
                'slopePlSizeSeedWeight'=seedPlSize,
                'slopePlDensSeedWeight'=avgSeedData$plDens_obs$mean,
                'slopeLbeeDistSeedWeight'=0 #centered log lbee dist
           )) %>% 
  getPreds(modSummaries_seed[[8]],parList = .,ENpar='lambdaSeedWeight',q=c(0.5,0.05,0.95)) %>% 
  transmute(plSize=exp(seedPlSize),mean,med,lwr,upr) 

(p6 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=lab)) %>% 
  ggplot(aes(x=plSize,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Plant size (g)',y='Seed size (mg)',fill=NULL,col=NULL)+
  scale_colour_manual(values=labCols)+scale_fill_manual(values=labCols)+
    scale_x_log10()+theme(legend.position = 'none'))
  

# (p <- ggarrange(p1,p2,p3,ncol=3,common.legend = TRUE,legend='bottom'))
# ggsave('allSeeds_plSize.png',p,bg='white',width = 12,height=4)

(p <- ggarrange(p1,p2,p3,p4,p5,p6,nrow=2,ncol=3,common.legend = FALSE))
ggsave('allSeeds.png',p,bg='white',width = 10,height=6)

# Total yield and seed size for commodity fields -----------------------------------------------------

# avgCommData$flDensMean$mean #Avg sqrt flw dens
# avgCommData$logHbeeDist$mean #avg log hbee dist
# avgCommData$logitFlwSurv$mean #avg logit flw surv

# debugonce(genCommYield)
# debugonce(getPreds)
# debugonce(rnormLim)
d <- list()
d_mean <- genCommYield(dat=d,simPar=TRUE)
d <- replicate(100,genCommYield(dat=d,simPar = TRUE),simplify = FALSE) %>% bind_rows()

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
{ #Takes about 10 s
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
  geom_line(data=d_mean,aes(x=hbeeDist,y=yield_tha,col=numHives))+
  labs(x='Hbee distance',y='Yield')
  
#Visitation effect

quantile(commData$hbeeVis/commData$totalTime,c(0.05,0.1,0.5,0.9,0.95))

datList <- list(
  hbeeVis=0:11,
  plDens=50 #Same plant density
) 

# debugonce(genCommYield)
# debugonce(getPreds)
d_mean <- genCommYield(datList) %>% as.data.frame() 

# debugonce(genCommYield)
# genCommYield(dat = datList, mods = modSummaries_commodity, simPar=TRUE)

{
  cl <- makeCluster(14)
  clusterExport(cl=cl,c('genCommYield'))
  d <- parLapply(cl,1:1000,function(i,d,m){
    source('../helperFunctions.R') #Helper functions
    genCommYield(dat = d, mods = m, simPar=TRUE)},d=datList,m=modSummaries_commodity) %>% 
    bind_rows(.id='rep')
  stopCluster(cl)
}

(p1 <- d %>% group_by(hbeeVis) %>% 
  summarize(med=median(seedWeight),lwr=quantile(seedWeight,0.05),upr=quantile(seedWeight,0.95)) %>% 
  ggplot(aes(x=hbeeVis*6))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(data=d_mean,aes(y=seedWeight))+
  labs(x='Honey bee visits per hour',y='Simulated seed size (mg/seed)'))

(p2 <- d %>% group_by(hbeeVis) %>% 
  summarize(med=median(yield_tha),lwr=quantile(yield_tha,0.05),upr=quantile(yield_tha,0.95)) %>% 
  ggplot(aes(x=hbeeVis*6))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(data=d_mean,aes(y=yield_tha))+
  labs(x='Honey bee visits per hour',y='Simulated yield (t/ha)'))

ggarrange(p1,p2,ncol=1,nrow=2)

# Total yield and seed size for seed fields -----------------------------------------------------

# avgSeedData$logHbeeDistCent #avg log hbee dist
# avgSeedData$logLbeeDistCent #avg logit flw surv
# avgSeedData$flDensMean$mean #Avg sqrt flw dens

d <- list()
d_mean <- genSeedYield(dat=d,simPar=FALSE)
# d <- replicate(100,genSeedYield(dat=d,simPar = TRUE),simplify = FALSE) %>% bind_rows()

#Distance effect
with(seedData,list(
  hbeeDist=c(hbee_dist,hbee_dist_extra),
  lbeeDist=c(lbee_dist,lbee_dist_extra))) %>% 
  lapply(.,function(x) quantile(x,c(0.05,0.1,0.5,0.9,0.95)))

summary(exp(seedData$plantSize))
summary(exp(seedData$plDens_obs))

datList <- expand.grid(hbeeDist=exp(seq(log(5),log(400),length.out=20)),
                       lbeeDist=exp(seq(log(4),log(60),length.out=20)),
                       plDens=38,plSize=25.2
                       ) %>% 
  lapply(.,function(x) x)

d_mean <- genSeedYield(datList,mods=modSummaries_seed) %>% as.data.frame()

b <- seq(3.2,3.6,0.05)
(p1 <- d_mean %>% 
    ggplot(aes(x=hbeeDist,y=lbeeDist,z=seedWeight))+
    geom_contour(col='black')+
    geom_label_contour(skip=0,breaks=b,fill='white',label.size=0)+
    labs(x='Honey bee distance',y='Leafcutter distance',fill='Simulated\nseed size\n(mg/seed)')+
    scale_fill_distiller())

b <- seq(2.5,3.7,0.1)
(p2 <- d_mean %>% 
    ggplot(aes(x=hbeeDist,y=lbeeDist,z=yield_tha))+
    geom_contour(col='black',breaks=b)+
    geom_label_contour(skip=1,breaks=b,fill='white',label.size=0)+
    labs(x='Honey bee distance',y='Leafcutter distance',fill='Simulated\nyield\n(t/ha)')
  )

ggarrange(p1,p2,ncol=2)

# #Faceted plot - can't control contours as easily
# d_mean %>% select(hbeeDist,lbeeDist,seedWeight,yield_tha) %>%
#   pivot_longer(-contains('Dist')) %>%
#   mutate(name=factor(name,labels=c('Seed size (mg/seed)','Yield (t/ha)'))) %>%
#   ggplot(aes(x=hbeeDist,y=lbeeDist,z=value))+
#   facet_wrap(~name,scales='free')+
#   geom_contour(col='black')+
#   geom_label_contour(skip=1,fill='white',label.size=0)+
#   labs(x='Honey bee distance',y='Leafcutter distance')

#Visitation effect alone
with(seedData,list( #Range of lbee/hbee visitations
  hbeeVis=(c(hbeeVis,hbeeVis_extra)/c(totalTime,totalTime_extra)),
  lbeeVis=(c(lbeeVis,lbeeVis_extra)/c(totalTime,totalTime_extra)))) %>% 
  lapply(.,function(x) quantile(x,c(0.05,0.1,0.5,0.9,0.95)))

datList <- expand.grid(hbeeVis=exp(seq(0,log(84),length.out=20)),
            lbeeVis=exp(seq(0,log(64),length.out=20)),
            plDens=38,plSize=25.2
            ) %>% lapply(.,function(x) x)

d_mean <- genSeedYield(datList,mods=modSummaries_seed) %>% as.data.frame()

b <- seq(3.41,3.44,0.005)
(p1 <- d_mean %>% 
  mutate(across(c(hbeeVis,lbeeVis),~.x*6)) %>% 
  ggplot(aes(x=hbeeVis,y=lbeeVis,z=seedWeight))+
  geom_contour(col='black',breaks=b)+
  geom_label_contour(skip=0,breaks=b,fill='white',label.size=0)+
  labs(x='Honey bee visits per hour',y='Leafcutter visits per hour')+
  scale_fill_distiller())

b <- seq(3.25,3.45,0.025)
(p2 <- d_mean %>% 
  mutate(across(c(hbeeVis,lbeeVis),~.x*6)) %>% 
  ggplot(aes(x=hbeeVis,y=lbeeVis,z=yield_tha))+
    geom_contour(col='black',breaks=b)+
    geom_label_contour(skip=0,breaks=b,fill='white',label.size=0)+
  labs(x='Honey bee visits per hour',y='Leafcutter visits per hour'))

ggarrange(p1,p2,ncol=1,nrow=2)


# (p1 <- d %>% group_by(hbeeVis) %>% 
#   summarize(med=median(seedWeight),lwr=quantile(seedWeight,0.05),upr=quantile(seedWeight,0.95)) %>% 
#   ggplot(aes(x=hbeeVis*6))+
#   geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
#   geom_line(aes(y=med))+# geom_line(data=d_mean,aes(y=yield_tha))+
#   labs(x='Honey bee visits per hour',y='Simulated seed size (mg/seed)'))
# 
# (p2 <- d %>% group_by(hbeeVis) %>% 
#   summarize(med=median(yield_tha),lwr=quantile(yield_tha,0.05),upr=quantile(yield_tha,0.95)) %>% 
#   ggplot(aes(x=hbeeVis*6))+
#   geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
#   # geom_line(data=d_mean,aes(y=yield_tha))+
#   geom_line(aes(y=med))+
#   labs(x='Honey bee visits per hour',y='Simulated yield (t/ha)'))

# Yield and seed size for both field types -------------------------------------

xaxis <- 'Distance from HB hives (m)'

#Commodity fields - distance and stocking effect
# summary(exp(commData$plantSize))
# summary(exp(commData$plDens_obs))

datList <- expand.grid(hbeeDist=exp(seq(log(1),log(400),length=30)),
                       plDens=48.5, #plSize=18,
                       numHives=c(0,40)
) %>% lapply(.,function(x) x)

d_mean <- genCommYield(datList) %>% as.data.frame() %>%
  mutate(numHives=factor(numHives))

{ #Takes about 10 s
  cl <- makeCluster(14) 
  clusterExport(cl=cl,c('genCommYield'))
  d <- parLapply(cl,1:1000,function(i,d,m){
    source('../helperFunctions.R') #Helper functions
    genCommYield(dat = d, mods = m, simPar=TRUE)},d=datList,m=modSummaries_commodity) %>% 
    bind_rows(.id='rep')
  stopCluster(cl)
  }

(p1 <- d %>% mutate(numHives=factor(numHives)) %>% group_by(hbeeDist,numHives) %>% 
  summarize(med=median(seedWeight),lwr=quantile(seedWeight,0.1),upr=quantile(seedWeight,0.9)) %>% 
  ggplot(aes(x=hbeeDist,y=med,group=numHives))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(data=d_mean,aes(x=hbeeDist,y=seedWeight))+
  annotate('text',x=c(300,300),y=c(2.85,2.77),label=c('40 hives','0 hives'))+
  labs(x=xaxis,y='Seed size (mg/seed)'))

(p2 <- d %>% mutate(numHives=factor(numHives)) %>% group_by(hbeeDist,numHives) %>% 
  summarize(med=median(yield_tha),lwr=quantile(yield_tha,0.1),upr=quantile(yield_tha,0.9)) %>% 
  ggplot(aes(x=hbeeDist,y=med,group=numHives))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(data=d_mean,aes(x=hbeeDist,y=yield_tha))+
  annotate('text',x=c(300,300),y=c(3.45,3.35),label=c('40 hives','0 hives'))+
  labs(x=xaxis,y='Yield (t/ha)'))

#Seed fields

#Distance effect
# with(seedData,list(
#   hbeeDist=c(hbee_dist,hbee_dist_extra),
#   lbeeDist=c(lbee_dist,lbee_dist_extra))) %>% 
#   lapply(.,function(x) quantile(x,c(0.05,0.1,0.5,0.9,0.95)))
# summary(exp(seedData$plantSize))
# summary(exp(seedData$plDens_obs))
yaxis <- "Distance from LCB shelters (m)"

datList <- expand.grid(hbeeDist=exp(seq(log(5),log(400),length.out=20)),
                       plDens=38,#plSize=25.2
                       lbeeDist=exp(seq(log(4),log(60),length.out=20))
                       
) %>% 
  lapply(.,function(x) x)

d_mean <- genSeedYield(datList,mods=modSummaries_seed) %>% as.data.frame()

b <- seq(3.2,3.6,0.05)
(p3 <- d_mean %>% 
    ggplot(aes(x=hbeeDist,y=lbeeDist,z=seedWeight))+
    geom_contour(col='black')+
    geom_label_contour(skip=0,breaks=b,fill='white',label.size=0)+
    labs(x=xaxis,y=yaxis,fill='Simulated\nseed size\n(mg/seed)')+
    scale_fill_distiller())

b <- seq(2.5,3.7,0.1)
(p4 <- d_mean %>% 
    ggplot(aes(x=hbeeDist,y=lbeeDist,z=yield_tha))+
    geom_contour(col='black',breaks=b)+
    geom_label_contour(skip=1,breaks=b,fill='white',label.size=0)+
    labs(x=xaxis,y=yaxis,fill='Simulated\nyield\n(t/ha)')
)

p <- ggarrange(p1,p2,p3,p4,labels = 'auto')
ggsave('allYield.png',p,height=10,width=10)

#Alternate version with biplots

datList <- list(plDens=38,#plSize=25.2
  hbeeDist=exp(seq(log(5),log(100),length.out=20)),
  lbeeDist=exp(avgSeedData$logLbeeDistCent))

{ #Takes about 10 s
  cl <- makeCluster(14) 
  clusterExport(cl=cl,c('genSeedYield'))
  d1 <- parLapply(cl,1:1000,function(i,d,m){
    source('../helperFunctions.R') #Helper functions
    genSeedYield(dat = d, mods = m, simPar=TRUE)},d=datList,m=modSummaries_seed) %>% 
    bind_rows(.id='rep')
  stopCluster(cl)
}

datList <- list(plDens=38,#plSize=25.2
                hbeeDist=exp(avgSeedData$logHbeeDistCent),
                lbeeDist=exp(seq(log(4),log(60),length.out=20)))

{ #Takes about 10 s
  cl <- makeCluster(14) 
  clusterExport(cl=cl,c('genSeedYield'))
  d2 <- parLapply(cl,1:1000,function(i,d,m){
    source('../helperFunctions.R') #Helper functions
    genSeedYield(dat = d, mods = m, simPar=TRUE)},d=datList,m=modSummaries_seed) %>% 
    bind_rows(.id='rep')
  stopCluster(cl)
}

p3 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=c('HB Dist','LCB Dist'))) %>% 
  mutate(dist=ifelse(type=='HB Dist',hbeeDist,lbeeDist)) %>% 
  group_by(type,dist) %>% 
  summarize(med=median(seedWeight),lwr=quantile(seedWeight,0.1),upr=quantile(seedWeight,0.9)) %>% 
  ggplot(aes(x=dist,y=med))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3,show.legend = FALSE)+
  geom_line(aes(col=type),show.legend = FALSE)+
  labs(x='Distance from apiary or shelter (m)',y='Seed size (mg/seed)',col='Visitor',fill='Visitor')+
  scale_colour_manual(values=c('darkorange','darkgreen'))+
  scale_fill_manual(values=c('darkorange','darkgreen'))

p4 <- bind_rows(d1,d2,.id = 'type') %>% 
  mutate(type=factor(type,labels=c('HB Dist','LCB Dist'))) %>% 
  mutate(dist=ifelse(type=='HB Dist',hbeeDist,lbeeDist)) %>% 
  group_by(type,dist) %>% 
  summarize(med=median(yield_tha),lwr=quantile(yield_tha,0.1),upr=quantile(yield_tha,0.9)) %>% 
  ggplot(aes(x=dist,y=med))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Distance from apiary or shelter (m)',y='Yield (t/ha)',col=NULL,fill=NULL)+
  scale_colour_manual(values=c('darkorange','darkgreen'))+
  scale_fill_manual(values=c('darkorange','darkgreen'))+
  theme(legend.position = c(0.8,0.8))
 
p <- ggarrange(p1,p2,p3,p4,labels = 'auto')
ggsave('allYield_alternate.png',p,height=10,width=10)

# General summary of variables for both models -----------------------

commDatSummary <- with(commData,list(
  'Number of hives'=numHives,'Distance to edge (m)'=dist,
  'HB visitation (hr$^{-1}$)'=6*hbeeVis/totalTime,
  'Flower density (m$^2$)'=(flDens+flDensMean)^2,
  'Pollen per stigma'=pollenCount,
  'Plant density (m$^2$)'=plDens_obs,
  'Plant vegetative mass (g)'=VegMass,
  'Plant seed mass (g)'=yield,
  'Flowers per plant'=flwCount,'Pods per plant'=podCount,'Seeds per pod'=seedCount,'Seed size (mg)'=seedMass
  )) %>% 
  lapply(.,function(x) data.frame(Mean=mean(x),Median=median(x),SD=sd(x),Min=min(x),Max=max(x))) %>% 
  bind_rows(.id='Variable') 

seedDatSummary <- with(seedData,list(
  'Distance to edge (m)'=c(hbee_dist,hbee_dist_extra),
  'Distance to LCB shelter (m)'=c(lbee_dist,lbee_dist_extra),
  'HB visitation (hr$^{-1}$)'=6*c(hbeeVis,hbeeVis_extra)/c(totalTime,totalTime_extra),
  'LCB visitation (hr$^{-1}$)'=6*c(lbeeVis,lbeeVis_extra)/c(totalTime,totalTime_extra),
  'Bay Edge/Centre'=c(isCent,isCent_extra),
  'Flower density (m$^2$)'=(c(flDens_obs,flDens_obs_extra)+flDensMean)^2,
  'Pollen per stigma'=pollenCount,
  'Plant density (m$^2$)'=plDens_obs,
  'Plant vegetative mass (g)'=exp(seedData$plantSize),
  'Plant seed mass (g)'=yield,
  'Flowers per plant'=flwCount,
  'Pods per plant'=podCount,
  'Seeds per pod'=seedCount_obs,'Seed size (mg)'=seedMass_obs
)) %>% 
  lapply(.,function(x) data.frame(Mean=mean(x),Median=median(x),SD=sd(x),Min=min(x),Max=max(x))) %>% 
  bind_rows(.id='Variable')
  
bind_rows(commDatSummary,seedDatSummary,.id='Field Type') %>% 
  mutate(`Field Type`=ifelse(`Field Type`==1,'Commodity','Seed')) %>% 
  xtable::xtable(.,digits=c(0,0,2,2,2,2,2,2)) %>% 
  print(.,include.rownames=FALSE,sanitize.text.function=identity)


# General summary of parameters for both models ---------------------------

lapply(modSummaries_commodity,function(x) x$summary) %>% #Get path coefficients
  bind_rows(.id='to') %>% 
  transmute(to,name=param,Z,pval) %>% 
  filter(to!='yield',!grepl('(int|sigma|lambda|Phi|phi|rho)',name)) %>% 
  mutate(to=case_when(
    to=='flDens' & grepl('PlDens$',name) ~ 'plDens',
    to=='flDens' & grepl('PlSize$',name) ~ 'plSize',
    TRUE ~ gsub('avg','seed',to)
  )) %>% 
  mutate(name=gsub('slope','',name)) %>% 
  filter(mapply(grepl,capFirst(to),name)) %>% 
  mutate(name=mapply(gsub,capFirst(to),'',name)) %>% 
  mutate(name=capFirst(name,TRUE)) 
  


