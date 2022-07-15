# MAKES PARTIAL EFFECTS PLOTS FOR SUB-MODELS

# Load everything ---------------------------------------------------------
library(tidyverse)
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
avgCommData <- lapply(commData,summaryDat) 
seedData$lbeeDistAll <- with(seedData,c(lbee_dist,lbee_dist_extra))
seedData$hbeeDistAll <- with(seedData,c(hbee_dist,hbee_dist_extra))

seedData$clogLbeeDist <- with(seedData,logCent(lbeeDistAll)) #Log-centered distances
seedData$clogHbeeDist <- with(seedData,logCent(hbeeDistAll))
avgSeedData <- lapply(seedData,summaryDat) 

avgSeedData$logHbeeDistCent <- mean(log(with(seedData,c(hbee_dist,hbee_dist_extra)))) #Center
avgSeedData$logLbeeDistCent <- mean(log(with(seedData,c(lbee_dist,lbee_dist_extra))))

load('../Models/modSummaries_commodity.Rdata') #Commodity data models
load('../Models/modSummaries_seed.Rdata') #Seed data models

# Visitation plots --------------------------------------------------------

timeOffset <- 6 #Use 1 hr offset (1 = 10 mins, 6 = 1 hr)

#Hbee stocking rate numbers
with(avgCommData,
           list('intHbeeVis'=1,
                'slopeHbeeDistHbeeVis'=log(1),
                'slopeNumHivesHbeeVis'= log(c(0,20,40)+1), 
                'slopeFlDensHbeeVis' = avgCommData$flDens$mean)) %>% 
  getPreds(modSummaries_commodity[[2]],parList = .,offset=timeOffset,
           trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(type='Commodity',numHives=c(0,20,40),mean,med,lwr,upr) 

#Hbee + Lbee visitation from edge of field
dists <- with(avgCommData,seq(min(log(dist)),max(log(400)),length=20)) #Distances

d1 <- with(avgCommData,
     list('intHbeeVis'=1,
     'slopeHbeeDistHbeeVis'=dists,
     'slopeNumHivesHbeeVis'= log(40+1), 
     'slopeFlDensHbeeVis' = avgCommData$flDens$mean)) %>% 
  getPreds(modSummaries_commodity[[2]],parList = .,offset=timeOffset,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(type='Commodity',dist=exp(dists),mean,med,lwr,upr)

d2 <- with(avgSeedData,
     list('intHbeeVis'=1,
          'slopeFlDensHbeeVis' = flDens_obs$mean,
          'slopeHbeeDistHbeeVis'=seq(min(clogHbeeDist),max(clogHbeeDist),length=20),
          'slopeLbeeDistHbeeVis'=clogLbeeDist$mean,
          'slopeCentHbeeVis'= 0
          )) %>% 
  getPreds(modSummaries_seed[[2]],parList = .,offset=timeOffset,ZIpar = 'thetaHbeeVis',trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(type='Seed',dist=exp(dists),mean,med,lwr,upr)

bind_rows(d1,d2) %>% 
  ggplot(aes(x=dist,y=mean))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=type),alpha=0.3)+
  geom_line(aes(col=type))+
  labs(x='Distance from apiary (m)',y='Visits per hour',
       fill='Field\nType',col='Field\nType')

#Hbee and Lbee visitation from shelters
  
dists <- log(1:30)-avgSeedData$logLbeeDistCent #Distances

d2 <- with(avgSeedData,
           list('intHbeeVis'=1,
                'slopeFlDensHbeeVis' = flDens_obs$mean,
                'slopeHbeeDistHbeeVis'=clogHbeeDist$mean,
                'slopeLbeeDistHbeeVis'=dists,
                'slopeCentHbeeVis'= 0
           )) %>% 
  getPreds(modSummaries_seed[[2]],parList = .,offset=timeOffset,ZIpar = 'thetaHbeeVis',trans='exp',q=c(0.5,0.05,0.95))

d2 <- with(avgSeedData,
           list('intHbeeVis'=1,
                'slopeFlDensHbeeVis' = flDens_obs$mean,
                'slopeHbeeDistHbeeVis'=seq(min(clogHbeeDist),max(clogHbeeDist),length=20),
                'slopeLbeeDistHbeeVis'=clogHbeeDist$mean,
                'slopeCentHbeeVis'= 0
           )) %>% 
  getPreds(modSummaries_seed[[2]],parList = .,offset=timeOffset,ZIpar = 'thetaHbeeVis',trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(type='Seed',dist=exp(dists),mean,med,lwr,upr)



  
