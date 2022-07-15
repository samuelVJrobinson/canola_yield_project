# MAKES PARTIAL EFFECTS PLOTS FOR SUB-MODELS

# Load everything ---------------------------------------------------------
library(tidyverse)
theme_set(theme_bw())
# setwd('~/Documents/canola_yield_project/Figures') #Galpern machine path
setwd('~/Projects/UofC/canola_yield_project/Models') #Multivac path

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

load('../Models/modSummaries_commodity.Rdata') #Commodity data models
load('../Models/modSummaries_seed.Rdata') #Seed data models

# Visitation plots --------------------------------------------------------

#Hbee + Lbee visitation from edge of field/

d1 <- with(avgCommData,
     list('intHbeeVis'=1,
     'slopeHbeeDistHbeeVis'=seq(min(log(dist)),max(log(400)),length=10),
     'slopeNumHivesHbeeVis'= log(40), 
     'slopeFlDensHbeeVis' = avgCommData$flDens$mean)) %>% 
  getPreds(modSummaries_commodity[[2]],parList = .,offset=1,trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(type='comm',dist=slopeHbeeDistHbeeVis,med,lwr,upr)

d2 <- with(avgSeedData,
     list('intHbeeVis'=1,
          'slopeFlDensHbeeVis' = flDens_obs$mean,
          'slopeHbeeDistHbeeVis'=seq(min(clogHbeeDist),max(clogHbeeDist),length=10),
          'slopeLbeeDistHbeeVis'=clogHbeeDist$mean,
          'slopeCentHbeeVis'= 0
          )) %>% 
  getPreds(modSummaries_seed[[2]],parList = .,offset=1,ZIpar = 'thetaHbeeVis',trans='exp',q=c(0.5,0.05,0.95)) %>% 
  transmute(type='seed',dist=with(avgSeedData,
                                  seq(min(clogHbeeDist),max(clogHbeeDist),length=10)),med,lwr,upr)

bind_rows(d1,d2) %>% 
  ggplot(aes(x=dist,y=med,col=type))+geom_line()


  #   ggplot(aes(x=exp(slopeHbeeDistHbeeVis),y=med))+
  # geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  # geom_line()+
  # labs(x='Distance from honeybee hives')

