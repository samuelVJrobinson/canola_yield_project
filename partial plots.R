# MAKES PARTIAL EFFECTS PLOTS FOR SUB-MODELS

# Load everything ---------------------------------------------------------
library(tidyverse)
theme_set(theme_bw())
setwd('~/Documents/canola_yield_project/Figures') #Galpern machine path

source('../helperFunctions.R') #Helper functions

load('../Models/datalist_commodity.Rdata') #Commodity data
commData <- datalist
load('../Models/datalist_seed.Rdata') #Seed data
seedData <- datalist; rm(datalist)

avgCommData <- lapply(commData,function(x) data.frame(min=min(x,na.rm=TRUE),med=median(x,na.rm=TRUE),max=max(x,na.rm=TRUE))) #Summary stats
avgSeedData <- lapply(seedData,function(x) data.frame(min=min(x,na.rm=TRUE),med=median(x,na.rm=TRUE),max=max(x,na.rm=TRUE))) #Summary stats

load('../Models/modSummaries_commodity.Rdata') #Commodity data models
load('../Models/modSummaries_seed.Rdata') #Seed data models

# Visitation plots --------------------------------------------------------

#Hbee + Lbee visitation from edge of field/

with(avgCommData,
     list('intHbeeVis'=1,
     'slopeHbeeDistHbeeVis'=seq(min(log(dist)),max(log(400)),length=10),
     'slopeNumHivesHbeeVis'= log(numHives$max), 'slopeFlDensHbeeVis' = avgCommData$flDens$med)) %>% 
  getPreds(modSummaries_commodity[[2]],parList = .,offset=1,ZIpar = 'phiHbeeVis',trans='exp',q=c(0.5,0.05,0.95)) 

  #   ggplot(aes(x=exp(slopeHbeeDistHbeeVis),y=med))+
  # geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  # geom_line()+
  # labs(x='Distance from honeybee hives')

