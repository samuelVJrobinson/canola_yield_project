library(ggplot2)
library(dplyr)
library(tidyr)

setwd('~/Projects/UofC/canola_yield_project/')

# Commodity field data -----------------

#Dryad commodity field visitation data - from other paper
tempSurvey <- read.csv('C://Users//Samuel//Desktop//doi_10.5061_dryad.vhhmgqnvj__v3//commodityVisitation.csv')

# data.frame(FieldNumber=tempSurvey$Field,FieldName=paste(surveyAllComm$Field,surveyAllComm$Year)) %>% 
#   unique %>% arrange(FieldName) #Numbered field IDs from Field name + Year

load("./Commodity field analysis/commodityfieldDataAll.RData")

rm(AICc,brix2mgul,deltaAIC,DIC,plotFixedTerms,predFixef,se,
   varComp,zeroCheck,conversion,visitorsAll,visitors2015)
fieldsAllComm <- fieldsAll; flowersAllComm <- flowersAll; plantsAllComm <- plantsAll;
seedsAllComm <- seedsAll; surveyAllComm <- surveyAll;
rm(fieldsAll,flowersAll,plantsAll,seedsAll,surveyAll)

#Unique field names
fieldNames <- unique(paste(surveyAllComm$Field,surveyAllComm$Year))

#Commodity plant data
seedTemp <- seedsAllComm %>% select(Field,Year,Distance,Plant,Pod:PodMass) %>% 
  mutate(Field=factor(paste(Field,Year),levels=fieldNames)) %>%
  mutate(Field=as.numeric(Field)) %>% 
  group_by(Field,Year,Distance,Plant,Pod) %>% 
  mutate(N=1:n()) %>% filter(N==1) %>% ungroup() %>% select(-N) %>% 
  pivot_longer(c(PodCount,PodMass)) %>% 
  pivot_wider(values_from=value,names_from = c(name,Pod),names_sep = '') 

plantsAllComm %>% 
  mutate(Field=factor(paste(Field,Year),levels=fieldNames)) %>%
  mutate(Field=as.numeric(Field)) %>% 
  mutate(SeedCount=ifelse(is.nan(SeedCount),NA,SeedCount)) %>% 
  mutate(AvPodCount=round(AvPodCount,1)) %>% 
  select(-contains('SEPod')) %>% 
  select(-BeeYard:-PropMissing,-Yield) %>% 
  select(-Area,-Variety,-Irrigated,-contains('AvPod')) %>% 
  left_join(seedTemp,by=c('Field','Year','Distance','Plant')) %>% 
  write.csv(.,'./DryadData/commodityPlants.csv',row.names = FALSE)

# Seed field data --------------

rm(list=ls())

#Seed field data
load("./Seed field analysis/seedfieldDataAll.RData")
fieldsAllSeed <- data.frame(allFields); plantsAllSeed <- data.frame(allPlants)
pollenAllSeed <- data.frame(allPollen); seedsAllSeed <- data.frame(allSeeds);
surveyAllSeed <- data.frame(allSurvey);
plantsAllSeed$Field <- gsub('Unrah','Unruh',plantsAllSeed$Field) #Fixes spelling error
seedsAllSeed$Field <- gsub('Unrah','Unruh',seedsAllSeed$Field)
seedsAllSeed$EdgeCent <- ifelse(seedsAllSeed$EdgeCent=='Cent','Center',seedsAllSeed$EdgeCent) #Fixes cent/center
rm(allFields,allPlants,allPollen,allSeeds,allSurvey,behav2015,visitors2016,nectar2016,folder)
#Set 'negative' missing pods (mistake in counting) to NA.
plantsAllSeed <- mutate(plantsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
seedsAllSeed <- mutate(seedsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
#Extra seed field data (from Riley W.)
rileyExtra <- read.csv('./Seed field analysis/rileyExtra.csv')

#Fix date structure in Riley's fields
rileyExtra$date <- as.character(rileyExtra$date)
rileyExtra$date[grepl('-',rileyExtra$date)] <- as.character(as.Date(rileyExtra$date,format='%d-%b-%y')[grepl('-',rileyExtra$date)])
rileyExtra$date[grepl('/',rileyExtra$date)] <- as.character(as.Date(rileyExtra$date,format='%m/%d/%Y')[grepl('/',rileyExtra$date)])
rileyExtra$date <- as.Date(rileyExtra$date,format='%F')

rileyExtra$site <- factor(rileyExtra$site)
samFields <- levels(fieldsAllSeed$Field) #Fields that I surveyed
rileyFields <- levels(rileyExtra$site)[!(levels(rileyExtra$site) %in% levels(fieldsAllSeed$Field))] #Removes fields that I also surveyed (only unique fields)

fieldsAllSeed$Field <- factor(fieldsAllSeed$Field,levels=c(samFields,rileyFields)) #Apply new ordering of fields to dataframes
rileyExtra$site <- factor(rileyExtra$site,levels=c(samFields,rileyFields))
surveyAllSeed$Field <- factor(surveyAllSeed$Field,levels=c(samFields,rileyFields))
pollenAllSeed$Field <- factor(pollenAllSeed$Field,levels=c(samFields,rileyFields))
plantsAllSeed$Field <- factor(plantsAllSeed$Field,levels=c(samFields,rileyFields))
seedsAllSeed$Field <- factor(seedsAllSeed$Field,levels=c(samFields,rileyFields))

rileyExtra <- rileyExtra %>%
  #Filter out NA hdist/ldist plots - hard to impute, and no downstream info
  filter(!is.na(hdist),!is.na(ldist)) %>%
  filter(hdist<=400) #Get rid of plots >400m away from bees (other side of the field)

rileyFieldsYear <- distinct(select(rileyExtra,Year,site,date))$Year[match(rileyFields,distinct(select(rileyExtra,Year,site,date))$site)]
rileyFieldsYear[is.na(rileyFieldsYear)] <- 2016 #Fixes one year

#Seed field visitation data
surveyAllSeed <- surveyAllSeed %>% 
  select(-Treatment,-EdgeDir,-Variety,-hbeePol,-hbeeNec) %>% 
  separate(StartTime,c('Date','StartTime'),sep=' ') %>% 
  mutate(Field=as.numeric(Field))
  
rileyExtra <- rileyExtra %>% select(-treatment) %>% 
  transmute(Field=as.numeric(site),Distance=hdist,minDist=ldist,Bay=factor(Bay,labels=c('F','M')),
            EdgeCent='Edge',Date=date,StartTime=NA,
            TotalTime=10,FlDens=flDens,lbee=lbee_vis,hbee=hbee_vis,otherBee=NA,hFly=NA,
            AirTemp=NA,WindSp=NA,RH=NA,Year,PlDens=NA) 
write.csv(rileyExtra,'./DryadData/seedVisitation_waytes.csv',row.names = FALSE)

#Survey data (+ pollen)
pollenAllSeed %>% mutate(Field=as.numeric(Field)) %>% 
  select(Field,Year,Distance,EdgeCent,Pollen) %>% 
  group_by(Field,Year,Distance,EdgeCent) %>% 
  mutate(p=paste0('Pollen',1:n())) %>% 
  pivot_wider(names_from=p,values_from=Pollen) %>%
  left_join(surveyAllSeed,.,by=c('Field','Year','Distance','EdgeCent')) %>% 
  mutate(across(contains('Pollen'),~ifelse(Bay=='M',NA,.x))) %>% 
  mutate(StartTime=gsub('\\:00$','',as.character(StartTime))) %>% 
  write.csv(.,'./DryadData/seedVisitation.csv',row.names = FALSE)

#Plant data

seedTemp <- seedsAllSeed %>% 
  mutate(Field=as.numeric(Field)) %>%
  select(Field,Year,Distance,EdgeCent,Plant,Pod:PodMass) %>%
  # group_by(Field,Year,Distance,EdgeCent,Plant,Pod) %>%
  # mutate(N=1:n()) %>% filter(N==1) %>% ungroup() %>% select(-N) %>%
  pivot_longer(c(PodCount,PodMass)) %>%
  pivot_wider(values_from=value,names_from = c(name,Pod),names_sep = '')

plantsAllSeed %>% 
  mutate(Field=as.numeric(Field)) %>%
  select(-contains('SEPod'),-EdgeDir,-Variety,-contains('AvPod')) %>%
  select(-lbee:-TotalTime) %>% 
  left_join(seedTemp,by=c('Field','Year','Distance','EdgeCent','Plant')) %>%
  group_by(Field,Distance,EdgeCent) %>% mutate(Plant=1:n()) %>% ungroup() %>%  #Reorder plant ID
  write.csv(.,'./DryadData/seedPlants.csv',row.names = FALSE)



