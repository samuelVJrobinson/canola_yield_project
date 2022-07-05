#Template for running claim sets
# wd <- getwd()
# setwd("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\canola_yield_project\\Models\\Seed model claims 3")

library(ggplot2)
library(dplyr)
library(tidyr)
library(beepr)
library(xtable)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 6)

#Big-text theme, no grid lines (used for Bayer 2016 presentation)
prestheme=theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15),
                axis.title=element_text(size=20),
                title=element_text(size=20),
                panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                panel.border=element_rect(size=1,colour='black'),
                strip.text=element_text(size=15))
theme_set(theme_bw()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

source("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\canola_yield_project\\helperFunctions.R")

load("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\canola_yield_project\\Models\\datalist_seed.Rdata")

#Get model names, and make list
modFiles <- dir(path='.',pattern = 'claim.+\\.stan',full.names = FALSE)
modFiles <- modFiles[!grepl('template',modFiles)]
modFiles <- modFiles[sapply(read.csv("claimsList_updated.csv")$Filename,function(x){
  l <- grep(x,modFiles)
  if(length(l)==0) 0 else l
})]

i <- 1

if(file.exists(modFiles[i])){
  print(paste0('Starting model ',modFiles[i]))
  #Run model
  # mod <- stan(file=modFiles[i],data=datalist,iter=3000,chains=4,control=list(adapt_delta=0.8),init=0)
  mod <- stan(file=modFiles[i],data=datalist,iter=400,chains=1,control=list(adapt_delta=0.8),init=0)
  temp <- parTable(mod) #Get parameter summaries
  #Save information to csv file
  modList <- read.csv('./Seed model claims 3/claimsList_updated.csv',sep=',',strip.white = TRUE)
  modList[i,match(names(temp),names(modList))] <- temp[grepl('claim',temp$param),]
  write.csv(modList,'./Seed model claims 3/claimsList_updated.csv',row.names = FALSE)
  
  #Check model
  parNam <- names(mod) #Model parameters
  
  #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
  useThese <- (!(grepl('(\\[[0-9]+,*[0-9]*\\]$|^lp|\\_miss)',parNam)) |
                 grepl('[sS]igma',parNam) 
  )
  parNam <- parNam[useThese] 
  
  #Removes pollen terms
  parNam <- parNam[!parNam %in% c('intPollen','slopeVisitPol','slopeHbeeDistPollen','sigmaPolField','pollenPhi')] 
  
  p <- traceplot(mod,pars=parNam,inc_warmup=FALSE)+labs(title=gsub('.*/','',modFiles[i]))
  ggsave(paste0('./',gsub('\\.stan','',modFiles[i]),'.png'),p,width=8,height=8)
  # print(mod,pars=parNam)
  # print(modList[i,match(names(temp),names(modList))])
  # fastPairs(mod,pars=parNam)
  # print(paste0('Model ',modFiles[i],' completed'))
  
} else print(paste0('Model ',i,' not found'))

setwd(wd)


