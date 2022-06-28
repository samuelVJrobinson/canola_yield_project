#Test of whether Riley's data, and then male bays from my data can be excluded without changing results

# Libraries and ggplot theme 
library(ggplot2)
theme_set(theme_classic())
library(dplyr)
library(tibble)
library(tidyr)
library(beepr)

library(xtable)

setwd('~/Projects/UofC/canola_yield_project') #Multivac path

source('helperFunctions.R')

# Load in data 

load('./Models/datalist_seed.Rdata')


#Same datalist, but no male bays

datalist2 <- datalist[!grepl('extra',names(datalist))]

datalist2 <- lapply(datalist2,function(x){
  isF <- !datalist$isMBay
  if(any(class(x)=='matrix')){
    return(x[isF,])
  } else if(length(x)==360){
    return(x[isF])
  } else {
    return(x)
  }
})

datalist2$Nplot <- sum(!datalist$isMBay)

#Get rid of flower and plant level measurements
datalist2 <- datalist2[1:(which(names(datalist2)=='Nflw')-1)]

#Fix flower density measurements
flDens <- rep(NA,360)
flDens[datalist$obsflDens_ind] <- datalist$flDens_obs
flDens <- flDens[!datalist$isMBay]
# str(datalist2)
datalist2$Nplot_flDensObs <- sum(!is.na(flDens))
datalist2$Nplot_flDensMiss <- sum(is.na(flDens))
datalist2$flDens_obs <- flDens[!is.na(flDens)]
datalist2$obsflDens_ind <- which(!is.na(flDens))
datalist2$missflDens_ind <- which(is.na(flDens))
rm(flDens)

library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE)
options(mc.cores = 6)

# Honey bee visitation ----------------------

#Models
modFiles <- c('seed_04hbeeVis.stan',
              'seed_04hbeeVis_dropMale.stan',
              'seed_04hbeeVis_noExtra.stan',
              'seed_04hbeeVis_noExtra_dropMale.stan',
              'seed_04hbeeVis_noMale.stan')
modList <- vector(mode = 'list',length = length(modFiles))
names(modList) <- c('original','dropMale','noExtra','noExtra_dropMale','noExtra_noMale')

#Original
modList[1] <- stan(file=modFiles[1],data=datalist,iter=2000,chains=4,init=0)
#Original - drop male term
modList[2] <- stan(file=modFiles[2],data=datalist,iter=2000,chains=4,init=0) 
#No extras
modList[3] <- stan(file=modFiles[3],data=datalist,iter=2000,chains=4,init=0)
#No extras - drop male term
modList[4] <- stan(file=modFiles[4],data=datalist,iter=2000,chains=4,init=0) 
#No extras - No male bays
modList[5] <- stan(file=modFiles[5],data=datalist2,iter=2000,chains=4,init=0) 

for(i in 1:length(modList)){ #Traceplots
  if(!is.null(modList[[i]])){
    n <- names(modList[[i]]) #Model parameters
    n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('[sS]igma',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
    n <- n[!grepl('_miss',n)] #Gets rid of imputed values
    p <- traceplot(modList[[i]],pars=n,inc_warmup=FALSE) #+ geom_hline(yintercept = 0) #Traceplots
    print(p)
    a <- readline('Press Return to continue: ')
    if(a!='') break()
  }
}

# PPplots(modList[[1]],c(datalist$hbeeVis,datalist$hbeeVis_extra),c('predHbeeVis_all','hbeeVis_resid','predHbeeVis_resid'),'Honeybee visits',jitterX=0.1)
# PPplots(modList[[2]],c(datalist$hbeeVis),c('predHbeeVis_all','hbeeVis_resid','predHbeeVis_resid'),'Honeybee visits',jitterX=0.1)
# PPplots(modList[[3]],c(datalist2$hbeeVis),c('predHbeeVis_all','hbeeVis_resid','predHbeeVis_resid'),'Honeybee visits',jitterX=0.1)

#Compare ranges of parameters
n <- names(modList[[1]]) #Model parameters
n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('[sS]igma',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
n <- n[!grepl('_miss',n)] #Gets rid of imputed values

(p <- lapply(modList,function(i){
  if(is.null(i)) return(NA)
  parnames <- n[n %in% names(i)]
  
  l <- t(sapply(extract(i,parnames),function(x) setNames(quantile(x,c(0.5,0.025,0.975)),c('median','lwr','upr'))))  
  l <- l %>% data.frame() %>% rownames_to_column('param')
}) %>% bind_rows(.id='mod') %>% 
  mutate(param=factor(param,levels=n),mod=factor(mod,levels=names(modList))) %>% 
  ggplot(aes(x=param,y=median,col=mod))+
  geom_pointrange(aes(ymax=upr,ymin=lwr),position=position_dodge(width=0.5))+
  facet_wrap(~param,scales='free')+
  geom_hline(yintercept = 0,linetype='dashed')+
  labs(x='Parameter',y='Value',col='Model',title='Honey bee visitation model'))+
  scale_colour_manual(values=c('black','blue','purple','darkorange','red'))

ggsave('../Figures/hbeeVisModels.png',p,width = 8,height=6)
  
#Leafcutter bee visitation ----------------------

#Models
modFiles <- c('seed_05lbeeVis.stan',
              'seed_05lbeeVis_dropMale.stan',
              'seed_05lbeeVis_noExtra.stan',
              'seed_05lbeeVis_noExtra_dropMale.stan',
              'seed_05lbeeVis_noMale.stan')

modList <- vector(mode = 'list',length = length(modFiles))
names(modList) <- c('original','dropMale','noExtra','noExtra_dropMale','noExtra_noMale')

#Original
modList[1] <- stan(file=modFiles[1],data=datalist,iter=2000,chains=4,init=0)
#Original - drop male term
modList[2] <- stan(file=modFiles[2],data=datalist,iter=2000,chains=4,init=0) 
#No extras
modList[3] <- stan(file=modFiles[3],data=datalist,iter=2000,chains=4,init=0)
#No extras - drop male term
modList[4] <- stan(file=modFiles[4],data=datalist,iter=2000,chains=4,init=0) 
#No extras - No male bays
modList[5] <- stan(file=modFiles[5],data=datalist2,iter=2000,chains=4,init=0) 

for(i in 1:length(modList)){
  if(!is.null(modList[[i]])){
    n <- names(modList[[i]]) #Model parameters
    n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('[sS]igma',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
    n <- n[!grepl('_miss',n)] #Gets rid of imputed values
    
    p <- traceplot(modList[[i]],pars=n,inc_warmup=FALSE) #+ geom_hline(yintercept = 0) #Traceplots
    
    print(p)
    a <- readline('Press Return to continue: ')
    if(a!='') break()
  }
}

#Compare ranges of parameters
n <- names(modList[[1]]) #Model parameters
n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('[sS]igma',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
n <- n[!grepl('_miss',n)] #Gets rid of imputed values

(p <- lapply(modList,function(i){
  if(is.null(i)) return(NA)
  parnames <- n[n %in% names(i)]
  
  l <- t(sapply(extract(i,parnames),function(x) setNames(quantile(x,c(0.5,0.025,0.975)),c('median','lwr','upr'))))  
  l <- l %>% data.frame() %>% rownames_to_column('param')
}) %>% bind_rows(.id='mod') %>% 
    mutate(param=factor(param,levels=n),mod=factor(mod,levels=names(modList))) %>% 
    ggplot(aes(x=param,y=median,col=mod))+
    geom_pointrange(aes(ymax=upr,ymin=lwr),position=position_dodge(width=0.5))+
    facet_wrap(~param,scales='free')+
    geom_hline(yintercept = 0,linetype='dashed')+
    labs(x='Parameter',y='Value',col='Model',title='Leafcutter visitation model')+
    scale_colour_manual(values=c('black','blue','purple','darkorange','red')))

ggsave('../Figures/lbeeVisModels.png',p,width = 8,height=6)
