#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN COMMODITY CANOLA FIELDS (2014+2015)

# Libraries and ggplot theme ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(beepr)
library(xtable)

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

# setwd('~/Projects/UofC/canola_yield_project') #Multivac path
setwd('~/Documents/canola_yield_project') #Galpern machine path

source('helperFunctions.R')

# Load in data ---------------------------------

load('./Models/datalist_commodity.Rdata')

# load("./Commodity field analysis/commodityfieldDataAll.RData")
# rm(AICc,brix2mgul,deltaAIC,DIC,plotFixedTerms,predFixef,se,varComp,zeroCheck,conversion,visitorsAll,visitors2015)
# fieldsAllComm <- fieldsAll; flowersAllComm <- flowersAll; plantsAllComm <- plantsAll;
# seedsAllComm <- seedsAll; surveyAllComm <- surveyAll;
# rm(fieldsAll,flowersAll,plantsAll,seedsAll,surveyAll)
# 
# #Set 'negative' missing pods (mistake in counting) to NA.
# plantsAllComm <- mutate(plantsAllComm,Missing=ifelse(Missing<0,NA,Missing))
# seedsAllComm <- mutate(seedsAllComm,Missing=ifelse(Missing<0,NA,Missing))
# setwd('~/Projects/UofC/canola_yield_project')
# 
# #List structure for Stan
# datalistField <- with(arrange(fieldsAllComm,as.factor(paste(Year,Field))),list( #Field-level measurements
#   Nfield=length(Year), #Number of fields
#   fieldName = paste(Year,Field), #Field name
#   numHives=NumHives, #Number of hives/field
#   is2015=Year==2015, #Is year 2015?
#   isGP=Area!='Lethbridge', #Is area Grand Prairie?
#   isIrrigated=Irrigated=='Irrigated', #Is field irrigated?
#   fieldSize=FieldSize_ha #Field size in hectares
# ))
# 
# datalistPlot <- with(arrange(surveyAllComm,factor(paste(Year,Field,Distance))),list( #Plot-level measurements
#   Nplot=length(Distance), #Number of plots
#   plotIndex=match(paste(Year,Field),datalistField$fieldName), #Index for field (which field?)
#   plotName=paste(Year,Field,Distance), #Name of plot location
#   dist=Distance, #Distance from edge
#   hbeeVis=Honeybee, #Visits by honeybees
#   flyVis=Fly, #Visits by flies
#   totalTime=TotalTime/10, #Total time (mins/10)
#   #Planting density
#   Nplot_densObs=sum(!is.na(PlDens)), #Number of plots with observed flower density
#   Nplot_densMiss=sum(is.na(PlDens)), #Plots with unobserved flower density
#   plDens_obs=log(PlDens[!is.na(PlDens)]), #(log) Planting density
#   obsPlDens_ind=which(!is.na(PlDens)), #Index for observed
#   missPlDens_ind=which(is.na(PlDens)), #Index for missing
#   flDens=sqrt(FlDens)-mean(sqrt(FlDens)), #(sqrt) Flower density - centered
#   flDensMean=mean(sqrt(FlDens))
# ))
# 
# datalistFlw <- flowersAllComm %>%
#   mutate(flowerIndex=as.numeric(factor(paste(Year,Field,Distance)))) %>%
#   arrange(flowerIndex) %>%
#   filter(!is.na(Pollen)) %>%
#   with(list(Nflw=sum(!is.na(Pollen)), #Number of pollen samples
#   flowerIndex=flowerIndex, #Index for flower (which plot?)
#   pollenCount=Pollen
# ))
# 
# datalistPlant <-
#   plantsAllComm %>% ungroup() %>%
#   filter(!is.na(Distance)) %>%
#   filter(SeedMass!=0) %>%
#   filter(!is.na(Pods),!is.na(Missing),!is.na(AvPodCount),!is.na(AvPodMass)) %>%
#   # filter(!is.na(Pods),!is.na(Missing)) %>%
#   mutate(plantIndex=match(paste(Year,Field,Distance),datalistPlot$plotName)) %>% #Index for plant (which plot?)
#   mutate(plantIndex_char=paste(Year,Field,Distance)) %>% #Index for plant (which plot?)
#   # mutate(plantIndex=match(paste(Year,Field),datalistField$fieldName)) %>% #Index for plant (which field?)
#   arrange(plantIndex) %>%
#   with(list(Nplant=length(VegMass), #Number of plant samples (some missing)
#             VegMass=VegMass,
#   # Nplant_obs=sum(!is.na(VegMass)), #Observed plants
#   # Nplant_miss=sum(is.na(VegMass)), #Missing plants
#   podCount=Pods, #Successful pods
#   flwCount=Pods+Missing, #Pods + Missing (total flw production)
#   plantIndex=plantIndex, #Index for plant (which plot?)
#   plantIndex_char=plantIndex_char,
#   plantSize=log(VegMass[!is.na(VegMass)]), #log weight of veg mass (g)
#   #Averaged seeds per pod and weight per seed
#   seedCount=AvPodCount,
#   seedMass=(AvPodMass/AvPodCount)*1000, #Weight per seed (mg)
#   yield=SeedMass[!is.na(SeedMass)] #observed weight of all seeds (g)
#   # obsPl_ind=which(!is.na(VegMass)),missPl_ind=which(is.na(VegMass)), #Indices for missing plant size
#   # obsYield_ind=which(!is.na(SeedMass)), missYield_ind=which(is.na(SeedMass)) #Indices for missing yield
# ))
# 
# # #Problem: plots exist at the pod level which do not exist at plant level
# # a <- plantsAllComm %>% ungroup() %>% filter(!is.na(Pods),!is.na(Missing),!is.na(AvPodCount),!is.na(AvPodMass)) %>%
# #   filter(SeedMass!=0) %>%  transmute(index=factor(paste(Year,Field,Distance,Plant))) %>%
# #   distinct() #Index for plant (which plot?)
# # b <- seedsAllComm %>% ungroup() %>% filter(!is.na(Plant)&!is.na(PodCount)&PodCount>0&!is.na(PodMass)) %>%
# #   filter(!is.na(Pods),!is.na(Missing)) %>% #Remove plants from plant level
# #   transmute(index=factor(paste(Year,Field,Distance,Plant))) %>% distinct()
# # keep <- which(a$index %in% b$index) #Plants to keep from seeds plants dataset
# #
# # datalistPod <- seedsAllComm %>% ungroup() %>%
# #   mutate(podIndex=as.numeric(factor(paste(Year,Field,Distance,Plant)))) %>%
# #   arrange(podIndex) %>%
# #   filter(podIndex %in% keep) %>%
# #   filter(!is.na(Plant)&!is.na(PodCount)&PodCount>0&!is.na(PodMass)) %>%
# #   filter(!is.na(Pods),!is.na(Missing)) %>% #Remove plants from plant level
# #   with(list(Npod=length(Distance), #Number of seeds measured
# #   seedCount=PodCount, #Number of seeds per pod
# #   seedMass=(PodMass/PodCount)*1000, #Weight per seed (mg)
# #   podIndex=podIndex) #Index for pod (which plant?)
# # )
# #Adds average flower count per plant, and average plant weight (plot-level), used in claim 6, 12, 15
# datalistPlot <- c(datalistPlot,with(datalistPlant,
#   list( #Both terms are missing from plots 248-271, so use the same N and indices
#     Nplot_flwCountObs = length(unique(plantIndex)),
#     Nplot_flwCountMiss = datalistPlot$Nplot-length(unique(plantIndex)),
#     flwCountPlot_obs = log(unname(tapply(flwCount,plantIndex,mean))), #Log-transformed average flower count
#     plSizePlot_obs = unname(tapply(plantSize,plantIndex,mean)), #Plant size
#     obsFlwCount_ind = unique(plantIndex),
#     missFlwCount_ind = (1:datalistPlot$Nplot)[!(1:datalistPlot$Nplot %in% unique(plantIndex))]
#   )
# ))
# 
# # datalistPod$seedMass[datalistPod$seedMass>8] <- with(datalistPod,seedMass[seedMass>8]/10) #Fixes weird outliers
# datalist <- c(datalistField,datalistPlot,datalistFlw,datalistPlant)
# rm(datalistField,datalistPlot,datalistFlw,datalistPlant,a,b,keep) #Cleanup
# str(datalist)
# 
# #Check for NAs
# if(any(sapply(datalist,function(x) sum(is.na(x)))>0)){
#   beep(1); print("NAs found in datalist")
# }
# 
# save(datalist,file = './Models/datalist_commodity.Rdata')

# Run main models ----------------------------

library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE) #Avoids recompilation
rstan_options(javascript = FALSE)
options(mc.cores = 6)

#Models split into separate stan files (faster)
modFiles <- dir(pattern = 'commodity.*\\.stan')
modFiles <- modFiles[!modFiles=='commodity_07flwSurv2.stan'] #Extra pod count model
modList <- vector(mode = 'list',length = length(modFiles))
names(modList) <- gsub('(commodity.*[0-9]{2}|\\.stan)','',modFiles)

modList[1] <- stan(file=modFiles[1],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #OK - Plant density, Plant size, Flower Density
modList[2] <- stan(file=modFiles[2],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.85),init=0) #OK - Visitation
modList[3] <- stan(file=modFiles[3],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #OK - Pollen
modList[4] <- stan(file=modFiles[4],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #OK - Flower count
modList[5] <- stan(file=modFiles[5],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0); beep(1) #OK - Pod Count (flw survival)
# modList[5] <- stan(file='commodity_07flwSurv2.stan',data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0); beep(1) #Pod count-only version - OK, but not as biologically interesting I think
modList[6] <- stan(file=modFiles[6],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #OK - Seed Count
modList[7] <- stan(file=modFiles[7],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #OK - Seed weight
modList[8] <- stan(file=modFiles[8],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #OK - Total yield
# allMod <- stan(file='visitation_pollen_model.stan',iter=100,chains=1,init=0,data=datalist) #Full model; avoid running

#Get model summaries into a list of tables
load('modSummaries_commodity.Rdata')
# modSummaries_commodity <- vector(mode = 'list',length = length(modList)) #Create new empty list
# names(modSummaries_commodity) <- names(modList)
#Update model summaries if needed
temp <- lapply(modList,parTable) #Get parameter summaries
for(i in 1:length(temp)){
  if(class(temp[[i]])=='data.frame'||length(temp[[i]])>1){ #If temp is not empty (model not run)
    modSummaries_commodity[[i]]$summary <- temp[[i]]$summary #Overwrite
    modSummaries_commodity[[i]]$covMat <- temp[[i]]$covMat
  }
}
parNames <- lapply(modSummaries_commodity,function(x) x$summary$param) #Clean up extra parameter names
for(i in 2:length(modSummaries_commodity)){
  if(!is.null(modSummaries_commodity[[i]])){
    chooseThese <- !parNames[[i]] %in% unlist(parNames[1:(i-1)])
    modSummaries_commodity[[i]]$summary <- modSummaries_commodity[[i]]$summary[chooseThese,]
    modSummaries_commodity[[i]]$covMat <- modSummaries_commodity[[i]]$covMat[chooseThese,chooseThese]
  }
}
save(modSummaries_commodity,file = 'modSummaries_commodity.Rdata')

#Traceplots
for(i in 1:length(modList)){
  if(!is.null(modList[[i]])){
    n <- names(modList[[i]]) #Model parameters
    n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('[sS]igma',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
    p <- traceplot(modList[[i]],pars=n,inc_warmup=FALSE)#+geom_hline(yintercept = 0) #Traceplots
    print(p)
    a <- readline('Press Return to continue: ')
    if(a!='') break()
  }
}



#Posterior predictive checks - OK 
PPplots(modList[[1]],datalist$plDens_obs,c('predPlDens','plDens_resid','predPlDens_resid'),
        index = datalist$obsPlDens_ind,main='Plant Density')
PPplots(modList[[1]],datalist$flDens,c('predFlDens','flDens_resid','predFlDens_resid'),'Flower density')
PPplots(modList[[1]],datalist$plantSize,c('predPlSize','plSize_resid','predPlSize_resid'),'Plant size')
PPplots(modList[[2]],datalist$hbeeVis,c('predHbeeVis','hbeeVis_resid','predHbeeVis_resid'),'Honeybee visits',jitterX=0.1) #Not great, tends to overpredict high visitation. 
PPplots(modList[[3]],datalist$pollenCount,c('predPollenCount','pollen_resid','predPollen_resid'),'Pollen')
PPplots(modList[[4]],datalist$flwCount,c('predFlwCount','flwCount_resid','predFlwCount_resid'),'Flowers per plant')
PPplots(modList[[5]],datalist$podCount,c('predPodCount','podCount_resid','predPodCount_resid'),'Pods per plant')
PPplots(modList[[6]],datalist$seedCount,c('predSeedCount','seedCount_resid','predSeedCount_resid'),'Seeds per pod')
PPplots(modList[[7]],datalist$seedMass,c('predSeedWeight','seedWeight_resid','predSeedWeight_resid'),'Seed size')
PPplots(modList[[8]],log(datalist$yield),c('predYield','yield_resid','predYield_resid'),'Seed mass per plant')

# Examine random intercept distributions
compareRE(modList[[1]],'intPlDens_field') 
compareRE(modList[[1]],'intPlSize_field')
# compareRE(modList[[1]],'intPlSize_plot') #Not great. Most overlap zero, and N_eff for sigma is low
compareRE(modList[[1]],'intFlDens_field') 
compareRE(modList[[2]],'intVisit_field') #skew-normal, but looks OK
compareRE(modList[[3]],'intPollen_field')
compareRE(modList[[4]],'intFlwCount_field')
compareRE(modList[[4]],'intFlwCount_plot')
compareRE(modList[[5]],'intFlwSurv_field')
compareRE(modList[[5]],'intFlwSurv_plot')
compareRE(modList[[6]],'intSeedCount_field')
# compareRE(modList[[6]],'intSeedCount_plot') #
compareRE(modList[[7]],'intSeedWeight_field')
# compareRE(modList[[7]],'intSeedWeight_plot')
compareRE(modList[[8]],'ranEffYield_field',1) #Intercepts
compareRE(modList[[8]],'ranEffYield_field',2) #Slopes
compareRE(modList[[8]],'ranEffYield_plot',1,0.3) #Intercepts
compareRE(modList[[8]],'ranEffYield_plot',2,0.3) #Slopes - right skew


# Run claims models ---------------------

library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE) #Avoids recompilation
rstan_options(javascript = FALSE)
options(mc.cores = 6)

# Dagitty claims list for commodity fields 

#Specify model
library(ggdag)
commDAG <- dagify(plDens ~ hbeeDist,
                  plSize ~ plDens + hbeeDist,
                  flDens ~ plSize + hbeeDist + plDens,
                  hbeeVis ~ hbeeDist + numHives + flDens,
                  pollen ~ hbeeVis + hbeeDist,
                  flwCount ~ plSize + flwSurv,
                  flwSurv ~ hbeeVis + plSize + pollen,
                  seedCount ~ hbeeVis + pollen + plSize + flwSurv + flwCount,
                  seedWeight ~ hbeeVis + pollen + plSize + seedCount + plDens + hbeeDist + flwSurv + flwCount
)

print(unlist(shipley.test(commDAG,TRUE))) #d-separation claims list

#Get model names, and make list
modFiles <- dir(path='./Commodity model claims 3',pattern = '*\\.stan',full.names = TRUE)
modFiles <- modFiles[!grepl('template',modFiles)]
modFiles <- modFiles[sapply(read.csv('./Commodity model claims 3/claimsList_updated2.csv')$Filename,function(x) grep(x,modFiles))]

# #Create storage list
# modList <- vector(mode = 'list',length = length(modFiles))
# names(modList) <- paste0('claim',formatC(1:length(modList),width=2,flag='0'))

runThese <- sapply(paste0('commodity_claim',formatC(c(2,3,6,9,10),width = 2,flag = 0),'.stan'),function(x) grep(x,modFiles))

for(i in runThese){
  # i <- 8
  overwrite <- TRUE
  if(file.exists(modFiles[i])){
    print(paste0('Starting model ',modFiles[i]))
    mod <- stan(file=modFiles[i],data=datalist,iter=3000,chains=4,control=list(adapt_delta=0.8),init=0)
    
    temp <- parTable(mod)$summary #Get parameter summaries
    #Save information to csv file
    modList <- read.csv('./Commodity model claims 3/claimsList_updated2.csv',sep=',',strip.white = TRUE)
    modList[i,match(names(temp),names(modList))] <- temp[grepl('claim',temp$param),]
    write.csv(modList,'./Commodity model claims 3/claimsList_updated2.csv',row.names = FALSE)

    #Check model
    parNam <- names(mod) #Model parameters
    
    #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
    useThese <- (!(grepl('(\\[[0-9]+,*[0-9]*\\]$|^lp)',parNam)) |
                   grepl('[sS]igma',parNam) 
                   )
    parNam <- parNam[useThese] 
    
    #Removes pollen terms
    parNam <- parNam[!parNam %in% c('intPollen','slopeVisitPol','slopeHbeeDistPollen','sigmaPolField','pollenPhi')] 
    
    print(traceplot(mod,pars=parNam,inc_warmup=FALSE)+labs(title=gsub('.*/','',modFiles[i])))
    print(mod,pars=parNam)
    # fastPairs(mod,pars=parNam)

    print(paste0('Model ',modFiles[i],' completed'))
    
  } else print(paste0('Model ',i,' not found'))
  rm(mod); gc() #Cleanup
}

#Calculate C-statistic

#Original model
read.csv('./Commodity model claims 3/claimsList_updated.csv',sep=',',strip.white = TRUE) %>% 
  mutate(pval=pnorm(-abs(Z))*2) %>% 
  shipley.dSep(.,pval,param)

#Updated
read.csv('./Commodity model claims 3/claimsList_updated2.csv',sep=',',strip.white = TRUE) %>% 
  mutate(pval=pnorm(-abs(Z))*2) %>% 
  shipley.dSep(.,pval,param)
