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

setwd('~/Projects/UofC/canola_yield_project') #Multivac path
# setwd('~/Documents/canola_yield_project') #Galpern machine path

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
#   isIrrigated=Irrigated=='Irrigated' #Is field irrigated?
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
# 
#   flDens=sqrt(FlDens) #(sqrt) Flower density
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
modList[2] <- stan(file=modFiles[2],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #OK - Visitation
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
    modSummaries_commodity[[i]] <- temp[[i]] #Overwrite
  }
}
parNames <- lapply(modSummaries_commodity,function(x) x$param) #Clean up extra parameter names
for(i in 2:length(modSummaries_commodity)){
  if(!is.null(modSummaries_commodity[[i]])){
    modSummaries_commodity[[i]] <- modSummaries_commodity[[i]][!parNames[[i]] %in% unlist(parNames[1:(i-1)]),]
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
PPplots(modList[[2]],datalist$hbeeVis,c('predHbeeVis','hbeeVis_resid','predHbeeVis_resid'),'Honeybee visits',jitterX=0.1) #Not great, tends to overpredict high visitation. Probably 
PPplots(modList[[3]],datalist$pollenCount,c('predPollenCount','pollen_resid','predPollen_resid'),'Pollen')
PPplots(modList[[4]],datalist$flwCount,c('predFlwCount','flwCount_resid','predFlwCount_resid'),'Flowers per plant')
PPplots(modList[[5]],datalist$podCount,c('predPodCount','podCount_resid','predPodCount_resid'),'Pods per plant')
PPplots(modList[[6]],datalist$seedCount,c('predSeedCount','seedCount_resid','predSeedCount_resid'),'Seeds per pod')
PPplots(modList[[7]],datalist$seedMass,c('predSeedWeight','seedWeight_resid','predSeedWeight_resid'),'Seed size')
PPplots(modList[[8]],log(datalist$yield),c('predYield','yield_resid','predYield_resid'),'Seed mass per plant')

# Examine random intercept distributions
compareRE(modList[[1]],'intPlDens_field') 
compareRE(modList[[1]],'intPlSize_field')
compareRE(modList[[1]],'intFlDens_field') 
compareRE(modList[[2]],'intVisit_field') #skew-normal, but looks OK
compareRE(modList[[3]],'intPollen_field')
compareRE(modList[[4]],'intFlwCount_field')
compareRE(modList[[4]],'intFlwCount_plot')
compareRE(modList[[5]],'intFlwSurv_field')
compareRE(modList[[5]],'intFlwSurv_plot')
compareRE(modList[[6]],'intSeedCount_field')
compareRE(modList[[6]],'intSeedCount_plot') #
compareRE(modList[[7]],'intSeedWeight_field')
compareRE(modList[[7]],'intSeedWeight_plot')
compareRE(modList[[8]],'ranEffYield_field',1) #Intercepts
compareRE(modList[[8]],'ranEffYield_field',2) #Slopes
compareRE(modList[[8]],'ranEffYield_plot',1,0.3) #Intercepts
compareRE(modList[[8]],'ranEffYield_plot',2,0.3) #Slopes - right skew


# Path diagram ------------------------------------------------------------

library(ggdag)
nodeCoords <- data.frame(name=c('numHives','hbeeDist','hbeeVis','pollen',
                                'plSize','plDens','flDens',
                                'flwCount','flwSurv','seedCount','seedWeight'),
                         labs=c('Number\nof Hives','Distance','Honey bee\nVisits','Pollen\nCount',
                                'Plant\nSize','Plant\nDensity','Flower\nDensity',
                                'Flowers\nper Plant','Pods\nper Plant','Seeds\nper Pod','Seed\nSize'),
                         x=c(0,1,0.5,0.5,
                             2.5,0,1,
                             3.5,4,2.5,3.5),
                         y=c(4,4,3,2,
                             0,1,1,
                             0,1,3,2))


#Specify model
commDAG <- dagify(plDens ~ hbeeDist,
                  plSize ~ plDens + hbeeDist,
                  flDens ~ plSize + hbeeDist + plDens,
                  hbeeVis ~ hbeeDist + numHives + flDens,
                  pollen ~ hbeeVis + hbeeDist,
                  flwCount ~ plSize + flwSurv,
                  flwSurv ~ hbeeVis + plSize + pollen,
                  seedCount ~ hbeeVis + pollen + plSize + flwSurv + flwCount,
                  seedWeight ~ hbeeVis + pollen + plSize + seedCount + plDens + hbeeDist + flwSurv + flwCount + flDens,
                  coords= list(x = setNames(nodeCoords$x,nodeCoords$name),
                               y = setNames(nodeCoords$y,nodeCoords$name)),
                  labels=setNames(nodeCoords$labs,nodeCoords$name)
)

# plot(commDAG)

#Get path coefficients (Z-scores) and match to edges

load('modSummaries_commodity.Rdata')

pathCoefs <- modSummaries_commodity %>% #Get path coefficients
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

commDAG <- commDAG %>% #Create tidy dagitty set
  tidy_dagitty() 

# ggplot(commDAG,aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_dag_edges() +
#   geom_dag_text(col='black') +
#   theme_dag_blank()

commDAG$data <- commDAG$data %>% #Match coefs to dagitty set
  rename(xstart=x,ystart=y) %>% 
  full_join(x=pathCoefs,y=.,by=c('to','name')) %>% 
  mutate(isNeg=Z<0,isSig=pval<0.05) %>% 
  mutate(L=sqrt(abs(Z)),1) %>% 
  mutate(edgeLab=ifelse(isSig,as.character(sign(Z)*round(L,1)),'')) %>% 
  mutate(C=ifelse(isNeg,'red','black')) %>% 
  mutate(A=ifelse(isSig,1,0.1))

ggplot(commDAG$data,aes(x = xstart, y = ystart, xend = xend, yend = yend))+
  # geom_dag_point(size=20,colour='black',shape='square')+
  # geom_dag_text(aes(label=label),col='white') +
  annotate('text',x=0.5,y=4.5+0.1,label='Plot Level',size=5)+
  annotate('rect',xmin=-0.5,ymin=0.5,xmax=1.5,ymax=4.5,fill=NA,col='black',
           linetype='dashed',linejoin='round',size=1)+
  
  annotate('text',x=3.5,y=3.5+0.1,label='Plant Level',size=5)+
  annotate('rect',xmin=2,ymin=-0.5,xmax=4.5,ymax=3.5,fill=NA,col='black',
           linetype='dashed',linejoin='round',size=1)+
  
  geom_dag_edges(aes(edge_width=L,edge_colour=C,edge_alpha=A),
                 arrow_directed=arrow(angle=20,type='open')) +
  geom_dag_edges(aes(label=edgeLab),label_pos=0.45,edge_alpha=0,label_size=6,fontface='bold',label_colour ='black') +
  geom_dag_edges(aes(label=edgeLab),label_pos=0.45,edge_alpha=0,label_size=6,fontface='plain',label_colour='white') +
  geom_dag_label_repel(aes(label=label),col='black',force=0) +
  theme_dag_blank()


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

# runThese <- 8:24 #1:length(modFiles)

# for(i in runThese){
{
  i <- 8
  overwrite <- TRUE
  if(file.exists(modFiles[i])){
    print(paste0('Starting model ',modFiles[i]))
    mod <- stan(file=modFiles[i],data=datalist,iter=3000,chains=4,control=list(adapt_delta=0.8),init=0)
    
    temp <- parTable(mod) #Get parameter summaries
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

modList <- read.csv('./Commodity model claims 3/claimsList_updated2.csv',sep=',',strip.white = TRUE)

modList %>% 
  mutate(pval=pnorm(-abs(Z))*2) %>% 
  shipley.dSep(.,pval,param)

# debugonce(shipley.dSep)

# Partial effects plots for commodity fields -----------------------------

#Plant density

list('intPlDens'=1,
           'slopeDistPlDens'=with(datalist,seq(min(log(dist)),max(log(dist)),length=10))) %>% 
  getPreds(modSummaries_commodity[[1]],parList = .,otherPars = c('sigmaPlDens','sigmaPlDens_field'))

list('intPlDens'=1,
     'slopeDistPlDens'=with(datalist,seq(min(log(dist)),max(log(dist)),length=10))) %>% 
  getPreds(modList[[1]],parList = .,otherPars = c('sigmaPlDens','sigmaPlDens_field'))


load('modPodcount.Rdata') #All extracted coefficients in list form (mod3)
mod3 <- extract(modPodcount) #Get coefficients from stan model
rm(modPodcount); gc()

#Plot-level data (271 rows)
plotDat <- with(datalist,data.frame( 
  field=plotIndex,year=is2015[plotIndex],irrigation=isIrrigated[plotIndex],
  dist=log(dist)-mean(log(dist)), #Centered on 3.412705
  GP=isGP[plotIndex],hives=log(numHives[plotIndex]+1), #Log-number of hives
  plDens=apply(mod3$plDens,2,mean),flDens=flDens,
  hbeeResid=log((hbeeVis+0.2)/totalTime)-apply(mod3$visitHbeeMu,2,median) #Log-residuals - #maybe not the best
))

#Model matrix for hbee visits
MM_hbee <- with(plotDat,data.frame(int=1,year,GP,yearGP=year*GP,dist,hives,flDens,irrigation)) %>% as.matrix()
#Coefficient matrix for hbee visits
coef_hbee <- with(mod3,data.frame(intVisit,slopeYearVis,slopeGpVis,slopeYearGpVis,
                                  slopeDistVis,slopeHiveVis,slopeFlDens,slopeIrrigVis))

#Partial effect of year/area, distance
MM_temp <- MM_hbee %>% as.data.frame() %>% mutate_at(vars(-year,-GP,-yearGP,-dist),mean) %>% 
   as.matrix()

#Partial effect of distance only
MM_temp2 <- MM_hbee %>% as.data.frame() %>% mutate_at(vars(-dist),mean) %>% 
  as.matrix()
MM_temp2 <- data.frame(MM_temp2,t(apply(MM_temp2 %*% t(coef_hbee),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% unique() %>% mutate_at(vars(pred:lwr),exp) %>% mutate(dist=exp(dist+3.41)) %>% 
  select(dist,pred) %>% unite(dat,dist,pred) %>% 
  with(.,expand.grid(dat=dat,year=c('2014','2015'),GP=c('Lethbridge','Grande Prairie'))) %>% 
  separate(dat,c('dist','pred'),sep='_',convert=T)

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_hbee),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred) %>% mutate_at(vars(pred:resid),exp) %>% 
  mutate(year=factor(year,labels=c('2014','2015')),GP=factor(GP,labels=c('Lethbridge','Grande Prairie'))) %>% 
  mutate(yearGP=factor(paste(GP,year)),dist=exp(dist+3.41)) %>%
  ggplot(aes(x=dist))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid))+
  geom_line(aes(y=pred),size=1)+facet_grid(year~GP)+ylim(0,5)+
  geom_line(data=MM_temp2,aes(x=dist,y=pred),col='red')+
  labs(x='Distance from edge(m)',y='Visits/10 mins')
ggsave('../Figures/Commodity/slopeYearGPVis.png',p1,width=8,height=6)

#Partial effect of year/area, using 5m plot only, using 20 hives only
MM_temp <- MM_hbee %>% as.data.frame() %>% mutate_at(vars(-year,-GP,-yearGP),mean) %>% 
  mutate(dist=(log(5)-3.41),hives=log(20+1)) %>% 
  as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_hbee),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred) %>% mutate_at(vars(pred:resid),exp) %>% 
  mutate(year=factor(year,labels=c('2014','2015')),GP=factor(GP,labels=c('Lethbridge','Grande\nPrairie'))) %>% 
  select(year,GP,pred,upr,lwr) %>%   
  distinct() %>% 
  ggplot(aes(x=GP,y=pred))+
  geom_pointrange(aes(ymax=upr,ymin=lwr,col=factor(year)),position=position_dodge(width=0.25))+
  labs(x=NULL,y='Honey bee visits/10 mins',col='Year')+
  scale_colour_manual(values=c('black','grey50'))


#Partial effect of stocking rates and distance
MM_temp <- MM_hbee %>% as.data.frame() %>% mutate_at(vars(-hives,-dist),mean) %>% 
  as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_hbee),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  mutate(dist=exp(dist+3.41),hives=round(exp(hives)-1)) %>% 
  filter(hives==0|hives==20|hives==40) %>% 
  # mutate(hives=factor(hives,levels=))
  ggplot(aes(x=dist,y=pred))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=factor(hives)),alpha=0.3,show.legend=F)+
  geom_line(aes(col=factor(hives)),size=1)+
  geom_jitter(aes(y=resid,col=factor(hives)),width=5,height=0)+
  labs(col='Number\nof Hives',y='Honey bee visits/10 mins',x='Distance (m)')+
  scale_colour_manual(values=c('black','blue','red'))+
  scale_fill_manual(values=c('black','blue','red'))
ggsave('../Figures/Commodity/slopeNumHivesDistVis.png',p1,width=8,height=6)


#Pollen data (1294 rows)
pollenDat <- with(datalist,data.frame(
  field=plotIndex[plotIndex[flowerIndex]],
  plot=plotIndex[flowerIndex],
  pol=pollenCount,
  visit=(hbeeVis/totalTime)[flowerIndex], #Visits/10 mins
  logvisit=log((hbeeVis/totalTime)+0.5)[flowerIndex],
  dist=dist[flowerIndex], #Distance
  logdist=(log(dist)-mean(log(dist)))[flowerIndex], 
  polResid=log(pollenCount+1)-apply(mod3$pollenMu,2,median)
))

#Model matrix for pollen
MM_pol <- with(pollenDat,data.frame(int=1,logvisit,logdist)) %>% as.matrix()

#Coefficient matrix for plant size
coef_pol <- with(mod3,data.frame(intPollen,slopeVisitPol,slopeHbeeDistPollen))

#Partial effect of hbee visitation
MM_temp <- MM_pol %>% as.data.frame() %>% mutate_at(vars(-logvisit),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_pol),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% mutate(resid=pollenDat$polResid+pred) %>%
  mutate_at(vars(upr,lwr,pred,resid), exp) %>% 
  mutate(visit=pollenDat$visit) %>% 
  ggplot(aes(x=visit))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid),alpha=0.7)+
  geom_line(aes(y=pred),size=1)+
  labs(x='Visits/10 mins',y='Pollen grains/stigma')+ylim(0,1000)
ggsave('../Figures/Commodity/slopeHbeeVisPol.png',p1,width=8,height=6)


#Plant-level data (792 rows)
plantDat <- with(datalist,data.frame(
  field=plotIndex[plantIndex],
  plot=plantIndex,
  visit=(hbeeVis/totalTime)[plotIndex[plantIndex]],
  logvisit=log((hbeeVis/totalTime)+0.5)[plotIndex[plantIndex]],
  pollen=apply(mod3$pollenPlot,2,median)[plotIndex[plantIndex]],
  plSize=plantSize, #Plant size
  plDens=apply(mod3$plDens,2,mean)[plantIndex], #Centered using 3.774427
  irrigation=isIrrigated[plotIndex[plantIndex]],
  year2015=is2015[plotIndex[plantIndex]],
  irrig2015=isIrrigated[plotIndex[plantIndex]]*is2015[plotIndex[plantIndex]],
  plSizeIrrig=plantSize*isIrrigated[plotIndex[plantIndex]],
  logdist=(log(dist)-mean(log(dist)))[plotIndex[plantIndex]], #Transformed distance
  dist=dist[plotIndex[plantIndex]], #Actual distance
  GP=isGP[plotIndex[plantIndex]],
  stocking=numHives[plotIndex[plantIndex]],
  distStocking=(log(dist)-mean(log(dist)))[plotIndex[plantIndex]]*numHives[plotIndex[plantIndex]],
  dist2015=(log(dist)-mean(log(dist)))[plotIndex[plantIndex]]*is2015[plotIndex[plantIndex]],
  plDensStocking=apply(mod3$plDens,2,mean)[plantIndex]*numHives[plotIndex[plantIndex]],
  irrigStocking=isIrrigated[plotIndex[plantIndex]]*numHives[plotIndex[plantIndex]],
  plSizeResid=apply(mod3$plSize_resid,2,median),
  flwCount=flwCount,
  flwCountResid=log(flwCount)-apply(mod3$flwCountMu,2,median), #Log-resid for flower count
  podCountResid=logit(podCount/flwCount)-apply(mod3$flwSurv,2,median) #Logit-resid for flower survival
  )) 

#Model matrix for plant size
MM_plSize <- with(plantDat,data.frame(int=1,plDens,logdist,GP,year2015,
                    irrigation)) %>% as.matrix()
#Coefficient matrix for plant size
coef_plSize <- with(mod3,data.frame(intPlSize,slopePlDensPlSize,slopeDistPlSize,
            slopeGpPlSize,slope2015PlSize,slopeIrrigPlSize))

#Partial effect of plant density
MM_temp <- MM_plSize %>% as.data.frame() %>% mutate_at(vars(-plDens),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_plSize),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plantDat$plSizeResid+pred) %>%
  mutate(plDens=exp(plDens+mean(log(surveyAllComm$PlDens),na.rm=T))) %>% #Untransform plant density
  mutate_at(vars(upr,lwr,pred,resid),function(x,y) exp(x+mean(log(plantsAllComm$VegMass),na.rm=T))) %>% 
  ggplot(aes(x=plDens))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_point(aes(y=resid),alpha=0.5)+
  geom_line(aes(y=pred),size=1)+
  labs(x=expression(paste('Plants per ',m^2)),y='Plant size (g)',title='Commodity')+ylim(0,50)
ggsave('../Figures/Commodity/slopePlDensPlSize.png',p1,width=6,height=6)

#Model matrix for flower production (per plant)
MM_flwCount <- with(plantDat,data.frame(int=1,plSize)) %>% as.matrix()

#Coefficient matrix for flower production
coef_flwCount <- with(mod3,data.frame(intFlwCount,slopePlSizeFlwCount))

#Effect of plant size on flower production
MM_temp <- MM_flwCount 

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
  mutate(resid=plantDat$flwCountResid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=exp(plantDat$plSize+2.63)))+ #Transform to real scale
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='red')+
  geom_jitter(aes(y=resid),alpha=0.5)+
  geom_line(aes(y=pred),size=1,col='red')+
  labs(x='Plant size (g)',y='Flowers per plant',title='Commodity')
ggsave('../Figures/Commodity/slopePlSizeFlwCount.png',p1,width=6,height=6)

#Model matrix for flower survival (pod count)
MM_flwSurv <- with(plantDat,data.frame(int=1,logvisit,pollen,plSize,
                                       plDens,irrigation,year2015)) %>% as.matrix()
#Coefficient matrix for flower survival
coef_flwSurv <- with(mod3,data.frame(intFlwSurv,slopeVisitSurv,slopePolSurv,slopePlSizeSurv,
                                     slopePlDensSurv,slopeIrrigSurv,slope2015Surv))

#Partial effects plot of visitation on flower survival
MM_temp <- MM_flwSurv %>% as.data.frame() %>% mutate_at(vars(-logvisit),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$podCountResid+pred) %>% 
  mutate_at(vars(pred:resid),invLogit) %>%
  ggplot(aes(x=plantDat$visit))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid),alpha=0.5)+ylim(0.5,0.9)+
  geom_line(aes(y=pred),size=1,col='red')+
  labs(x='Visits/10 mins',y='Pod survival')
  
#Partial effects of pollen deposition on flower survival
MM_temp <- MM_flwSurv %>% as.data.frame() %>% mutate_at(vars(-pollen),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$podCountResid+pred) %>% mutate_at(vars(pred:resid),invLogit) %>%
  ggplot(aes(x=exp(pollen+median(mod3$intPollen))))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid),alpha=0.5)+#ylim(0.5,0.9)+
  geom_line(aes(y=pred),size=1,col='red')+
  labs(x='Pollen per stigma',y='Pod survival')

#Partial effects plot for plant size
MM_temp <- MM_flwSurv %>% as.data.frame() %>% mutate_at(vars(-plSize),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$podCountResid+pred) %>% mutate_at(vars(pred:resid),invLogit) %>%
  ggplot(aes(x=exp(plSize+2.63)))+
  geom_ribbon(aes(ymax=upr*100,ymin=lwr*100),alpha=0.3)+
  geom_jitter(aes(y=resid*100),alpha=0.5)+
  ylim(45,90)+
  # xlim(-2,2)+
  geom_line(aes(y=pred*100),size=1)+
  labs(x='Plant size (g)',y='Fruit set (%)',title='Commodity')
ggsave('../Figures/Commodity/slopePlSizeFlwSurv.png',p1,width=6,height=6)

#Partial effects plot for plant density
MM_temp <- MM_flwSurv %>% as.data.frame() %>% mutate_at(vars(-plDens),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$podCountResid+pred) %>% mutate_at(vars(pred:resid),invLogit) %>%
  ggplot(aes(x=plDens))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid),alpha=0.5)+
  geom_line(aes(y=pred),size=1,col='red')+
  labs(x='Plant density',y='Pod survival')

#Partial effects plot of plant size on flower and pod number
MM_temp <- MM_flwSurv %>% as.data.frame() %>% mutate(plSize=cut(plSize,breaks=c(-5,-1,1,5))) %>%
  mutate(plSize=as.numeric(as.character(factor(plSize,labels=c('-1','0','1'))))) %>% 
  mutate_at(vars(-plSize),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(plSize=factor(plSize,labels=c('Low','Med','High'))) %>% 
  mutate(resid=plantDat$podCountResid+pred,flwCount=plantDat$flwCount) %>% 
  mutate_at(vars(pred:resid),invLogit) %>%
  arrange(plSize,flwCount) %>% 
  ggplot(aes(x=flwCount,y=resid*flwCount))+
  geom_point(alpha=0.5,col='black')+
  geom_line(aes(y=pred*flwCount,col=plSize),size=1)+
  geom_ribbon(aes(ymax=upr*flwCount,ymin=lwr*flwCount,fill=plSize),alpha=0.3)+
  geom_abline(intercept=0,slope=1)+
  labs(x='Flower number',y='Pod number',col='Plant Size',fill='Plant Size')+
  theme(legend.position=c(0.1,0.82),legend.background=element_rect(fill='white',colour='black',linetype='solid'))+
  # scale_x_log10()+scale_y_log10()+
  scale_colour_manual(values=c('blue','purple','red'))+
  scale_fill_manual(values=c('blue','purple','red'))



#Pod-level data (3872 rows)
podDat <- with(datalist,data.frame(
  field=plotIndex[plantIndex[podIndex]],
  plot=plantIndex[podIndex],
  plant=podIndex,
  visit=log((hbeeVis/totalTime)+0.5)[plotIndex[plantIndex[podIndex]]],
  pollen=apply(mod3$pollenPlot,2,median)[plotIndex[plantIndex[podIndex]]],
  seedCount=seedCount, plSize=plantSize[podIndex], #Plant size
  irrigation=isIrrigated[plotIndex[plantIndex[podIndex]]],
  year2015=is2015[plotIndex[plantIndex[podIndex]]],
  irrig2015=isIrrigated[plotIndex[plantIndex[podIndex]]]*is2015[plotIndex[plantIndex[podIndex]]],
  seedCountPlSize=seedCount*plantSize[podIndex],
  plSizeIrrig=plantSize[podIndex]*isIrrigated[plotIndex[plantIndex[podIndex]]],
  seedCount2015=seedCount*is2015[plotIndex[plantIndex[podIndex]]],
  plSize2015=plantSize[podIndex]*is2015[plotIndex[plantIndex[podIndex]]],
  plSizeIrrig2015=plantSize[podIndex]*isIrrigated[plotIndex[plantIndex[podIndex]]]*is2015[plotIndex[plantIndex[podIndex]]],
  seedWeight=seedMass,seedWeightResid=apply(mod3$seedMass_resid,2,median), #Median resid for seedMass
  # seedCountResid=apply(mod3$seedCount_resid,2,median) #Residuals on real scale
  seedCountResid=apply(log(matrix(rep(datalist$seedCount,nrow(mod3$seedCountMu)),ncol=nrow(mod3$seedCountMu)))-t(mod3$seedCountMu),1,mean)
  )) 

#Model matrix for seedCount
MM_seedCount <- with(podDat,data.frame(int=1,visit,pollen,plSize,year2015)) %>% as.matrix()
#Coefficient matrix for seedCount
coef_seedCount <- with(mod3,data.frame(intSeedCount,slopeVisitSeedCount,slopePolSeedCount,slopePlSizeCount,
            slope2015SeedCount))

#Partial effects plot for year
MM_temp <- MM_seedCount %>% as.data.frame() %>% mutate_at(vars(-year2015),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=exp(podDat$seedCountResid+pred),pred=exp(pred),lwr=exp(lwr),upr=exp(upr)) %>% 
  mutate(year2015=factor(year2015,labels=c('2014','2015'))) %>% 
  ggplot()+geom_jitter(aes(x=year2015,y=resid),width=0.25,alpha=0.3)+
  geom_pointrange(aes(x=year2015,y=pred,ymax=upr,ymin=lwr),col='red',size=1)+
  labs(x='Year',y='Seeds per pod')+ylim(0,40)

#Model matrix for seedMass 
MM_seedMass <- with(podDat,data.frame(int=1,visit,pollen,seedCount=seedCount,plSize,irrigation,year2015,irrig2015)) %>% as.matrix()

#Coefficent matrix for seedMass
coef_seedMass <- with(mod3,data.frame(intSeedWeight,slopeVisitSeedWeight,slopePolSeedWeight,slopeSeedCount,
                      slopePlSizeWeight,slopeIrrigSeedWeight,slope2015SeedWeight,slope2015IrrigSeedWeight)) %>% as.matrix()

#Partial effect of year
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-year2015),mean) %>% 
  as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=podDat$seedWeightResid+pred) %>% 
  mutate(year2015=factor(year2015,labels=c('2014','2015'))) %>% 
  ggplot()+geom_jitter(aes(x=year2015,y=resid),width=0.25,alpha=0.3)+
  geom_pointrange(aes(x=year2015,y=pred,ymax=upr,ymin=lwr),col='red',size=1)+
  labs(x='Year',y='Weight per seed')

  
#Partial effect of seed count
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-seedCount),mean) %>% 
  as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=podDat$seedWeightResid+pred) %>% 
  ggplot(aes(x=seedCount))+
  # geom_ribbon(data=temp2,aes(x=x,ymax=upr,ymin=lwr),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_line(aes(y=pred),col='red',size=1)+
  geom_jitter(aes(y=resid),width=0.2,alpha=0.3)+xlim(0,40)+
  labs(x='Seeds per pod',y='1000 seed weight (g)')
ggsave('../Figures/Commodity/slopePodCountPodWeight.png',p1,width=6,height=4)

#Partial effect of seed count - seed size relationship, with 3 levels of plant size (0.1,0.5,0.9 percentiles)
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-seedCount,-plSize),mean) %>%
  mutate(plSize=cut(plSize,breaks=c(-5,-0.9,0.9,5))) %>%
  mutate(plSize=as.numeric(as.character(factor(plSize,labels=c('-0.9','0','0.9'))))) %>%
  as.matrix()
 
p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=podDat$seedWeightResid+pred,plSize=factor(plSize,labels=c('Low','Med','High'))) %>% #arrange(seedCount) %>%
  ggplot()+
  geom_point(aes(seedCount,resid),col='black',alpha=0.2)+
  geom_ribbon(aes(x=seedCount,ymax=upr,ymin=lwr,fill=plSize),alpha=0.3)+
  geom_line(aes(seedCount,pred,col=plSize),size=1)+
  labs(x='Seeds per pod',y='1000 seed weight (g)')+
  scale_fill_manual(values=c('blue','purple','red'))+
  scale_colour_manual(values=c('blue','purple','red'))+
  labs(col='Plant Size',fill='Plant Size')+
  theme(legend.position=c(0.9,0.82),legend.background=element_rect(fill='white',colour='black',linetype='solid'),
        legend.title=element_text(size=12),legend.text=element_text(size=12))
ggsave('../Figures/Commodity/slopePodCountPodWeightPlSize.png',p1,width=6,height=4)
  

#Partial effect of irrigation
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-irrigation),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=pred+podDat$seedWeightResid) %>%
  mutate(irrigation=factor(irrigation,labels=c('Dryland','Irrigated'))) %>% 
  ggplot()+
  # geom_jitter(aes(x=irrigation,y=resid),width=0.25,alpha=0.3)+
  geom_violin(aes(x=irrigation,y=resid))+
  geom_pointrange(aes(x=irrigation,y=pred,ymax=upr,ymin=lwr),col='red',size=1)+
  labs(x='Irrigation',y='Weight per seed')

# Total effects for commodity fields --------------------------------------

#Idea 1: take all plots/plants, simulate effect of honeybees/no honeybees, see what yield would have been ("what if")
#Idea 2: simulate everything from the ground up at an "average" field, using honeybees/no honeybees

#Trying idea 2

# load('modPodcount.Rdata') #Seed count, size, and yield
# mod1 <- extract(modPodcount)
# load('modPodcount2.Rdata') #All other coefficients (plot/plant level) - 22 mins for 1000 iter
# mod2 <- extract(modPodcount)
# mod3 <- c(mod1,mod2)
# mod3 <- mod3[unique(c(names(mod1),names(mod2)))] #Removes doubles
# mod3$lp__ <- NULL #Removes LP
# rm(modPodcount,mod1,mod2); gc();
# save(mod3,file='modPodcount3.Rdata')
setwd('~/Projects/UofC/canola_yield_project/Models')
load('modPodcount.Rdata')
modPodcount <- extract(modPodcount)

3.412705 #Centering value for distance
3.774427 #Centering value for plant density
2.626568 #Centering value for plant size
21 #centering value for flower density: sqrt(FlDens)-21
1.0653 #centering value for logit-survival

#Checks argument lengths and returns correct sample size
getN <- function(args){
  N <- sapply(args,length) #Lengths of arguments
  #Stops function if arguments mismatch
  if(length(unique(N[N>1]))>1) stop('Arguments have mismatching lengths') 
  return(max(N))
}

#Generate plant densities from distances
genPlDens <- function(dist,is2015,isIrrig,isGP,int,slopeDist,slope2015,slopeIrrig,
                      slopeIrrig2015,slopeGP,stDev,transform=F){
  N <- getN(list(dist,is2015,isIrrig,isGP))
  #dist= distances, intercept,slope of distance, stDev, and back-transform
  transdist <- log(dist)-3.412705
  mu <- int+slopeDist*transdist+is2015*slope2015+isIrrig*slopeIrrig+is2015*isIrrig*slopeIrrig2015+isGP*slopeGP
  a <- rnorm(N,mu,stDev)
  if(transform) return(exp(a+3.774427)) else return(a)
}

#Generate plant sizes from plant density/distance
genPlSize <- function(plDens,dist,isGP,is2015,isIrrig,int,slopePlDens,slopeDist,slopeGP,
                      slope2015,slopeIrrig,stDev,transform=F){
  N <- getN(list(plDens,dist,isGP,is2015,isIrrig))
  #plDens=densities, dist= distances, intercept,slope of density, slope of distance, stDev, and back-transform
  transdist <- log(dist)-3.412705
  transdens <- log(plDens)-3.774427
  mu <- int+slopePlDens*transdens+slopeDist*transdist+slopeGP*isGP+slope2015*is2015+slopeIrrig*isIrrig
  a <- rnorm(plDens,mu,stDev)
  if(transform) return(exp(a+2.626568)) else return(a)
}

#Generate flower density for plot
genFlwDens <- function(avgPlSize,dist,int,slopePlSize,slopeDist,stDev,transform=F){
  N <- getN(list(avgPlSize,dist)) 
  #average log(plSize)-2.626568 plant size for plot, distance, intercept, slopes, stDev
  transdist <- log(dist)-3.412705
  a <- rnorm(N,int + avgPlSize*slopePlSize + transdist*slopeDist,stDev)
  if(transform) return((a+21)^2) else return(a)
}

#Generate hbee visits (per 10 mins)
genHbeeVis <- function(dist,numHives,flDens,is2015,isGP,isIrrig,int,slopeDist,
                       slopeHives,slopeFlDens,slope2015,slopeGP,slope2015GP,slopeIrrig,phi,
                       addVar=T){
  #addVar: should use rbinom to simulate variance, or return expected value?
  N <- getN(list(dist,numHives,flDens,is2015,isGP,isIrrig))
  transdist <- log(dist)-3.412705
  transfldens <- sqrt(flDens)-21
  transNumHives <- log(numHives+1)
  lambda <- int+slopeDist*transdist+slopeHives*transNumHives+slopeFlDens*transfldens+slope2015*is2015+
    slopeGP*isGP+slope2015GP*is2015*isGP+slopeIrrig*isIrrig
  if(addVar){
    a <- rnbinom(length(transdist),mu=exp(lambda),size=phi)
    return(a)
  } else return(exp(lambda))
}

#Generate avg pollen counts - use stDev from plot level
genAvgPol <- function(hbee,dist,int,slopeHbee,slopeDist,stDev,center=T,transform=F){
  N <- getN(list(hbee,dist))
  #hbee counts, distance, intercept, slopes, center data, transform
  #Note: transform doesn't make sense unless center=F (Jensen's inequality)
  transhbee <- log(hbee+0.5)
  transdist <- log(dist)-3.412705
  a <- rnorm(N,ifelse(center,0,int)+slopeHbee*transhbee+slopeDist*transdist,stDev)
  if(transform) return(exp(a)) else return(a)
}

#Generate flower survival rate
genFlwSurv <- function(hbee,pol,plSize,plDens,isIrrig,is2015,
                       int,slopeHbee,slopePol,slopePlSize,slopePlDens,slopeIrrig,slope2015,
                       intPhi,slopePlSizePhi,addVar=T){ 
  N <- getN(list(hbee,pol,plSize,plDens,isIrrig,is2015))
  transhbee <- log(hbee+0.5)
  transsize <- log(plSize)-2.626568
  transdens <- log(plDens)-3.774427
  mu <- invLogit(int+slopeHbee*transhbee+slopePol*pol+slopePlSize*transsize+ #Expected value
                   slopePlDens*transdens+slopeIrrig*isIrrig+slope2015*is2015)
  if(!addVar) return(mu)
  #Generate dispersion (phi)
  phi <- exp(intPhi+slopePlSizePhi*transsize) #Dispersion
  theta <- rbeta(N,mu*phi,(1-mu)*phi) #Generate from beta
  return(theta)
  # a <- rbinom(N,flwCount,theta)
}

#Generate flower counts from plant sizes
genFlwCount <- function(plSize,surv,is2015,int,slopePlSize,slopeSurv,slope2015,
                        intPhi,slopePlSizePhi,addVar=T){
  N <- getN(list(plSize,surv))
  #addVar: should use rbinom to simulate variance, or return expected value?
  #plSize=plant size, intercept,slope of plant size, phi
  transsize <- log(plSize)-2.626568
  transsurv <- logit(surv)-1.0653 #Centered logit of survival
  mu <- exp(int+slopePlSize*transsize+slopeSurv*transsurv+slope2015*is2015) #Expected value
  if(addVar){
    phi <- exp(intPhi+slopePlSizePhi*transsize) #Dispersion
    a <- rnbinom(N,mu=mu,size=phi)
    return(a)  
  } else (return(round(mu)))
}

#Generate pod counts from flower counts, flower survival rates
genPodCount <- function(flwCount,flwSurv){
  N <- getN(list(flwCount,flwSurv))
  a <- rbinom(N,flwCount,flwSurv)
  return(a)
}

#Generate avg seed counts per plant (seed count) - use stDev from plant level
genSeedCount <- function(hbee,pol,plSize,is2015,int,slopeHbee,slopePol,slopePlSize,
                         slope2015,stDev,transform=T){
  N <- getN(list(hbee,pol,plSize,is2015))
  transhbee <- log(hbee+0.5)
  transsize <- log(plSize)-2.626568
  mu <- rnorm(max(N),int + transhbee*slopeHbee + pol*slopePol + transsize*slopePlSize + is2015*slope2015,stDev)
  if(transform) return(exp(mu)) else return(mu)
}

#Generate avg seed size per plant (seed weight, mg) - use stDev from plant level
genSeedWeight <- function(hbee,pol,avgSeedCount,plSize,isIrrig,is2015,int,slopeHbee,
                          slopePol,slopeSeedCount,slopePlSize,slopeIrrig,
                          slope2015,lambda,stDev){
  N <- getN(list(hbee,pol,avgSeedCount,plSize,isIrrig,is2015))
  transhbee <- log(hbee+0.5) #Transform
  transsize <- log(plSize)-2.626568
  mu <- int + transhbee*slopeHbee + pol*slopePol + avgSeedCount*slopeSeedCount + transsize*slopePlSize + 
    isIrrig*slopeIrrig + is2015*slope2015 + (1/lambda)
  a <- rnorm(max(N),mu,stDev)
  return(a)
}

#Generate weight of all seeds at a given plot
genYield <- function(pods,seedCount,seedWeight,int,slopeYield,plotLev=T){
  calcYield <- pods*seedCount*seedWeight/1000 #Seed weight per plant (g)
  a <- exp(int+log(calcYield)*slopeYield)
  if(plotLev) return(sum(a)) else return(a) #Return plot-level yield (default)
}

#Simulation
simCommodity <- function(dist,is2015,isIrrig,isGP,numHives,dat,returnAll=F,useMean=F,plotVar=F,
                         other=NA){
  #returnAll: return all simulation parameters, or just yield?
  #useMean: use mean of distribution, or use random draw from it?
  #plotVar: simulate plot-level variance?
  #other: list of endogenous variables to fix at a specific value (e.g. fix plant density/size). Works only for plDens right now
  # rm(dist,is2015,isIrrig,isGP,numHives,dat,returnAll,useMean) #Cleanup
  # rm(simPlDens,simPlSize,simFlwDens,simHbeeVis,simAvgPol,simFlwCount,simPodCount,simSeedCount,simSeedWeight,simYield)
  # 
  # dist=1:400
  # is2015=0
  # isIrrig=0
  # isGP=0
  # numHives=40
  # dat=modPodcount
  # returnAll=F
  # useMean=F
  # plotVar=T
  
  if(useMean) dat <- lapply(dat,mean)
  
  
  #Simulated plant density
  if(!is.na(other)&grepl('simPlDens',names(other))){
    simPlDens <- other$simPlDens
  } else {
  simPlDens <- round(with(dat,genPlDens(dist=dist,is2015=is2015,isIrrig=isIrrig,isGP=isGP,
                      int=sample(intPlDens,1),slopeDist=sample(slopeDistPlDens,1),
                      slope2015=sample(slope2015PlDens,1),slopeIrrig=sample(slopeIrrigPlDens,1),
                      slopeIrrig2015=sample(slope2015IrrigPlDens,1),slopeGP=sample(slopeGPPlDens,1),
                      stDev=ifelse(plotVar,sample(sigmaPlDens,1),0),T)))
  }
  
  # #Compare to actual - looks OK
  # par(mfrow=c(2,1))
  # hist(surveyAllComm$PlDens,breaks=seq(0,200,10),xlab='Actual Plant Density',main=NULL)
  # hist(simPlDens,breaks=seq(0,200,10),xlab='Simulated Plant Density',main=NULL)
  
  #Simulated plant sizes within plots
  simPlSize <- with(dat,mapply(genPlSize,plDens=simPlDens,dist=dist,isGP=isGP,is2015=is2015,
                   isIrrig=isIrrig,int=sample(intPlDens,1),slopePlDens=sample(slopePlDensPlSize,1),
                   slopeDist=sample(slopeDistPlSize,1),slopeGP=sample(slopeGpPlSize,1),
                   slope2015=sample(slope2015PlSize,1),slopeIrrig=sample(slopeIrrigPlSize,1),
                   stDev=ifelse(plotVar,sample(sigmaPlSize,1),0),transform=T,SIMPLIFY=F))
  
  # #Compare to actual - looks OK
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simPlSize),plantsAllComm$VegMass,na.rm=T)+10
  # hist(plantsAllComm$VegMass,xlab='Actual Plant Size',main=NULL,breaks=seq(0,upr,10))
  # hist(unlist(simPlSize),xlab='Sim Plant Size',main=NULL,breaks=seq(0,upr,10))
  
  #Simulate flower density 
  simFlwDens <- with(dat,genFlwDens(avgPlSize=sapply(simPlSize,mean),
                       dist=dist,int=sample(intFlDens,1),slopePlSize=sample(slopePlSizeFlDens,1),
                       slopeDist=sample(slopeHbeeDistFlDens,1),
                       stDev=ifelse(plotVar,sample(sigmaFlDens,1),0),transform=T))
  # #Compare to actual - looks OK, but tends to vary quite a bit
  # par(mfrow=c(2,1))
  # hist(surveyAllComm$FlDens,breaks=seq(0,max(c(simFlwDens,surveyAllComm$FlDens))+100,100),xlab='Actual Flw Dens',main=NULL)
  # hist(simFlwDens,breaks=seq(0,max(c(simFlwDens,surveyAllComm$FlDens))+100,100),xlab='Sim Flw Dens',main=NULL)
  
  #Simulate hbee visits
  simHbeeVis <- with(dat,genHbeeVis(dist=dist,numHives=numHives,flDens=simFlwDens,is2015=is2015,
                      isGP=isGP,isIrrig=isIrrig,int=sample(intVisit,1),slopeDist=sample(slopeDistVis,1),
                      slopeHives=sample(slopeHiveVis,1),slopeFlDens=sample(slopeFlDens,1),
                      slope2015=sample(slopeYearVis,1),slopeGP=sample(slopeGpVis,1),
                      slope2015GP=sample(slopeYearGpVis,1),slopeIrrig=sample(slopeIrrigVis,1),
                      phi=sample(visitHbeePhi,1),addVar=plotVar))
  
  # par(mfrow=c(2,1))
  # upr=max(c(surveyAllComm$Honeybee,simHbeeVis))
  # with(surveyAllComm,hist(Honeybee/(TotalTime/10),breaks=seq(0,upr+2,2),xlab='Actual Visits/10mins',main=NULL))
  # hist(simHbeeVis,breaks=seq(0,upr+2,2),xlab='Sim Visits/10mins',main=NULL)
  
  #Simulate (average) pollen deposition
  simAvgPol <- with(dat,genAvgPol(hbee=simHbeeVis,dist=dist,int=sample(intPollen,1),
                      slopeHbee=sample(slopeVisitPol,1),slopeDist=sample(slopeHbeeDistPollen,1),
                      stDev=ifelse(plotVar,sample(sigmaPolPlot,1),0),center=T,transform=F))
  
  # #Compare to actual - looks OK
  # par(mfrow=c(2,1))
  # hist(with(dat,genAvgPol(hbee=simHbeeVis,dist=dist,int=sample(intPollen,1),slopeHbee=sample(slopeVisitPol,1),
  #     slopeDist=sample(slopeHbeeDistPollen,1),stDev=sample(sigmaPolPlot,1),center=F,transform=T)),
  #     xlab='Sim Avg pollen/plot',main=NULL,breaks=seq(0,3000,100))
  # hist(with(flowersAllComm,tapply(Pollen,paste(Field,Year,Distance),mean)),xlab='Actual Avg Pollen/plot',main=NULL,breaks=seq(0,3000,100))
  
  #Simulate flower survival rate
  simFlwSurv <- with(dat,mapply(genFlwSurv,hbee=simHbeeVis,pol=simAvgPol,plSize=simPlSize,
                                plDens=simPlDens,isIrrig=isIrrig,is2015=is2015,
                                int=sample(intFlwSurv,1),slopeHbee=sample(slopeVisitSurv,1),
                                slopePol=sample(slopePolSurv,1),slopePlSize=sample(slopePlSizeSurv,1),
                                slopePlDens=sample(slopePlDensSurv,1),slopeIrrig=sample(slopeIrrigSurv,1),
                                slope2015=sample(slope2015Surv,1),intPhi=sample(intPhiFlwSurv,1),
                                slopePlSizePhi=sample(slopePlSizePhiFlwSurv,1),addVar=plotVar,SIMPLIFY=F))
  
  # #Flower survivorship (proportion) - looks OK
  # par(mfrow=c(2,1))
  # hist(unlist(simFlwSurv),xlab='Sim flower survival',main=NULL,breaks=seq(0,1,0.05))
  # hist(with(plantsAllComm,Pods/(Pods+Missing)),xlab='Actual flower survival',main=NULL,breaks=seq(0,1,0.05))
  
  #Simulate flower count per plant
  simFlwCount <- with(dat,mapply(genFlwCount,plSize=simPlSize,surv=simFlwSurv,is2015=is2015,
                                 int=sample(intFlwCount,1),slopePlSize=sample(slopePlSizeFlwCount,1),
                                 slopeSurv=sample(slopeSurvFlwCount,1),slope2015=sample(slope2015FlwCount,1),
                                 intPhi=sample(intPhiFlwCount,1),
                                 slopePlSizePhi=sample(slopePlSizePhiFlwCount,1),addVar=plotVar,SIMPLIFY=F))
  
  # #Compare to actual - looks OK
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simFlwCount),with(plantsAllComm,Pods+Missing),na.rm=T)
  # hist(unlist(simFlwCount),xlab='Sim Flowers/plant',main=NULL,breaks=seq(0,upr+50,50))
  # abline(v=mean(unlist(simFlwCount)),col='red');abline(v=median(unlist(simFlwCount)),col='red',lty='dashed')
  # hist(with(plantsAllComm,Pods+Missing),xlab='Actual Flowers/plant',main=NULL,breaks=seq(0,upr+50,50))
  # abline(v=mean(with(plantsAllComm,Pods+Missing),na.rm=T),col='red');abline(v=median(with(plantsAllComm,Pods+Missing),na.rm=T),col='red',lty='dashed')

  simPodCount <- with(dat,mapply(genPodCount,simFlwCount,simFlwSurv,SIMPLIFY=F))
  
  # #Compare to actual - looks OK
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simPodCount),with(plantsAllComm,Pods),na.rm=T)
  # hist(unlist(simPodCount),xlab='Sim Pods/plant',main=NULL,breaks=seq(0,upr+50,50))
  # hist(with(plantsAllComm,Pods),xlab='Actual Pods/plant',main=NULL,breaks=seq(0,upr+50,50))
  
  #Simulate avg seeds per pod
  simSeedCount <- with(dat,mapply(genSeedCount,hbee=simHbeeVis,pol=simAvgPol,plSize=simPlSize,is2015=is2015,
                     int=sample(intSeedCount,1),slopeHbee=sample(slopeVisitSeedCount,1),
                     slopePol=sample(slopePolSeedCount,1),slopePlSize=sample(slopePlSizeCount,1),
                     slope2015=sample(slope2015SeedCount,1),
                     stDev=ifelse(plotVar,sample(sigmaSeedCount_plant,1),0),SIMPLIFY=F))
  
  # #OK, but spread seems low - maybe needs a t-dist at the plant level?
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simSeedCount),plantsAllComm$AvPodCount,na.rm=T)+2
  # hist(unlist(simSeedCount),breaks=seq(0,upr,2),main=NULL,xlab='Sim Avg Seed Count')
  # hist(plantsAllComm$AvPodCount,breaks=seq(0,upr,2),main=NULL,xlab='Actual Avg Seed Count')
  
  #Simulate avg seed weight (mg)
  simSeedWeight <- with(dat,mapply(genSeedWeight,hbee=simHbeeVis,pol=simAvgPol,avgSeedCount=simSeedCount,
                      plSize=simPlSize,isIrrig=isIrrig,is2015=is2015,int=sample(intSeedWeight,1),
                      slopeHbee=sample(slopeVisitSeedWeight,1),
                      slopePol=sample(slopePolSeedWeight,1),slopeSeedCount=sample(slopeSeedCount,1),
                      slopePlSize=sample(slopePlSizeWeight,1),slopeIrrig=sample(slopeIrrigSeedWeight,1),
                      slope2015=sample(slope2015SeedWeight,1),lambda=sample(lambdaSeedWeight,1),
                      stDev=ifelse(plotVar,sample(sigmaSeedWeight_plant,1),0),SIMPLIFY=F))
  
  # #Looks good
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simSeedWeight),(plantsAllComm$AvPodMass*1000)/plantsAllComm$AvPodCount,na.rm=T)+1
  # hist(unlist(simSeedWeight),xlab='Sim weight per seed (mg)',main=NULL,breaks=seq(0,upr,0.5))
  # hist(with(plantsAllComm,AvPodMass*1000/AvPodCount),xlab='Actual weight per seed (mg)',
  #      main=NULL,breaks=seq(0,upr,0.5))
  
  
  #Simulate total yield for each plot
  simYield <- with(dat,mapply(genYield,pods=simPodCount,seedCount=simSeedCount,seedWeight=simSeedWeight,
                     int=sample(intYield,1),slopeYield=sample(slopeYield,1),SIMPLIFY=T))
  
  # #Looks good. Tails aren't as long, but this is probably due to tails from sim seed count
  # par(mfrow=c(2,1))
  # hist(simYield,xlab='Sim yield/m2',main=NULL,breaks=seq(0,1500,50))
  # hist(with(plantsAllComm,{
  #  a <- SeedMass*PlDens
  #  b <- paste(Field,Distance)
  #  return(tapply(a,b,mean))
  #  }),xlab='Actual yield/m2',main=NULL,breaks=seq(0,1500,50))
  
  if(!returnAll){
    return(simYield)
  } else {
    a <- list(simPlDens,simPlSize,simFlwDens,simHbeeVis,simAvgPol,simFlwCount,simPodCount,simSeedCount,simSeedWeight,simYield)
    return(lapply(a,function(x) lapply(x,mean))) #Mean of each distance
  }
}


#Try out generic fields, using distance/irrig/year data from data
results <- with(surveyAllComm,replicate(500,simCommodity(dist=Distance,is2015=Year==2015,
                    isIrrig=Irrigated=='Irrigated',isGP=Area=='Grand Prairie',
                    numHives=NumHives,dat=mod3)))

results2 <- select(surveyAllComm,Field,Distance,Year) %>% #Simulated yield per plot
  unite(ID,Field:Year,sep='_',remove=F) %>% 
  bind_cols(data.frame(t(apply(results,1,function(x) quantile(x,c(0.5,0.05,0.95),na.rm=T))))) %>%
  rename('pred'='X50.','lwr'='X5.','upr'='X95.')
  
results3 <- ungroup(plantsAllComm) %>% select(Field,Distance,Year,SeedMass,PlDens) %>% #Actual yield per plot
  unite(ID,Field:Year) %>% group_by(ID) %>% 
  summarize(meanSeeds=mean(SeedMass,na.rm=T),plDens=first(PlDens)) %>% 
  transmute(ID,Yield=meanSeeds*plDens) %>% ungroup()
  
left_join(results2,results3,by='ID') %>% #Range looks OK
  ggplot(aes(x=Yield,y=pred))+
  # geom_point()+
  geom_pointrange(aes(ymax=upr,ymin=lwr),alpha=0.4)+
  geom_abline(intercept=0,slope=1,col='red')+
  labs(x='Actual Yield',y='Predicted Yield')


#Simulate pollination effects at generic fields
dists <- rep(seq(1,501,10),3)
hives <- rep(c(0,40,400),each=length(seq(1,501,10)))
results <- replicate(1000,simCommodity(dist=dists,is2015=1,isIrrig=0,
                     isGP=1,numHives=hives,dat=mod3,plotVar=F))

#Simulated yield distribution
results <- data.frame(dist=dists,numHives=hives,t(apply(results,1,function(x) quantile(x,c(0.5,0.05,0.95),na.rm=T)))) %>%
  rename('pred'='X50.','lwr'='X5.','upr'='X95.')

# #Results in g/m2
# results %>% mutate(numHives=factor(numHives,labels=c('NoHives','40Hives','400Hives'))) %>% 
#   ggplot(aes(x=dist,y=pred,col=numHives,fill=numHives))+
#   geom_ribbon(aes(ymax=upr,ymin=lwr,col=NULL),alpha=0.3,show.legend=F)+
#   geom_line(size=1)+
#   labs(y='Predicted yield (g/m2)',x='Distance',col='Stocking')

#Results in bu/acre
results %>% mutate(numHives=factor(numHives,labels=c('No Hives','40 Hives','400 Hives'))) %>% 
  mutate_at(vars(pred,upr,lwr),g2bushels) %>% 
  ggplot(aes(x=dist,y=pred,col=numHives,fill=numHives))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,col=NULL),alpha=0.3,show.legend=F)+
  geom_line(size=1)+
  labs(y='Predicted yield (bu/acre)',x='Distance',col='Stocking')


#Actual yield distribution
ggplot(aes(x=Distance,y=Yield))+geom_jitter(width=5,alpha=0.3)+
  geom_pointrange(data=allResults,aes(x=dist,y=pred,ymax=upr,ymin=lwr,col=scenario),position=position_dodge(width=1))


#Simulate plant density/size effects at generic fields
test <- replicate(100,simCommodity(dist=50,is2015=0,isIrrig=0,isGP=0,numHives=20,dat=modPodcount,returnAll=F,useMean=F,plotVar=F,other=list(simPlDens=c(1:200))))

plot(apply(test,1,mean))


# Test commodity field analysis, using piecewiseSEM --------------------------
library(piecewiseSEM)
library(lme4)

#This model is coarse, because quantities are all modeled at plot level, but this gets the general sense of the model, and allows d-sep tests to be made

#Datalist taken from analysis in Stan above
str(datalist)

#Plot-level 
plotlev <- with(datalist,data.frame(plot=1:Nplot,field=plotIndex[1:Nplot],
                                logDist=log(dist),hbeeVis=hbeeVis,totalTime=totalTime,
                                is2015=is2015[plotIndex[1:Nplot]],isIrrig=isIrrigated[plotIndex[1:Nplot]],
                                isGP=isGP[plotIndex[1:Nplot]],numHives=numHives[plotIndex[1:Nplot]],
                                pollen=NA,plantDens=NA,flDens=flDens,plSize=NA,flwCount=NA,
                                avgSurv=NA,avgSeedWeight=NA,avgSeedCount=NA))

plotlev$pollen[sort(unique(datalist$flowerIndex))] <- #log pollen
  with(datalist,tapply(pollenCount,flowerIndex,function(x) mean(log(x+1),na.rm=T)))
plotlev$plantDens[datalist$obsPlDens_ind] <- datalist$plDens_obs #log plant density
plotlev$plSize[sort(unique(datalist$plantIndex[datalist$obsPl_ind]))] <-  #Plant size
  with(datalist,tapply(plantSize_obs,plantIndex[obsPl_ind],function(x) mean(x,na.rm=T)))
#(log) Avg flower count
plotlev$flwCount[sort(unique(datalist$plantIndex))] <- 
  log(with(datalist,tapply(flwCount,plantIndex,mean,na.rm=T)))
#Avg (logit) flower survival
plotlev$avgSurv[sort(unique(datalist$plantIndex))] <- 
  logit(with(datalist,tapply(podCount/flwCount,plantIndex,mean,na.rm=T)))
#Avg weight per seed
plotlev$avgSeedWeight[sort(unique(datalist$plantIndex[datalist$podIndex]))] <- 
  with(datalist,tapply(seedMass*1000,plantIndex[podIndex],mean,na.rm=T))
#Avg (log) seed count
plotlev$avgSeedCount[sort(unique(datalist$plantIndex[datalist$podIndex]))] <- 
  with(datalist,tapply(seedCount,plantIndex[podIndex],function(x) mean(log(x),na.rm=T)))
plotlev <- plotlev %>%  drop_na() %>% mutate_at(vars(plot,field),factor) %>% 
  mutate(logVisRate=log(hbeeVis+0.5/totalTime)) %>% 
  mutate(pollen=scale(pollen))

#Plant size model
mod1 <- lmer(plSize~plantDens+logDist+numHives+
               isIrrig*plantDens+
               plantDens:numHives+
               logDist:numHives+
               # logDist:is2015+
               (1|field),data=plotlev,REML=F)
summary(mod1)
drop1(mod1,test='Chisq')

plantsAllComm %>% 
  mutate(logVeg=log(VegMass),logDens=scale(log(PlDens)),logDist=log(Distance),NumHives=scale(NumHives)) %>% 
  mutate(Year=Year==2015,Irrigated=Irrigated=='Irrigated') %>% 
  mutate(Plant=paste0(Field,Plant)) %>% 
  lmer(logVeg~logDens+NumHives*logDist+Irrigated+Year*logDist+
         (1|Field/Plant),data=.,REML=T) %>% 
  summary()


  drop1(.,test='Chisq')
summary(mod1)

pred1 <- expand.grid(plantDens=seq(-0.60,0.60,0.05),is2015=c(T,F),logDist=c(0,3,6),numHives=c(0,10,40),field=28)
pred1 <- cbind(pred1,pred=predict(mod1,newdata=pred1)) %>% 
  mutate(logDist=factor(logDist,labels=c('Dist=Near','Dist=Mid','Dist=Far'))) %>% 
  mutate(is2015=factor(is2015,labels=c('2014','2015')))
ggplot(pred1,aes(plantDens,pred,col=factor(numHives)))+geom_line(size=1)+
  facet_grid(logDist~is2015)+labs(col='NumHives')+labs(y='PlantSize',col='Hive\nNumber')

#Plant density model
mod2 <- lmer(plantDens~is2015*isIrrig+isGP+logDist+(logDist|field),data=plotlev,REML=F)
summary(mod2)
drop1(mod2,test='Chisq')  
summary(lm(plotlev$plantDens~cbind(predict(mod2)))) #R2=0.69

#Seed size model
plantsAllComm %>% ungroup() %>% 
  filter(Pods!=0,PropMissing>0,PropMissing<1) %>% 
  mutate(logVeg=log(VegMass),logDens=scale(log(PlDens)),logDist=log(Distance),NumHives=scale(NumHives)) %>% 
  mutate(Year=Year==2015,Irrigated=Irrigated=='Irrigated',Plot=paste0(Field,Distance)) %>% 
  lmer(logit(PropMissing)~PlDens+VegMass*Irrigated*Year-VegMass:Irrigated:Year-VegMass:Year+
         (1|Field/Plot),data=.,REML=F) %>% summary()



mod1 <- psem(lmer(plantDens~logDist+(1|field),data=plotlev),
             lmer(plSize~plantDens+isIrrig+isGP+is2015+(1|field),data=plotlev),
             lmer(flDens~plSize+isIrrig+isGP+is2015+(1|field),data=plotlev),
             lmer(flwCount~plSize+(1|field),data=plotlev),
             # glmer.nb(hbeeVis~offset(log(totalTime))+logDist+is2015+isGP+(1|field),data=plotlev),
             lmer(logVisRate~logDist+is2015+isGP+(1|field),data=plotlev),
             lmer(pollen~logVisRate+(1|field),data=plotlev),
             lmer(avgSurv~pollen+plSize+(1|field),data=plotlev),
             lmer(avgSeedCount~pollen+plSize+(1|field),data=plotlev),
             lmer(avgSeedWeight~pollen+plSize+(1|field),data=plotlev)
)

#Summary function doesn't work, so using lapply
summary(mod1)
lapply(mod1,summary)
dSep(mod1)
fisherC(mod1,conserve=T) #603.2, df=76, p=0

#Missing paths between plant size and hbeeVis, flDens and logDist, avgSeedCount and is2015
mod2 <- psem(lmer(plantDens~isGP+logDist+(1|field),data=plotlev),
             lmer(plSize~plantDens+isIrrig+isGP+is2015+(1|field),data=plotlev),
             lmer(flDens~plSize+logDist+isIrrig+isGP+is2015+(1|field),data=plotlev),
             lmer(flwCount~plSize+is2015+(1|field),data=plotlev),
             # glmer.nb(hbeeVis~offset(log(totalTime))+plSize+logDist+is2015+isGP+(1|field),data=plotlev),
             lmer(logVisRate~flDens+logDist+is2015*isGP+(1|field),data=plotlev),
             lmer(pollen~logVisRate+logDist+(1|field),data=plotlev),
             lmer(avgSurv~flwCount+pollen+logVisRate+plSize+(1|field),data=plotlev),
             lmer(avgSeedCount~flwCount+is2015+pollen+logVisRate+plSize+(1|field),data=plotlev),
             lmer(avgSeedWeight~flwCount+avgSurv+isIrrig+pollen+logVisRate+plSize+(1|field),data=plotlev)
)
summary(mod2[[8]])
lapply(mod2,summary)
dSep(mod2,conditioning=F)
fisherC(mod2,conserve=T) 

#This has a really weird set of dependencies at the end. Almost all metrics of yield are related to flower count

#Plant-level
plantlev <- with(datalist,data.frame(plant=1:Nplant,plot=plantIndex[1:Nplant],
                           field=plotIndex[plantIndex[1:Nplant]],plantSize=NA,podCount,flwCount,
                           avgSeedCount=NA,avgSeedWeight=NA,visit=log(0.5+hbeeVis/totalTime)[plantIndex],
                           distance=dist[plantIndex[1:Nplant]],pollen=NA))
#Plant size
plantlev$plantSize[datalist$obsPl_ind] <- datalist$plantSize_obs
#Average seed number
plantlev$avgSeedCount[unique(datalist$podIndex)] <- tapply(datalist$seedCount,datalist$podIndex,mean)
#Average weight per seed
plantlev$avgSeedWeight[unique(datalist$podIndex)] <- tapply(datalist$seedMass,datalist$podIndex,mean)
#Average pollen per plot
plantlev$pollen <-tapply(datalist$pollen,datalist$flowerIndex,function(x) log(mean(x)))[match(plantlev$plot,as.numeric(names(tapply(datalist$pollen,datalist$flowerIndex,mean))))]

library(lme4) #Model of pod success at plant level. No real effect of visitation, pollen deposition, or plant size
mod1 <- glmer(cbind(podCount,flwCount)~visit+scale(pollen)+plantSize+
                (1|field)+(1|plot)+(1|plant),data=plantlev,family='binomial')
summary(mod1)

#Pod-level
podlev <- with(datalist,data.frame(pod=1:Npod,plant=podIndex[1:Npod],
                           plot=plantIndex[podIndex[1:Npod]],
                           field=plotIndex[plantIndex[podIndex[1:Npod]]],
                           plantSize=NA,avgFlwSurv=NA,
                           seedCount=seedCount,seedMass=seedMass*1000,
                           visit=log(0.5+hbeeVis/totalTime)[plantIndex[podIndex[1:Npod]]],
                           distance=dist[plantIndex[podIndex[1:Npod]]],pollen=NA))
#Plant size
podlev$plantSize <- with(datalist,plantSize_obs[obsPl_ind[podIndex[podlev$pod]]])
#Pollen
podlev$pollen <- with(datalist,tapply(pollenCount,flowerIndex,function(x) log(mean(x))))[match(podlev$plot,sort(unique(datalist$flowerIndex)))]
#Pod survival (logit)
podlev$avgFlwSurv <- logit(with(datalist,(podCount/flwCount))[podlev$plant])
podlev <- mutate_at(podlev,vars(pod,plant,plot,field),factor) %>% drop_na() #Convert pod/plant/plot/field index to factor

mod2 <- lmer(sqrt(seedMass)~scale(seedCount)+scale(pollen)+visit+plantSize+(1|field)+(1|plot)+(1|plant),data=podlev)
plot(mod2)
qqnorm(resid(mod2)); qqline(resid(mod2));
par(mfrow=c(2,1));hist(resid(mod2),xlim=c(-1,1));curve(dnorm(x,0,0.31),-1,1);par(mfrow=c(1,1))
summary(mod2); 
drop1(mod2,test='Chisq') #Strong effect of seedcount


# Field-level visitation, nectar, and pollen deposition model (JAGS) -------------

datalistField <- with(fieldsAllComm,list(
  NFields=nrow(fieldsAllComm),
  NVar=length(unique(sapply(strsplit(as.character(Variety),' '),function(x) x[1]))),
  repOnes=rep(1,length(unique(sapply(strsplit(as.character(Variety),' '),function(x) x[1])))),
  Year=as.numeric(factor(Year)),
  Field=as.numeric(factor(paste(Year,Field))), #Fields
  NumHives=NumHives, #Number of hives/field
  BeeYard=as.numeric(BeeYard=='Stocked'), #Stocked/unstocked
  Variety=as.numeric(factor(sapply(strsplit(as.character(Variety),' '),function(x) x[1]))), #Variety
  Location=as.numeric(Area), #Location (GP or Lethbridge)
  Irrigated=as.numeric(Irrigated), #Irrigated
  AirTemp=Temp-mean(Temp),
  muVarietyNect=rep(0,length(unique(sapply(strsplit(as.character(Variety),' '),function(x) x[1])))),
  precVarietyNect=diag(0.1,length(unique(sapply(strsplit(as.character(Variety),' '),function(x) x[1]))))
))

datalistPlot <- with(filter(surveyAllComm,Distance!=1),list(
  NPlots=length(Distance),
  FieldPlot=as.numeric(factor(paste(Year,Field))), #Index for field
  Distance=Distance,
  centLogDist=log(Distance)-mean(log(c(5,20,100,500))), #Centered log of distance
  FlDens=FlDens, #Floral density (/m2)
  TotalTime=TotalTime, #Total time (mins)
  HbeeCount=Honeybee,
  Soilwater=Soilwater,
  HbeeVisRate=Honeybee/(FlDens*TotalTime/60) #Interaction rate (visits/flw*hr)
))

datalistFlw <- with(filter(flowersAllComm,Distance!=1),list(
  NFlowers=length(Distance),
  FieldFlw=as.numeric(factor(paste(Year,Field))), #Index for field
  PlotFlw=as.numeric(factor(paste(Year,Field,Distance))),
  PollenCount=Pollen,
  NectVol=Vol,
  Sugar=Sugar
))

datalist<-c(datalistField,datalistPlot,datalistFlw)

library(jagsUI)
library(coda)

startlist <- function() list(
  intNumHives=-0.15, #Hbee Count
  intNumHivesShape=-0.21,
  slopeAirtempHbee=0,
  slopeNumHives=1,
  r=.1,
  intMaxNect=0, #Max nectar
  slopeAirTempMaxNect=0,
  slopeIrrigMaxNect=0,
  intNectarProd=-2.7, #Nectar production
  slopeAirTempNect=-0.4,
  slopeIrrigNect=3.1,
  shapeNectVol=1.4 #shape factor
)

mod1<- jags(data=datalist,inits=startlist, #nonlinear version
            parameters.to.save=c(
              'intNumHives','intNumHivesShape', #Hbee visits
              'slopeAirtempHbee','slopeNumHives','slopeNumHivesShape','r',
              'intMaxNect','slopeAirTempMaxNect','slopeIrrigMaxNect', #Nectar levels
              'intNectarProd','slopeAirTempNect','slopeIrrigNect', 
              'shapeNectVol',
              'intPollen','slopePollen', #Pollen counts
              'hbeeCountFit','hbeeCountFitNew', #Posterior checks
              'nectVolFit','nectVolFitNew',
              'pollenCountFit','pollenCountFitNew'
              ),
               model.file='field_flower_level.txt',
               n.chains=1,n.adapt=8000,n.iter=5500,n.burnin=500,n.thin=5,parallel=F)
beepr::beep(2)
summary(mod1)
xyplot(mod1) #Mixing
densityplot(mod1) #Density plots

pp.check(mod1,'hbeeCountFit','hbeeCountFitNew') #Fit is OK, but strange.
pp.check(mod1,'nectVolFit','nectVolFitNew') #Good fit. Gamma version is better than normal version, generally.
pp.check(mod1,'pollenCountFit','pollenCountFitNew')

samples <- mod1$samples[[1]]
# autocorr(mcmc(samples))
pairs(matrix(samples,ncol=ncol(samples),dimnames=dimnames(samples)),upper.panel=function(x,y) text(mean(range(x)),mean(range(y)),round(cor(x,y),2),cex=0.75+1.5*abs(cor(x,y)))) #Pairplot

#Predict hbee counts - based on distance and air temperature
predFun<-function(samples,centLogDist,centAirTemp,Stocked){
  #powfun <- function(x,b0,LRC) b0+(-10-b0)*(1-exp(-exp(LRC)*x)) #nonlinear version
  temp<-matrix(NA,nrow=nrow(samples),ncol=length(centLogDist))
  for(i in 1:nrow(samples)){
    maxCount<-samples[i,'intNumHives']+samples[i,'slopeNumHives']*Stocked+samples[i,'slopeAirtempHbee']*centAirTemp
    beeCountShape <- samples[i,'intNumHivesShape']
    #temp[i,]<- exp(powfun(Distance,maxCount,beeCountShape))
    temp[i,]<-exp(maxCount+beeCountShape*centLogDist)
  }
  return(data.frame(fit=apply(temp,2,mean),
              upr=apply(temp,2,function(x) quantile(x,0.975)),
              lwr=apply(temp,2,function(x) quantile(x,0.025))))
}
mean(fieldsAllComm$Temp) #24 C is mean temperature

#Predictions at minimum temperature (16.9 C)
predMinTemp<-data.frame(Distance=seq(5,500,5),centLogDist=log(seq(5,500,5))-mean(log(c(5,20,100,500))),
                rbind(predFun(samples,log(seq(5,500,5))-mean(log(c(5,20,100,500))),min(datalist$AirTemp),0)),Temp=16.9)

predMeanTemp<-data.frame(Distance=seq(5,500,5),centLogDist=log(seq(5,500,5))-mean(log(c(5,20,100,500))),
                rbind(predFun(samples,log(seq(5,500,5))-mean(log(c(5,20,100,500))),0,0)),Temp=24)

#Max temperature (30.8 C)
predMaxTemp<-data.frame(Distance=seq(5,500,5),centLogDist=log(seq(5,500,5))-mean(log(c(5,20,100,500))),
                rbind(predFun(samples,log(seq(5,500,5))-mean(log(c(5,20,100,500))),max(datalist$AirTemp),0)),Temp=30.8)

pred<-rbind(predMinTemp,predMeanTemp,predMaxTemp) %>%
  # mutate(Temp=factor(Temp)) %>% 
  # str()
  mutate(Temp=factor(Temp,labels=c('Low (17 C)','Med (24 C)','High (31 C)')))

(p1<-ggplot(pred,aes(Distance,fit,col=Temp))+ #Fits for honeybee counts ~ distance
  geom_line(size=1)+
  geom_jitter(data=filter(surveyAllComm,Distance!=1),aes(Distance,Honeybee),col='black',size=1)+
  #geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  ylim(0,15)+
  labs(x='Distance (m)',y='Honeybee visits/10 mins',col='Air\nTemperature')+
  scale_colour_manual(values=c('blue','black','red')))
ggsave('../Figures/hbeeVisDistTemp.png',p1,width=9,height=6)


#Predict nectar volume
predFun<-function(samples,visRate,airtemp,irrig){ #Predicts nectar volume | visRate,airtemp,irrig
  possfun <- function(visitsHrFls,maxNect,nectProd) (maxNect/((visitsHrFls/nectProd)*maxNect+1))
  temp<-matrix(NA,nrow=nrow(samples),ncol=length(visRate))
  for(i in 1:nrow(samples)){
    nectProd<-exp(samples[i,'intNectarProd']+
                    samples[i,'slopeAirTempNect']*airtemp+
                    samples[i,'slopeIrrigNect']*irrig)
    maxNect<-exp(samples[i,'intMaxNect']+
                   samples[i,'slopeAirTempMaxNect']*airtemp+ #Effect of air temperature
                   samples[i,'slopeIrrigMaxNect']*irrig)
    temp[i,]<- possfun(visRate,maxNect,nectProd)
  }
  return(data.frame(fit=apply(temp,2,median),
                    upr=apply(temp,2,function(x) quantile(x,0.975)),
                    lwr=apply(temp,2,function(x) quantile(x,0.025))))
}

#Predictions at different temperatures
predMinTemp<-data.frame(visitsHrFls=seq(0,0.7,0.005),
                        rbind(predFun(samples,seq(0,0.7,0.005),min(datalist$AirTemp),0)),
                        Temp=16.9)
predMeanTemp<-data.frame(visitsHrFls=seq(0,0.7,0.005),
                        rbind(predFun(samples,seq(0,0.7,0.005),mean(datalist$AirTemp),0)),
                        Temp=24)
predMaxTemp<-data.frame(visitsHrFls=seq(0,0.7,0.005),
                        rbind(predFun(samples,seq(0,0.7,0.005),max(datalist$AirTemp),0)),
                        Temp=30.8)
pred<-rbind(predMinTemp,predMeanTemp,predMaxTemp) %>% 
  mutate(Temp=factor(Temp,labels=c('Low (17 C)','Med (24 C)','High (31 C)')))

(p2<-ggplot(pred,aes(visitsHrFls,fit,col=Temp))+
  geom_line(size=1)+
  geom_point(data=flowersAllComm,aes(x=Honeybee/(FlDens*Time/60),y=Vol,col=NULL))+
  ylim(0,2)+labs(y='Nectar volume (uL)',x='Visits per Hr per Flower',col='Air\nTemperature')+
  scale_colour_manual(values=c('blue','black','red')))

ggsave('../Figures/nectVolVisTemp.png',p2,width=9,height=6)

#Predictions at different irrigations

predControl<-data.frame(visitsHrFls=seq(0,0.7,0.005),Irrig='Dryland',
                        rbind(predFun(samples,seq(0,0.7,0.005),mean(datalist$AirTemp),0)))
predIrrig<-data.frame(visitsHrFls=seq(0,0.7,0.005),Irrig='Irrigated',
                         rbind(predFun(samples,seq(0,0.7,0.005),mean(datalist$AirTemp),1)))

pred<-rbind(predControl,predIrrig) %>% 
  mutate(Irrig=factor(Irrig))

(p3<-ggplot(pred,aes(visitsHrFls,fit,col=Irrig))+
  geom_line(size=1)+
  geom_point(data=flowersAllComm,aes(x=Honeybee/(FlDens*Time/60),y=Vol,col=NULL))+
  ylim(0,3)+labs(y='Nectar volume (uL)',x='Visits per Hr per Flower',col='Irrigation')+
  scale_colour_manual(values=c('darkorange','blue')))
ggsave('../Figures/nectVolVisIrrig.png',p3,width=9,height=6)

#Pollen count predictions

predFun<-function(samples,HbeeCount){
  temp<-matrix(NA,nrow=nrow(samples),ncol=length(HbeeCount))
  for(i in 1:nrow(samples)){
    temp[i,]<- exp(samples[i,'intPollen']+samples[i,'slopePollen']*HbeeCount)
  }
  return(data.frame(fit=apply(temp,2,mean),
                    upr=apply(temp,2,function(x) quantile(x,0.975)),
                    lwr=apply(temp,2,function(x) quantile(x,0.025))))
}

pred<-data.frame(HbeeCount=0:40,rbind(predFun(samples,0:40)))

(p4<-ggplot(pred,aes(HbeeCount,fit))+
  geom_line(size=1)+
  geom_point(data=flowersAllComm,aes(x=Honeybee,y=Pollen))+
  labs(y='Pollen count per stigma',x='Honeybee visits/10mins')+
  ylim(0,1500))
ggsave('../Figures/pollenVisits.png',p4,width=9,height=6)

#Nectar production
predFun<-function(samples,airtemp,irrig){ #Predicts nectar production | airtemp,irrig
  temp<-matrix(NA,nrow=nrow(samples),ncol=length(airtemp))
  for(i in 1:nrow(samples)){
    temp[i,]<-exp(samples[i,'intNectarProd']+
                    samples[i,'slopeAirTempNect']*airtemp+
                    samples[i,'slopeIrrigNect']*irrig)
    
  }
  return(data.frame(fit=apply(temp,2,mean),
                    upr=apply(temp,2,function(x) quantile(x,0.975)),
                    lwr=apply(temp,2,function(x) quantile(x,0.025))))
}

predControl<-data.frame(airtemp=24+c(-7:7),Irrig='Control',
                 rbind(predFun(samples,c(-7:7),0)))
predIrrig<-data.frame(airtemp=24+c(-7:7),Irrig='Irrigated',
                        rbind(predFun(samples,c(-7:7),1)))
pred<-rbind(predControl,predIrrig) %>% 
  mutate(Irrig=factor(Irrig))

p5<-ggplot(pred,aes(airtemp,fit,col=Irrig))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,col=NULL,fill=Irrig),alpha=0.3,show.legend=F)+
  geom_line(size=1)+
  scale_colour_manual(values=c('darkorange','blue'))+
  scale_fill_manual(values=c('darkorange','blue'))+
  labs(y=expression(paste('Nectar production (',mu,'L/hr)')),x='Air temperature (C)',col='Irrigation')+xlim(24,31)+
  ylim(0,1.5)


#ML functions for experimentation
library(MASS)
ggplot(surveyAllComm,aes(Distance,Honeybee,col=BeeYard))+geom_jitter()+
  geom_smooth(method='glm.nb',formula=y~log(x),se=F)+
  scale_colour_manual(values=c('orange','green'))+
  ylim(0,20)

filter(surveyAllComm,Distance!=1,BeeYard=='Stocked'|Honeybee<20) %>% 
  ggplot(aes(Distance,Honeybee))+geom_jitter()+
  facet_wrap(~BeeYard,ncol=1)+
  geom_smooth(method='glm.nb',formula=y~log(x),se=F)+
  scale_colour_manual(values=c('orange','green'))#+
  #ylim(0,20)

mod1ML<-glm.nb(Honeybee~log(Distance)*NumHives,data=surveyAllComm)
summary(mod1ML)

# Unused code (mixture of commodity and seed field data) -------------------------------------------------------------

# Try simulating data for seed count per pod
set.seed(1)
p1<-0.5 #Transition prob. from N1 to N2
N1 <- rep(100,1000) #Number of individuals at start
N2 <- rbinom(1000,N1,p1) #Number of individuals surviving first step

logit_p2<-7+N2*(-0.13) #Transition prob. from N2 to N3 depends on size of N2
p2<-exp(logit_p2)/(1+exp(logit_p2)) #Inverse-logit transform

N3 <- rbinom(1000,N2,p2) #Number of individuals surviving seconds step

par(mfrow=c(3,1))
hist(N1,xlim=c(0,100),breaks=50)
hist(N2,xlim=c(0,100),breaks=50)
hist(N3,xlim=c(0,100),breaks=50)
par(mfrow=c(1,1))
# 
# datalist <- list(
#   N1=100,
#   N3=N3,
#   N3_2=N3 #Hacky approach
# )
# 
# startlist <-function() list( #Starting values
#   a1=rnorm(1,0,1),
#   a2=rnorm(1,0,1),
#   b2=rnorm(1,0,1)
# )
# 
# writeLines("
# model{ 
#           #Priors
#            a1 ~ dnorm(0,0.1)
#            a2 ~ dnorm(0,0.1) 
#            b2 ~ dnorm(0,0.1)
#            
#              for(i in 1:length(N3)){ #Likelihood
#              
#              logit(p1[i]) <- a1 #Prob. of surviving to N2
#              N2[i] ~ dbin(p1[i],N1) T(N3_2[i],)  #Number of N1 that survive to N2
#              
#              logit(p2[i]) <- a2 + b2*N2[i] #Prob. of surviving to N3
#              N3[i] ~ dbin(p2[i],N2[i])  #Number of N2 that survive to N3		
#              }	
#            }
#            ", con="sim_pod_count_model.txt")
# 
# model <- jags(data=datalist,inits=startlist,
#               parameters.to.save=c('a1','a2','b2'),
#               model.file='model.txt',
#               n.chains=3,n.adapt=500,n.iter=2500,n.burnin=500,n.thin=2,parallel=T)
# 
# summary(model) #Doesn't converge. This seems to be indicative of a deeper problem.
# xyplot(model)


#Stan version

# datalistPod=with(temp[!is.na(temp$PodCount),],list(
#   Npod=length(PodCount),
#   # FertCount=rep(NA,length(PodCount)), #Trying this to see if I can trick Stan into doing discrete variables
#   #This doesn't work. See: https://stackoverflow.com/questions/35243449/how-to-deal-with-the-missing-data-in-stan
#   #Also, chapter 11 of Stan manual
#   SeedCount=PodCount
# ))
# 
# library("rstan")
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# mod1_stan = stan(file = 'pod_level.stan',data=datalistPod)
# summary(mod1_stan)
# pairs(mod1_stan)
#       
# stan_plot(mod1_stan,pars=c('int_Pol'))
# stan_dens(mod1_stan,pars=c('int_Pol'))
# 
# traceplot(mod1_stan)
# 
# #Pod count ~ weight model
# temp=filter(seedsAllComm,Year==2014,PodMass>0)
# 
# datalistPod=with(temp[!is.na(temp$PodCount),],list(
#   Npod=length(PodCount),
#   #SeedCount=PodCount,
#   PodMass=PodMass
# ))
# detach("package:rstan", unload=TRUE)

powfun <- function(x,b0,b1,LRC) b0+(b1-b0)*(1-exp(-exp(LRC)*x))
plot(exp(powfun(c(0:500),3,-5,-5)),type='l',ylab='Count',ylim=c(0,20))
lines(exp(powfun(c(0:500),3,-5,-2)),col='green')
lines(exp(powfun(c(0:500),3,-5,-8)),col='blue')

# #Playing around with fits for honeybee visitation
# library(lme4)
# temp <- with(datalist,{
#   data.frame(hbeeVis=c(hbeeVis,hbeeVis_extra),hbeeDist=c(hbee_dist,hbee_dist_extra),
#              lbeeDist=c(lbee_dist,lbee_dist_extra),Field=factor(c(plotIndex,plotIndex_extra)),
#              cent=c(isCent,isCent_extra),Fbay=c(isFBay,isFBay_extra),
#              lbeeStocking=c(lbeeStocking,lbeeStocking_extra),
#              totalTime=c(totalTime,totalTime_extra))}) %>%
#   mutate(hbeeDist=log(hbeeDist)-mean(log(hbeeDist)),lbeeDist=log(lbeeDist)-mean(log(lbeeDist)))
# hbeeMod <-glmer.nb(hbeeVis~offset(log(totalTime))+hbeeDist*lbeeDist+cent+Fbay+lbeeStocking+(1|Field),data=temp)
# summary(hbeeMod)
# AIC(hbeeMod)
# #Marginal plots
# mm1 <- data.frame(model.matrix(hbeeMod))
# data.frame(model.matrix(hbeeMod)) %>% #Lbee dist
#   mutate_at(vars(-lbeeDist),funs(mean)) %>%
#   transmute(pred=as.vector(as.matrix(.) %*% unname(fixef(hbeeMod))),res=resid(hbeeMod),
#             lbeeDist=with(datalist,c(lbee_dist,lbee_dist_extra))) %>%
#   arrange(lbeeDist) %>%
#   ggplot(aes(x=lbeeDist))+geom_point(aes(y=exp(pred+res)))+geom_line(aes(y=exp(pred)),col='red')+
#   labs(x='lbee Distance',y='Visitation rate')+xlim(0,50)
# 
# data.frame(model.matrix(hbeeMod)) %>% #hbee dist/bay center
#   mutate_at(vars(-hbeeDist,-centTRUE),funs(mean)) %>%
#   transmute(pred=as.vector(as.matrix(.) %*% unname(fixef(hbeeMod))),res=resid(hbeeMod),
#             hbeeDist=with(datalist,c(hbee_dist,hbee_dist_extra)),
#             cent=with(datalist,c(isCent,isCent_extra))) %>%
#   mutate(cent=factor(cent,labels=c('Edge','Center'))) %>%
#   arrange(cent,hbeeDist) %>%
#   ggplot(aes(x=hbeeDist,col=cent))+geom_point(aes(y=exp(pred+res)))+geom_line(aes(y=exp(pred)))+
#   labs(x='hbee Distance',y='Visitation rate')+xlim(0,400)+
#   scale_colour_manual(values=c('red','blue'))
# 
# #Same for leafcutter visitation
# temp <- with(datalist,{
#   data.frame(lbeeVis=c(lbeeVis,lbeeVis_extra),hbeeDist=c(hbee_dist,hbee_dist_extra),
#              lbeeDist=c(lbee_dist,lbee_dist_extra),Field=factor(c(plotIndex,plotIndex_extra)),
#              cent=c(isCent,isCent_extra),Fbay=c(isFBay,isFBay_extra),
#              lbeeStocking=c(lbeeStocking,lbeeStocking_extra),
#              totalTime=c(totalTime,totalTime_extra))}) %>%
#   mutate(hbeeDist=log(hbeeDist)-mean(log(hbeeDist)),lbeeDist=log(lbeeDist)-mean(log(lbeeDist)))
# lbeeMod <- glmer.nb(lbeeVis~offset(log(totalTime))+hbeeDist*cent+lbeeDist+lbeeStocking+lbeeStocking:hbeeDist+Fbay+(1|Field),data=temp)
# summary(lbeeMod)
# par(mfrow=c(2,1));plot(fitted(lbeeMod),resid(lbeeMod));qqnorm(resid(lbeeMod));qqline(resid(lbeeMod));par(mfrow=c(1,1));
# AIC(lbeeMod)
# 
# with(datalist,plot(c(hbee_dist,hbee_dist_extra),log(c(hbeeVis/totalTime,hbeeVis_extra/totalTime_extra)+0.5),pch=19,cex=0.5,xlab='Distance from hives (m)',ylab='Hbee visits/10 mins')) #Lots of zeros
# lines(seq(0,800,10),predict(hbeeMod,newdata=data.frame(hbeeDist=seq(0,800,10)),interval='none'))
# 
# # lbeeMod <- with(datalist,{
# #   data.frame(lbeeVis=c(lbeeVis/totalTime,lbeeVis_extra/totalTime_extra),
# #              lbeeDist=c(lbee_dist,lbee_dist_extra),Field=factor(c(plotIndex,plotIndex_extra)),
# #              hbeeDist=c(hbee_dist,hbee_dist_extra),
# #              cent=c(isCent,isCent_extra),Fbay=c(isFBay,isFBay_extra),
# #              totalTime=c(totalTime,totalTime_extra))
# # }) %>%  mutate(lbeeVis=log((lbeeVis/totalTime)+0.5),hbeeDist=log(hbeeDist)-mean(log(hbeeDist))) %>%
# #   nlmer(lbeeVis~SSasymp(lbeeDist,Asym,R0,lrc)~ hbeeDist + (Asym|Field) ,data=.,
# #         start=c(Asym=0.7493 ,R0=4.9948,lrc=-2.9177))
# # logLik(lbeeMod)
# # plot(lbeeMod)
# # qqnorm(resid(lbeeMod)); qqline(resid(lbeeMod));
# 
# #Marginal plots
# mm1 <- data.frame(model.matrix(lbeeMod))
# data.frame(model.matrix(lbeeMod)) %>% #Lbee dist
#   mutate_at(vars(-lbeeDist),funs(mean)) %>%
#   transmute(pred=as.vector(as.matrix(.) %*% unname(fixef(lbeeMod))),res=resid(lbeeMod),
#             lbeeDist=with(datalist,c(lbee_dist,lbee_dist_extra))) %>%
#   arrange(lbeeDist) %>%
#   ggplot(aes(x=lbeeDist))+geom_point(aes(y=exp(pred+res)))+geom_line(aes(y=exp(pred)),col='red')+
#   labs(x='lbee Distance',y='Visitation rate')+xlim(0,50)
# 
# data.frame(model.matrix(lbeeMod)) %>% #hbee dist/bay center
#   mutate_at(vars(-hbeeDist,-centTRUE,-hbeeDist.centTRUE),funs(mean)) %>%
#   transmute(pred=as.vector(as.matrix(.) %*% unname(fixef(lbeeMod))),res=resid(lbeeMod),
#             hbeeDist=with(datalist,c(hbee_dist,hbee_dist_extra)),
#             cent=with(datalist,c(isCent,isCent_extra))) %>%
#   mutate(cent=factor(cent,labels=c('Edge','Center'))) %>%
#   arrange(hbeeDist) %>%
#   ggplot(aes(x=hbeeDist,col=cent))+geom_point(aes(y=exp(pred+res)))+geom_line(aes(y=exp(pred)))+
#   labs(x='hbee Distance',y='Visitation rate')+xlim(0,400)+ylim(0,30)+
#   scale_colour_manual(values=c('red','blue'))
# 
# data.frame(model.matrix(lbeeMod)) %>% #hbee dist/stocking
#   mutate_at(vars(-hbeeDist,-lbeeStockingTRUE,-hbeeDist.lbeeStockingTRUE),funs(mean)) %>%
#   transmute(pred=as.vector(as.matrix(.) %*% unname(fixef(lbeeMod))),res=resid(lbeeMod),
#             hbeeDist=with(datalist,c(hbee_dist,hbee_dist_extra)),
#             stocking=with(datalist,c(lbeeStocking,lbeeStocking_extra))) %>%
#   mutate(stocking=factor(stocking,labels=c('Regular','Half'))) %>%
#   arrange(stocking,hbeeDist) %>%
#   ggplot(aes(x=hbeeDist,col=stocking))+geom_point(aes(y=exp(pred+res)))+geom_line(aes(y=exp(pred)))+
#   labs(x='hbee Distance',y='Visitation rate')+xlim(0,400)+ylim(0,30)+
#   scale_colour_manual(values=c('red','blue'))

#' # Data from Wang et al 2011 - ovule counts - unused --------------------------------
#' ovDat <- read.csv('wang2011ovData.csv') %>% 
#'   mutate(Freq=ifelse(Freq<0,0,Freq))
#' 
#' ovDatNew <- expand.grid(Freq=NA,Count=20:50,Day=unique(ovDat$Day))
#' for(i in unique(ovDat$Day)){
#'   ovDatNew[ovDatNew$Day==i,'Freq'] <- with(ovDat,approx(x=Count[Day==i],y=Freq[Day==i],xout=c(20:50)))$y #Interpolate
#'   ovDatNew$Freq[is.na(ovDatNew$Freq[ovDatNew$Day==i])] <- 0 #Turn NAs to zeros
#'   ovDatNew[ovDatNew$Day==i,'Freq'] <- ovDatNew[ovDatNew$Day==i,'Freq']/sum(ovDatNew[ovDatNew$Day==i,'Freq']) #Normalize
#' }
#' ggplot(ovDatNew,aes(Count,Freq,col=Day))+geom_line()
#' 
#' #Negative binomial/poisson from data has too much variance. Trying normal dist
#' par(mfrow=c(5,1))
#' for(i in unique(ovDatNew$Day)){
#'   pars <- data.frame(mean=rep(NA,1000),sd=rep(NA,1000))
#'   for(rep in 1:1000){
#'     #Simulated count data
#'     simCount <- with(ovDatNew[ovDatNew$Day==i,],sample(ovDatNew$Count,100,prob=ovDatNew$Freq,replace=T))
#'     simPar <- optim(c(0.5,1),function(x,count) -sum(dnorm(x=count,mean=x[1],sd=x[2],log=T)),count=simCount)$par
#'     pars$mean[rep] <- simPar[1]
#'     pars$sd[rep] <- simPar[2]
#'   }
#'   plot(c(20:45),dnorm(c(20:45),mean=median(pars$mean),sd=median(pars$sd)),ylab='Density',type='l',ylim=c(0,max(ovDatNew$Freq)),
#'        main=i)
#'   with(ovDatNew[ovDatNew$Day==i,],lines(Count,Freq,col='red'))
#' }
#' 
#' #overall
#' par(mfrow=c(1,1))
#' pars <- data.frame(mean=rep(NA,1000),sd=rep(NA,1000))
#' for(rep in 1:1000){
#'   #Simulated count data
#'   simCount <- with(ovDatNew,sample(ovDatNew$Count,100,prob=ovDatNew$Freq,replace=T))
#'   simPar <- optim(c(0.5,1),function(x,count) -sum(dnorm(x=count,mean=x[1],sd=x[2],log=T)),count=simCount)$par
#'   pars$mean[rep] <- simPar[1]
#'   pars$sd[rep] <- simPar[2]
#' }
#' plot(c(20:45),dnorm(c(20:45),mean=median(pars$mean),sd=median(pars$sd)),ylab='Density',type='l',ylim=c(0,0.15))
#' ovDatNew %>% group_by(Count) %>% summarize(Freq=sum(Freq)) %>% ungroup() %>% mutate(Freq=Freq/sum(Freq)) %>% 
#'   with(.,lines(Count,Freq,col='red'))
#' points(c(20:45),dpois(c(20:45),32.15),col='blue') #Poisson is too wide
#' 
#' #Mean=32.15, SD=2.66 for overall ovule distribution. Discrete version:
#' plot(10:55,log(pnorm(11:56,32.15,2.66)-pnorm(10:55,32.15,2.66)),ylab='LogLik')
#' 
#' 
#' 
#' # Pod weight-count modelling (JAGS) - unused -----------------------------
#' library(jagsUI)
#' setwd('./Models')
#' #Pod count model
#' datalist <- with(temp[!is.na(temp$PodCount)&temp$PodCount>0,],list( #Strips NAs and 0 counts
#'   Npod=length(PodCount), #Number of observed seed counts
#'   Nunique=length(unique(PodCount)), #Number of unique seed counts
#'   uniquePodCount=sort(unique(PodCount)),
#'   NuniquePodCount=as.vector(table(PodCount)),
#'   uniqueMatches=match(PodCount,sort(unique(PodCount))), #Matching index
#'   SeedCount=PodCount
#' ))
#' 
#' datalist <- c(datalist,with(flowersAllComm[!is.na(flowersAllComm$Pollen),],{list(
#'   Npollen=length(Pollen),
#'   PollenCount=Pollen
#' )})) #Append pollen data to list
#' 
#' datalist <- c(datalist,with(plantsAllComm[!is.na(plantsAllComm$Pods)&!is.na(plantsAllComm$Missing)&plantsAllComm$Missing>0,],{list(
#'   Nplants=length(Pods),
#'   Pods=Pods,PodsMissing=Missing
#' )})) #Append flower success data to list
#' 
#' datalist$SeedCount2 <- datalist$SeedCount
#' 
#' startlist <- function() list(
#'   #lambda=5.6, r=0.6,
#'   intPol=-1.15,
#'   slopePol1=-3.1,
#'   intPod=1.16,
#'   slopePod=-0.001,
#'   simPollenCount=datalist$SeedCount+10,
#'   simOvCount=datalist$SeedCount+10
#' )
#' 
#' params <- c(#'lambda','r', #Params for negbin pollination process
#'   # 'simPollenCount',
#'   # 'simOvCount',
#'   'simPollenZero','simFertFailPod','simPodFailPod','totalPodFail',
#'   #'LLpod','LLseed','lik',
#'   'simSeedCount','fitSeedCount','fitSimSeedCount',
#'   'intPol','slopePol'
#'   #'intPod','slopePod'
#' )
#' 
#' tick <- Sys.time()
#' mod2 <- jags(data=datalist,
#'              inits=startlist,
#'              parameters.to.save=params,
#'              model.file='pod_count_weight.txt',
#'              n.chains=3,n.adapt=200,n.iter=300,n.burnin=100,n.thin=1,parallel=T)
#' tock <- Sys.time()
#' tock-tick 
#' beep(1)
#' # save(mod2,file='./podCountResults.RData')
#' # load('podCountResults.Rdata') 
#' 
#' print(mod2)
#' traceplot(mod2,parameters=c('intPol','slopePol'))
#' # traceplot(mod2,parameters='totalPodFail')
#' pp.check(mod2,actual='fitSeedCount',new='fitSimSeedCount') #Variance = higher in actual dataset. Poisson performs OK, but is weirdly above 1:1 line
#' 
#' mod2samp <- mod2$sims.list #Samples from mod2
#' str(mod2samp)
#' 
#' plot(datalist$SeedCount,apply(mod2samp$LLpod,2,median),pch=19) #LL decreases with SeedCount. This should be in the opposite direction
#' points(datalist$SeedCount,apply(mod2samp$LLpod,2,max),pch=3) 
#' points(datalist$SeedCount,apply(mod2samp$LLpod,2,min),pch=3) 
#' 
#' hist(apply(mod2samp$LLseed,2,median))
#' 
#' 
#' 
#' #Predicted (median) seed counts vs actual
#' data.frame(med=apply(mod2samp$simSeedCount,2,median),
#'            upr=apply(mod2samp$simSeedCount,2,function(x) quantile(x,0.95)),
#'            lwr=apply(mod2samp$simSeedCount,2,function(x) quantile(x,0.05)),
#'            actual=datalist$SeedCount) %>%
#'   ggplot(aes(actual,med,ymax=upr,ymin=lwr))+geom_pointrange()+
#'   geom_abline(intercept=0,slope=1)+labs(x='Actual',y='Predicted')+
#'   coord_cartesian(xlim=c(0,55),ylim=c(0,70))
#' 
#' #Simulated seed count
#' hist(mod2samp$simSeedCount[10,])
#' table(mod2samp$simSeedCount[10,])
#' 
#' with(mod2samp,{ #Shape of likelihood functions
#'   lambda <- 5.685
#'   r <- 0.615
#'   p <- r/(r+exp(lambda))
#'   meanNegBin <- (r*(1-p))/p
#'   
#'   par(mfrow=c(3,1))
#'   # curve(invLogit(median(intPol)+(median(slopePol1)*((x-meanNegBin)/1000))+(median(slopePol2)*((x-meanNegBin)/1000)^2)),1,2500,
#'   #       ylab='p(pollenSurv)')
#'   # plot(invLogit(median(intPol)+(median(slopePol1)*((c(1:2500)-meanNegBin)/1000))+(median(slopePol2)*((c(1:2500)-meanNegBin)/1000)^2))*
#'   #         dnbinom(c(1:2500),prob=p,size=r),ylab='p(pollenSurv)*p(pollen)',pch=19)
#'   curve(invLogit(bpoint(x/1000,median(intPol),median(slopePol1),median(polBp),median(slopePol2))),1,2500,ylab='p(pollenSurv)')
#'   plot(invLogit(bpoint(1:2500/1000,median(intPol),median(slopePol1),median(polBp),median(slopePol2)))*dnbinom(c(1:2500),
#'                                                                                                               prob=p,size=r),ylab='p(pollenSurv)*p(pollen)',pch=19)
#'   
#'   plot(1:2500,((c(1:2500)-meanNegBin)/1000),pch=19,type='l')
#'   print(((c(1)-meanNegBin)/1000))
#'   par(mfrow=c(1,1))
#' })
#' 
#' #Pairplot of mean seeds vs other params
#' with(mod2samp,data.frame(naSeeds=apply(simSeedCount,1,function(x) sum(is.na(x))),
#'                          meanSeeds=apply(simSeedCount,1,function(x) mean(x,na.rm=T)),
#'                          #meanOv=apply(simOvCount,1,mean),
#'                          deviance,
#'                          intPol,slopePol1,slopePol2,polBp)) %>% pairs()
#' 
#' # #Ovule count compared to seed count
#' # data.frame(SeedCount=datalist$SeedCount,simOvMed=apply(mod2samp$simOvCount,2,median),
#' #            upr=apply(mod2samp$simOvCount,2,max),
#' #            lwr=apply(mod2samp$simOvCount,2,min)) %>% 
#' #   ggplot(aes(SeedCount,simOvMed))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
#' #   geom_abline(intercept=0,slope=1)
#' # 
#' # #Pollen count compared to seed count
#' # data.frame(SeedCount=datalist$SeedCount,simPolMed=apply(mod2samp$simPollenCount,2,median),
#' #            upr=apply(mod2samp$simPollenCount,2,max),
#' #            lwr=apply(mod2samp$simPollenCount,2,min)) %>%
#' #   ggplot(aes(SeedCount,simPolMed))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
#' #   geom_abline(intercept=0,slope=1)
#' # 
#' # data.frame(isna=is.na(mod2samp$simSeedCount[5,]),simSeed=mod2samp$simSeedCount[5,],
#' #            simPol=mod2samp$simPollenCount[5,],#simOv=mod2samp$simOvCount[5,],
#' #            seedCount=datalist$SeedCount) %>% 
#' #   ggplot(aes(simOv,simSeed))+geom_point()+
#' #   geom_point(aes(simOv,seedCount),col='red')
#' 
#' 
#' hist(mod2samp$simSeedCount[100,])
#' 
#' with(mod2samp,{
#'   par(mfrow=c(2,1))
#'   iter <- 51
#'   hist(simPollenCount[iter,],main='Simulated pollen',xlab='Count') #Simulated pollen (50th iteration)
#'   hist(simOvCount[iter,],main='Simulated ovules',xlab='Count') #Simulated fertilized ovules
#'   # hist(simFertSeed[iter,],main='Simulated fertilized seeds',xlab='Count',breaks=seq(0,52,1)) #Simulated fertilized ovules 
#'   hist(datalist$SeedCount,main='Acutual seed count',xlab='Count',breaks=seq(1,52,1))
#'   print(c(intPol[iter],slopePol1[iter],slopePol2[iter],polBp[iter]))
#'   par(mfrow=c(1,1))
#' })
#' 
#' 
#' 
#' #Test
#' par(mfrow=c(2,1))
#' simPolCount <- 
#'   simFertCount <- rpois(datalistPod$Npod,120*exp(mod2$q50$int.Pol)) #Simulate fert ov count
#' simSeedCount <- rpois(datalistPod$Npod,(exp(mod2$q50$slope.Fert*simFertCount)*simFertCount)+0.01) #Simulate seed count
#' hist(simFertCount,main=NULL,xlab='Fert Ov Count',xlim=c(0,60),breaks=30) 
#' hist(simSeedCount,main=NULL,xlab='Sim Seed Count',xlim=c(0,60),breaks=30)
#' hist(datalistPod$SeedCount,main=NULL,xlab='Actual Seed Count',xlim=c(0,60),breaks=30)
#' 
#' pp.check(mod2,'fit.SeedCount','fit.SeedCount.new') #Good
#' pp.check(mod2,'fit.PodMass','fit.PodMass.new') #Weird. Very close to 1:1, but offset by constant value. Does the variance change with podcount too?
#' 
#' par(mfrow=c(2,1))
#' hist(temp$PodMass,breaks=50,xlim=c(0,0.7),main='Actual Pod Mass',xlab='g')
#' hist(rlnorm(1000,-3.04,1/sqrt(1.964)),breaks=100,xlim=c(0,0.7),main='Simulated Pod Mass',xlab='g')
#' 
#' detach("package:jagsUI", unload=TRUE)
#' 
#' 
#' # Pod weight-count modeling - custom likelihood (Stan) - unused ------------------------------------
#' 
#' library(rstan)
#' setwd('./Models')
#' rstan_options(auto_write = TRUE)
#' options(mc.cores = 3)
#' 
#' #Rows without missing seed counts/weights
#' # noMiss <- with(temp,rowSums(is.na(cbind(PodCount,PodMass)))==0 & PodCount!=0)
#' 
#' #List structure for Stan
#' 
#' datalist <- with(temp[!is.na(temp$PodCount)&temp$PodCount>0,],list( #Strips NAs and 0 counts
#'   Npod=length(PodCount), #Number of observed seed counts
#'   Nunique=length(unique(PodCount)), #Number of unique seed counts
#'   uniquePodCount=sort(unique(PodCount)),
#'   NuniquePodCount=as.vector(table(PodCount)),
#'   uniqueMatches=match(PodCount,sort(unique(PodCount))), #Matching index
#'   SeedCount=PodCount
#' ))
#' 
#' datalist <- c(datalist,with(flowersAllComm[!is.na(flowersAllComm$Pollen),],{list(
#'   Npollen=length(Pollen),
#'   PollenCount=Pollen
#' )})) #Append pollen data to list
#' 
#' datalist <- c(datalist,with(plantsAllComm[!is.na(plantsAllComm$Pods)&!is.na(plantsAllComm$Missing)&plantsAllComm$Missing>0,],{list(
#'   Nplants=length(Pods),
#'   Pods=Pods,PodsMissing=Missing
#' )})) #Append flower success data to list
#' 
#' modPodcount = stan(file='pod_level.stan',data=datalist,iter=100,chains=1,control=list(adapt_delta=0.8),
#'                    init=function() list(intPolSurv=0,slopePolSurv=-3,intPodSurv=-2,slopePodSurv=3))
#' print(modPodcount)
#' traceplot(modPodcount) #Very poor traces, esp for pollen mu and phi. Likely what is happening is the mu and phi terms are being "influenced" by seed counts. Not sure how to separate these two, aside from heirarchical terms.
#' #Try running model without pollen included, and use results from below for pollen input.
#' #Same results occur even when using pollen. Try to run this in JAGS instead?
#' 
#' modPolcount <- stan(file='pollenMod.stan',data=datalist[c('Npollen','PollenCount')],iter=1000,chains=3)
#' print(modPolcount) #mu=293.75, phi=0.61
#' traceplot(modPolcount)
#' 
#' # podMod <- stan_model(file='pod_level.stan')
#' # optFit <- optimizing(podMod,data=datalist,init=list(
#' #   intPolSurv=-1,slopePolSurv=-3.5,intPodSurv=-2.9,slopePodSurv=3.5))
#' # round(optFit$par,4) #Weird. Chooses super-low survival probs.
#' 
#' #Check output
#' pars <- c('intPol','intFert','slopeFert')
#' print(modPodcount,pars=pars)
#' traceplot(modPodcount,pars=pars)
#' pairs(modPodcount,pars=pars)
#' 
#' # #Simulate seed count
#' # simcounts <- lapply(extract(modPodcount,pars=pars),median)
#' # invLogit <- function(x) exp(x)/(1+exp(x))
#' # par(mfrow=c(3,1))
#' # simFertCount <- rbinom(datalist$Npod,120,invLogit(simcounts[[1]])) #Simulate fert ov count
#' # simSeedCount <- rbinom(datalist$Npod,simFertCount,invLogit(simcounts[[2]]+simcounts[[3]]*simFertCount)) 
#' # hist(simFertCount,main=NULL,xlab='Fert Ov Count',xlim=c(0,100)) 
#' # hist(simSeedCount,main=NULL,xlab='Sim Seed Count',xlim=c(0,100),breaks=length(unique(simSeedCount)))
#' # hist(datalist$SeedCount,main=NULL,xlab='Actual Seed Count',xlim=c(0,100),breaks=length(unique(datalist$SeedCount)))
#' 
#' #Doesn't look very good. Perhaps try an arrival model of pollen?
#' # #Fit distributions
#' # poisFun <- function(x,dat) -sum(dpois(dat,x,log=T)) #Poisson
#' # logLikPois <- optimize(f=poisFun,lower=0,upper=1000,dat=flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)])
#' # geomFun <- function(x,dat) -sum(dgeom(dat,x,log=T)) #Geometric
#' # logLikGeom <- optimize(f=geomFun,c(0,1),dat=flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)])
#' # negbinFun <- function(x,dat) -sum(dnbinom(dat,mu=x[1],size=x[2],log=T)) #Neg Bin
#' # logLikNB <- optim(c(1,1),negbinFun,dat=flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)])
#' # 
#' # 2*2+(2*logLikNB$value) #AIC for NB - best
#' # 2*1+(2*logLikGeom$objective) #AIC for geom - worse
#' # 2*1+(2*logLikPois$objective) #AIC for poisson - much worse
#' # 
#' # par(mfrow=c(4,1)) #Plot results
#' # hist(rnbinom(datalist$Npod,mu=logLikNB[[1]][1],size=logLikNB[[1]][2]),main=NULL,xlab='Neg Bin',xlim=c(0,4000),breaks=seq(0,4000,20))
#' # hist(rgeom(datalist$Npod,logLikGeom[[1]][1]),main=NULL,xlab='Geom',xlim=c(0,4000),breaks=seq(0,4000,20))
#' # hist(rpois(datalist$Npod,logLikPois[[1]][1]),main=NULL,xlab='Pois',xlim=c(0,4000),breaks=seq(0,4000,20))
#' # hist(flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)],xlim=c(0,4000),breaks=seq(0,4000,20),main=NULL,xlab='Actual pollen counts')
#' 
#' 
#' # Test model to simulate seed counts - unused --------------------------------------
#' 
#' #Likelihood method: (non-vectorized version)
#' 
#' #p(Seed)*p(Ovule)*p(Pollen)
#' # prob <- array(NA,c(maxPolNum,maxOvNum,length(seedCount)))
#' # # prob <- matrix(NA,nrow=maxPolNum,ncol=maxOvNum)
#' # for(seed in 1:length(seedCount)){ #For each seed count
#' #   for(polCount in seedCount[seed]:maxPolNum){ #For each potential pollen count
#' #     for(ovCount in seedCount[seed]:min(maxOvNum,polCount)){ #For each potential ovule count
#' #       prob[polCount,ovCount,seed] <- dnbinom(polCount,mu=292.9317,size=0.6136358,log=T)+
#' #         dbinom(ovCount,polCount,0.05,log=T)+
#' #         dbinom(seedCount[seed],ovCount,0.01,log=T)
#' #     }
#' #   }
#' # }
#' 
#' #Values for pollination process:
#' #mu=292.9316584
#' #size=0.6136358
#' 
#' margProb <- function(coefs,dat,maxOvNum,maxPolNum,mu,size,lambda){
#'   
#'   seedCount <- c(1:maxOvNum) #Number of seeds to try
#'   
#'   #Replicate all levels of seed/ovule/pollen counts
#'   prob <- expand.grid(seedCount=seedCount,ovCount=1:maxOvNum,polCount=1:maxPolNum)
#'   #Remove impossible categories
#'   prob <- prob[(prob$ovCount>=prob$seedCount) & (prob$polCount>=prob$seedCount),]
#'   #Pollen survial prob
#'   prob$probPol <- invLogit(coefs[1]+coefs[2]*(prob$polCount/1000))
#'   #Pod survival prob
#'   # prob$probSeed <- invLogit(coefs[3]+coefs[4]*prob$seedCount)
#'   
#'   lpPollen <- dnbinom(1:maxPolNum,mu=mu,size=size,log=T) #Pollen LP
#'   lpOv <- dpois(1:maxOvNum,lambda=lambda,log=T) #Ovule LP
#'   # lpPolSurv <- matrix(NA,nrow=maxPolNum,ncol=maxPolNum) #Pollen survival LP
#'   # for(i in 1:maxPolNum){ #For each potential pollen receipt
#'   #   lpPolSurv[i,1:i] <- dbinom(1:i,i,invLogit(coefs[1]+coefs[2]*i)) #Pollen survival lp
#'   # }
#'   lpPodSurv <- dbinom(1,1,invLogit(coefs[3]+coefs[4]*(seedCount/10)),log=T) #Pod survival prob
#'   
#'   #Calculate log-prob for each
#'   prob$lp <- lpPollen[prob$polCount] +
#'     lpOv[prob$ovCount]+
#'     ifelse(prob$seedCount==prob$ovCount, #Fert prob
#'            #if seeds==ovules, successful pollen count could have been anything between #seeds and #pollen
#'            pbinom(prob$seedCount-1,prob$polCount,prob$probPol,log=T,lower.tail=F),
#'            #if seeds<ovules, only one way to get seed count
#'            # lpPolSurv[prob$seedCount,prob$polCount]
#'            dbinom(prob$seedCount,prob$polCount,prob$probPol,log=T)
#'     )+
#'     lpPodSurv[prob$seedCount]
#'   #Sum by seed counts
#'   seedLp <- with(prob,tapply(lp,seedCount,function(x) log(sum(exp(x))))) 
#'   # barplot(c(1-sum(exp(seedLp)),exp(seedLp)),pch=19,ylab='p(x)',xlab='Seed count')
#'   
#'   
#'   #Match -lp to actual counts, sum
#'   lp <- sum(seedLp[match(dat$SeedCount,names(seedLp))])+ #Seed count lp
#'     sum(dbinom(dat$FailFlw,dat$AllFlw,1-sum(exp(seedLp)),log=T)) #Flw abortion lp
#'   
#'   return(-lp) 
#' }
#' 
#' margProb(c(-2,-1,-3.0506205,4.580386),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32)
#' 
#' library(doSNOW)
#' cl <- makeCluster(3,type='SOCK')
#' clusterExport(cl,list('margProb','invLogit','datalist'))
#' 
#' #Try across various pollen survival probs (1st arg)
#' ll3 <- parSapply(cl,seq(-3,3,length.out=15),function(x) 
#'   margProb(c(x,-3.5249114,-2.8723591,0.3460556),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))
#' 
#' #Pollen survival slopes (2nd arg)
#' ll4 <- parSapply(cl,seq(-6,-2,length.out=15),function(x) 
#'   margProb(c(-1.0314562,x,-2.8723591,0.3460556),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))
#' 
#' #Try across various pod abortion intercepts (3rd arg)
#' ll2 <- parSapply(cl,seq(-6,0,length.out=15),function(x) 
#'   margProb(c(-1.0314562,-3.5249114,x,0.3460556),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))
#' 
#' #Try across various pod abortion slopes (4th arg)
#' ll <- parSapply(cl,seq(-3,3,length.out=15),function(x)
#'   margProb(c(-1.0314562,-3.5249114,-2.8723591,x),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))
#' 
#' par(mfrow=c(2,2))
#' plot(seq(-3,3,length.out=15),ll3,xlab='Pollen surv int | Pod surv Intercept,Slope',ylab='-Log Lik',pch=19,main='Arg1: pol surv int')
#' plot(seq(-6,-2,length.out=15),ll4,xlab='Pollen surv int | Pod surv Intercept,Slope',ylab='-Log Lik',pch=19,main='Arg2: pol surv slope') 
#' plot(seq(-6,0,length.out=15),ll2,xlab='Intercept | Slope values, Pollen surv',ylab='-Log Lik',pch=19,main='Arg3: pod surv int')
#' plot(seq(-3,3,length.out=15),ll,xlab='Pod surv slope | Pod surv int, Pol surv int, slope',ylab='-Log Lik',pch=19,main='Arg4: pod surv slope')
#' par(mfrow=c(1,1))
#' beep(1)
#' 
#' stopCluster(cl)
#' 
#' library(optimParallel)
#' cl <- makeCluster(3,type='SOCK')
#' clusterExport(cl,list('margProb','invLogit','datalist'))
#' setDefaultCluster(cl=cl)
#' 
#' # -1.031621,-3.524373,-2.871942,3.460029; LL:23946.96
#' optFit <- optimParallel(c(-1.031,-3.524,-2.872,3.46),margProb,
#'                         dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),
#'                         maxOvNum=52,maxPolNum=1500,mu=292.9316584,size=0.6136358,lambda=32,
#'                         method='L-BFGS-B',lower=c(-10,-10,-20,-1),upper=c(10,10,10,13),parallel=cl)
#' # optFit <- c(-1.031621,-3.524373,-2.871942,3.460029)
#' 
#' stopCluster(cl)
#' detach('package:optimParallel')
#' 
#' simSeeds <- function(maxNovules,coefs,plotResults=T,actualSeeds=NA,propDiff=F){
#'   #Idea: 
#'   #Pollination is a Negative Binomial process with mu=292.9,size=0.61
#'   mu <- coefs[1]
#'   size <- coefs[2]
#'   #Ovule number per flower is a Poisson process with lambda =~32. Wang et al 2011 found that there are cross-season difference, but they aren't that big (~30-35 ovules).
#'   lambda <- 32
#'   #Fertilization is a Binomial process, where Fert Ovules ~ Bin(Pollen,p)
#'   #If more ovules are fertilized than available, fert ovules = total ovules (integrate all outcomes < Novules)
#'   
#'   #Pod abortion is a Bernoulli process where Abortion ~ Bern(invLogit(int1+slope1*Nfert))
#'   int1 <- coefs[3]
#'   slope1 <- coefs[4]
#'   
#'   int2 <- coefs[5]
#'   slope2 <- coefs[6]
#'   
#'   set.seed(1)
#'   #Generate pollen counts for 1000 flowers
#'   simPol <- rnbinom(5000,mu=mu,size=size)
#'   #Generate ovule counts for flowers
#'   simOv <- rpois(5000,32)
#'   simOv[simOv>maxNovules] <- maxNovules #Set all ovules > maxNovules to maxNovules
#'   
#'   #Ovules are fertilized by pollen
#'   probFert <- invLogit(int1+slope1*simPol/1000)
#'   simFert <- rbinom(length(simPol),simPol,probFert) #Fixed perc of pollen makes it to ovules
#'   simFert[simFert>simOv] <- simOv[simFert>simOv] #If Npollen>Novules, reverts to ov number "first come first served"
#'   
#'   #Fertilized ovules become seeds
#'   simSeed <- simFert[simFert>0] #Fertilized (suriving pollen>1) ovules become seeds
#'   probSeed <- invLogit(int2+slope2*simSeed/10) #Prob of pod abortion ~ # of fertilized ovules
#'   simSeed <- simSeed[rbinom(length(simSeed),1,probSeed)==1] #Plant aborts pods depending on number of ovules
#'   
#'   if(plotResults==T & (length(actualSeeds)!=1 & !is.na(actualSeeds[1]))){
#'     #Plots
#'     par(mfrow=c(5,1),mar=c(5,5,2,5))
#'     #Pollen
#'     hist(simPol,main=NULL,xlab='Pollen grains',breaks=50,xlim=c(0,max(simPol)))
#'     par(new=T) #Pollen survival
#'     curve(invLogit(int1+slope1*x/1000),0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
#'     axis(side=4,cex.axis=1,col='red',col.axis='red') 
#'     mtext(side=4,line=3,'Pollen survival ',col='red')
#'     #Ovules
#'     hist(simOv,xlab='Ovules',xlim=c(0,max(simOv)),breaks=seq(0,max(simOv),1),main=NULL)
#'     #Fertilized ovules
#'     hist(simFert,main=NULL,xlab='Fertilized Ovules',breaks=seq(0,max(simFert),1),xlim=c(0,max(simOv)))
#'     par(new=T) #Survival prob
#'     curve(invLogit(int2+slope2*x/10),0,max(simFert),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
#'     axis(side=4,cex.axis=1,col='red',col.axis='red') 
#'     mtext(side=4,line=3,'Pod survival',col='red')
#'     #Seed counts
#'     hist(simSeed,main=NULL,xlab='Seeds/pod',xlim=c(0,max(simOv)),breaks=seq(0,max(simOv),1))
#'     #Actual seeds
#'     hist(actualSeeds,main=NULL,xlab='Actual Seed Count',xlim=c(0,max(simOv)),breaks=length(unique(datalist$SeedCount)))
#'     par(mfrow=c(1,1),mar= c(5, 4, 4, 2) + 0.1)
#'   } else if(propDiff==T) {
#'     uniqueCounts <- sort(unique(c(actualSeeds,simSeed))) #Unique counts (from sim or actual)
#'     props <- matrix(0,nrow=length(uniqueCounts),ncol=2) #Matrix to store proportions
#'     propsActual <- table(actualSeeds)/length(actualSeeds)
#'     propsSim <- table(simSeed)/length(simSeed)
#'     props[match(names(propsActual),uniqueCounts),1] <- propsActual
#'     props[match(names(propsSim),uniqueCounts),2] <- propsSim
#'     return(sum(abs(props[,1]-props[,2])))
#'   } else {
#'     return(list(simPol=simPol,simFert=simFert,simSeed=simSeed))
#'   }
#' }
#' 
#' simSeeds(maxNovules=52,coefs=c(292.9316584,0.6136358,optFit$par),plotResults=T,actualSeeds=datalist$SeedCount)
#' #Looks pretty good, actually
#' 
#' #Proportion missing seeds. Should be related to 
#' plantsAllComm %>% mutate(PropMis=Missing/(Pods+Missing)) %>% filter(!is.na(PropMis),PropMis>0) %>% 
#'   ggplot(aes(PropMis))+geom_histogram()
#' 

