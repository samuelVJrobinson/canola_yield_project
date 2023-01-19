#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN SEED CANOLA FIELDS (2014+2015)

# Libraries, functions, and seed DAG  ---------------------------------------------------------
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

# setwd('~/Projects/UofC/canola_yield_project/Models') #Multivac path
setwd('~/Documents/canola_yield_project/Models') #Galpern machine path

source('../helperFunctions.R')

library(ggdag)
nodeCoords <- data.frame(name=c('numHives','hbeeDist','lbeeDist','hbeeVis','lbeeVis','pollen',
                                'plSize','plDens','flDens','cent',
                                'flwCount','flwSurv','seedCount','seedWeight'),
                         labs=c('Number\nof Hives','Honey bee\nDistance','Leafcutter\nDistance',
                                'Honey bee\nVisits','Leafcutter\nVisits','Pollen\nCount',
                                'Plant\nSize','Plant\nDensity','Flower\nDensity','Bay Centre',
                                'Flowers\nper Plant','Pods\nper Plant','Seeds\nper Pod','Seed\nSize'),
                         x=c(0,0,1,0,1,0.5,
                             2.5,0,1,0.5,
                             3.5,4,2.5,3.5),
                         y=c(4,4,4,3,3,2,
                             0,1,1,0.5,
                             0,1,3,2))

#Specify model
seedDAG <- dagify(plDens ~ hbeeDist,
                  plSize ~ plDens + hbeeDist,
                  flDens ~ hbeeDist,
                  hbeeVis ~ flDens + hbeeDist + lbeeDist + cent,
                  lbeeVis ~ flDens + hbeeDist + lbeeDist + cent,
                  pollen ~ hbeeVis  + lbeeVis + cent + hbeeDist + lbeeDist + flDens,
                  flwCount ~ plSize + cent + flwSurv,
                  flwSurv ~ pollen + plSize + cent + hbeeDist + lbeeDist + flDens,
                  seedCount ~ pollen + plSize + cent + hbeeDist + flDens + flwSurv + flwCount,
                  seedWeight ~ pollen + seedCount + plSize + plDens + lbeeDist,
                  coords= list(x = setNames(nodeCoords$x,nodeCoords$name),
                               y = setNames(nodeCoords$y,nodeCoords$name)),
                  labels=setNames(nodeCoords$labs,nodeCoords$name)
)

# Load in data ----------------------

load('./datalist_seed.Rdata')

# #Seed field data
# load("../Seed field analysis/seedfieldDataAll.RData")
# fieldsAllSeed <- data.frame(allFields); plantsAllSeed <- data.frame(allPlants)
# pollenAllSeed <- data.frame(allPollen); seedsAllSeed <- data.frame(allSeeds);
# surveyAllSeed <- data.frame(allSurvey);
# plantsAllSeed$Field <- gsub('Unrah','Unruh',plantsAllSeed$Field) #Fixes spelling error
# seedsAllSeed$Field <- gsub('Unrah','Unruh',seedsAllSeed$Field)
# seedsAllSeed$EdgeCent <- ifelse(seedsAllSeed$EdgeCent=='Cent','Center',seedsAllSeed$EdgeCent) #Fixes cent/center
# rm(allFields,allPlants,allPollen,allSeeds,allSurvey,behav2015,visitors2016,nectar2016,folder,survey2015,survey2016)
# #Set 'negative' missing pods (mistake in counting) to NA.
# plantsAllSeed <- mutate(plantsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
# seedsAllSeed <- mutate(seedsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
# #Extra seed field data (from Riley W.)
# rileyExtra <- read.csv('../Seed field analysis/rileyExtra.csv')
# setwd('~/Projects/UofC/canola_yield_project')
# 
# # Filter for yearly model runs, if needed
# {
#   yfilt <- 2016
#   rileyExtra <- filter(rileyExtra,Year==yfilt)
#   fieldsAllSeed <- filter(fieldsAllSeed,Year==yfilt) %>% droplevels()
#   surveyAllSeed <- filter(surveyAllSeed,Year==yfilt) %>% droplevels()
#   pollenAllSeed <- filter(pollenAllSeed,Year==yfilt) %>% droplevels()
#   plantsAllSeed <- filter(plantsAllSeed,Year==yfilt) %>% droplevels()
#   seedsAllSeed <- filter(seedsAllSeed,Year==yfilt) %>% droplevels()
# }
# 
# #Organize names of extra sites from Riley
# rileyExtra$site <- factor(rileyExtra$site)
# rileyFields <- levels(rileyExtra$site)[!(levels(rileyExtra$site) %in% levels(fieldsAllSeed$Field))]
# samFields <- levels(fieldsAllSeed$Field)
# fieldsAllSeed$Field <- factor(fieldsAllSeed$Field,levels=c(samFields,rileyFields)) #Apply new ordering of fields to dataframes
# rileyExtra$site <- factor(rileyExtra$site,levels=c(samFields,rileyFields))
# surveyAllSeed$Field <- factor(surveyAllSeed$Field,levels=c(samFields,rileyFields))
# pollenAllSeed$Field <- factor(pollenAllSeed$Field,levels=c(samFields,rileyFields))
# plantsAllSeed$Field <- factor(plantsAllSeed$Field,levels=c(samFields,rileyFields))
# seedsAllSeed$Field <- factor(seedsAllSeed$Field,levels=c(samFields,rileyFields))
# plantsAllSeed <- plantsAllSeed %>% group_by(Field,Distance,EdgeCent) %>% mutate(Plant=1:n()) %>% ungroup()
# #Fix date structure in Riley's fields
# rileyExtra$date <- as.character(rileyExtra$date)
# rileyExtra$date[grepl('-',rileyExtra$date)] <- as.character(as.Date(rileyExtra$date,format='%d-%b-%y')[grepl('-',rileyExtra$date)])
# rileyExtra$date[grepl('/',rileyExtra$date)] <- as.character(as.Date(rileyExtra$date,format='%m/%d/%Y')[grepl('/',rileyExtra$date)])
# rileyExtra$date <- as.Date(rileyExtra$date,format='%F')
# 
# rileyExtra <- rileyExtra %>%
#   #Filter out NA hdist/ldist plots - hard to impute, and no downstream info
#   filter(!is.na(hdist),!is.na(ldist)) %>%
#   filter(hdist<=400) #Get rid of plots >400m away from bees (other side of the field)
# 
# # #Joins plant ID to seedsAllSeed
# # temp <- select(plantsAllSeed,Year,Field,Distance,EdgeCent,Branch,Pods,Missing,Plant) %>%
# #   unite(ID,Year:Missing)
# # seedsAllSeed <- seedsAllSeed %>%
# #   unite(ID,Year,Field,Distance,EdgeCent,Branch,Pods,Missing,remove=F) %>%
# #   left_join(temp,by='ID') %>%  select(-ID)
# # rm(temp)
# # #Add average pollen per plot
# # surveyAllSeed <- surveyAllSeed %>%
# #   unite(plot,Field,Distance,EdgeCent,remove=F) %>%
# #   left_join(summarize(group_by(unite(pollenAllSeed,plot,Field,Distance,EdgeCent),plot),
# #                       polCountPlot=log(mean(Pollen))),by='plot')
# 
# datalistField <- list( #Field-level measurements
#   Nfield=length(fieldsAllSeed$Year), #Number of fields
#   Nfield_extra=length(rileyFields) #Number of extra fields from Riley
# )
# 
# datalistPlot <- with(surveyAllSeed,list( #Plot-level measurements
#   Nplot=length(Distance), #Number of plots
#   plotIndex=as.numeric(Field), #Index for field (which field is plot from?)
#   lbeeStocking=Treatment=='Double tent', #Half leafcutter stocking?
#   lbeeStocking2=as.matrix(model.matrix(~Treatment)[,c(2:3)]), #Model matrix for 3 levels of stocking (2x tent or 2x tent + bees)
#   is2016=Year==2016, #Is field from 2016?
#   hbee_dist=Distance, #Distance from honeybee hives
#   hbeeVis=hbee, #Visits by honeybees
#   lbee_dist=minDist, #Centered distance from leafcutters
#   lbeeVis=lbee, #Leafcutter bee visits
#   isCent=EdgeCent=='Center', #Is plot from bay center?
#   isMBay=Bay=='M', #Is plot from M bay?
#   totalTime=TotalTime/10, #Total time (mins/10)
#   plotList=paste(Field,Distance,Bay,EdgeCent),
# 
#   # Flower density
#   Nplot_flDensObs = sum(!is.na(FlDens)), #Number of plots with observed flDens
#   Nplot_flDensMiss = sum(is.na(FlDens)), #Number of plots missing flDens
#   flDens_obs=sqrt(FlDens*4)[!is.na(FlDens)], #(sqrt) Flower density - these were recorded in a 0.25m2 plot
#   obsflDens_ind = which(!is.na(FlDens)), #Index for observed
#   missflDens_ind = which(is.na(FlDens)), #Index for missing
# 
#   #Female-only plots (plant density)
#   Nplot_F=sum(Bay=='F'), #Number of female plots
#   plotIndex_F=match(1:length(Distance),which(Bay=='F')), #Index for female-only plots (which female plot j does plot i belong to?)
#   plotIndex_F2=match(which(Bay=='F'),1:length(Distance)), #Reverse index (which plot i does female plot j belong to?)
# 
#   #Plant density
#   Nplot_plDensObs = sum(Bay=='F' & !is.na(PlDens)), #Number of F plots with observed plDens
#   Nplot_plDensMiss = sum(Bay == 'F' & is.na(PlDens)), #Number of F plots missing plDens
#   plDens_obs=log(PlDens)[Bay=='F' & !is.na(PlDens)], # (log) Plant density in F plots
#   obsPlDens_ind = which(Bay=='F' & !is.na(PlDens)), #Index for observed
#   missPlDens_ind = which(Bay=='F' & is.na(PlDens)) #Index for missing
# 
# ))
# datalistPlot$plotIndex_F[is.na(datalistPlot$plotIndex_F)] <- 0 #Set male plot indices to zero
# datalistPlot$totalTime[is.na(datalistPlot$totalTime)] <- 0.5 #Fix one missing time point
# 
# #Join in extra data from Riley
# datalistPlot_extra <- with(rileyExtra,
#                             list(
#                               Nplot=length(ldist), #Number of extra plots
#                               plotIndex=as.numeric(site), #Index for field (which field?)
#                               is2016=Year==2016,
#                               lbeeStocking=treatment=='Double tent',
#                               lbeeStocking2=as.matrix(cbind(as.numeric(treatment=='Double tent'),as.numeric(treatment=='Double tent and bees'))),
#                               hbee_dist=hdist,
#                               hbeeVis=hbee_vis,
#                               lbee_dist=ldist,
#                               lbeeVis=lbee_vis,
#                               isCent=rep(FALSE,length(ldist)), #Riley's plots were at edge of bay
#                               isMBay=Bay=='Male', #Is plot from M bay?
#                               totalTime=rep(10,length(ldist))/10, #Riley used 10 mins for everything
#                               #(sqrt) Flower density
#                               Nplot_flDensObs = sum(!is.na(flDens)), #Number of observed plots
#                               Nplot_flDensMiss = sum(is.na(flDens)), #Number of missing plots
#                               flDens_obs=sqrt(flDens)[!is.na(flDens)], #Observed - these were recorded in a 1m2 plot
#                               obsflDens_ind = which(!is.na(flDens)), #Index for observed
#                               missflDens_ind = which(is.na(flDens)) #Index for missing
#                             ))
# 
# names(datalistPlot_extra) <- paste0(names(datalistPlot_extra),'_extra') #Append "extra"
# 
# datalistPlot <- c(datalistPlot,datalistPlot_extra); rm(datalistPlot_extra)
# 
# datalistFlw <- with(pollenAllSeed,list( #Pollen samples - no NAs
#   Nflw=length(Pollen), #Number of pollen samples
#   flowerIndex=match(paste(Field,Distance,'F',EdgeCent),datalistPlot$plotList),#Index for flower (which plot?)
#   pollenCount=Pollen
# ))
# 
# datalistPlant <- plantsAllSeed %>%
#   select(Year:EdgeCent,VegMass,SeedMass,Pods:Plant) %>%
#   filter(!is.na(Pods),!is.na(Missing)) %>% #Filter out plants with missing pods/flw counts - only 6 plants, and hard to impute
#   mutate(AvgSeedMass=1000*AvPodMass/AvPodCount) %>% #Average weight per seed
#   with(.,list(
#     Nplant=length(Distance), #Number of plant samples
#     plantIndex=match(paste(Field,Distance,'F',EdgeCent),datalistPlot$plotList), #Index for plant (which F plot?)
#     plantSize=log(VegMass[!is.na(VegMass)]), #(log) weight of veg mass (g)
# 
#     podCount=Pods, #Successful pods
#     flwCount=(Pods+Missing), #Pods + Missing (total flw production)
#     logitFlwSurv = logit(Pods/(Pods+Missing)), #Logit flw survival
# 
#     #Average seeds per pod
#     Nplant_seedCountObs = sum(!is.na(AvPodCount)), #Number of observed
#     Nplant_seedCountMiss = sum(is.na(AvPodCount)), #Number of missing
#     seedCount_obs=AvPodCount[!is.na(AvPodCount)], #Avg seeds per pod - NAs present
#     obsSeedCount_ind = which(!is.na(AvPodCount)), #Observed
#     missSeedCount_ind = which(is.na(AvPodCount)), #Missing
# 
#     #Average mass per seed
#     Nplant_seedMassObs = sum(!is.na(AvgSeedMass)), #Number of observed
#     Nplant_seedMassMiss = sum(is.na(AvgSeedMass)), #Number of missing
#     seedMass_obs = AvgSeedMass[!is.na(AvgSeedMass)], #Avg weight per seed (g/1000 seeds)
#     obsSeedMass_ind = which(!is.na(AvgSeedMass)), #Index for observed
#     missSeedMass_ind = which(is.na(AvgSeedMass)), #Index for missing
# 
#     yield=SeedMass[!is.na(SeedMass)], #Weight of all seeds (g)
#     calcYield = AvPodCount*AvgSeedMass*Pods/1000
#     # plantList=paste(Field,Distance,'F',EdgeCent,Plant) #Name of plot (character)
#   ))
# datalistPlant$calcYield <- ifelse(is.na(datalistPlant$calcYield),-1,datalistPlant$calcYield) #Replaces NAs with -1s
# 
# # datalistPod <- seedsAllSeed %>% ungroup() %>%
# #   mutate(plantList=paste(Field,Distance,'F',EdgeCent,Plant)) %>% select(-Field,-Distance,-EdgeCent,-Plant) %>%
# #   mutate(inPlantList=plantList %in% datalistPlant$plantList) %>%
# #   filter(inPlantList) %>% #Filter out pods where plants had missing flower counts
# #   filter(!is.na(PodCount),!is.na(PodMass)) %>% #Filter out pods with missing seed counts
# #   with(.,list(
# #   Npod=length(PodCount), #Number of seeds measured
# #   seedCount=PodCount, #Number of seeds per pod
# #   seedMass=1000*PodMass/PodCount, #Weight per seed (mg)
# #   #Index for pod (which plant?)
# #   podIndex=match(plantList,datalistPlant$plantList)
# # ))
# datalistPlant$plantList <- datalistPlot$plotList <-  NULL
# datalist <- c(datalistField,datalistPlot,datalistFlw,datalistPlant)
# 
# flDensMean <- mean(with(datalist,c(flDens_obs,flDens_obs_extra))) #Center flower density
# datalist$flDens_obs <- datalist$flDens_obs-flDensMean
# datalist$flDens_obs_extra <- datalist$flDens_obs_extra-flDensMean
# datalist$flDensMean <- flDensMean
# 
# if(any(sapply(datalist,function(x) sum(is.na(x)))!=0)){ beep(1); print("NAs found in datalist")}
# 
# rm(datalistField,datalistPlot,datalistFlw,datalistPlant,rileyFields,samFields,flDensMean) #Cleanup
# str(datalist)
# 
# save(datalist,file = './datalist_seed.Rdata')

# Run models ----------------------------

# library(shinystan)

#Models split into separate stan files (faster)
modFiles <- dir(pattern = 'seed_.*\\.stan')
modList <- vector(mode = 'list',length = length(modFiles))
names(modList) <- gsub('(seed_.*[0-9]{2}|\\.stan)','',modFiles)

modList[1] <- stan(file=modFiles[1],data=datalist,iter=2000,chains=4,init=0.1) #Plant density, Plant size, Flower Density - OK
modList[2] <- stan(file=modFiles[2],data=datalist,iter=2000,chains=4,init=0) #Hbee visitation - OK, but traces for RE aren't great
modList[3] <- stan(file=modFiles[3],data=datalist,iter=2000,chains=4,init=0) #Lbee visitation - OK
modList[4] <- stan(file=modFiles[4],data=datalist,iter=2000,chains=4,init=0) #Pollen - OK
modList[5] <- stan(file=modFiles[5],data=datalist,iter=2000,chains=4,init=0) #Flower count per plant - OK
modList[6] <- stan(file=modFiles[6],data=datalist,iter=2000,chains=4,init=0) #Pod count (flw survival) per plant
modList[7] <- stan(file=modFiles[7],data=datalist,iter=2000,chains=4,init=0) #Seeds per pod
modList[8] <- stan(file=modFiles[8],data=datalist,iter=2000,chains=4,init=0) #Weight per seed
modList[9] <- stan(file=modFiles[9],data=datalist,iter=3000,chains=4,init=0,control=list(adapt_delta=0.9)) #Yield
beepr::beep(1)

apply(extract(modList[[4]],pars='pollenMu_plot')[[1]],2,mean) %>% range

#Traceplots
for(i in 1:length(modList)){
  if(!is.null(modList[[i]])){
    n <- names(modList[[i]]) #Model parameters
    n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('([sS]igma|slope)',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term) or "slope"
    n <- n[!grepl('_miss',n)] #Gets rid of imputed values
    
    p <- traceplot(modList[[i]],pars=n,inc_warmup=FALSE) #+ geom_hline(yintercept = 0) #Traceplots
    fastPairs(modList[[i]],pars=n)
    print(modList[[i]],pars=n)
    # p <- stan_plot(modList[[i]],pars=n) #Pointrange plot
    
    print(p)
    a <- readline('Press Return to continue: ')
    if(a!='') break()
  }
}

# Get model summaries into a list of tables
# modSummaries_seed <- vector(mode = 'list',length = length(modList)) #Create new empty list
# names(modSummaries_seed) <- names(modList)
load('modSummaries_seed.Rdata') #Load existing list
#Update model summaries if needed
temp <- lapply(modList,parTable) #Get parameter summaries
for(i in 1:length(temp)){
  if(class(temp[[i]])=='data.frame'||length(temp[[i]])>1){ #If temp is not empty (model not run)
    modSummaries_seed[[i]]$summary <- temp[[i]]$summary #Overwrite
    modSummaries_seed[[i]]$covMat <- temp[[i]]$covMat
  }
}
parNames <- lapply(modSummaries_seed,function(x) x$summary$param) #Clean up extra parameter names
for(i in 2:length(modSummaries_seed)){
  if(!is.null(modSummaries_seed[[i]])){
    chooseThese <- !parNames[[i]] %in% unlist(parNames[1:(i-1)])
    modSummaries_seed[[i]]$summary <- modSummaries_seed[[i]]$summary[chooseThese,]
    modSummaries_seed[[i]]$covMat <- modSummaries_seed[[i]]$covMat[chooseThese,chooseThese]
  }
}
save(modSummaries_seed,file = 'modSummaries_seed.Rdata'); rm(temp)

#Posterior predictive checks
PPplots(modList[[1]],datalist$plDens_obs,c('predPlDens','plDens_resid','predPlDens_resid'),
        index = datalist$obsPlDens_ind,main='Plant Density') #OK, but bimodal peaks
PPplots(modList[[1]],datalist$flDens_obs,c('predFlDens','flDens_resid','predFlDens_resid'),
        index=datalist$obsflDens_ind,main='Flower density') #OK, but skewed distribution. Not improved by MBay or Year terms
PPplots(modList[[1]],datalist$plantSize,c('predPlSize','plSize_resid','predPlSize_resid'),
        'Plant size') #OK
PPplots(modList[[2]],c(datalist$hbeeVis,datalist$hbeeVis_extra),
        c('predHbeeVis_all','hbeeVis_resid','predHbeeVis_resid'),'Honeybee visits',
        ZIpar = 'thetaHbeeVis') #Not good
PPplots(modList[[3]],c(datalist$lbeeVis,datalist$lbeeVis_extra),
        c('predLbeeVis_all','lbeeVis_resid','predLbeeVis_resid'),'Leafcutter visits',
        ZIpar='thetaLbeeVis') #OK, but overpredicts at low values
PPplots(modList[[4]],datalist$pollenCount,c('predPollenCount','pollen_resid','predPollen_resid'),
        'Pollen') #OK
PPplots(modList[[5]],datalist$flwCount,c('predFlwCount','flwCount_resid','predFlwCount_resid'),
        'Flowers per plant') #PP not great, but everything else is fine - probably due to log-scaling issue
PPplots(modList[[6]],datalist$podCount,c('predPodCount','podCount_resid','predPodCount_resid'),'Pods per plant') #OK
PPplots(modList[[7]],datalist$seedCount_obs,c('predSeedCount','seedCount_resid','predSeedCount_resid'),
        index=datalist$obsSeedCount_ind,'Seeds per pod') #OK, but distribution is weird, even with exp-normal. 
PPplots(modList[[8]],datalist$seedMass_obs,c('predSeedMass','seedMass_resid','predSeedMass_resid'),
        index=datalist$obsSeedMass_ind,'Seed Mass')
PPplots(modList[[9]],log(datalist$yield),c('predYield','yield_resid','predYield_resid'),'Seed mass per plant')

# Examine random intercept distributions
compareRE(modList[[1]],'intPlDens_field') 
compareRE(modList[[1]],'intPlSize_field')
compareRE(modList[[1]],'intFlDens_field') #heavy tails
compareRE(modList[[2]],'intHbeeVis_field') #right skew
compareRE(modList[[3]],'intLbeeVis_field')
compareRE(modList[[4]],'intPollen_field')
compareRE(modList[[4]],'intPollen_plot')
compareRE(modList[[5]],'intFlwCount_field')
compareRE(modList[[5]],'intFlwCount_plot')
compareRE(modList[[6]],'intFlwSurv_field')
compareRE(modList[[6]],'intFlwSurv_plot')
compareRE(modList[[7]],'intSeedCount_field')
# compareRE(modList[[7]],'intSeedCount_plot') #
compareRE(modList[[8]],'intSeedWeight_field')
# compareRE(modList[[8]],'intSeedWeight_plot')
compareRE(modList[[9]],'ranEffYield_field',1) #Intercepts
compareRE(modList[[9]],'ranEffYield_field',2) #Slopes
compareRE(modList[[9]],'ranEffYield_plot',1,0.3) #Intercepts
compareRE(modList[[9]],'ranEffYield_plot',2,0.3) #Slopes

# Dagitty claims list for seed fields -------------------------------------

unlist(shipley.test(seedDAG,TRUE))

# test <- stan(file = "C:\\Users\\Samuel\\Documents\\Projects\\UofC\\canola_yield_project\\Models\\Seed model claims 3\\seed_claim04.stan",
#              data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0)
# n <- names(test) #Model parameters
# n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('([sS]igma|slope)',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term) or "slope"
# n <- n[!grepl('_miss',n)] #Gets rid of imputed values
# traceplot(test,pars=n,inc_warmup=FALSE) #+ geom_hline(yintercept = 0) #Traceplots


#Get model names, and make list
modFiles <- dir(path='./Seed model claims 3',pattern = '*\\.stan',full.names = TRUE)
modFiles <- modFiles[!grepl('template',modFiles)]
modFiles <- modFiles[sapply(read.csv('./Seed model claims 3/claimsList_updated2.csv')$Filename,function(x){
    l <- grep(x,modFiles)
    if(length(l)==0) 0 else l
  })]


for(i in 14){
  overwrite <- TRUE
  if(file.exists(modFiles[i])){
    print(paste0('Starting model ',modFiles[i]))
    #Run model
    mod <- stan(file=modFiles[i],data=datalist,iter=3000,chains=4,control=list(adapt_delta=0.8),init=0)
    # mod <- stan(file=modFiles[i],data=datalist,iter=400,chains=1,control=list(adapt_delta=0.8),init=0)
    temp <- parTable(mod) #Get parameter summaries
    #Save information to csv file
    modList <- read.csv('./Seed model claims 3/claimsList_updated2.csv',sep=',',strip.white = TRUE)
    modList[i,match(names(temp),names(modList))] <- temp[grepl('claim',temp$param),]
    write.csv(modList,'./Seed model claims 3/claimsList_updated2.csv',row.names = FALSE)
    
    #Check model
    parNam <- names(mod) #Model parameters
    
    #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
    useThese <- (!(grepl('(\\[[0-9]+,*[0-9]*\\]$|^lp|\\_miss)',parNam)) |
                   grepl('[sS]igma',parNam) 
    )
    parNam <- parNam[useThese] 
    
    #Removes pollen terms
    parNam <- parNam[!parNam %in% c('intPollen','slopeVisitPol','slopeHbeeDistPollen','sigmaPolField','pollenPhi')] 
    
    print(traceplot(mod,pars=parNam,inc_warmup=FALSE)+labs(title=gsub('.*/','',modFiles[i])))
    print(mod,pars=parNam)
    print(modList[i,match(names(temp),names(modList))])
    # fastPairs(mod,pars=parNam)
    print(paste0('Model ',modFiles[i],' completed'))
    
  } else print(paste0('Model ',i,' not found'))
}

#Calculate C-statistic

#Original model
read.csv('./Seed model claims 3/claimsList_updated.csv',sep=',',strip.white = TRUE) %>% 
  mutate(pval=pnorm(-abs(Z))*2) %>% 
  shipley.dSep(.,pval,param)

read.csv('./Seed model claims 3/claimsList_updated2.csv',sep=',',strip.white = TRUE) %>% 
  mutate(pval=pnorm(-abs(Z))*2) %>% 
  shipley.dSep(.,pval,param)


# Marginal plots ---------------------

load('modSummaries_seed.Rdata') #Load existing list

#Input parameters get turned into a model matrix, multiplied, compiled
pl <- list('intPlDens'=1,
           'slopeHbeeDistPlDens'=with(datalist,seq(min(log(hbee_dist)),max(log(hbee_dist)),length=30)))

getPreds(modSummaries_seed[[1]],pl) %>% 
  ggplot(aes(x=slopeHbeeDistPlDens,y=med))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.2)+
  geom_line()



# #Something going wrong with the visitation sub-model. Trying this with glmmTMB instead
# temp <- with(datalist,data.frame(
#   totalTime = c(totalTime,totalTime_extra),
#   hbeeDist = c(hbee_dist,hbee_dist_extra),
#   lbeeDist = c(lbee_dist,lbee_dist_extra),
#   dTent = c(lbeeStocking2[,1],lbeeStocking2_extra[,1])==1,
#   dTentBees = c(lbeeStocking2[,2],lbeeStocking2_extra[,2])==1,
#   isCent = c(isCent,isCent_extra),
#   isMBay=c(isMBay,isMBay_extra),
#   is2016 = c(is2016,is2016_extra),
#   flDens=NA,
#   hbeeVis=c(hbeeVis,hbeeVis_extra),
#   lbeeVis=c(lbeeVis,lbeeVis_extra),
#   field = factor(c(plotIndex,plotIndex_extra))
# ))
# #Flower density
# temp$flDens[datalist$obsflDens_ind] <- datalist$flDens_obs
# temp$flDens[datalist$obsflDens_ind_extra+datalist$Nplot] <- datalist$flDens_obs_extra
# temp$flDens <- temp$flDens^2
# temp <- na.omit(temp)
# 
# library(glmmTMB)
# library(DHARMa)
# 
# m1 <- glmmTMB(hbeeVis ~ offset(log(totalTime)) + flDens + log(hbeeDist)*log(lbeeDist)+
#                 isMBay + isCent + lbeeVis + 
#                 (1|field), 
#               family='nbinom1',ziformula = ~1,
#               data=temp)
# 
# data.frame(pred=predict(m1),actual=temp$hbeeVis) %>% 
#   ggplot(aes(pred,actual))+
#   geom_point(position=position_jitter(width=0.1))+
#   geom_abline(intercept = 0, slope = 1)
# summary(m1)
# m1Sim <- simulateResiduals(m1,n=1000)
# testDispersion(m1Sim)
# plot(m1Sim)
# qqnorm(residuals(m1Sim))
	

# #Full model
# modPodcount_seed <- stan(file='visitation_pollen_model_seed.stan',data=datalist,
#                          iter=1000,chains=3,control=list(adapt_delta=0.8),init=inits)
# beep(1)
# # Setting max_treedepth=15 takes about 2-3x as long to run model. Use with care.

# save(modPodcount_seed,file='modPodcount_seed.Rdata') #Flw count, seed, yield model - 4hrs for 2000 iter
# save(modPodcount_seed,file='modPodcount_seed2.Rdata') #All other params (visits,plSize,etc)
# save(modPodcount_seed,file='modPodcount_seed3.Rdata') #No bee visits, no seed size/count, only flw count/survival

load('modPodcount_seed3.Rdata')
# csvFiles <- paste("./From cedar/",dir("./From cedar",pattern='[0-9].csv'),sep='')
# modPodcount_seed <- read_stan_csv(csvFiles[3])  #model from cedar

# #ML version, no variance estimates
# modPodcount_seed2 <- optimizing(stan_model(file='visitation_pollen_model_seed.stan'),data=datalist,init=inits)
# str(modPodcount_seed2)
# modPodcount_seed2$par[names(modPodcount_seed2$par) %in% pars] #Get optim results

pars <- c('intPlDens','slopeHbeeDistPlDens',#'slopeHbeeDistSqPlDens', #Planting density
          'sigmaPlDens','sigmaPlDens_field') 
pars <- c('intPlSize','slopePlDensPlSize','slopeDistPlSize','sigmaPlSize_field','sigmaPlSize','nuPlSize') #Plant size
pars <- c('intFlDens','slopePlSizeFlDens',
          'slope2016FlDens','slopeDistFlDens',
          'sigmaFlDens','sigmaFlDens_field','nuFlDens') #Flower density
pars <- c('intVisitLbee','slopeHbeeDistLbee','slopeLbeeDistLbee','slopeCentLbee','slopeMBayLbee', #Lbee vis
          'slopeStockingLbee','slope2016Lbee','slopeCentHbeeDistLbee','slopeStockingHbeeDistLbee',
          #       'slopeStockingLbeeDistLbee',
          'slopeFlDensLbee','sigmaLbeeVisField','visitLbeePhi','zeroVisLbeeTheta')
pars <- c('intVisitHbee','slopeHbeeDistHbee','slopeLbeeDistHbee','slopeLbeeHbeeDistHbee',
          'slopeLbeeVisHbee', #Hbee vis 
          'slopeCentHbee','slopeFlDensHbee','slopeMBayHbee', 
          'visitHbeePhi','zeroVisHbeeTheta') 
pars <- c('intPol','slopeHbeePol','slopeLbeePol','slopeCentPol','slopeHbeeDistPol','slopeFlDensPol', #Pollen
          'pollenPhi','sigmaPolField','sigmaPolPlot')
pars <- c('intFlwCount','slopePlSizeFlwCount', #Flower count per plant
          'slopeCentFlwCount',#'slopePolFlwCount','slopeFlDensFlwCount',
          'slopeFlwSurvFlwCount','sigmaFlwCount_field',
          'intPhiFlwCount','slopePlSizePhiFlwCount','sigmaPhiFlwCount_field')
pars <- c('intFlwSurv','slopePolSurv','slopePlSizeSurv', #Flower survival
          'slopeEdgeCentSurv','slopeHbeeDistSurv','slopeLbeeDistSurv',
          'slopeFlwDensSurv','sigmaFlwSurv_field','sigmaFlwSurv_plot',
          'intPhiFlwSurv','slopePlSizePhiFlwSurv','sigmaPhiFlwSurv_field')
pars <- c('intSeedCount','slopePolSeedCount','slopePlSizeCount', #Seeds per pod
          'slopeEdgeCentSeedCount','slopeHbeeDistSeedCount','slopeFlDensSeedCount',
          'slopeSurvSeedCount','sigmaSeedCount_field','sigmaSeedCount_plot',
          'sigmaSeedCount_plant','seedCountPhi')
pars <- c('intSeedWeight','slopePolSeedWeight','slopeSeedCount', #Weight per seed
          'slopePlSizeSeedWeight',
          'slope2016SeedWeight','slopeLbeeDistSeedWeight','slopePlDensSeedWeight',
          'slopeStockingSeedWeight','slopePlDensPlSizeSeedWeight',
          'sigmaSeedWeight','sigmaSeedWeight_field', #Random effects for plot don't converge well
          'sigmaSeedWeight_plant','lambdaSeedWeight')

# pars <- c('intYield','slopeYield','sigmaYield','sigmaYield_field','sigmaYield_plot',
#           'L_field','L_plot')
stan_hist(modPodcount_seed,pars=pars)+geom_vline(xintercept=0,linetype='dashed')
traceplot(modPodcount_seed,pars=pars,inc_warmup=F)
print(modPodcount_seed,pars=pars)

#Check model fit:
mod3 <- extract(modPodcount_seed)

# fastPairs(mod3[pars])

# # #Coefficients
# (mod3coefs <- data.frame(par='Weight per seed',parname=rownames(coefs(mod3[pars])),coefs(mod3[pars]),row.names=NULL))
# print(xtable(mod3coefs,digits=c(0,0,0,3,3,3,3,3,3,0,4)),include.rownames=F)

# #Distribution of random effects intercepts
# t(apply(mod3$intPlDens_field,2,function(x) quantile(x,c(0.5,0.025,0.975)))) %>% 
#   as.data.frame() %>% rename(med='50%',lwr='2.5%',upr='97.5%') %>% 
#   arrange(med) %>% mutate(plot=1:n()) %>% 
#   ggplot()+geom_pointrange(aes(x=plot,y=med,ymax=upr,ymin=lwr),alpha=0.5)+
#   geom_hline(yintercept=0,col='red')

#planting density - good
with(mod3,PPplots(apply(plDens_resid,1,function(x) sum(abs(x))),
                  apply(predPlDens_resid,1,function(x) sum(abs(x))),
                  datalist$plDens_obs,apply(mod3$predPlDens,2,median)[datalist$obsPlDens_ind],
                  main='Plant density'))
r2calc(mod3,'plDensMu','sigmaPlDens_field','sigmaPlDens') 

#plant size - gets way worse after getting rid of flower count term. Why is this?
with(mod3,PPplots(resid=apply(plSize_resid,1,function(x) sum(abs(x))),
                  predResid=apply(predPlSize_resid,1,function(x) sum(abs(x))),
                  actual=datalist$plantSize,pred=apply(plSizeMu,2,median),
                  main='Plant size'))
r2calc(mod3,'plSizeMu',c('sigmaPlSize_field'),'sigmaPlSize') 

#flower density - good
with(mod3,PPplots(apply(flDens_resid,1,function(x) sum(abs(x))),
                  apply(predFlDens_resid,1,function(x) sum(abs(x))),
                  apply(flDens,2,median),apply(predFlDens,2,median),
                  main='Flower density'))
r2calc(mod3,'flDensMu','sigmaFlDens_field','sigmaFlDens') 

#hbee visits - not the best, but OK
with(mod3,PPplots(apply(hbeeVis_resid,1,function(x) sum(abs(x))),
                  apply(predHbeeVis_resid,1,function(x) sum(abs(x))),
                  with(datalist,c(hbeeVis,hbeeVis_extra)),apply(mod3$predHbeeVis_all,2,median),
                  main='Hbee visits')) 

#lbee visits - not the best, but OK
# Negbin (ZI or regular) is not good, but Poisson (ZI or regular) is far worse, and traces for intercept are bad. Sticking with regular NB for now.
with(mod3,PPplots(apply(lbeeVis_resid,1,function(x) sum(abs(x))),
                  apply(predLbeeVis_resid,1,function(x) sum(abs(x))),
                  with(datalist,c(lbeeVis,lbeeVis_extra)),apply(mod3$predLbeeVis_all,2,median),
                  main='Lbee visits')) 

#pollen - good
with(mod3,PPplots(apply(pollen_resid,1,function(x) sum(abs(x))),
                  apply(predPollen_resid,1,function(x) sum(abs(x))),
                  datalist$pollenCount,apply(mod3$predPollenCount,2,median),
                  main='Pollen counts')) 
r2calc(mod3,'pollenMu',c('sigmaPolField','sigmaPolPlot'),'pollenPhi') 

# #Flower count - OK
# with(mod3,PPplots(apply(flwCount_resid,1,function(x) sum(abs(x))),
#                apply(predFlwCount_resid,1,function(x) sum(abs(x))),
#                datalist$flwCount,apply(mod3$predFlwCount,2,median),
#                main='Flower count per plant'))

#flower survival - OK
with(mod3,PPplots(resid=apply(podCount_resid,1,function(x) sum(abs(x))),
                  predResid=apply(predPodCount_resid,1,function(x) sum(abs(x))),
                  actual=datalist$podCount,pred=apply(predPodCount,2,median),main='Pods per plant'))
r2calc(mod3,'flwSurv',c('sigmaFlwSurv_field','sigmaFlwSurv_plot'),'intPhiFlwSurv') 


#Posterior predictive check plots
PPplots <- function(resid,predResid,actual,pred,main=NULL){
  par(mfrow=c(2,1))
  plot(resid~predResid,ylab='Sum actual residuals',xlab='Sum simulated residuals',main=main)
  x <- sum(resid<predResid)/length(resid)
  legend('topleft',paste('p =',round(min(x,1-x),3)))
  abline(0,1,col='red') #PP plot
  
  plot(actual~pred,type='p', #Predicted vs Actual
       xlab=paste('Predicted',main),ylab=paste('Actual',main),pch=19,col=rgb(0,0,0,0.3),xlim=range(mod3$predPodCount)) 
  
  for(i in 1:ncol(mod3$predPodCount)){
    lines(x=range(mod3$predPodCount[,i]),rep(actual[i],2),col=rgb(0,0,0,0.3))
  }
  dim(mod3$predPodCount)
  
  points(pred,actual, #Predicted vs Actual
         col='red',pch=19) 
  abline(0,1,col='red')
  abline(lm(actual~pred),col='red',lty=2)
  par(mfrow=c(1,1))
}


#seeds per pod - not good; weird distribution of seeds
with(mod3,PPplots(apply(seedCount_resid,1,function(x) sum(abs(x))),
                  apply(predSeedCount_resid,1,function(x) sum(abs(x))),
                  jitter(datalist$seedCount),jitter(apply(mod3$predSeedCount,2,median)),
                  main='Seeds per pod')) 

#weight per seed - good
with(mod3,PPplots(apply(seedMass_resid,1,function(x) sum(abs(x))),
                  apply(predSeedMass_resid,1,function(x) sum(abs(x))),
                  datalist$seedMass,apply(mod3$predSeedMass,2,median),
                  main='Weight per seed'))

# Run yearly version of main models --------------------------------------------

library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE) #Avoids recompilation
rstan_options(javascript = FALSE)
options(mc.cores = 6)

load('./datalist2015_seed.Rdata')
datalist2015 <- datalist
load('./datalist2016_seed.Rdata')
datalist2016 <- datalist; rm(datalist)

# 2015 models
modFiles <- dir(pattern = 'seed_.*\\.stan')
modList2015 <- vector(mode = 'list',length = length(modFiles))
names(modList2015) <- gsub('(seed_.*[0-9]{2}|\\.stan)','',modFiles)

modList2015[1] <- stan(file=modFiles[1],data=datalist2015,iter=4000,warmup=3000,
                       control=list(adapt_delta=0.9),chains=4,init=0.1) #Plant density, Plant size, Flower Density - OK
modList2015[2] <- stan(file=modFiles[2],data=datalist2015,iter=2000,chains=4,init=0) #Hbee visitation - OK, but traces for RE aren't great
modList2015[3] <- stan(file=modFiles[3],data=datalist2015,iter=2000,chains=4,init=0) #Lbee visitation - OK
modList2015[4] <- stan(file=modFiles[4],data=datalist2015,iter=2000,chains=4,init=0) #Pollen - OK
modList2015[5] <- stan(file=modFiles[5],data=datalist2015,iter=2000,chains=4,init=0) #Flower count per plant - OK
modList2015[6] <- stan(file=modFiles[6],data=datalist2015,iter=2000,chains=4,init=0) #Pod count (flw survival) per plant
modList2015[7] <- stan(file=modFiles[7],data=datalist2015,iter=2000,chains=4,init=0) #Seeds per pod
modList2015[8] <- stan(file=modFiles[8],data=datalist2015,iter=2000,chains=4,init=0) #Weight per seed
modList2015[9] <- stan(file=modFiles[9],data=datalist2015,iter=3000,chains=4,init=0,control=list(adapt_delta=0.9)) #Yield

# #Traceplots
# for(i in 1:length(modList2015)){
#   if(!is.null(modList2015[[i]])){
#     n <- names(modList2015[[i]]) #Model parameters
#     n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('[sS]igma',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
#     p <- traceplot(modList2015[[i]],pars=n,inc_warmup=FALSE)#+geom_hline(yintercept = 0) #Traceplots
#     print(p); rm(p)
#     a <- readline('Press Return to continue: ')
#     if(a!='') break()
#   }
# }

# Get model summaries into a list of tables
# modSummaries2015_seed <- vector(mode = 'list',length = length(modList2015)) #Create new empty list
# names(modSummaries2015_seed) <- names(modList2015)
load('modSummaries2015_seed.Rdata') #Load existing list
#Update model summaries if needed
temp <- lapply(modList2015,parTable) #Get parameter summaries
for(i in 1:length(temp)){
  if(class(temp[[i]])=='data.frame'||length(temp[[i]])>1){ #If temp is not empty (model not run)
    modSummaries2015_seed[[i]]$summary <- temp[[i]]$summary #Overwrite
    modSummaries2015_seed[[i]]$covMat <- temp[[i]]$covMat
  }
}
parNames <- lapply(modSummaries2015_seed,function(x) x$summary$param) #Clean up extra parameter names
for(i in 2:length(modSummaries2015_seed)){
  if(!is.null(modSummaries2015_seed[[i]])){
    chooseThese <- !parNames[[i]] %in% unlist(parNames[1:(i-1)])
    modSummaries2015_seed[[i]]$summary <- modSummaries2015_seed[[i]]$summary[chooseThese,]
    modSummaries2015_seed[[i]]$covMat <- modSummaries2015_seed[[i]]$covMat[chooseThese,chooseThese]
  }
}
save(modSummaries2015_seed,file = 'modSummaries2015_seed.Rdata'); rm(temp)

#2016 models
datalist2016$missPlDens_ind <- 0

modFiles <- dir(pattern = 'seed_.*\\.stan')
modList2016 <- vector(mode = 'list',length = length(modFiles))
names(modList2016) <- gsub('(seed_.*[0-9]{2}|\\.stan)','',modFiles)

modList2016[1] <- stan(file=modFiles[1],data=datalist2016,iter=4000,warmup=3000,
                       control=list(adapt_delta=0.9),chains=4,init=0.1) #Plant density, Plant size, Flower Density - OK
modList2016[2] <- stan(file=modFiles[2],data=datalist2016,iter=2000,chains=4,init=0) #Hbee visitation - OK, but traces for RE aren't great
modList2016[3] <- stan(file=modFiles[3],data=datalist2016,iter=2000,chains=4,init=0) #Lbee visitation - OK
modList2016[4] <- stan(file=modFiles[4],data=datalist2016,iter=2000,chains=4,init=0) #Pollen - OK
modList2016[5] <- stan(file=modFiles[5],data=datalist2016,iter=2000,chains=4,init=0) #Flower count per plant - OK
modList2016[6] <- stan(file=modFiles[6],data=datalist2016,iter=2000,chains=4,init=0) #Pod count (flw survival) per plant
modList2016[7] <- stan(file=modFiles[7],data=datalist2016,iter=2000,chains=4,init=0) #Seeds per pod
modList2016[8] <- stan(file=modFiles[8],data=datalist2016,iter=2000,chains=4,init=0) #Weight per seed
modList2016[9] <- stan(file=modFiles[9],data=datalist2016,iter=3000,chains=4,init=0,control=list(adapt_delta=0.9)) #Yield

#Traceplots
for(i in 1:length(modList2016)){
  if(!is.null(modList2016[[i]])){
    n <- names(modList2016[[i]]) #Model parameters
    n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('[sS]igma',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term)
    p <- traceplot(modList2016[[i]],pars=n,inc_warmup=FALSE)#+geom_hline(yintercept = 0) #Traceplots
    print(p); rm(p)
    a <- readline('Press Return to continue: ')
    if(a!='') break()
  }
}

# Get model summaries into a list of tables
# modSummaries2016_seed <- vector(mode = 'list',length = length(modList2016)) #Create new empty list
# names(modSummaries2016_seed) <- names(modList2016)
load('modSummaries2016_seed.Rdata') #Load existing list
#Update model summaries if needed
temp <- lapply(modList2016,parTable) #Get parameter summaries
for(i in 1:length(temp)){
  if(class(temp[[i]])=='data.frame'||length(temp[[i]])>1){ #If temp is not empty (model not run)
    modSummaries2016_seed[[i]]$summary <- temp[[i]]$summary #Overwrite
    modSummaries2016_seed[[i]]$covMat <- temp[[i]]$covMat
  }
}
parNames <- lapply(modSummaries2016_seed,function(x) x$summary$param) #Clean up extra parameter names
for(i in 2:length(modSummaries2016_seed)){
  if(!is.null(modSummaries2016_seed[[i]])){
    chooseThese <- !parNames[[i]] %in% unlist(parNames[1:(i-1)])
    modSummaries2016_seed[[i]]$summary <- modSummaries2016_seed[[i]]$summary[chooseThese,]
    modSummaries2016_seed[[i]]$covMat <- modSummaries2016_seed[[i]]$covMat[chooseThese,chooseThese]
  }
}
save(modSummaries2016_seed,file = 'modSummaries2016_seed.Rdata'); rm(temp)
