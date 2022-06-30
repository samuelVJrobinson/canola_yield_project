#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN SEED CANOLA FIELDS (2014+2015)

# Libraries, functions, and seed DAG  ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(beepr)
library(xtable)

library(rstan)
setwd('./Models')
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

# setwd('~/Projects/UofC/canola_yield_project') #Multivac path
# setwd('~/Documents/canola_yield_project') #Galpern machine path

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
                  pollen ~ hbeeVis  + lbeeVis + cent + hbeeDist + flDens,
                  flwCount ~ plSize + cent + flwSurv,
                  flwSurv ~ pollen + plSize + cent + hbeeDist + lbeeDist + flDens,
                  seedCount ~ pollen + plSize + cent + hbeeDist + flDens + flwSurv,
                  seedWeight ~ pollen + seedCount + plSize + plDens + lbeeDist,
                  coords= list(x = setNames(nodeCoords$x,nodeCoords$name),
                               y = setNames(nodeCoords$y,nodeCoords$name)),
                  labels=setNames(nodeCoords$labs,nodeCoords$name)
)

# Load in data ----------------------

load('./datalist_seed.Rdata')

# #Seed field data
# load("./Seed field analysis/seedfieldDataAll.RData")
# fieldsAllSeed <- data.frame(allFields); plantsAllSeed <- data.frame(allPlants)
# pollenAllSeed <- data.frame(allPollen); seedsAllSeed <- data.frame(allSeeds);
# surveyAllSeed <- data.frame(allSurvey); 
# plantsAllSeed$Field <- gsub('Unrah','Unruh',plantsAllSeed$Field) #Fixes spelling error
# seedsAllSeed$Field <- gsub('Unrah','Unruh',seedsAllSeed$Field)
# seedsAllSeed$EdgeCent <- ifelse(seedsAllSeed$EdgeCent=='Cent','Center',seedsAllSeed$EdgeCent) #Fixes cent/center
# rm(allFields,allPlants,allPollen,allSeeds,allSurvey,behav2015,visitors2016,nectar2016,folder)
# #Set 'negative' missing pods (mistake in counting) to NA.
# plantsAllSeed <- mutate(plantsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
# seedsAllSeed <- mutate(seedsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
# #Extra seed field data (from Riley W.)
# rileyExtra <- read.csv('./Seed field analysis/rileyExtra.csv')
# setwd('~/Projects/UofC/canola_yield_project')
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
# 
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
# 
# datalist <- c(datalistField,datalistPlot,datalistFlw,datalistPlant)
# 
# if(any(sapply(datalist,function(x) sum(is.na(x)))!=0)){ beep(1); print("NAs found in datalist")}
# 
# rm(datalistField,datalistPlot,datalistFlw,datalistPlant,rileyFields,samFields) #Cleanup
# str(datalist)
# 
# save(datalist,file = './Models/datalist_seed.Rdata')


# Run models ----------------------------

# library(shinystan)

#Models split into separate stan files (faster)
modFiles <- dir(pattern = 'seed_.*\\.stan')
modList <- vector(mode = 'list',length = length(modFiles))
names(modList) <- gsub('(seed_.*[0-9]{2}|\\.stan)','',modFiles)

modList[1] <- stan(file=modFiles[1],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0.1) #Plant density, Plant size, Flower Density - OK
modList[2] <- stan(file=modFiles[2],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #Hbee visitation - OK
modList[3] <- stan(file=modFiles[3],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #Lbee visitation - OK
modList[4] <- stan(file=modFiles[4],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #Pollen - OK
modList[5] <- stan(file=modFiles[5],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #Flower count per plant - OK
modList[6] <- stan(file=modFiles[6],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #Pod count (flw survival) per plant
modList[7] <- stan(file=modFiles[7],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #Seeds per pod
modList[8] <- stan(file=modFiles[8],data=datalist,iter=2000,chains=4,control=list(adapt_delta=0.8),init=0) #Weight per seed
modList[9] <- stan(file=modFiles[9],data=datalist,iter=3000,chains=4,control=list(adapt_delta=0.8),init=0) #Yield
beepr::beep(1)

#Traceplots
for(i in 1:length(modList)){
  if(!is.null(modList[[i]])){
    n <- names(modList[[i]]) #Model parameters
    n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('([sS]igma|slope)',n)] #Gets rid of parameter vectors, unless it contains "sigma" (variance term) or "slope"
    n <- n[!grepl('_miss',n)] #Gets rid of imputed values
    
    p <- traceplot(modList[[i]],pars=n,inc_warmup=FALSE) #+ geom_hline(yintercept = 0) #Traceplots
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
    modSummaries_seed[[i]] <- temp[[i]] #Overwrite
  }
}
parNames <- lapply(modSummaries_seed,function(x) x$param) #Clean up extra parameter names
for(i in 2:length(modSummaries_seed)){
  if(!is.null(modSummaries_seed[[i]])){
    modSummaries_seed[[i]] <- modSummaries_seed[[i]][!parNames[[i]] %in% unlist(parNames[1:(i-1)]),]
  }
}
save(modSummaries_seed,file = 'modSummaries_seed.Rdata'); rm(temp)

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
compareRE(modList[[7]],'intSeedCount_plot') #
compareRE(modList[[8]],'intSeedWeight_field')
compareRE(modList[[8]],'intSeedWeight_plot')
compareRE(modList[[9]],'ranEffYield_field',1) #Intercepts
compareRE(modList[[9]],'ranEffYield_field',2) #Slopes
compareRE(modList[[9]],'ranEffYield_plot',1,0.3) #Intercepts
compareRE(modList[[9]],'ranEffYield_plot',2,0.3) #Slopes

#Posterior predictive checks
PPplots(modList[[1]],datalist$plDens_obs,c('predPlDens','plDens_resid','predPlDens_resid'),
        index = datalist$obsPlDens_ind,main='Plant Density') #OK, but bimodal peaks
PPplots(modList[[1]],datalist$flDens_obs,c('predFlDens','flDens_resid','predFlDens_resid'),
        index=datalist$obsflDens_ind,main='Flower density') #OK, but skewed distribution. Not improved by MBay or Year terms
PPplots(modList[[1]],datalist$plantSize,c('predPlSize','plSize_resid','predPlSize_resid'),
        'Plant size') #OK
PPplots(modList[[2]],c(datalist$hbeeVis,datalist$hbeeVis_extra),
        c('predHbeeVis_all','hbeeVis_resid','predHbeeVis_resid'),'Honeybee visits') #Not good
PPplots(modList[[3]],c(datalist$lbeeVis,datalist$lbeeVis_extra),
        c('predLbeeVis_all','lbeeVis_resid','predLbeeVis_resid'),'Leafcutter visits') #OK, but overpredicts at low values

PPplots(modList[[4]],datalist$pollenCount,c('predPollenCount','pollen_resid','predPollen_resid'),
        'Pollen') #OK
PPplots(modList[[5]],datalist$flwCount,c('predFlwCount','flwCount_resid','predFlwCount_resid'),
        'Flowers per plant') #OK
PPplots(modList[[6]],datalist$podCount,c('predPodCount','podCount_resid','predPodCount_resid'),'Pods per plant') #OK
PPplots(modList[[7]],datalist$seedCount_obs,c('predSeedCount','seedCount_resid','predSeedCount_resid'),
        index=datalist$obsSeedCount_ind,'Seeds per pod') #OK, but distribution is weird, even with exp-normal. 
PPplots(modList[[8]],datalist$seedMass_obs,c('predSeedMass','seedMass_resid','predSeedMass_resid'),
        index=datalist$obsSeedMass_ind,'Seed Mass')
PPplots(modList[[9]],log(datalist$yield),c('predYield','yield_resid','predYield_resid'),'Seed mass per plant')

# Path diagram ------------------------------------------------------------

load('modSummaries_seed.Rdata')

pathCoefs <- modSummaries_seed %>% #Get path coefficients
  bind_rows(.id='to') %>% 
  transmute(to,name=param,Z,pval) %>% 
  filter(to!='yield',!grepl('(int|sigma|lambda|Phi|phi|rho|nu|theta)',name)) %>% 
  mutate(to=case_when(
    to=='flDens' & grepl('PlDens$',name) ~ 'plDens',
    to=='flDens' & grepl('PlSize$',name) ~ 'plSize',
    TRUE ~ gsub('avg','seed',to)
  )) %>% 
  mutate(name=gsub('slope','',name)) %>% 
  filter(mapply(grepl,capFirst(to),name)) %>% 
  mutate(name=mapply(gsub,capFirst(to),'',name)) %>% 
  mutate(name=capFirst(name,TRUE)) 

tidySeedDAG <- seedDAG %>% #Create tidy dagitty set
  tidy_dagitty() 

# ggplot(tidySeedDAG,aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_dag_edges() +
#   geom_dag_text(col='black') +
#   theme_dag_blank()

tidySeedDAG$data <- tidySeedDAG$data %>% #Match coefs to dagitty set
  rename(xstart=x,ystart=y) %>% 
  full_join(x=pathCoefs,y=.,by=c('to','name')) %>% 
  mutate(isNeg=Z<0,isSig=pval<0.05) %>% 
  mutate(L=sqrt(abs(Z)),1) %>% 
  mutate(edgeLab=ifelse(isSig,as.character(sign(Z)*round(L,1)),'')) %>% 
  mutate(C=ifelse(isNeg,'red','black')) %>% 
  mutate(A=ifelse(isSig,1,0.1))

ggplot(tidySeedDAG$data,aes(x = xstart, y = ystart, xend = xend, yend = yend))+
  annotate('text',x=0.5,y=4.5+0.1,label='Plot Level',size=5)+
  annotate('rect',xmin=-0.5,ymin=0,xmax=1.5,ymax=4.5,fill=NA,col='black',
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
modFiles <- modFiles[sapply(read.csv('./Seed model claims 3/claimsList_updated.csv')$Filename,function(x){
    l <- grep(x,modFiles)
    if(length(l)==0) 0 else l
  })]


for(i in 13){
  overwrite <- TRUE
  if(file.exists(modFiles[i])){
    print(paste0('Starting model ',modFiles[i]))
    #Run model
    mod <- stan(file=modFiles[i],data=datalist,iter=3000,chains=4,control=list(adapt_delta=0.8),init=0)
    # mod <- stan(file=modFiles[i],data=datalist,iter=500,chains=1,control=list(adapt_delta=0.8),init=0)
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
    
    print(traceplot(mod,pars=parNam,inc_warmup=FALSE)+labs(title=gsub('.*/','',modFiles[i])))
    print(mod,pars=parNam)
    print(modList[i,match(names(temp),names(modList))])
    # fastPairs(mod,pars=parNam)
    print(paste0('Model ',modFiles[i],' completed'))
    
  } else print(paste0('Model ',i,' not found'))
}



#Marginal plots ---------------------


#Input parameters get turned into a model matrix, multiplied, compiled
pl <- list('intPlDens'=1,
           'slopeHbeeDistPlDens'=with(datalist,seq(min(log(hbee_dist)),max(log(hbee_dist)),length=10)))

getPreds(modList[[1]],pl) %>% 
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

# Partial effects plots for seed fields -----------------------------------
load('modPodcount_seed3.Rdata')
mod3 <- extract(modPodcount_seed)
rm(modPodcount_seed); gc()

#Plot-level data (647 rows) - visitation data
plotDat <- with(datalist,data.frame(
  plot=1:(Nplot+Nplot_extra),
  hbeeDist=c(hbee_dist,hbee_dist_extra),
  lbeeDist=c(lbee_dist,lbee_dist_extra), #Centered on 3.113714
  lbeeVis=c(lbeeVis,lbeeVis_extra),isMBay=c(isMBay,isMBay_extra), #Centered on 4.641845
  flDens=apply(mod3$flDens,2,mean),isCent=c(isCent,isCent_extra), #Sqrt tranformed - 22
  halfStock=c(lbeeStocking,lbeeStocking_extra), #half-stocking
  year=c(is2016,is2016_extra),
  logLbeeVis=log(1+c(lbeeVis,lbeeVis_extra)/c(totalTime,totalTime_extra)),
  logHbeeVis=log(1+c(hbeeVis,hbeeVis_extra)/c(totalTime,totalTime_extra)),
  hbeeResid=log((c(hbeeVis,hbeeVis_extra)+1)/c(totalTime,totalTime_extra))-apply(mod3$visitMu_hbee,2,median), #Log-residuals
  lbeeResid=log((c(lbeeVis,lbeeVis_extra)+1)/c(totalTime,totalTime_extra))-apply(mod3$visitMu_lbee,2,median)
)) %>% 
  mutate(logHbeeDist=log(hbeeDist)-mean(log(hbeeDist)),logLbeeDist=log(lbeeDist)-mean(log(lbeeDist)))

#Model matrix for lbees
MM_lbeeVis <- with(plotDat,data.frame(int=1,logLbeeDist,logHbeeDist,isMBay,isCent,halfStock,year,flDens,
                                      centHbeeDist=isCent*logHbeeDist,stockHbeeDist=halfStock*logHbeeDist)) %>% as.matrix()

#Coefficent matrix for lbee visits
coef_lbeeVis <- with(mod3,data.frame(intVisitLbee,slopeLbeeDistLbee,slopeHbeeDistLbee,slopeMBayLbee,slopeCentLbee,                                     slopeStockingLbee=slopeStockingLbee[,2],slope2016Lbee,slopeFlDensLbee,slopeCentHbeeDistLbee,
                                     slopeStockingHbeeDistLbee=slopeStockingHbeeDistLbee[,2])) %>% as.matrix()

#Partial effect of lbee distance
MM_temp <- MM_lbeeVis %>% as.data.frame() %>% 
  mutate_at(vars(-logLbeeDist),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_lbeeVis),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
  mutate(resid=plotDat$hbeeResid+pred) %>%  
  mutate_at(vars(pred:resid),exp) %>% 
  mutate(lbeeDist=exp(logLbeeDist+3.113714)) %>% 
  ggplot(aes(x=lbeeDist))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  # geom_jitter(aes(y=resid),alpha=0.5,size=1,width=0.5)+
  geom_point(aes(y=resid),alpha=0.5,size=1)+
  lims(x=c(0,50))+
  labs(x='Distance from shelter',y='Visits/10 mins')
ggsave('../Figures/Seed/slopeLbeeDistLbee.png',p1,width=8,height=4)

#Partial effect of hbee distance/stocking
MM_temp <- MM_lbeeVis %>% as.data.frame() %>% 
  mutate_at(vars(-logHbeeDist,-halfStock,-stockHbeeDist),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_lbeeVis),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
  mutate(resid=plotDat$hbeeResid+pred) %>% 
  mutate(halfStock=factor(halfStock,labels=c('Full','Half'))) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=exp(logHbeeDist+4.641845)))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=halfStock),alpha=0.3,show.legend=F)+
  geom_line(aes(y=pred,col=halfStock),size=1)+
  geom_jitter(aes(y=resid,col=halfStock),alpha=0.5,size=1,width=5)+
  xlim(0,400)+ylim(0,50)+
  labs(x='Distance from edge',y='Visits/10 mins',col='Stocking')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  theme(legend.position=c(0.8,0.85),legend.background=element_rect(fill='white',colour='black',linetype='solid'),
        legend.title=element_text(size=18))
ggsave('../Figures/Seed/slopeStockingHbeeDistLbee.png',p1,width=8,height=4.5)

#Partial effect of edge/center of bay
MM_temp <- MM_lbeeVis %>% as.data.frame() %>% 
  mutate_at(vars(-isCent),mean) %>% as.matrix()

##Plot of bay centre effect on lbees
# p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_lbeeVis),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
#   rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
#   mutate(resid=plotDat$hbeeResid+pred) %>% 
#   mutate(edgeCent=factor(isCent,labels=c('Edge','Center'))) %>% 
#   mutate_at(vars(pred:resid),exp) %>% 
#   ggplot(aes(x=edgeCent))+
#   geom_jitter(aes(y=resid),alpha=0.5,size=1,width=0.2,col='forestgreen')+
#   geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),size=1)+
#   lims(y=c(0,50))+labs(x='Bay position',y='Visits/10 mins',title='Leafcutters')
# ggsave('../Figures/Seed/slopeEdgeCentLbee.png',p1,width=4,height=4)

#Data for combined lbee-hbee plot
lbeeEdgeCent <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_lbeeVis),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred) %>%
  mutate(edgeCent=factor(isCent,labels=c('Edge','Center'))) %>%
  mutate_at(vars(pred:resid),exp) %>% mutate(Spp='lbee')

#Model matrix for hbees 
MM_hbeeVis <- with(plotDat,data.frame(int=1,flDens,logHbeeDist,logLbeeDist,lbeeHbeeDist=logHbeeDist*logLbeeDist,lbeeVis=logLbeeVis,isMBay,isCent)) %>% as.matrix()

#Coefficent matrix for hbee visits
coef_hbeeVis <- with(mod3,data.frame(intVisitHbee,slopeFlDensHbee,slopeHbeeDistHbee,
                                     slopeLbeeDistHbee,slopeLbeeHbeeDistHbee,slopeLbeeVisHbee,slopeMBayHbee,
                                     slopeCentHbee,zeroVisHbeeTheta)) %>% as.matrix()

#Partial effect of hbee distance
MM_temp <- MM_hbeeVis %>% as.data.frame() %>% 
  mutate_at(vars(-logHbeeDist),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_hbeeVis[,-which(colnames(coef_hbeeVis)=='zeroVisHbeeTheta')]),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred) %>% 
  # mutate_at(vars(pred:resid),function(x) exp(x)*(1-mean(coef_hbeeVis[,'zeroVisHbeeTheta']))) %>%
  mutate_at(vars(pred:resid),exp) %>%
  ggplot(aes(x=exp(logHbeeDist+4.641845)))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,show.legend=F)+
  geom_line(aes(y=pred),size=1)+
  geom_jitter(aes(y=resid),alpha=0.5,size=1,width=5)+
  lims(x=c(0,400),y=c(0,100))+
  labs(x='Distance from edge',y='Visits/10 mins')
ggsave('../Figures/Seed/slopeHbeeDistHbee.png',p1,width=8,height=4.5)

#Partial effect of edge/center
MM_temp <- MM_hbeeVis %>% as.data.frame() %>% 
  mutate_at(vars(-isCent),mean) %>% as.matrix()

# p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_hbeeVis[,-which(colnames(coef_hbeeVis)=='zeroVisHbeeTheta')]),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
#   rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
#   mutate(resid=plotDat$hbeeResid+pred,edgeCent=factor(isCent,labels=c('Edge','Center'))) %>% 
#   mutate_at(vars(pred:resid),function(x) exp(x)*(1-mean(coef_hbeeVis[,'zeroVisHbeeTheta']))) %>%
#   ggplot(aes(x=edgeCent))+
#   geom_jitter(aes(y=resid),alpha=0.5,size=1,width=0.2,col='darkorange')+
#   geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),size=1)+
#   lims(y=c(0,50))+
#   labs(x='Bay position',y='Visits/10 mins',title='Honeybees')
# ggsave('../Figures/Seed/slopeEdgeCentHbee.png',p1,width=4,height=4)

#Data for combined lbee-hbee plot
hbeeEdgeCent <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_hbeeVis[,-which(colnames(coef_hbeeVis)=='zeroVisHbeeTheta')]),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred,edgeCent=factor(isCent,labels=c('Edge','Center'))) %>%
  mutate_at(vars(pred:resid),function(x) exp(x)*(1-mean(coef_hbeeVis[,'zeroVisHbeeTheta']))) %>% 
  mutate(Spp='hbee')

#Combined lbee-hbee plot, showing effect of bay center
p1 <- bind_rows(hbeeEdgeCent,lbeeEdgeCent) %>% 
  mutate(Spp=ifelse(Spp=='hbee','Honey bee','Leafcutter')) %>% 
  mutate(edgeCent=ifelse(edgeCent=='Center','Centre','Edge')) %>% 
  mutate(edgeCent=factor(edgeCent,levels=c('Edge','Centre'))) %>% 
  ggplot(aes(x=edgeCent,col=Spp))+
  geom_point(aes(y=resid),alpha=0.5,size=1,position=position_jitterdodge(jitter.width=0.4,dodge.width=0.75))+
  geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),position=position_dodge(width=0.75,preserve='total'),size=1,show.legend=F)+
  lims(y=c(0,50))+
  labs(x='Bay position',y='Visits/10 mins',col='Visitor')+
  scale_colour_manual(values=c('darkorange','forestgreen'))+
  theme(legend.position=c(0.88,0.88),
        legend.background=element_rect(fill='white',colour='black',linetype='solid'),
        legend.title=element_text(size=18))
ggsave('../Figures/Seed/slopeEdgeCent_both.png',p1,width=8,height=6)


#Partial effect of lbee distance on both hbees and lbees
#Model matrix for lbees
MM_lbeeVis <- with(plotDat,data.frame(int=1,logLbeeDist,logHbeeDist,isMBay,isCent,halfStock,year,flDens,
                                      centHbeeDist=isCent*logHbeeDist,stockHbeeDist=halfStock*logHbeeDist)) %>% as.matrix()

#Coef matrix for lbees
coef_lbeeVis <- with(mod3,data.frame(intVisitLbee,slopeLbeeDistLbee,slopeHbeeDistLbee,slopeMBayLbee,slopeCentLbee,
                                     slopeStockingLbee,slope2016Lbee,slopeFlDensLbee,slopeCentHbeeDistLbee,
                                     slopeStockingHbeeDistLbee)) %>% as.matrix()

#Partial effect of lbee distance
MM_temp_lbees <- MM_lbeeVis %>% as.data.frame() %>% 
  mutate_at(vars(-logLbeeDist),mean) %>% as.matrix()

lbeePartial <- data.frame(MM_temp_lbees,t(apply(MM_temp_lbees %*% t(coef_lbeeVis),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
  mutate(resid=plotDat$lbeeResid+pred) %>%  
  mutate_at(vars(pred:resid),exp) %>% 
  mutate(lbeeDist=exp(logLbeeDist+3.113714)) %>% mutate(Type='Leafcutter')

#Model matrix for hbees
MM_hbeeVis <- with(plotDat,data.frame(int=1,flDens,logHbeeDist,logLbeeDist,lbeeHbeeDist=logHbeeDist*logLbeeDist,lbeeVis=logLbeeVis,isMBay,isCent)) %>% as.matrix()

#Coefficent matrix for hbee visits
coef_hbeeVis <- with(mod3,data.frame(intVisitHbee,slopeFlDensHbee,slopeHbeeDistHbee,
                                     slopeLbeeDistHbee,slopeLbeeHbeeDistHbee,slopeLbeeVisHbee,slopeMBayHbee,
                                     slopeCentHbee,zeroVisHbeeTheta)) %>% as.matrix()

MM_temp_hbees <- MM_hbeeVis %>% as.data.frame() %>% 
  mutate_at(vars(-logLbeeDist),mean) %>% as.matrix()

hbeePartial <- data.frame(MM_temp_hbees,t(apply(MM_temp_hbees %*% t(coef_hbeeVis[,-which(colnames(coef_hbeeVis)=='zeroVisHbeeTheta')]),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred) %>% 
  mutate_at(vars(pred:resid),function(x) exp(x)*(1-mean(coef_hbeeVis[,'zeroVisHbeeTheta']))) %>% 
  mutate(lbeeDist=exp(logLbeeDist+3.113714)) %>% mutate(Type='Honeybee')

p1 <- bind_rows(lbeePartial,hbeePartial) %>% 
  ggplot(aes(x=lbeeDist))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=Type),alpha=0.3,show.legend=F)+
  geom_line(aes(y=pred,col=Type),size=1)+
  geom_jitter(aes(y=resid,col=Type),alpha=0.5,size=1,width=0.5)+
  lims(x=c(0,50),y=c(0,100))+
  labs(x='Distance from shelter',y='Visits/10 mins',col='Visitor')+
  scale_colour_manual(values=c('darkorange','forestgreen'))+scale_fill_manual(values=c('darkorange','forestgreen'))+
  theme(legend.position=c(0.8,0.85),legend.background=element_rect(fill='white',colour='black',linetype='solid'),
        legend.title=element_text(size=18))
ggsave('../Figures/Seed/slopeLbeeDistBothBees.png',p1,width=8,height=4.5)

#Flower-level data (pollen)
polDat <- with(datalist,data.frame(
  flw=1:Nflw,
  hbeeVis=log(1+hbeeVis/totalTime)[flowerIndex], #hbee Visitation
  lbeeVis=log(1+lbeeVis/totalTime)[flowerIndex], #lbee visitation
  edgeCent=isCent[flowerIndex], #Edge
  hbeeDist=(log(hbee_dist)-4.641845)[flowerIndex], #Hbee dist, centered on 4.641845
  flDens=apply(mod3$flDens,2,median)[flowerIndex],
  pol_resid=apply(mod3$pollenMu,2,median)-log(pollenCount+1) #Log-residuals
))

#Model matrix
MM_pollen <- with(polDat,data.frame(int=1,hbeeVis,lbeeVis,edgeCent,hbeeDist,flDens)) %>% 
  as.matrix()

#Coef matrix for pollen
coef_pollen <- with(mod3,data.frame(intPol,slopeHbeePol,slopeLbeePol,slopeCentPol,slopeHbeeDistPol,slopeFlDensPol)) %>% 
  as.matrix()

#Partial effect of hbee visits
MM_temp <- MM_pollen %>% as.data.frame() %>% 
  mutate_at(vars(-hbeeVis),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_pollen),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=polDat$pol_resid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=exp(hbeeVis))) +
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='darkorange',show.legend=F)+
  geom_line(aes(y=pred),col='darkorange',size=1)+
  geom_jitter(aes(y=resid),alpha=0.5,size=1,width=0,col='darkorange')+
  lims(y=c(0,75),x=c(0,50))+
  labs(x='Visits/10 mins',y='Pollen grains/stigma',title='Honeybee')
ggsave('../Figures/Seed/slopeHbeeVisPol.png',p1,width=4,height=4)

#Partial effect of lbee visits
MM_temp <- MM_pollen %>% as.data.frame() %>% 
  mutate_at(vars(-lbeeVis),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_pollen),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=polDat$pol_resid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>%
  ggplot(aes(x=exp(lbeeVis))) +
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='forestgreen',show.legend=F)+
  geom_line(aes(y=pred),col='forestgreen',size=1)+
  geom_jitter(aes(y=resid),alpha=0.5,size=1,width=0,col='forestgreen')+
  lims(y=c(0,75),x=c(0,50))+
  labs(x='Visits/10 mins',y='Pollen grains/stigma',title='Leafcutter')
ggsave('../Figures/Seed/slopeLbeeVisPol.png',p1,width=4,height=4)

#Partial effect of distance and bay center
MM_temp <- MM_pollen %>% as.data.frame() %>% 
  mutate_at(vars(-hbeeDist,-edgeCent),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_pollen),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=polDat$pol_resid+pred,edgeCent=factor(edgeCent,labels=c('Edge','Center'))) %>% 
  mutate_at(vars(pred:resid),exp) %>%
  ggplot(aes(x=exp(hbeeDist+4.641845)))+
  # geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr,col=edgeCent),size=1,alpha=0.3,show.legend=F)+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=edgeCent),alpha=0.3,show.legend=F)+
  geom_jitter(aes(y=resid,col=edgeCent),alpha=0.5,size=1,width=5)+
  geom_line(aes(y=pred,col=edgeCent),size=1)+
  lims(y=c(0,75),x=c(0,400))+
  labs(x='Distance from edge',y='Pollen grains/stigma',col='Bay position')+
  scale_colour_manual(values=c('blue','red'))+scale_fill_manual(values=c('blue','red'))+
  theme(legend.position=c(0.8,0.85),legend.background=element_rect(fill='white',colour='black',linetype='solid'),legend.title=element_text(size=18))
ggsave('../Figures/Seed/slopeDistCentPol.png',p1,width=8,height=4.5)


#Plant-level data (582 rows)
plantDat <- with(datalist,data.frame(
  plant=1:Nplant,
  plDens=apply(mod3$plDens,2,median)[plotIndex_F[plantIndex]], #Centered around 3.573633
  plSize=plantSize, #Centered around 3.197572
  hbeeDist=(log(hbee_dist)-mean(log(c(hbee_dist,hbee_dist_extra))))[plotIndex_F[plantIndex[1:Nplant]]], #Centered around 4.641845
  hbeeDist_trans=hbee_dist[plotIndex_F[plantIndex[1:Nplant]]], #Back-transformed hbee distance
  lbeeDist=(log(lbee_dist)-mean(log(c(lbee_dist,lbee_dist_extra))))[plotIndex_F[plantIndex[1:Nplant]]], #Centered around 3.113714
  edgeCent=isCent[plotIndex_F[plantIndex[1:Nplant]]],
  pol=apply(mod3$pollenMu_plot,2,median)[plotIndex_F[plantIndex]], #Centered on 2.523904 (mean(mod3$intPol))
  lbeeVis=log(1+hbeeVis/totalTime)[plotIndex_F[plantIndex[1:Nplant]]],
  hbeeVis=log(1+lbeeVis/totalTime)[plotIndex_F[plantIndex[1:Nplant]]],
  flDens=apply(mod3$flDens,2,median)[plotIndex_F[plantIndex[1:Nplant]]],
  flwCount=flwCount,podCount=podCount,
  plSize_resid=apply(mod3$plSize_resid,2,median),
  flwCount_resid=apply(mod3$flwCountMu,2,median)-log(flwCount), #Log-residuals for flower count
  surv=podCount/flwCount, #Propotion Survival
  centSurvLogit=logit(podCount/flwCount)-mean(logit(podCount/flwCount)),
  # surv_resid= podCount/flwCount-invLogit(apply(mod3$flwSurv,2,median))
  surv_resid= logit(podCount/flwCount)-apply(mod3$flwSurv,2,median) #Logit-resid for flower surv
))

#Model matrix for plant size
MM_plSize <- with(plantDat,data.frame(int=1,plDens,hbeeDist)) %>% as.matrix()
#Coefficient matrix for plant size
coef_plSize <- with(mod3,data.frame(intPlSize,slopePlDensPlSize,slopeDistPlSize)) %>% as.matrix()

#Partial effect of density on plant size
MM_temp <- MM_plSize %>% as.data.frame() %>% mutate_at(vars(-plDens),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_plSize),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$plSize_resid+pred,plDens=exp(plDens+3.573633)) %>% arrange(plDens) %>% 
  mutate_at(vars(pred:resid),function(x) exp(x+3.197572)) %>%
  ggplot(aes(x=plDens))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='black')+
  geom_point(aes(y=resid),alpha=0.5)+#ylim(250,650)+
  geom_line(aes(y=pred),size=1,col='black')+
  labs(x=expression(paste('Plants per ',m^2)),y='Plant size (g)',title='Seed')
ggsave('../Figures/Seed/slopePlDensPlSize.png',p1,width=6,height=6)

#Model matrix for flower count 
MM_flwCount <- with(plantDat,data.frame(int=1,plSize,edgeCent,centSurvLogit)) %>% 
  as.matrix()
#Coefficent matrix for flower count
coef_flwCount <- with(mod3,data.frame(intFlwCount,slopePlSizeFlwCount,slopeCentFlwCount,slopeFlwSurvFlwCount)) %>% as.matrix()

#Partial effect of plant size
MM_temp <- MM_flwCount %>% as.data.frame() %>% mutate_at(vars(-plSize),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$flwCount_resid+pred,plSize=exp(plSize+3.204388)) %>% arrange(plSize) %>% 
  mutate_at(vars(pred:resid),exp) %>%
  filter(!(resid>1900&plSize<50&plSize>20)) %>% #Get rid of weird outlier
  ggplot(aes(x=plSize))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3,fill='red')+
  geom_point(aes(y=resid),alpha=0.5)+#ylim(250,650)+
  geom_line(aes(y=pred),size=1,col='red')+
  labs(x='Plant size (g)',y='Flowers per plant',title='Seed')
ggsave('../Figures/Seed/slopePlSizeFlwCount.png',p1,width=6,height=6)


#Partial effect of bay center
MM_temp <- MM_flwCount %>% as.data.frame() %>% mutate_at(vars(-edgeCent),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$flwCount_resid+pred) %>% arrange(edgeCent) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  mutate(edgeCent=factor(edgeCent,labels=c('Bay Edge','Bay Center'))) %>% 
  ggplot(aes(x=edgeCent))+
  geom_jitter(aes(y=resid),alpha=0.5)+
  geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),col='red',size=1)+
  ylim(250,650)+
  labs(x=NULL,y='Flowers per plant')

#Model matrix for flower survivorship
MM_flwSurv <- with(plantDat,data.frame(int=1,pol,plSize,edgeCent,hbeeDist,lbeeDist,flDens)) %>% 
  as.matrix()
#Coefficent matrix for survivorship
coef_flwSurv <- with(mod3,data.frame(intFlwSurv,slopePolSurv,slopePlSizeSurv,
                                     slopeEdgeCentSurv,slopeHbeeDistSurv,slopeLbeeDistSurv,slopeFlwDensSurv)) %>% as.matrix()

#Partial effect of plant size
MM_temp <- MM_flwSurv %>% as.data.frame() %>% 
  mutate_at(vars(-plSize),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$surv_resid+pred) %>% 
  mutate_at(vars(pred:resid),function(x) invLogit(x)*100) %>% 
  ggplot(aes(x=exp(plSize+3.197572)))+
  geom_jitter(aes(y=resid),alpha=0.5,width=0)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  labs(x='Plant size (g)',y='Flower survival (%)',title='Seed')+
  lims(y=c(45,90))
ggsave('../Figures/Seed/slopePlsizeSurv.png',p1,width=6,height=6)

#Partial effect of pollen amount
MM_temp <- MM_flwSurv %>% as.data.frame() %>% 
  mutate_at(vars(-pol),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$surv_resid+pred) %>% 
  mutate_at(vars(pred:resid),function(x) invLogit(x)*100) %>% 
  ggplot(aes(x=exp(pol+2.523904)))+
  geom_jitter(aes(y=resid),alpha=0.5,width=0)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  labs(x='Pollen/stigma',y='Pod survival (%)')
ggsave('../Figures/Seed/slopePolSurv.png',p1,width=4,height=4)

#Partial effect of shelter distance
MM_temp <- MM_flwSurv %>% as.data.frame() %>% 
  mutate_at(vars(-lbeeDist),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$surv_resid+pred) %>% 
  mutate_at(vars(pred:resid),function(x) invLogit(x)*100) %>% 
  ggplot(aes(x=exp(lbeeDist+3.113714)))+
  geom_jitter(aes(y=resid),alpha=0.5,width=0)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  labs(x='Distance from Shelter',y='Pod survival (%)')
ggsave('../Figures/Seed/slopeLbeeDistSurv.png',p1,width=4,height=4)

#Partial effect of hbee distance
MM_temp <- MM_flwSurv %>% as.data.frame() %>% 
  mutate_at(vars(-hbeeDist,-edgeCent),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$surv_resid+pred,edgeCent=factor(edgeCent,labels=c('Edge','Center'))) %>% 
  mutate_at(vars(pred:resid),function(x) invLogit(x)*100) %>% 
  ggplot(aes(x=exp(hbeeDist+4.641845)))+
  geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr,col=edgeCent),size=1,alpha=0.3,show.legend=F)+
  geom_jitter(aes(y=resid,col=edgeCent),alpha=0.5,width=)+
  geom_line(aes(y=pred,col=edgeCent),size=1)+
  labs(x='Distance from edge of field',y='Pod survival (%)',col='Bay position')+
  scale_colour_manual(values=c('blue','red'))+scale_fill_manual(values=c('blue','red'))+
  theme(legend.position=c(0.75,0.8),legend.background=element_rect(fill='white',colour='black',linetype='solid'),legend.title=element_text(size=15))
ggsave('../Figures/Seed/slopeHbeeDistSurv.png',p1,width=4,height=4)


#Pod-level data (2885 rows)
podDat <- with(datalist,data.frame(
  field=plotIndex[plantIndex[podIndex[1:Npod]]],plot=plantIndex[podIndex[1:Npod]],plant=podIndex[1:Npod],
  pol=tapply(apply(mod3$pollenMu-matrix(rep(mod3$intPol,1050),ncol=1050),2,median),flowerIndex,mean)[
    plotIndex_F[plantIndex[podIndex]]],
  plSize=plantSize[podIndex], #Plant size - centered on 3.197572
  surv=(logit(podCount/flwCount)-mean(logit(podCount/flwCount)))[podIndex], #Pod survival - centered on 0.7384344
  edgeCent=isCent[plantIndex[podIndex]], 
  lbeeDist=(log(c(lbee_dist,lbee_dist_extra))-mean(log(c(lbee_dist,lbee_dist_extra))))[plantIndex[podIndex]],
  hbeeDist=(log(c(hbee_dist,hbee_dist_extra))-mean(log(c(hbee_dist,hbee_dist_extra))))[plantIndex[podIndex]],
  flDens=apply(mod3$flDens,2,mean)[plantIndex[podIndex]], 
  plDens=apply(mod3$plDens,2,mean)[plantIndex[podIndex]], #centered on 3.573633
  stocking=lbeeStocking[plantIndex[podIndex]],surv=logit(podCount/flwCount)[plantIndex[podIndex]],
  year2016=is2016[plantIndex[podIndex]],seedCount=seedCount,seedWeight=seedMass,
  seedCountResid=log(seedCount)-apply(mod3$seedCountMu,2,median), #Median log-resid for seedCount
  seedMassResid=apply(mod3$seedMass_resid,2,median) #Median resid for seedMass
))

#Model matrix for seed count
MM_seedCount <- with(podDat,data.frame(int=1,pol,plSize,edgeCent,lbeeDist,flDens,surv))
#Coef matrix for seed count
coef_seedCount <- with(mod3,data.frame(intSeedCount,slopePolSeedCount,slopePlSizeCount,
                                       slopeEdgeCentSeedCount,slopeHbeeDistSeedCount,slopeFlDensSeedCount,slopeSurvSeedCount)) %>% as.matrix()

#Partial effect of year
MM_temp <- MM_seedCount %>% as.data.frame() %>% mutate_at(vars(-year2016),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedCountResid+pred,year2016=factor(year2016,labels=c('2015','2016'))) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=year2016))+
  geom_jitter(aes(y=resid),alpha=0.3,width=0.3)+
  geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),col='red')+
  labs(x='Year',y='Seeds per pod')


#Partial effects of plant size
MM_temp <- MM_seedCount %>% as.data.frame() %>% 
  mutate_at(vars(-plSize),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedCount),1,function(x) 
  quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedCountResid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=exp(plSize+3.197572)))+
  geom_jitter(aes(y=resid),alpha=0.3,width=1)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_line(aes(y=pred),col='red',size=1)+
  lims(y=c(0,40))+
  labs(x='Plant size (g)',y='Seeds per pod')
ggsave('../Figures/Seed/slopePlSizeSeedCount.png',p1,width=4,height=4)

#Partial effects of pollen
MM_temp <- MM_seedCount %>% as.data.frame() %>% 
  mutate_at(vars(-pol),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedCount),1,function(x) 
  quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedCountResid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=exp(pol+2.523904)))+
  geom_jitter(aes(y=resid),alpha=0.3,width=1)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_line(aes(y=pred),col='red',size=1)+
  coord_cartesian(ylim=c(0,40))+
  labs(x='Pollen/stigma',y='Seeds per pod')
ggsave('../Figures/Seed/slopePolSeedCount.png',p1,width=4,height=4)

#Partial effects of survival
MM_temp <- MM_seedCount %>% as.data.frame() %>% 
  mutate_at(vars(-surv),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedCount),1,function(x) 
  quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedCountResid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=invLogit(surv+0.7384344)*100))+
  geom_jitter(aes(y=resid),alpha=0.3,width=0)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_line(aes(y=pred),col='red',size=1)+
  lims(y=c(0,40),x=c(25,100))+
  labs(x='Pod survival (%)',y='Seeds per pod')
ggsave('../Figures/Seed/slopeSurvSeedCount.png',p1,width=4,height=4)

#Partial effects of edge/cent
MM_temp <- MM_seedCount %>% as.data.frame() %>% 
  mutate_at(vars(-edgeCent),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedCount),1,function(x) 
  quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedCountResid+pred,edgeCent=factor(edgeCent,labels=c('Edge','Center'))) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=edgeCent))+
  geom_jitter(aes(y=resid),alpha=0.3,width=0.2)+
  geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),col='red',size=1)+
  lims(y=c(0,40))+labs(x='Bay Position',y='Seeds per pod')
ggsave('../Figures/Seed/slopeCentSeedCount.png',p1,width=4,height=4)

#Model matrix for seedMass 
MM_seedMass <- with(podDat,data.frame(int=1,pol,seedCount,plSize,plDens,year2016,
                                      lbeeDist,stocking,plDensPlSize=plSize*plDens)) %>% 
  as.matrix()


#Coefficent matrix for seedMass
coef_seedMass <- with(mod3,data.frame(intSeedWeight,slopePolSeedWeight,slopeSeedCount,
                                      slopePlSizeSeedWeight,slopePlDensSeedWeight,slope2016SeedWeight,
                                      slopeLbeeDistSeedWeight,slopeStockingSeedWeight,slopePlDensPlSizeSeedWeight)) %>% 
  as.matrix()

#Partial effect of year
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-year2016),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedCountResid+pred,year2016=factor(year2016,labels=c('2015','2016'))) %>% 
  ggplot(aes(x=year2016))+
  geom_jitter(aes(y=resid),alpha=0.3,width=0.3)+
  geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),col='red')+
  labs(x='Year',y='Weight per seed (mg)')

#Partial effect of seed count
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-seedCount),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred) %>% #arrange(seedCount) %>% 
  ggplot()+
  geom_point(aes(seedCount,resid),alpha=0.3)+
  geom_ribbon(aes(x=seedCount,ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_line(aes(seedCount,pred),col='red',size=1)+
  labs(x='Seeds per pod',y='1000 seed weight (g)')
ggsave('../Figures/Seed/slopePodCountPodWeight.png',p1,width=4,height=4)

#Partial effect of plant size
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-plSize),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) 
  quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred) %>% arrange(plSize) %>% 
  ggplot(aes(x=exp(plSize+3.197572)))+
  geom_point(aes(y=resid),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_line(aes(y=pred),size=1,col='red')+
  ylim(0,7)+labs(x='Plant size (g)',y='1000 seed weight (g)')
ggsave('../Figures/Seed/slopePlSizePodWeight.png',p1,width=4,height=4)

#Partial effect of plant density
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-plDens),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred) %>% arrange(plSize) %>% 
  ggplot(aes(x=exp(plDens+3.573633)))+
  geom_point(aes(y=resid),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_line(aes(y=pred),size=1,col='red')+
  ylim(0,7)+labs(x=expression(paste('Plants/',m^2)),y='1000 seed weight (g)')
ggsave('../Figures/Seed/slopePlDensPodWeight.png',p1,width=4,height=4)

#Partial plot of half-stocking
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-stocking),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) 
  quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred,stocking=factor(stocking,labels=c('Full','Half'))) %>%
  ggplot(aes(x=stocking))+geom_jitter(aes(y=resid),alpha=0.3,width=0.2)+
  geom_pointrange(aes(y=pred,ymax=upr,ymin=lwr),size=1,col='red')+
  ylim(0,7)+labs(x='Tent stocking',y='1000 seed weight (g)')
ggsave('../Figures/Seed/slopeStockingPodWeight.png',p1,width=4,height=4)


#Partial effect of seed count - seed size relationship, with 3 levels of plant size (0.1,0.5,0.9 percentiles)
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-seedCount,-plSize),mean) %>% 
  mutate(plSize=cut(plSize,breaks=c(-5,-0.9,0.9,5))) %>% 
  mutate(plSize=as.numeric(as.character(factor(plSize,labels=c('-0.9','0','0.9'))))) %>% 
  as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred,plSize=factor(plSize,labels=c('Low','Med','High'))) %>% #arrange(seedCount) %>% 
  ggplot()+
  geom_point(aes(seedCount,resid),col='black',alpha=0.2)+
  geom_ribbon(aes(x=seedCount,ymax=upr,ymin=lwr,fill=plSize),alpha=0.3)+
  geom_line(aes(seedCount,pred,col=plSize),size=1)+
  labs(x='Seeds per pod',y='1000 seed weight (g)')+
  scale_fill_manual(values=c('blue','purple','red'))+
  scale_colour_manual(values=c('blue','purple','red'))+
  labs(col='Plant Size',fill='Plant Size',main='Seed canola')+
  theme(legend.position=c(0.9,0.82),legend.background=element_rect(fill='white',colour='black',linetype='solid'),
        legend.title=element_text(size=12),legend.text=element_text(size=12))
ggsave('../Figures/Seed/slopePodCountPodWeightPlSize.png',p1,width=6,height=4)


# Total effects for seed fields -------------------------------------------

#Load all coefficients from models
setwd('./Models')
library(rstan)
modPodcount_seed <- read_stan_csv("./From cedar/visitation_pollen_model_seed4.csv")  #model from cedar
mod3 <- extract(modPodcount_seed)
# #Fixes Fbay-Mbay names
# names(mod3)[grepl('FBay',names(mod3))] <- sub('FBay','MBay',names(mod3)[grepl('FBay',names(mod3))])

# transdens <- log(plDens)-3.573633 #Plant density
# transPlSize <- log(plSize)-3.197572 #Plant size
# transHbeeDist <- log(hbeeDist)-4.641845 #Hbee distance
# transLbeeDist <- log(lbeeDist)-3.113714 #Lbee distance
# transFlDens <- sqrt(flDens)-22 #Flw dens
# transLbeeVis <- log(1+lbeeVis)
# transSurv <- logit(surv) - 0.7384344 #Survival

#Checks argument lengths and returns correct sample size
getN <- function(args){
  N <- sapply(args,length) #Lengths of arguments
  #Stops function if arguments mismatch
  if(length(unique(N[N>1]))>1) stop('Arguments have mismatching lengths') 
  return(max(N))
}

#Generate plant densities from distances
genPlDens <- function(hbeeDist,int,slopeDist,stDev,transform=F){
  N <- getN(list(hbeeDist))
  if(stDev<0) stop('Negative sigma')
  #dist= distances, intercept,slope of distance, stDev, and back-transform
  transHbeeDist <- log(hbeeDist)-4.641845
  mu <- int+slopeDist*transHbeeDist
  a <- rnorm(N,mu,stDev)
  if(any(is.na(a))) stop('NAs from rnorm')
  if(transform) return(exp(a+3.197572)) else return(a)
}

#Generate plant sizes from plant density/distance
genPlSize <- function(plDens,hbeeDist,int,slopePlDens,slopeDist,stDev,transform=F){
  N <- getN(list(plDens,hbeeDist))
  if(stDev<0) stop('Negative sigma')
  transHbeeDist <- log(hbeeDist)-4.641845
  transPlDens <- log(plDens)-3.573633
  mu <- int+slopePlDens*transPlDens+slopeDist*transHbeeDist
  a <- rnorm(plDens,mu,stDev)
  if(any(is.na(a))) stop('NAs from rnorm')
  if(transform) return(exp(a+3.197572)) else return(a)
}

#Generate flower density for plot - pl size must be centered/transformed before
genFlwDens <- function(avgPlSize,isMbay=0,is2016,hbeeDist,int,slopePlSize,slopeMbay,slope2016,
                       slopeDist,stDev,transform=F){
  N <- getN(list(avgPlSize,isMbay,is2016,hbeeDist))
  if(stDev<0) stop('Negative sigma')
  # transPlSize <- log(avgPlSize)-3.197572 #Plant size - already transformed, so no need
  transHbeeDist <- log(hbeeDist)-4.641845
  mu <- int + avgPlSize*slopePlSize + isMbay*slopeMbay + is2016*slope2016 + transHbeeDist*slopeDist
  a <- rnorm(N,mu,stDev)
  if(any(is.na(a))) stop('NAs from rnorm')
  if(transform) return((a+22)^2) else return(a)
}

#Generate lbee visits (per 10 mins)
genLbeeVis <- function(lbeeDist,hbeeDist,isMbay=0,isCent,isHalfStock,is2016,flDens,
                       int,slopeLbeeDist,slopeHbeeDist,slopeMbay,slopeCent,slopeStock,
                       slope2016,slopeFlDens,slopeCentHbeeDist,slopeStockHbeeDist,phi,addVar=T){
  N <- getN(list(lbeeDist,hbeeDist,isMbay,isCent,isHalfStock,is2016,flDens))
  if(phi<0) stop('Negative phi')
  transHbeeDist <- log(hbeeDist)-4.641845 #Hbee distance
  transLbeeDist <- log(lbeeDist)-3.113714 #Lbee distance
  transFlDens <- sqrt(flDens)-22 #Flw dens
  lambda <- int+slopeHbeeDist*transHbeeDist+slopeLbeeDist*transLbeeDist+    
    slopeMbay*isMbay+slopeCent*isCent+slopeStock*isHalfStock+ slope2016*is2016+
    slopeFlDens*transFlDens+slopeCentHbeeDist*isCent*transHbeeDist+
    slopeStockHbeeDist*isHalfStock*transHbeeDist
  if(addVar){
    a <- rnbinom(N,mu=exp(lambda),size=phi)
    return(a)
  } else return(exp(lambda))
  
}

#Generate hbee visits (per 10 mins)
genHbeeVis <- function(flDens,hbeeDist,lbeeDist,lbeeVis,isMbay=0,isCent,int,
                       slopeFlDens,slopeHbeeDist,slopeLbeeDist,slopeHbeeLbeeDist,
                       slopeLbeeVis,slopeMbay,slopeCent,phiNB,thetaZI,addVar=T){
  N <- getN(list(flDens,hbeeDist,lbeeDist,lbeeVis,isMbay,isCent))
  if(phiNB<0) stop('Negative phiNB')
  if(thetaZI<0|thetaZI>1) stop('ThetaZI not between 0 and 1')
  if(any(!is.finite(lbeeDist)|lbeeDist<0|!is.finite(hbeeDist)|hbeeDist<0)) stop('Hbee dist or Lbee dist either negative, infinite, or NA')
  
  transHbeeDist <- log(hbeeDist)-4.641845 #Hbee distance
  transLbeeDist <- log(lbeeDist)-3.113714 #Lbee distance
  transLbeeVis <- log(1+lbeeVis) #Lbee visits
  transFlDens <- sqrt(flDens)-22 #Flw dens
  
  lambda <- int+slopeFlDens*transFlDens+slopeHbeeDist*transHbeeDist+
    slopeLbeeDist*transLbeeDist+slopeHbeeLbeeDist*transHbeeDist*transLbeeDist+
    slopeLbeeVis*transLbeeVis+slopeMbay*isMbay+slopeCent*isCent
  if(addVar){
    a <- rep(NA,N)
    a[rbinom(N,1,thetaZI)==1] <- 0 #ZI process generates extra zeros
    a[is.na(a)] <- rnbinom(sum(is.na(a)),mu=exp(lambda)[is.na(a)],size=phiNB) #Other slots filled by NB process
    return(a)
  } else {
    a <- exp(lambda)*(1-thetaZI)
    return(a)
  }
}

#Generate avg pollen counts - use stDev from plot level
genAvgPol <- function(hbeeVis,lbeeVis,isCent,hbeeDist,flDens,
                      int,slopeHbeeVis,slopeLbeeVis,slopeCent,slopeHbeeDist,slopeFlDens,
                      stDev,center=T,transform=F){
  N <- getN(list(hbeeVis,lbeeVis,isCent,hbeeDist,flDens))
  if(stDev<0) stop('Negative sigma')
  #hbee counts, distance, intercept, slopes, center data, transform
  #Note: transform doesn't make sense unless center=F (Jensen's inequality)
  transHbeeVis <- log(hbeeVis+1)
  transLbeeVis <- log(lbeeVis+1)
  transHbeeDist <- log(hbeeDist)-4.641845
  transFlDens <- sqrt(flDens)-22
  mu <- ifelse(center,0,int)+slopeHbeeVis*transHbeeVis+slopeLbeeVis*transLbeeVis+
    slopeCent*isCent+slopeHbeeDist*transHbeeDist+slopeFlDens*transFlDens
  a <- rnorm(N,mu,stDev)
  if(any(is.na(a))) stop('NAs from rnorm. stDev:',stDev,' slopeHbeeVis:',slopeHbeeVis,
                         ' slopeLbeeVis:',slopeLbeeVis,' slopeCent:',slopeCent,' slopeHbeeDist:',slopeHbeeDist,
                         ' slopeFlDens:',slopeFlDens,' center:',center,' int:',int,'\nmeanHbeeVis:',mean(hbeeVis),
                         ' meanLbeeVis:',mean(lbeeVis),' meanIsCent:',mean(isCent),' meanHbeeDist:',mean(hbeeDist),
                         ' meanFlDens:',mean(flDens),' transform:',transform,
                         '\nmeanTransHbee:',mean(transHbeeVis),' meanTransLbeeVis:',mean(transLbeeVis),' meanTransHbee:',
                         mean(transHbeeDist),' meanTransFlDens:',mean(transFlDens))
  if(transform) return(exp(a)) else return(a) 
}

#Generate flower survival (proportion from beta-binomial) 
genPodCount <- function(pol,plSize,isCent,hbeeDist,lbeeDist,flwDens,
                        int,slopePol,slopePlSize,slopeCent,slopeHbeeDist,slopeLbeeDist,slopeFlwDens,
                        intPhi,slopePlSizePhi,addVar=T){
  #noBeta removes beta-binomial step
  N <- getN(list(pol,plSize,isCent,hbeeDist,lbeeDist,flwDens))
  transPlSize <- log(plSize)-3.197572 #Plant size
  transHbeeDist <- log(hbeeDist)-4.641845 #Hbee distance
  transLbeeDist <- log(lbeeDist)-3.113714 #Lbee distance
  transFlDens <- sqrt(flwDens)-22 #Flw dens
  #Expected value
  mu <- invLogit(int+slopePol*pol+slopePlSize*transPlSize+slopeCent*isCent+
                   slopeHbeeDist*transHbeeDist+slopeLbeeDist*transLbeeDist+
                   slopeFlwDens*transFlDens)
  if(!addVar) return(mu) #Without dispersion
  phi <- exp(intPhi+slopePlSizePhi*transPlSize) #Dispersion
  theta <- rbeta(N,mu*phi,(1-mu)*phi) #Adds beta
  return(theta)
}

#Generate flower counts from plant sizes, edgeCent, and survival rate
genFlwCount <- function(plSize,isCent,surv,
                        int,slopePlSize,slopeCent,slopeSurv,
                        intPhi,slopePlSizePhi,addVar=T){
  N <- getN(list(plSize,isCent,surv))
  
  #plSize=plant size, intercept,slope of plant size, phi
  transPlSize <- log(plSize)-3.197572 #Plant size
  transSurv <- logit(surv) - 0.7384344 #Survival
  mu <- int+slopePlSize*transPlSize+slopeCent*isCent+
    slopeSurv*transSurv #Expected value
  if(!addVar) return(exp(mu)) #No variance
  phi <- exp(intPhi+slopePlSizePhi*transPlSize) #Dispersion
  a <- round(exp(rnorm(N,mu,phi)))
  if(any(is.na(a))) stop('NAs from rnorm')
  return(a)
}

#Generate avg seed counts per plant (seed count) - use stDev from plant level
genSeedCount <- function(pol,plSize,isCent,hbeeDist,flDens,flwSurv,
                         int,slopePol,slopePlSize,slopeCent,slopeHbeeDist,slopeFlDens,slopeSurv,
                         stDev,transform=T){
  N <- getN(list(pol,plSize,isCent,hbeeDist,flDens,flwSurv))
  if(stDev<0) stop('Negative sigma')
  transPlSize <- log(plSize)-3.197572 #Plant size
  transHbeeDist <- log(hbeeDist)-4.641845 #Hbee distance
  transFlDens <- sqrt(flDens)-22 #Flw dens
  transFlSurv <- logit(flwSurv) - 0.7384344 #Transformed flw survival
  
  #Expected value
  mu <- int+slopePol*pol+transPlSize*slopePlSize+slopeCent*isCent+slopeHbeeDist*transHbeeDist+
    slopeFlDens*transFlDens+slopeSurv*transFlSurv
  a <- rnorm(N,mu,stDev)
  if(any(is.na(a))) stop('NAs from rnorm')
  if(transform) return(exp(a)) else return(a)
}

#Generate avg seed size per plant (seed weight, mg) - use stDev from plant level
genSeedWeight <- function(pol,avgSeedCount,plSize,plDens,is2016,lbeeDist,isHalfStock,
                          int,slopePol,slopeSeedCount,slopePlSize,slopePlDens,slope2016,
                          slopeLbeeDist,slopeStocking,slopePlDensPlSize,lambda,stDev,addVar=T){
  N <- getN(list(pol,avgSeedCount,plSize,plDens,is2016,lbeeDist,isHalfStock))
  if(stDev<0) stop('Negative sigma')
  transPlDens <- log(plDens)-3.573633 #Plant density
  transPlSize <- log(plSize)-3.197572 #Plant size
  transLbeeDist <- log(lbeeDist)-3.113714 #Lbee distance
  
  mu <- int + slopePol*pol+slopeSeedCount*avgSeedCount+slopePlSize*transPlSize+
    slopePlDens*transPlDens+slope2016*is2016+slopeLbeeDist*transLbeeDist+slopeStocking*isHalfStock+
    slopePlDensPlSize*transPlDens*transPlSize
  if(!addVar) {
    a <- mu+(1/lambda)
    a[a<0] <- 0.22 #Set negative seed weights to 0.22 (lowest observed)
    return(a)
  } else {
    a <- rnorm(N,mu,stDev)+rexp(N,lambda)
    if(any(is.na(a))) stop('NAs from rnorm')
    a[a<0] <- 0.22 
    return(a)
  }
}

#Generate weight of all seeds at a given plot
genYield <- function(pods,seedCount,seedWeight,int,slopeYield,plotLev=T){
  calcYield <- pods*seedCount*seedWeight/1000 #Seed weight per plant (g)
  a <- exp(int+log(calcYield)*slopeYield)
  if(plotLev) return(sum(a)) else return(a) #Return plot-level yield (default)
}

#Simulation
simSeed <- function(hbeeDist,lbeeDist,isCent,is2016,isHalfStock,dat,returnAll=F,useMean=F,plotVar=T){
  #returnAll: return all simulation parameters, or just yield?
  #useMean: use mean of distribution, or add plot-level variance?
  
  #Test values
  hbeeDist=seq(1,401,50)
  lbeeDist=10
  isCent=0
  is2016=0
  isHalfStock=0
  dat=mod3
  returnAll=F
  useMean=F
  plotVar=T
  rm(hbeeDist,lbeeDist,isCent,is2016,returnAll,useMean,isHalfStock) #Cleanup
  rm(simPlDens,simPlSize,simFlwDens,simLbeeVis,simHbeeVis,simAvgPol,simFlwCount,
     simPodCount,simPodSurv,simSeedCount,simSeedWeight,simYield)
  
  if(useMean) dat <- lapply(dat,mean)
  
  #Simulated plant density
  simPlDens <- round(with(dat,genPlDens(hbeeDist=hbeeDist,int=sample(intPlDens,1),
                                        slopeDist=sample(slopeHbeeDistPlDens,1),
                                        stDev=ifelse(plotVar,sample(sigmaPlDens,1),0),transform=T)))
  
  # # Compare to actual - looks OK, but a bit too sparse at the high end. Mainly caused by random effect
  # par(mfrow=c(2,1))
  # hist(surveyAllSeed$PlDens,breaks=seq(0,150,5),xlab='Actual Plant Density',main=NULL,freq=F)
  # hist(simPlDens,breaks=seq(0,150,5),xlab='Sim Plant Density',main=NULL,freq=F)
  
  #Simulated plant sizes within plots
  simPlSize <- with(dat,mapply(genPlSize,plDens=simPlDens,hbeeDist=hbeeDist,
                               int=sample(intPlDens,1),slopePlDens=sample(slopePlDensPlSize,1),
                               slopeDist=sample(slopeDistPlSize,1),
                               stDev=ifelse(plotVar,sample(sigmaPlSize,1),0),transform=T,SIMPLIFY=F))
  
  # #Compare to actual - looks OK
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simPlSize),plantsAllSeed$VegMass,na.rm=T)+10
  # hist(plantsAllSeed$VegMass,xlab='Actual Plant Size',main=NULL,breaks=seq(0,upr,10))
  # hist(unlist(simPlSize),xlab='Sim Plant Size',main=NULL,breaks=seq(0,upr,10))
  
  #Simulate flower density 
  simFlwDens <- with(dat,genFlwDens(avgPlSize=sapply(simPlSize,function(x) mean(log(x)-3.197572)),
                                    is2016=is2016,hbeeDist=hbeeDist,int=sample(intFlDens,1),
                                    slopePlSize=sample(slopePlSizeFlDens,1),slopeMbay=sample(slopeMBayFlDens,1),
                                    slope2016=sample(slope2016FlDens,1),slopeDist=sample(slopeDistFlDens,1),
                                    stDev=ifelse(plotVar,sample(sigmaFlDens,1),0),transform=T))
  # #Compare to actual - OK
  # par(mfrow=c(2,1))
  # upr <- max(surveyAllSeed$FlDens*4,simFlwDens,na.rm=T)+50
  # hist(surveyAllSeed$FlDens*4,breaks=seq(0,upr,50),xlab='Actual Flw Dens',main=NULL)
  # hist(simFlwDens,breaks=seq(0,upr,50),xlab='Sim Flw Dens',main=NULL)
  
  #Simulate lbee visits
  simLbeeVis <- with(dat,genLbeeVis(lbeeDist=lbeeDist,hbeeDist=hbeeDist,isCent=isCent,
                                    isHalfStock=isHalfStock,is2016=is2016,flDens=simFlwDens,
                                    int=sample(intVisitLbee,1),slopeLbeeDist=sample(slopeLbeeDistHbee,1),
                                    slopeHbeeDist=sample(slopeHbeeDistLbee,1),slopeMbay=sample(slopeMBayLbee,1),
                                    slopeCent=sample(slopeCentLbee,1),slopeStock=sample(slopeStockingLbee,1),
                                    slope2016=sample(slope2016Lbee,1),slopeFlDens=sample(slopeFlDensLbee,1),
                                    slopeCentHbeeDist=sample(slopeCentHbeeDistLbee,1),
                                    slopeStockHbeeDist=sample(slopeStockingHbeeDistLbee,1),phi=sample(visitLbeePhi,1),
                                    addVar=plotVar))
  
  # #Looks a bit spread out
  # par(mfrow=c(2,1))
  # upr <- max(c(with(surveyAllSeed,lbee/(TotalTime/10)),simLbeeVis),na.rm=T)
  # with(surveyAllSeed,hist(lbee/(TotalTime/10),breaks=seq(0,upr+5,5),xlab='Actual Visits/10mins',main=NULL))
  # hist(simLbeeVis,breaks=seq(0,upr+5,5),xlab='Sim Visits/10mins',main=NULL)
  
  #Simulate hbee visits 
  simHbeeVis <- with(dat,genHbeeVis(flDens=simFlwDens,hbeeDist=hbeeDist,lbeeDist=lbeeDist,lbeeVis=simLbeeVis,
                                    isCent=isCent,int=sample(intVisitHbee,1),slopeFlDens=sample(slopeFlDensHbee,1),
                                    slopeHbeeDist=sample(slopeHbeeDistHbee,1),slopeLbeeDist=sample(slopeLbeeDistHbee,1),
                                    slopeHbeeLbeeDist=sample(slopeLbeeHbeeDistHbee,1),slopeLbeeVis=sample(slopeLbeeVisHbee,1),
                                    slopeMbay=sample(slopeMBayHbee,1),slopeCent=sample(slopeCentHbee,1),
                                    phiNB=sample(visitHbeePhi,1),thetaZI=sample(zeroVisHbeeTheta,1),addVar=plotVar))
  
  # #Looks OK
  # par(mfrow=c(2,1))
  # upr=max(c(with(surveyAllSeed,hbee/(TotalTime/10)),simHbeeVis),na.rm=T)
  # with(surveyAllSeed,hist(hbee/(TotalTime/10),breaks=seq(0,upr+2,2),xlab='Actual Visits/10mins',main=NULL))
  # hist(simHbeeVis,breaks=seq(0,upr+2,2),xlab='Sim Visits/10mins',main=NULL)
  
  #Simulate (average) pollen deposition
  simAvgPol <- with(dat,genAvgPol(hbeeVis=simHbeeVis,lbeeVis=simLbeeVis,isCent=isCent,hbeeDist=hbeeDist,
                                  flDens=simFlwDens,int=sample(intPol,1),slopeHbeeVis=sample(slopeHbeePol,1),
                                  slopeLbeeVis=sample(slopeLbeePol,1),slopeCent=sample(slopeCentPol,1),
                                  slopeHbeeDist=sample(slopeHbeeDistPol,1),slopeFlDens=sample(slopeFlDensPol,1),
                                  stDev=ifelse(plotVar,sample(sigmaPolPlot,1),0),center=T,transform=F))
  
  # #Compare to actual - looks OK
  # par(mfrow=c(2,1))
  # hist(with(pollenAllSeed,tapply(Pollen,paste(Field,Year,Distance),mean)),
  #      xlab='Actual Avg Pollen/plot',main=NULL,breaks=seq(0,200,5))
  # hist(with(dat,genAvgPol(hbeeVis=simHbeeVis,lbeeVis=simLbeeVis,isCent=isCent,hbeeDist=hbeeDist,
  #           flDens=simFlwDens,int=sample(intPol,1),slopeHbeeVis=sample(slopeHbeePol,1),
  #           slopeLbeeVis=sample(slopeLbeePol,1),
  #           slopeCent=sample(slopeCentPol,1),slopeHbeeDist=sample(slopeHbeeDistPol,1),
  #           slopeFlDens=sample(slopeFlDensPol,1),
  #           stDev=sample(sigmaPolPlot,1),center=F,transform=T)),xlab='Sim Avg pollen/plot',
  #      main=NULL,breaks=seq(0,200,5))
  
  
  #Simulate flower survival to pod
  simPodSurv <- with(dat,mapply(genPodCount,pol=simAvgPol,plSize=simPlSize,isCent=isCent,
                                hbeeDist=hbeeDist,lbeeDist=lbeeDist,flwDens=simFlwDens,int=sample(intFlwSurv,1),
                                slopePol=sample(slopePolSurv,1),slopePlSize=sample(slopePlSizeSurv,1),
                                slopeCent=sample(slopeEdgeCentSurv,1),slopeHbeeDist=sample(slopeHbeeDistSurv,1),
                                slopeLbeeDist=sample(slopeLbeeDistSurv,1),slopeFlwDens=sample(slopeFlwDensSurv,1),
                                intPhi=sample(intPhiFlwSurv,1),slopePlSizePhi=sample(slopePlSizePhiFlwSurv,1),
                                addVar=plotVar,SIMPLIFY=F))
  
  # #Flower survivorship (proportion) - looks OK
  # par(mfrow=c(2,1))
  # hist(unlist(simPodSurv),xlab='Sim flower survival',main=NULL,breaks=seq(0,1,0.05))
  # hist(with(plantsAllSeed,Pods/(Pods+Missing)),xlab='Actual flower survival',main=NULL,breaks=seq(0,1,0.05))
  
  
  #Simulate flower count per plant
  simFlwCount <- with(dat,mapply(genFlwCount,plSize=simPlSize,isCent=isCent,surv=simPodSurv,
                                 int=sample(intFlwCount,1),slopePlSize=sample(slopePlSizeFlwCount,1),
                                 slopeCent=sample(slopeCentFlwCount,1),slopeSurv=sample(slopeFlwSurvFlwCount,1),
                                 intPhi=sample(intPhiFlwCount,1),slopePlSizePhi=sample(slopePlSizePhiFlwCount,1),
                                 addVar=plotVar,SIMPLIFY=F))
  
  # #Compare to actual - looks good
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simFlwCount),with(plantsAllSeed,Pods+Missing),na.rm=T)
  # hist(unlist(simFlwCount),xlab='Sim Flowers/plant',main=NULL,breaks=seq(0,upr+50,50))
  # hist(with(plantsAllSeed,Pods+Missing),xlab='Actual Flowers/plant',main=NULL,breaks=seq(0,upr+50,50))
  
  #Generate pod count
  simPodCount <- mapply(function(x,y) round(x*y),x=simPodSurv,y=simFlwCount,SIMPLIFY=F) 
  
  # #Compare to actual - looks good
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simPodSurv)*unlist(simFlwCount),with(plantsAllSeed,Pods),na.rm=T)
  # hist(round(unlist(simPodSurv)*unlist(simFlwCount)),xlab='Sim Pods/plant',main=NULL,breaks=seq(0,upr+50,50))
  # hist(with(plantsAllSeed,Pods),xlab='Actual Pods/plant',main=NULL,breaks=seq(0,upr+50,50))
  
  #Simulate avg seeds per pod
  simSeedCount <- with(dat,mapply(genSeedCount,pol=simAvgPol,plSize=simPlSize,isCent=isCent,hbeeDist=hbeeDist,
                                  flDens=simFlwDens,flwSurv=simPodSurv,int=sample(intSeedCount,1),
                                  slopePol=sample(slopePolSeedCount,1),slopePlSize=sample(slopePlSizeCount,1),
                                  slopeCent=sample(slopeEdgeCentSeedCount,1),slopeHbeeDist=sample(slopeHbeeDistSeedCount,1),
                                  slopeFlDens=sample(slopeFlDensSeedCount,1),slopeSurv=sample(slopeSurvSeedCount,1),
                                  stDev=ifelse(plotVar,sample(sigmaSeedCount_plant,1),0),transform=T,SIMPLIFY=F))
  
  # #OK, but spread seems too small - perhaps a t-dist would be better?
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simSeedCount),plantsAllSeed$AvPodCount,na.rm=T)+2
  # hist(unlist(simSeedCount),breaks=seq(0,upr,2),main=NULL,xlab='Sim Avg Seed Count')
  # hist(plantsAllSeed$AvPodCount,breaks=seq(0,upr,2),main=NULL,xlab='Actual Avg Seed Count')
  
  #Simulate avg seed weight (mg)
  simSeedWeight <- with(dat,mapply(genSeedWeight,pol=simAvgPol,avgSeedCount=simSeedCount,plSize=simPlSize,
                                   plDens=simPlDens,is2016=is2016,lbeeDist=lbeeDist,isHalfStock=isHalfStock,
                                   int=sample(intSeedWeight,1),slopePol=sample(slopePolSeedWeight,1),
                                   slopeSeedCount=sample(slopeSeedCount,1),slopePlSize=sample(slopePlSizeSeedWeight,1),
                                   slopePlDens=sample(slopePlDensSeedWeight,1),slope2016=sample(slope2016SeedWeight,1),
                                   slopeLbeeDist=sample(slopeLbeeDistSeedWeight,1),slopeStocking=sample(slopeStockingSeedWeight,1),
                                   slopePlDensPlSize=sample(slopePlDensPlSizeSeedWeight,1),
                                   lambda=sample(lambdaSeedWeight,1),stDev=sample(sigmaSeedWeight_plant,1),addVar=plotVar,
                                   SIMPLIFY=F))
  
  # #Looks good
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simSeedWeight),with(plantsAllSeed,AvPodMass*1000/AvPodCount),na.rm=T)+1
  # min <- min(unlist(simSeedWeight),with(plantsAllSeed,AvPodMass*1000/AvPodCount),na.rm=T)-1
  # hist(unlist(simSeedWeight),xlab='Sim weight per seed (mg)',main=NULL,breaks=seq(min,upr,0.2))
  # hist(with(plantsAllComm,AvPodMass*1000/AvPodCount),xlab='Actual weight per seed (mg)',main=NULL,breaks=seq(min,upr,0.2))
  
  #Simulate total yield for each plot
  simYield <- with(dat,mapply(genYield,pods=simPodCount,seedCount=simSeedCount,seedWeight=simSeedWeight,
                              int=sample(intYield,1),slopeYield=sample(slopeYield,1)))
  
  #Looks OK. Tails aren't as long, but this is probably due to tails from sim seed count
  par(mfrow=c(2,1))
  hist(simYield,xlab='Sim yield/m2',main=NULL,breaks=seq(0,1000,20))
  hist(within(datalist,tapply(yield*exp(plDens+3.573633)[plantIndex],plantIndex,mean)),
       xlab='Actual yield/m2',main=NULL,breaks=seq(0,1000,20))
  
  # if(any(!is.finite(simYield))) stop('Non-finite yield estimate:',yield[!is.finite(simYield)])
  
  if(!returnAll){
    return(simYield)
  } else {
    a <- list(simPlDens,simPlSize,simFlwDens,simHbeeVis,simAvgPol,simFlwCount,
              simPodCount,simSeedCount,simSeedWeight,simYield)
    return(lapply(a,function(x) lapply(x,mean))) #Mean of each distance
  }
}

# #Try out generic fields
# results <- with(filter(surveyAllSeed,Bay=='F'),replicate(50,simSeed(hbeeDist=Distance,lbeeDist=minDist,
#                 isCent=EdgeCent=='Center',is2016=Year==2016,
#                 isHalfStock=Treatment=='Double tent',dat=mod3)))
# 
# results2 <- filter(surveyAllSeed,Bay=='F') %>%
#   select(Field,Distance,EdgeCent,PlDens) %>% #Simulated yield per plot
#   unite(ID,Field:EdgeCent,sep='_') %>%
#   bind_cols(data.frame(t(apply(results,1,function(x) quantile(x,c(0.5,0.05,0.95),na.rm=T))))) %>%
#   rename('pred'='X50.','lwr'='X5.','upr'='X95.')
# 
# #Looks OK, at least for average field
# ungroup(plantsAllSeed) %>% select(Field,Distance,EdgeCent,SeedMass) %>% #Actual yield per plot
#   unite(ID,Field:EdgeCent) %>% group_by(ID) %>% summarize(totalSeeds=mean(SeedMass,na.rm=T)) %>%
#   left_join(results2,by='ID') %>% mutate(actual=totalSeeds*PlDens) %>% select(-totalSeeds,-PlDens) %>%
#   ggplot(aes(x=actual,y=pred))+
#   # geom_point()+
#   geom_pointrange(aes(ymax=upr,ymin=lwr),alpha=0.4)+
#   geom_abline(intercept=0,slope=1,col='red')+
#   labs(x='Actual Yield',y='Predicted Yield')

#Simulate pollination effects away from shelters
scenario <- expand.grid(hDist=c(1),lDist=seq(1,51,2),cent=c(0,1),halfStock=c(0))
results <- replicate(500,with(scenario,simSeed(hbeeDist=hDist,lbeeDist=lDist,isCent=cent,
                                               is2016=0,isHalfStock=halfStock,dat=mod3,plotVar=F)))
# beep()
#Simulated yield distribution
data.frame(scenario,t(apply(results,1,function(x) quantile(x,c(0.5,0.05,0.95),na.rm=T)))) %>%
  rename('pred'='X50.','lwr'='X5.','upr'='X95.') %>% 
  mutate_at(vars(pred,upr,lwr),g2bushels) %>% #Results in bu/acre
  mutate(edgeCent=factor(cent,labels=c('Edge','Center'))) %>% 
  ggplot(aes(x=lDist,y=pred))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=edgeCent),alpha=0.2,show.legend=F)+
  geom_line(aes(col=edgeCent),size=1)+
  labs(y='Predicted yield (bu/acre)',x='Shelter Distance',col='Bay\nPosition')+
  scale_colour_manual(values=c('blue','red'))+scale_fill_manual(values=c('blue','red'))

#Simulate pollination effects away from edge of field
scenario <- expand.grid(hDist=seq(1,401,10),lDist=c(10),cent=c(0,1),halfStock=c(0))
results <- replicate(500,with(scenario,simSeed(hbeeDist=hDist,lbeeDist=lDist,isCent=cent,
                                               is2016=0,isHalfStock=halfStock,dat=mod3,plotVar=F)))
# beep()
#Simulated yield distribution
p1 <- data.frame(scenario,t(apply(results,1,function(x) quantile(x,c(0.5,0.05,0.95),na.rm=T)))) %>%
  rename('pred'='X50.','lwr'='X5.','upr'='X95.') %>% 
  mutate(edgeCent=factor(cent,labels=c('Edge','Center'))) %>% 
  mutate_at(vars(pred,upr,lwr),g2bushels) %>% #Results in bu/acre
  ggplot(aes(x=hDist,y=pred))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=edgeCent),alpha=0.2,show.legend=F)+
  geom_line(aes(col=edgeCent),size=1)+
  labs(y='Predicted yield (bu/acre)',x='Distance into field',col='Bay\nPosition')+
  scale_colour_manual(values=c('blue','red'))+scale_fill_manual(values=c('blue','red'))
ggsave('../Figures/Seed/simHbeeDistYield.png',p1,width=8,height=5)


#Simulate regular bay scenario

#Create standard field (7m bays * 86 bays across = 602m)
scal <- (800/602) #Scaling factor
bayType <- matrix(rep(rep(c(0,1,1,1,1,1,1),86),602),ncol=602,byrow=T) #1 = F, 0 = M
centType <- matrix(rep(rep(c(NA,0,0.5,1,1,0.5,0),86),602),ncol=602,byrow=T) #1 = cent, 0 = edge, 0.5 = half-cent, NA = Mbay
xloc <- matrix(rep(c(1:602),602),ncol=602,byrow=T)*scal - (scal/2) #Dist across
yloc <- matrix(rep(c(1:602),each=602),ncol=602,byrow=T)*scal - (scal/2) #Dist down
dists <- matrix(NA,ncol=602,nrow=602) #Distance matrix for distance from edge
for(row in 1:602){
  for(col in 1:602){
    # #Distance from hives - takes ~10 sec
    # dists[row,col] <- min(dist(matrix(c(xloc[row,col],yloc[row,col],0,0,0,800,800,0,800,800),ncol=2,byrow=T))[1:4])
    dists[row,col] <- min(xloc[row,col],abs(xloc[row,col]-800),yloc[row,col],abs(yloc[row,col]-800))
  }
}

#Lbee spacing (72 tents total)
#1st col= rows(y), 2nd col= cols(x)
lbeeLocs <- cbind(round(rep(c((602/8)*(1:8)-(602/8)/2,(602/7)*(0:7)),length.out=8*9)),
                  rep(2+(7*9)*(1:9),each=8))*scal
lbeeLocs[lbeeLocs[,1]==0,1] <- 1*scal #Makes sure all tents within field
plot(lbeeLocs) #OK

#Functions to extract lbee distances
getLbeeDist <- function(x,y) dist(rbind(x,y))
minDist <- function(x,y,lbeeLocs) {
  startLim <- 30 
  foundTents <- F
  while(!foundTents){
    xlims <- x + startLim*c(-1,1); ylims <- y + startLim*c(-1,1) #Search limits 
    lookHere <- lbeeLocs[,1]>xlims[1]  & lbeeLocs[,1]<xlims[2] & lbeeLocs[,2]>ylims[1]  & lbeeLocs[,2]<ylims[2]
    if(any(lookHere)) foundTents <- T
    if(startLim>500) stop('No tents within 500 m')
    startLim <- startLim+50 #Look 50m further
  }
  
  if(sum(lookHere)==1) { #If only 1 tent found
    return(dist(rbind(lbeeLocs[lookHere,],c(x,y))))
  } else { 
    return(min(apply(lbeeLocs[lookHere,],1,getLbeeDist,y=c(x,y))))
  }
}

lbeeDists <- matrix(NA,ncol=602,nrow=602) #Takes a minute to generate
for(i in 1:602){
  lbeeDists[,i] <- sapply(1:602*scal,minDist,y=i*scal,lbeeLocs=lbeeLocs)
}
lbeeDists[lbeeDists==0] <- 0.5 #Set minimum distance to half a meter

#Looks OK
par(mfrow=c(2,2))
image(1:602*scal,1:602*scal,t(lbeeDists),xlab='',ylab='',main='Lbee Distances')
points(lbeeLocs[,c(2:1)],pch=19,col='black')
image(1:602*scal,1:602*scal,t(dists),main='Hbee distances',xlab='',ylab='')
image(1:602*scal,1:602*scal,t(bayType),xlab='',ylab='',main='Bays')
image(1:602*scal,1:602*scal,t(centType),xlab='',ylab='',main='Cent/Edge')
par(mfrow=c(1,1))

scenario <- data.frame(x=rep(1:602,each=602),y=rep(1:602,602),hbeeDist=as.vector(dists),lbeeDist=as.vector(lbeeDists)) %>% 
  mutate(bay=as.vector(bayType),cent=as.vector(centType)) %>% filter(bay!=0) #Get rid of male bays
scenario$yield <- with(scenario,simSeed(hbeeDist=hbeeDist,lbeeDist=lbeeDist,isCent=cent,
                                        is2016=0,isHalfStock=0,dat=mod3,plotVar=F))
beep(1)
mean(g2bushels(scenario$yield)) #~36 bushels per acre on average

scenario %>% 
  ggplot(aes(x=x*scal,y=y*scal))+geom_raster(aes(fill=g2bushels(yield)))+
  labs(x=NULL,y=NULL,fill='Yield\n(bu/\nacre)')+
  geom_point(data=data.frame(lbeeLocs),aes(X2,X1),col='black')+
  # scale_fill_gradient(low='yellow',high='blue')+
  coord_cartesian(xlim=c(200,600),ylim=c(200,600))

# Test seed field model using piecewiseSEM -------------------------------------------
library(piecewiseSEM)
library(lme4)

#This model is very coarse, because quantities are all modeled at plot level, but this gets the general sense of the model, and allows d-sep tests to be made

#Model all sub-plot quantities at the plot level:
#Pollen
polMod <- glmer.nb(Pollen~(1|Field/FieldPlot),data=pollenAllSeed)
pollenAllSeed$pred <- predict(polMod)
surveyAllSeed <- pollenAllSeed %>% #Joins plot-level pollen
  group_by(Field,Distance,EdgeCent) %>% 
  summarize(predPol=first(pred)) %>% 
  unite(ID,Field,Distance,EdgeCent,sep='_') %>% 
  left_join(unite(surveyAllSeed,ID,Field,Distance,EdgeCent,sep='_',remove=F),.,by='ID') %>% 
  dplyr::select(-ID)
#Plant size
plantsAllSeed <- plantsAllSeed %>% 
  unite(FieldPlot,Field,Distance,EdgeCent,remove=F) #Create FieldPlot column
plantMod <- lmer(log(VegMass)~(1|Field/FieldPlot),data=plantsAllSeed)
plantsAllSeed$predPlSize <- predict(plantMod,newdata=plantsAllSeed,allow.new.levels=T) #Guesses at population level for NAs
surveyAllSeed <- plantsAllSeed %>% #Joins plot-level plant size
  group_by(Field,Distance,EdgeCent) %>% 
  summarize(predPlSize=first(predPlSize)) %>% 
  unite(ID,Field,Distance,EdgeCent,sep='_') %>% 
  left_join(unite(surveyAllSeed,ID,Field,Distance,EdgeCent,sep='_',remove=F),.,by='ID') %>% 
  dplyr::select(-ID)
#Flowers per plant 
plantsAllSeed <- mutate(plantsAllSeed,flwCount=Pods+Missing)
flwCountMod <- glmer.nb(flwCount~(1|Field/FieldPlot),data=plantsAllSeed)
plantsAllSeed$predFlw <- predict(flwCountMod,newdata=plantsAllSeed,allow.new.levels=T) #Guesses at population level for NAs
surveyAllSeed <- plantsAllSeed %>% #Joins plot-level plant size
  group_by(Field,Distance,EdgeCent) %>% 
  summarize(predFlwCount=first(predFlw)) %>% 
  unite(ID,Field,Distance,EdgeCent,sep='_') %>% 
  left_join(unite(surveyAllSeed,ID,Field,Distance,EdgeCent,sep='_',remove=F),.,by='ID') %>% 
  dplyr::select(-ID)
#Seeds per pod
seedsAllSeed <- seedsAllSeed %>% 
  unite(FieldPlot,Field,Distance,EdgeCent,remove=F) #Create FieldPlot column
seedCountMod <- lmer(sqrt(AvPodCount)~(1|Field/FieldPlot),data=plantsAllSeed) #Uses avg seed size, but OK for now
plantsAllSeed$predSeedCount <- predict(seedCountMod,na.action=na.exclude) #NAs included
surveyAllSeed <- plantsAllSeed %>% #Joins plot-level plant size
  group_by(Field,Distance,EdgeCent) %>% 
  summarize(predSeedCount=first(predSeedCount)) %>% 
  unite(ID,Field,Distance,EdgeCent,sep='_') %>% 
  left_join(unite(surveyAllSeed,ID,Field,Distance,EdgeCent,sep='_',remove=F),.,by='ID') %>% 
  dplyr::select(-ID)
#Seed size (mg)
seedSizeMod <- lmer((AvPodMass*1000/AvPodCount)~(1|Field/FieldPlot),data=plantsAllSeed) #Uses avg seed size, but OK for now
plantsAllSeed$predSeedSize <- predict(seedSizeMod,na.action=na.exclude) #NAs included
surveyAllSeed <- plantsAllSeed %>% #Joins plot-level plant size
  group_by(Field,Distance,EdgeCent) %>% 
  summarize(predSeedSize=first(predSeedSize)) %>% 
  unite(ID,Field,Distance,EdgeCent,sep='_') %>% 
  left_join(unite(surveyAllSeed,ID,Field,Distance,EdgeCent,sep='_',remove=F),.,by='ID') %>% 
  dplyr::select(-ID)
#Pod success
podSuccessMod <- glmer(cbind(Pods,Missing)~(1|Field/FieldPlot),family=binomial,data=plantsAllSeed)
plantsAllSeed$predPodSuccess <- predict(podSuccessMod,na.action=na.exclude) #NAs included
surveyAllSeed <- plantsAllSeed %>% #Joins plot-level plant size
  group_by(Field,Distance,EdgeCent) %>% 
  summarize(predPodSuccess=mean(predPodSuccess,na.rm=T)) %>% 
  unite(ID,Field,Distance,EdgeCent,sep='_') %>% 
  left_join(unite(surveyAllSeed,ID,Field,Distance,EdgeCent,sep='_',remove=F),.,by='ID') %>% 
  dplyr::select(-ID)

#Scale data and choose only columns of interest
temp <- surveyAllSeed %>% ungroup() %>% 
  mutate(Stocking=ifelse(Treatment=='Double tent','Half','Full')) %>%  
  dplyr::select(Field,Year,Distance,minDist,Bay,EdgeCent,Stocking,TotalTime,FlDens,PlDens,lbee,hbee,predPol:predPodSuccess) %>% 
  mutate(TotalTime=ifelse(is.na(TotalTime),5,TotalTime),is2016=Year==2016) %>% 
  mutate(logDist=scale(log(Distance)),logMinDist=scale(log(minDist)),logTime=log(TotalTime)) %>%  #Scales predictors
  mutate(sqrtFlDens=sqrt(FlDens)) %>% 
  mutate_at(vars(sqrtFlDens,predPol:predPodSuccess),funs('scale')) %>% #Scales plot-level terms
  filter(!is.na(predPlSize),!is.na(predSeedCount),!is.na(FlDens)) %>% 
  mutate(lbeeVisRate=log(1+lbee/(TotalTime/10)),hbeeVisRate=log(1+hbee/(TotalTime/10)))

temp %>% mutate(minDist=log(minDist)-mean(log(minDist)),Distance=log(Distance)-mean(log(Distance))) %>%
  mutate(FlDens=sqrt(FlDens)-mean(sqrt(FlDens))) %>% 
  lmer(lbeeVisRate~minDist+Distance*EdgeCent+Stocking+Bay+Distance:Stocking+factor(Year)*Distance+FlDens+(Distance|Field),data=.) %>% 
  summary()

#Fill in missing values for PlDens (M bay, mainly)
temp$PlDens[is.na(temp$PlDens)] <- predict(lmer(PlDens~logDist+(1|Field),data=temp),newdata=temp,allow.new.levels=T)[is.na(temp$PlDens)] 
visMod1 <- psem( #Original model
  lmer(PlDens~logDist+(1|Field),data=temp), #Good
  lmer(predPlSize~is2016+PlDens+(1|Field),data=temp), #Good
  lmer(sqrtFlDens~predPlSize+(1|Field),data=temp), #Good
  # glmer.nb(hbee~sqrtFlDens+logMinDist+logDist+Bay+EdgeCent+offset(logTime)+(1|Field), #Good
  #          data=temp,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))), #hbee visits
  # glmer.nb(lbee~sqrtFlDens+logMinDist+Bay+EdgeCent+Stocking+offset(logTime)+(1|Field), #Good
  #          data=temp,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))), #lbee visits
  lmer(hbeeVisRate~sqrtFlDens+logMinDist+logDist+Bay+EdgeCent+(1|Field),data=temp), #hbee visits - good
  lmer(lbeeVisRate~sqrtFlDens+logMinDist+Bay+EdgeCent+Stocking+(1|Field),data=temp), #lbee visits - good
  lmer(predPol~lbeeVisRate+hbeeVisRate+(1|Field),data=temp), #log Pollen counts - good
  lmer(predFlwCount~predPlSize+(1|Field),data=temp), #Flower count per plant - good
  lmer(predPodSuccess~predPol+predPlSize+(1|Field),data=temp), #logit pod succes
  lmer(predSeedCount~predPol+predPlSize+(1|Field),data=temp), #log Seed counts
  lmer(predSeedSize~predSeedCount+predPol+predPlSize+(1|Field),data=temp) #log Seed size
)

dSep(visMod1,direction=c('lbee -> hbee'),conditioning=F)
fisherC(visMod1,conserve=T) 
summary(visMod1)

# lbeeVis~hbeeDist*cent+lbeeDist+lbeeStocking*hbeeDist+Fbay

visMod2 <- psem( #Updated model
  lmer(PlDens~logDist+(1|Field),data=temp),
  lmer(predPlSize~logDist+PlDens+is2016+(1|Field),data=temp),
  lmer(sqrtFlDens~is2016+predPlSize+logDist+logMinDist+(1|Field),data=temp),
  glmer.nb(lbee~sqrtFlDens+logDist*EdgeCent+logMinDist+Stocking*logDist+Bay+offset(logTime)+(1|Field),
           data=temp,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))), #lbee visits
  glmer.nb(hbee~is2016+logMinDist+logDist+Bay+EdgeCent+offset(logTime)+(1|Field),
           data=temp,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))), #hbee visits
  lmer(predPol~logMinDist+logDist+EdgeCent+lbee+hbee+(1|Field),data=temp), #log Pollen counts
  lmer(predSeedCount~predPodSuccess+predPol+predPlSize+EdgeCent+logDist+(1|Field),data=temp), #log Seed counts
  lmer(predSeedSize~is2016+predPol+predPlSize+predSeedCount+(1|Field),data=temp), #log Seed size
  lmer(predPodSuccess~predPol+predPlSize+logMinDist+logDist+EdgeCent+(1|Field),data=temp) #logit pod success
)

summary(visMod2)

dSep(visMod2,direction=c('lbee -> hbee'),conditioning=F)
fisherC(visMod2,conserve=T) #31.37, p=0.77

#dSep test points to relationship b/w pollen and distance into field, independent of visitation rate. This suggests that jump lengths of pollinators matter.
#Alternatively, it could be that the overall rate is important, and not the observed visitation rate -> pollen related to distance, but not lbee or hbee counts.
#Seedcount model has strong effect of hbee dist, pollen, and edge/center
#Seedsize model has only seed count
#Podsuccess model has strong effect of seedCount, lbee and hbee dist, pollen, and edge/center. Seems to have very similar results to seedcount model, so perhaps this is actually a reflection of the same process? Need to think about this a bit more.

#Next steps: 1) write "original" model in Stan, get dSep tests from PSEM, write/run dSep tests in Stan, 
#2) write "final" model in Stan, get dSep tests from PSEM, write/run dSep tests in Stan

#Summary function doesn't work, so using lapply
lapply(visMod2[1:6],summary)

dSep(visMod2,direction=c('lbee -> hbee'),conditioning=F)
fisherC(visMod2,conserve=T) #31.37, p=0.77

temp %>% dplyr::select(predPol:predPodSuccess) %>% pairs()

ggplot(temp,aes(exp(logDist),hbee/exp(logTime)))+geom_point()+
  labs(x='Distance',y='Visits/min')

hbeeMod <- glmer.nb(hbee~logDist+logMinDist+Bay+EdgeCent+offset(logTime)+(1|Field),data=temp)
summary(hbeeMod)
plot(hbeeMod)



