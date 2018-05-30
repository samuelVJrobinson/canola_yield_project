#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN SEED CANOLA FIELDS (2015)
# Libraries and ggplot theme ----------------------------------------------
library(ggplot2)
library(MCMCglmm)
library(nlme)
library(mgcv)
library(MASS)
library(dplyr)
library(tidyr)
library(reshape2)

#Big-text theme, no grid lines (used for Bayer 2016 presentation)
prestheme=theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15),
                axis.title=element_text(size=20),
                title=element_text(size=20),
                axis.line.x=element_line(colour='black'),
                axis.line.y=element_line(colour='black'),
                #panel.grid.major=element_line(size=0.5,colour='black',linetype='dotted'),
                #panel.grid.minor=element_blank,
                panel.border=element_blank())
theme_set(theme_classic()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

# Load and organize data  ---------------------------------------------------------

load("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\seedfieldData2015.RData")

# fields=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\field info 2015.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
#
# survey=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\survey data 2015.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
#
# plants=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\plant info 2015 seed.csv",stringsAsFactors=T,strip.white=T,na.strings = c("",'NA'))
#
# se=function(x) sd(x)/sqrt(length(x)) #Convenience SE function
#
# AICc=function(model,...){ #Second-order bias corrected AIC
#   #require(plyr)
#   corAIC=function(x){
#     k=length(coef(x))
#     n=length(residuals(x))
#     return(AIC(x)+(2*k*k-1)/(n-k-1))
#   }
#   if(missing(...)){
#     return(corAIC(model))
#   }
#   else{
#     mlist=list(model,...)
#     len=plyr::laply(mlist,function(x) length(residuals(x)))
#     if(length(len)!=sum(len %in% len[1]))
#       warning('Models use different # data points. Comparing AICc values is invalid')
#     forms=sapply(mlist, function(x) paste(formula(x)))
#     data.frame(row.names=paste(forms[2,],forms[1,],forms[3,]),AICc=plyr::laply(mlist,corAIC))
#   }
# }
#
# deltaAIC=function(m1,m2,...,correct=TRUE){
#   if(correct) a=AICc(m1,m2,...) else a=AIC(m1,m2,...)
#   a$delta_AICc=round(a$AICc-min(a$AICc),4)
#   a$w_i=round(exp(-0.5*a$delta_AICc)/sum(exp(-0.5*a$delta_AICc)),4) #Model Probabilities (Akaike weights)
#   a
# }
#
# DIC=function(m1,...){ #Convenience function to extract DIC values from MCMCglmm models
#   modlist=list(m1,...)
#   dicvals=unlist(lapply(modlist,'[[',c('DIC')))
#   fixeffs=as.character(lapply((lapply(modlist,'[[','Fixed')),'[[','formula'))
#   raneffs=as.character(lapply((lapply(modlist,'[[','Random')),'[[','formula'))
#   data.frame(Fixed=fixeffs,Random=raneffs,DIC=dicvals,deltaDIC=dicvals-min(dicvals))
# }
#
# plotFixedTerms=function(m1){ #Convenience function to plot 95% distributions for MCMCglmm fixed effects
#   dat=as.data.frame(summary(m1)$solutions)
#   dat$params=as.factor(rownames(dat))
#   colnames(dat)[c(2,3)]=c('lwr','upr')
#   ggplot(dat,aes(x=params,y=post.mean,ymax=upr,ymin=lwr))+geom_pointrange()+labs(x='Fixed Effects',y='Posterior Mean')+geom_hline(yintercept=0,colour='red')
# }
#
# zeroCheck=function(mod,dat) { #Posterior predictive checks for zero-inflation (MCMCglmm course notes pg. 104)
#   require(MCMCglmm)
#   require(ggplot2)
#   try(if(is.null(mod$X)) stop("Needs fixed effect design matrix (use saveX=T in MCMCglmm arguments)"))
#   nz=1:nrow(mod$Sol) #Number of samples to hold
#   oz=sum(dat[as.character(mod$Fixed$formula[2])]==0) #Actual number of 0s in the dataset
#   for (i in 1:nrow(mod$Sol)) {
#     pred.l <- rnorm(nrow(dat), (mod$X %*% mod$Sol[i,])@x, sqrt(mod$VCV[i])) #Predicted mean
#     nz[i] <- sum(rpois(nrow(dat), exp(pred.l)) == 0)} #Expected number under poisson
#   qplot(nz,xlab='Expected number of zeros')+geom_vline(aes(xintercept=oz))
#   #Expected zeros, given the overdispersed poisson model
#   }
#
# varComp=function(mod,units=T) { #Tabulates variance components of random effects for MCMCglmm models
#   VCV=mod$VCV
#   if(units==F) VCV=VCV[,c(1:ncol(VCV)-1)] #Strips out "units" column (for binomial models)
#   round(data.frame(mode=posterior.mode(VCV/rowSums(VCV)),HPDinterval(VCV/rowSums(VCV))),3)
# }
#
# #Load data
#
# #field-level data
# fields=filter(fields,Type=='Seed') %>% #Strips out commodity field data
#   select(Field,Deg.N,Deg.W,Hectares,Treatment,NumHives=Number.of.Hives, #Selects columns of interest
#          Surveyed,Start.Time,End.Time,BowlStart,BowlEnd,
#          AirTemp=Temperature,Wind=Wind.Speed.km.hr.,RH=Humidity,Harvested) %>%
#   mutate(Start.Time=as.POSIXct(paste(Surveyed,Start.Time),format='%b %d, %Y %H:%M'), #Converts dates to single format
#          End.Time=as.POSIXct(paste(Surveyed,End.Time),format='%b %d, %Y %H:%M'),
#          BowlStart=as.POSIXct(paste(Surveyed,BowlStart),format='%b %d, %Y %H:%M'),
#          BowlEnd=as.POSIXct(paste(Surveyed,BowlEnd),format='%b %d, %Y %H:%M'),
#          SurveyTime=End.Time-Start.Time, #Time spent on survey
#          BowlTime=BowlEnd-BowlStart, #Hours that bowl traps were used
#          SurveyStart=Start.Time) %>%
#   select(-BowlEnd,-Start.Time,-End.Time) #Removes end times for bowls and survey
#
# #plot-level data
#
# #Bag numbers and associated data to match with plants
# survey=filter(survey,Type=='Seed') %>% #Strips out commodity field data
#   select(Field,Distance,Bay:EndTime,FlDens=Fls.50cm.2,HoneyTP=HoneybeeTopPollen,HoneyTN=HoneybeeTopNectar, #Selects columns of interest
#          HoneySP=HoneybeeSidePollen,HoneySN=HoneybeeSideNectar,NumHTP:NumHSN,
#          Leafbee=Leafcutterbee,NumLeafbee=NumLeafcutterbee,Otherbee:NumFly,Len1:Len5,Microcap.volume,
#          Soilwater:MountingOK,PlDens=Plants.m.2,Bags=BagNumbers)
#
# #Temporary plot-level dataframe (saved for converting to flower-level data)
# flowers=survey
#
# survey=transmute(rowwise(survey),Field,Distance,Bay,EdgeCent=Edge.Cent,EdgeDir,StartTime,EndTime,
#                  FlDens,HoneyTP,HoneyTN,HoneySP,HoneySN,NumHTP,NumHTN,NumHSP,NumHSN,
#                  Leafbee,NumLeafbee,Otherbee,NumOtherbee,Hoverfly,NumHoverfly,Fly,NumFly,
#                  AvLen=mean(Len1,Len2,Len3,Len4,Len5,na.rm=T), #Avg Nectar(uL)
#                  hasNectar=sum(cbind(Len1,Len2,Len3,Len4,Len5)>0), #Number of nectar samples >0
#                  Soilwater,minShelter=min(LeafShelter1,LeafShelter2,LeafShelter3,LeafShelter4,na.rm=T), #Distance to closest leafcutter shelter
#                  AvgPol=mean(as.numeric(Pollen1,Pollen2,Pollen3,Pollen4,Pollen5)), #Average pollen per stigma -
#                  #WEIRD PROBLEM: previous line doesn't work unless casted to numeric
#                  #http://stackoverflow.com/questions/29224719/dplyr-error-strange-issue-when-combining-group-by-mutate-and-ifelse-is-it-a-b
#                  hasPol=rowSums(cbind(Pollen1,Pollen2,Pollen3,Pollen4,Pollen5)>0), #Number of stigmas with pollen >0
#                  PlDens,Bags,FieldPlot=paste(Field,Distance,Bay,Edge.Cent,sep=',')) %>%  #Plant density, and unique ID
#   left_join(select(fields,Field,Hectares:Surveyed,AirTemp:RH),by='Field') %>% #Joins field-level data
#   mutate(StartTime=as.POSIXct(paste(Surveyed,StartTime),format='%b %d, %Y %H:%M'),
#          EndTime=as.POSIXct(paste(Surveyed,EndTime),format='%b %d, %Y %H:%M'),TotalTime=EndTime-StartTime) %>%
#   select(-EndTime)
#
# #Honeybee-only dataframe to examine behaviour (top/side-working & nectar/pollen)
# hbees=select(survey,Field:EdgeCent,FieldPlot,HoneyTP:NumHSN) %>%
#   unite(Top_Pollen,HoneyTP,NumHTP) %>% unite(Top_Nectar,HoneyTN,NumHTN) %>%  #Unite 8 honeybee columns into 4
#   unite(Side_Pollen,HoneySP,NumHSP) %>% unite(Side_Nectar,HoneySN,NumHSN) %>%
#   gather(Behaviour,Count,-Field:-FieldPlot) %>% # Gather into 1 long column
#   separate(Count,c('Visits','Count'),convert=T) %>% #Separate into visits and counts
#   separate(Behaviour,c('Approach','Foraging'),remove=F) %>% #Separates Top/Side and Nectar/Pollen
#   mutate(Approach=as.factor(Approach),Foraging=as.factor(Foraging)) %>% #Converts new columns to factor
#   left_join(select(survey,FieldPlot,StartTime,FlDens,Leafbee,NumLeafbee,minShelter,AvgPol,hasPol,Hectares:RH),
#             by='FieldPlot') %>% #Joins in leafcutter visits (could add other visitors later), along with other variables
#   select(-FieldPlot)
#
# survey=survey %>%
#   mutate(Honeybee=sum(HoneyTP,HoneyTN,HoneySP,HoneySN),NumHoneybee=sum(NumHTP,NumHTN,NumHSP,NumHSN)) %>%
#   select(-HoneyTP:-NumHSN) #Cleans up main survey dataframe (removes honeybee behaviour info)
#
# #flower-level data - Pollen and nectar dataframes could be avoided using the strategy from hbees (i.e. unite, gather, separate)
# flowers=rowwise(flowers)%>%
#   mutate(StartTime=as.POSIXct(StartTime,format='%H:%M'), #Start time
#                TotalTime=as.POSIXct(EndTime,format='%H:%M')-as.POSIXct(StartTime,format='%H:%M'), #Total time
#                minShelter=min(LeafShelter1,LeafShelter2,LeafShelter3,LeafShelter4,na.rm=T)) %>% #Closest leafcutter shelter
#   select(Field:EdgeDir,minShelter,FlDens:Soilwater,Pollen1:Pollen5,PlDens) %>% #Cuts out lat/lon, individual leafcutter distances
#   rename(EdgeCent=Edge.Cent) %>%
#   unite(FieldPlot,Field,Distance,Bay,EdgeCent,sep=',',remove=F) %>% #Unique ID for join with flower-level data
#   unite(LenPol1,Len1,Pollen1) %>% unite(LenPol2,Len2,Pollen2) %>% #Unite 10 pollen/nectar cols into 5
#   unite(LenPol3,Len3,Pollen3) %>% unite(LenPol4,Len4,Pollen4) %>% unite(LenPol5,Len5,Pollen5) %>%
#   gather(Flower,Count,LenPol1:LenPol5) %>%
#   separate(Count,c('Length','Pollen'),sep='_',convert=T) %>% #Separate into visits and counts
#   mutate(Flower=as.numeric(Flower)) %>%
#   left_join(select(survey,FieldPlot:NumHoneybee),by='FieldPlot') %>% #Joins plot-level data
#   select(-HoneyTP:-NumHSN) #Removes unique ID, and honeybee behavioural data
#
# #Makes honeybee-only dataframe
# hbees=filter(hbees,Behaviour!='Side_Pollen') %>% #There were no side-working pollen foragers anywhere.
#   mutate(Field=as.factor(Field),Bay=as.factor(Bay),EdgeCent=as.factor(EdgeCent))
#
# #Makes wide-form dataframe for comparing proportions of visitors
# hbeesWide=unite(hbees,Visits_Count,Visits,Count) %>%
#   select(-Approach,-Foraging) %>%
#   spread(Behaviour,Visits_Count) %>%
#   separate(Top_Pollen,c('TopPolVis','TopPolCount'),sep='_',convert=T) %>% #Top working pollen, Visits and Counts
#   separate(Top_Nectar,c('TopNecVis','TopNecCount'),sep='_',convert=T) %>% #Top working nectar, Visits and Counts
#   separate(Side_Nectar,c('SideNecVis','SideNecCount'),sep='_',convert=T) %>% #Side working nectar, Visits and Counts
#   rowwise() %>%
#   mutate(hbeeVis=TopPolVis+TopNecVis+SideNecVis,hbeeCount=TopPolCount+TopNecCount+SideNecCount,
#          VisCount=hbeeVis/hbeeCount) #Totals of everything
#
# #Leafcutter bees only
# lbees=select(hbeesWide,Field:RH,hbeeVis:hbeeCount,VisCount) %>%
#   rename(hbeeVisCount=VisCount,lbeeVis=Leafbee,lbeeCount=NumLeafbee) %>%
#   ungroup()
#
# #Plant data
# bagdat=filter(survey,Bay=='F') %>%
#   select(Bags,Field,Distance,EdgeCent,EdgeDir,Leafbee:NumFly,minShelter,AvgPol,PlDens,Hectares,Treatment,NumHives,Honeybee,NumHoneybee) %>% #Data to merge from survey dataframe
#   separate(Bags,c('Bag1','Bag2','Bag3','Bag4','Bag5'),convert=T) %>%
#   gather(BagNum,Bag,1:5) %>%
#   select(-BagNum) %>%
#   rename(BagNum=Bag)
#
# plants=rename(plants,BagNum=Bag.) %>%
#   arrange(FieldName,BagNum) #Re-orders plant measurements
#
# seeds=transmute(plants,Field=FieldName,BagNum,VegMass=TotalMass-SeedMass,SeedMass,Branch,
#                 Pods=Pods+Bag+BagTscar,Missing=Missing-Bag,
#                 Pod1,Pod2,Pod3,Pod4,Pod5,
#                 Weigh1,Weigh2,Weigh3,Weigh4,Weigh5,SeedCount) %>% #Selects correct columns
#   gather('Pod','Value',Pod1:Weigh5) %>% #Melts pod count and pod weight
#   separate(Pod,into=c('Param','PodNum'),sep=-2) %>%
#   mutate(Param=as.factor(Param),PodNum=as.numeric(PodNum)) %>%
#   spread(Param,Value) %>% #Casts pod count and pod weight back
#   rename(PodCount=Pod,PodMass=Weigh,Pod=PodNum) %>%
#   group_by(Field,BagNum,Pod)
#
# plants=summarise_each(group_by(seeds,Field,BagNum),funs(mean,se))%>% #Mean and SE of metrics for each plant (all mean values except for pod measurements will equal plant-level metrics)
#   select(Field,BagNum,VegMass=VegMass_mean,SeedMass=SeedMass_mean,
#          Branch=Branch_mean,Pods=Pods_mean,Missing=Missing_mean,SeedCount=SeedCount_mean,
#          AvPodCount=PodCount_mean,SEPodCount=PodCount_se,AvPodMass=PodMass_mean,SEPodMass=PodMass_se)
#
# #Merges Distances and other data into seed and plant dataframes
# seeds=left_join(seeds,select(bagdat,-Field),by='BagNum') %>%
#   rename(Plant=BagNum) %>%
#   arrange(Field,Distance,Plant) %>%
#   ungroup()%>%
#   unite(FieldPlot,Field,Distance,EdgeCent,remove=F)%>%
#   mutate(FieldPlot=as.factor(FieldPlot))%>%
#   filter(!is.na(VegMass))
#
# plants=left_join(plants,select(bagdat,-Field),by='BagNum') %>%
#   rename(Plant=BagNum) %>%
#   arrange(Field,Distance,Plant) %>%
#   ungroup() %>%
#   unite(FieldPlot,Field,Distance,EdgeCent,remove=F) %>%
#   mutate(FieldPlot=as.factor(FieldPlot)) %>%
#   filter(!is.na(VegMass))
#
# pollen=filter(flowers,Bay=='F') %>%
#   select(Field,Distance,EdgeCent,minShelter,Leafbee:NumLeafbee,Flower,Pollen,Treatment,Honeybee:NumHoneybee) %>%
#   mutate(FieldPlot=paste(Field,Distance,sep='_')) #Unique plot names (for random effect)
#
# rm(bagdat) #Remove distance
#
# #Folder to store figures/graphs:
# figs="C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\Figures"
#
# # Save workspace
# save.image("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\seedfieldData2015.RData")

# Honeybee visitation ---------------------------------------

#PREDICTIONS:
#1) Visits decrease with distance from the colony (TRUE)
  #b) Counts decrease at the same rate (i.e. visits/individual doesn't change) (TRUE)
#2) Visits are lower in male bays (Visit~Distance intercept changes for male/female bays) (FALSE)
#3) If distance-based competition operates differently for pollen/nectar
  #a) Visits~Distance slope would change for male/female bays (FALSE)
  #b) Prop(pollen foragers)~Distance (FALSE - difficult to model because of zeros)
#4) If side-working is somehow more efficient, or related to visit density (FALSE - again, difficult to model)

#Bee visits with distance - not much difference b/w bays
ggplot(hbeesWide,aes(Distance,hbeeVis*12,colour=Bay))+geom_point(position=position_jitter(width=10))+geom_smooth(method='glm.nb',se=T)+labs(x='Distance (m)',y='Honey bee visits/hr')+
  scale_colour_manual(values=c('red','blue'))
#Number of bees with distance - possible difference?
ggplot(hbeesWide,aes(Distance,hbeeCount*12,colour=Bay))+
  geom_point(position=position_jitter(width=10,height=.5))+
  geom_smooth(method='glm.nb',se=T)+labs(x='Distance (m)',y='Honey bee individuals/hr')+
  scale_colour_manual(values=c('red','blue'))

#Honeybee visits with distance from leafcutter shelter -flat for the most part. Possible interaction with beehive distance at 20m, but probably not...
ggplot(hbeesWide,aes(minShelter,hbeeVis,colour=Bay))+
  geom_point(position=position_jitter(width=10,height=.5))+
  geom_smooth(method='glm.nb',se=F)+facet_wrap(~Distance,ncol=1)+
  labs(x='Leafcutter distance (m)',y='Honey bee visits/hr')+
  scale_colour_manual(values=c('red','blue'))
#Again, possible interaction with distance at 20m?
ggplot(hbeesWide,aes(minShelter,hbeeCount,colour=Bay))+
  geom_point(position=position_jitter(width=10,height=.5))+
  geom_smooth(method='glm.nb',se=F)+
  facet_wrap(~Distance,ncol=1)+
  labs(x='Leafcutter distance (m)',y='Honey bee visits/hr')+
  scale_colour_manual(values=c('red','blue'))

#VISITATION WITH DISTANCE
#Total honeybee visits ~ Distance : confirmed in ordinary glm and ZI-glmm
#No real effect of bay, bay:dist, pollination treatments, hive density, or leafcutter shelter distance. Tried regressions at individual beehive distances for hbeeVis~minShelter, but poor fitting occurs


#Using a zero-altered poisson model
pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        #G = list(G1 = list(V = diag(2), nu = 0.002))) #Old prior (unused)
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
          G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process
#Mixed Prior: uses different distributions for Poisson and ZI processes (see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018769.html)

#Distance model is best (deltaDIC=-25.19); Bay (0.68), Bay:Distance (1.72), and leafcutter shelter distance (1.31) don't add anything more. Treatment and Distance:Treatment also don't add anything
#Included ZI terms because according to the Course Notes, both ZI and Poisson terms should be included
hbeeVisDist1=MCMCglmm(hbeeVis~trait*(Distance),
                      rcov=~idh(trait):units,
                      random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                      family='zapoisson',prior=pr,
                      nitt=50000+5000,thin=50,burnin=5000,verbose=F,
                      data=hbeesWide)

pred1=unique(with(hbeesWide,data.frame(Distance,predict(hbeeVisDist1,interval='confidence'))))
p1=ggplot(pred1,aes(x=Distance,y=fit*12))+
  geom_ribbon(aes(ymax=upr*12,ymin=lwr*12),alpha=0.3)+
  geom_line(size=1)+
  geom_point(data=hbeesWide,aes(Distance,hbeeVis*12),position=position_jitter(width=5,height=5),size=2)+
  ylim(0,600)+ylab('Honeybee visits/hr')+xlab('Distance to honeybees (m)')
ggsave(paste(figs,'hbeeVisDist.png',sep='\\'),p1,width=8,height=6)

#Plot of Dist:Treatment

pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        #G = list(G1 = list(V = diag(2), nu = 0.002))) #Old prior (unused)
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1000, alpha.mu=0, alpha.V=1))) #Prior for ZI process

hbeeVisDist2=MCMCglmm(hbeeVis~trait*(Distance+Treatment),
                      rcov=~idh(trait):units,
                      random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                      family='zapoisson',prior=pr,
                      nitt=50000+5000,thin=50,burnin=5000,verbose=F,
                      data=hbeesWide)

pred1a=unique(with(hbeesWide,data.frame(Distance,Treatment,predict(hbeeVisDist2,interval='confidence'))))
p1a=ggplot(pred1a,aes(Distance,fit*12))+
  geom_ribbon(aes(ymax=upr*12,ymin=lwr*12,fill=Treatment),alpha=0.3)+
  geom_line(aes(col=Treatment),size=1)+
  geom_point(data=hbeesWide,aes(Distance,hbeeVis*12,col=Treatment),position=position_jitter(width=5,height=5),size=2)+
  scale_colour_manual(values=c('blue','darkorange','darkgreen'))+scale_fill_manual(values=c('blue','orange','green'))+
  ylim(0,600)+ylab('Honeybee visits/hr')+xlab('Distance (m)')

#Total number of honeybees ~ Distance
#No real effect of Bay, or Treatment. Tried regressions at individual beehive distances for hbeeCount~minShelter, but null model is always best.

pr=list(R = list(V = diag(1), nu = 0.002),
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for poisson process

temp=select(hbeesWide,hbeeCount,Distance,Bay,Field,Treatment,minShelter) %>%
  na.omit() #Strips out NAs
#Again, Count ~ Distance model is the best (deltaDIC=-10.70). Bay is marginally important (-0.47), and Dist:Bay is "significant" (p=0.05), but doesn't improve DIC (1.73). Treatment (0.66) and minShelter (-0.16) don't appear to add anything .
#It may be that the NA values (usually meant we lost count during the survey) could be filled in with estimated counts based on count:visit ratios. Look into this later...
hbeeCountDist1=MCMCglmm(hbeeCount~Distance*Bay,
                        family='poisson',prior=pr,
                        random=~idh(1):Field,
                        burnin=5000,nitt=150000+5000,thin=150,verbose=F,
                        data=temp)

pred2=unique(with(temp,data.frame(Distance,Bay,predict(hbeeCountDist1,interval='confidence'))))
p2=ggplot(pred2,aes(Distance,fit*12,col=Bay))+
  geom_ribbon(aes(ymax=upr*12,ymin=lwr*12,col=NULL,fill=Bay),alpha=0.3)+geom_line(size=1)+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  geom_point(data=temp,aes(Distance,hbeeCount*12),position=position_jitter(width=5,height=3),size=2)+
  ylab('Honeybee count/hr')+xlab('Distance (m)')
ggsave(paste(figs,'hbeeCountDist.png',sep='\\'),p2,width=8,height=6)

#Visits/bee ~ 1
#No real effect of distance, bay, or treatment
temp=select(hbeesWide,hbeeCount,hbeeVis,VisCount,Distance,Bay,Field,Treatment) %>%
  na.omit()

c1=glm.nb(hbeeVis~offset(hbeeCount)+Bay*Distance,data=temp,na.action=na.omit) #None of these terms appear important

pr=list(R = list(V = diag(1), nu = 0.002),
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for poisson process

# pr_var=diag(15)*1e+6
# pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
# pr=list(B=list(mu=matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),15),V=pr_var),
#         R=list(V=1,nu=0.02),
#         G=list(G1=list(V = 1, nu = 1,alpha.mu=0,alpha.V=2000)))

#Distance (deltaDIC=-0.03), bay (1.70), and treatment(1.01) appear to have no effect on Visits/worker
hbeeVisCount1=MCMCglmm(log(VisCount)~Distance*Bay,random=~idh(1):Field,prior=pr,verbose=F,
                       nitt=50000+5000,burnin=5000,thin=50,
                       data=mutate(hbeesWide,VisCount=hbeeVis/hbeeCount))

pred3=unique(with(hbeesWide,data.frame(Distance,Bay,predict(hbeeVisCount1,interval='confidence')))) %>%
  mutate(fit=exp(fit),lwr=exp(lwr),upr=exp(upr))
p3=ggplot(pred3,aes(Distance,fit,col=Bay))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,col=NULL,fill=Bay),alpha=0.3)+geom_line()+
  scale_colour_manual(values=c('red','blue'))+
  geom_point(data=temp,aes(Distance,VisCount),position=position_jitter(width=5))+
  ylab('Visits/honeybee')+xlab('Distance (m)')
ggsave(paste(figs,'hbeeVisCount.png',sep='\\'),p3,width=8,height=6)


hbeeVisCount1a=MCMCglmm(log(VisCount)~Distance*Treatment,random=~idh(1):Field,prior=pr,verbose=F,
                       nitt=50000+5000,burnin=5000,thin=50,
                       data=mutate(hbeesWide,VisCount=hbeeVis/hbeeCount))

pred3a=unique(with(hbeesWide,data.frame(Distance,Treatment,predict(hbeeVisCount1a,interval='confidence')))) %>%
  mutate(fit=exp(fit),lwr=exp(lwr),upr=exp(upr))
p3=ggplot(pred3a,aes(Distance,fit,col=Treatment))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,col=NULL,fill=Treatment),alpha=0.3)+geom_line()+
  scale_colour_manual(values=c('blue','darkorange','darkgreen'))+scale_fill_manual(values=c('blue','orange','green'))+
  geom_point(data=temp,aes(Distance,VisCount),position=position_jitter(width=5))+
  ylab('Visits/honeybee')+xlab('Distance (m)')
ggsave(paste(figs,'hbeeVisCount.png',sep='\\'),p3,width=8,height=6)


#VISITATION WITH DISTANCE AND BEHAVIOUR TYPE
#Visits ~ Distance
#For top-nectar, no difference between visitation in male/female(deltaDIC=-0.29)
#Non-estimable in other behaviour types because of number of zeros

#Plots of visitation types
temp=group_by(hbees,Behaviour,Bay) %>%
  summarize(Visits=sum(Visits,na.rm=T))%>%
  ungroup()%>%
  mutate(Behaviour=factor(Behaviour,labels=c('Side \n Nectar','Top \n Nectar','Top \n Pollen'))) %>%
  mutate(Visits=100*Visits/sum(Visits))

#Number of visits/counts, broken down by bay and behaviour type
p4a=ggplot(temp,aes(Behaviour,Visits,fill=Bay))+geom_bar(stat='identity',position="dodge")+
  scale_fill_manual(values=c('red','blue'))+
  theme(axis.title.x=element_blank())+
  ylab('Percent visits')
ggsave(paste(figs,'hbeePercVis.png',sep='\\'),p4a,width=8,height=6)

#Only estimable for Top-working nectar foragers, and Side nectar foragers in M bay
temp=filter(hbees,Behaviour=='Top_Nectar')
pr=list(R = list(V = diag(2), nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
        G2=list(V=1, nu=1000, alpha.mu=0, alpha.V=1))) #Prior for ZI process

hbeeVisBehav1=MCMCglmm(Visits~trait*(Distance+Bay),
                      rcov=~idh(trait):units,
                      random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                      family='zapoisson',prior=pr,
                      nitt=150000+5000,thin=150,burnin=5000,verbose=F,
                      data=temp)
pred4=unique(with(temp,data.frame(Distance,Bay,predict(hbeeVisBehav1,interval='confidence')*12,Behaviour='Top_Nectar'))) #Predictions for visits/hr
p4=ggplot(hbees,aes(Distance,Visits*12,colour=Bay))+
  geom_point(position=position_jitter(width=5))+
  geom_smooth(data=filter(hbees,Behaviour!='Top_Nectar'),method='glm.nb',se=F)+
  facet_wrap(~Behaviour,ncol=3)+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  ylab('Honeybee visits/hr')+
  geom_ribbon(data=pred4,aes(x=Distance,y=NULL,ymax=upr,ymin=lwr,fill=Bay,col=NULL),alpha=0.3)+
  geom_line(data=pred4,aes(x=Distance,y=fit),size=1)
ggsave(paste(figs,'hbeeVisBehav1.png',sep='\\'),p4,width=8,height=6)

#VISITATION WITH DISTANCE AND BAY POSITION (F BAYS ONLY at 5 and 400m)
temp=group_by(hbees,Behaviour,EdgeCent,Bay) %>%
  filter(Bay=='F',Distance!=20,Distance!=100) %>%
  mutate(Distance=factor(Distance,labels=c('Near','Far')))

#Number of visits/counts, broken down by bay position and behaviour type -  Nectar foragers common across F bay
# Side-workers and pollen foragers likely "spillover" from M bay
group_by(hbees,Behaviour,EdgeCent,Bay) %>%
  filter(Bay=='F',Distance!=20,Distance!=100) %>%
  mutate(Distance=factor(Distance,labels=c('Near','Far'))) %>%
  summarize(Visits=sum(Visits,na.rm=T),Counts=sum(Count,na.rm=T)) %>%
  gather("type","n",Visits:Counts) %>%
  ggplot(aes(Behaviour,n,fill=EdgeCent))+geom_bar(stat='identity',position="dodge")+
  facet_wrap(~type,ncol=1,scale='free_y')+
  scale_fill_manual(values=c('purple','green'))+
  theme(axis.title.x=element_blank())+
  labs(y='Number',x=NULL,fill='Position')

#Number of visits, broken down by bay position across distance. ZA model used because of weird numerical problems in poisson models.
#Visits ~ Distance + EdgeCent (kind of). No effect of EdgeCent(deltaDIC=2.53) or EdgeCent:Distance(4.17), but ZA term is significant for Distance and Poisson term is significant for EdgeCent.
#My interpretation: less likely to see honeybees at the center of the field, but when you do, they visit the same number of times as the edge. No less likely to see honeybees at the edge of the bay, but they visit in lower numbers.
#Difficult to display p-values because ZI terms have different significance than poisson terms.

temp=ungroup(hbeesWide) %>%
  filter(Bay=='F',Distance!=20,Distance!=100) %>%
  mutate(Distance=factor(Distance,labels=c('Near','Far')))

#ZA version
pr=list(R = list(V = 1, nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process


hbeeVisFbay1=MCMCglmm(hbeeVis~trait*(Distance+EdgeCent),
          rcov=~trait:units,
          random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
          family='zapoisson',prior=pr,
          nitt=150000+5000,thin=150,burnin=5000,verbose=F,
          data=temp)
pred=unique(with(temp,data.frame(Distance,EdgeCent,predict(hbeeVisFbay1,interval='confidence')*12))) #Predictions for visits/hr

p5=ggplot(temp,aes(Distance,hbeeVis*12,col=EdgeCent))+
  geom_point(position=position_jitterdodge(jitter.width=0.3,dodge.width=0.5),size=2)+
  scale_colour_manual(values=c('purple','forestgreen'))+
  geom_pointrange(data=pred,aes(x=Distance,y=fit,ymax=upr,ymin=lwr,col=NULL,group=EdgeCent),
                  position=position_dodge(width=.5),size=1)+
  labs(col='Position',fill='Position',y='Honeybee visits/hr')
ggsave(paste(figs,'hbeeVisFbay.png',sep='\\'),p5,width=8,height=6)

#Proportional model (i.e. P(Behaviour)~Bay) - not the best, considering the number of 0 visits.
# # Removes 0-visit observations (can't calculate prop with 0 as denominator)
# temp=select(rowwise(hbeesWide),TopPolVis,TopNecVis,SideNecVis,Field,Distance,EdgeCent,Bay) %>%
#   mutate(total=sum(TopPolVis,TopNecVis,SideNecVis)) %>%
#   filter(total>0) %>%
#   mutate(propTP=TopPolVis/total,propTN=TopNecVis/total,propSN=SideNecVis/total) %>%
#   select(Field:Bay,propTP:propSN) %>%
#   gather(behaviour,visits,-Field:-Bay)
#
# # Prior (as recommended by Course Notes/Tutorial for multinomial models)
# k=3 # categorical outcomes (trinomial)
# I=diag(k-1) #Identity matrix
# J=matrix(rep(1, (k-1)^2), c(k-1, k-1)) #Matrix of 1s
# IJ=.5*(I+J) #Identity matrix with 0.5 on off-diagonals
# pr = list(R = list(V=IJ, nu=1,fix=1), #Fixes residual variance (see pg.9 of MCMCglmm Tutorial)
#           G = list(G1=list(V=1,nu=4)))
#
# #Appears that full model (trait*Bay-1) is best,
# behavDist1=MCMCglmm(cbind(TopPolVis,TopNecVis,SideNecVis)~trait*Bay-1,
#                     rcov=~us(trait):units,
#                     nitt=2e5+3000,thin=400,burnin=3000,
#                     prior=pr,random=~Field,verbose=F,
#                     family='multinomial3',
#                     data=as.data.frame(hbeesWide))
#
# Delta=cbind(c(-1, 1, 0), c(-1, 0, 1)) #Contrast matrix
# c2=((16*sqrt(3))/(15*pi))^2 #Constant for logit transforms (see pg. 52 of Course Notes)
# D=ginv(Delta %*% t(Delta)) %*% Delta #Rescaled contrast matrix
# Int=t(apply(behavDist1$Sol[, 1:2], 1, function(x) {D %*% (x/sqrt(1 + c2 * diag(IJ)))})) #Applies correction to each run
# colnames(Int)=cbind('TopPolVis','TopNecVis','SideNecVis') #Changes column names
# Int=mcmc(exp(Int)/rowSums(exp(Int))) #Transforms to probability scale, and converts to MCMC object
# #Male bays only:
# Int_M=t(apply(with(behavDist1,cbind(rowSums(Sol[,c(1,3)]),rowSums(Sol[,c(2,3,4)]))), 1, function(x) {D %*% (x/sqrt(1 + c2 * diag(IJ)))})) #Applies correction to each run
# colnames(Int_M)=cbind('TopPolVis','TopNecVis','SideNecVis') #Changes column names
# Int_M=mcmc(exp(Int_M)/rowSums(exp(Int_M))) #Transforms to probability scale, and converts to MCMC object
#
# #NOTE: THIS DOESN'T APPEAR TO BE ASSIGNING VALUES CORRECTLY (WRONG COLUMN NAMING?)
# pred=data.frame(rbind(HPDinterval(Int),HPDinterval(Int_M)),mode=c(posterior.mode(Int),posterior.mode(Int_M)),
#            behav=rep(varnames(Int),2),bay=rep(c('F','M'),each=3),row.names=NULL)
#
# ggplot(pred,aes(bay,mode,colour=behav))+
#   geom_pointrange(aes(ymax=upper,ymin=lower),position=position_dodge(width=0.3))+
#   labs(y='Probability',x='Bay',colour='Behaviour')+
#   scale_colour_manual(labels=c('Side-Nectar','Top-Nectar','Top-Pollen'),
#                         values=c("red", "blue", "darkorange"))

# Leafcutter visitation -------------------------------------

#Predictions:
#1) Leafcutter visits decrease with distance from shelters (using nearest shelter, but possibly could use inverse distance or negative exponential) - FALSE
#2a) Leafcutter visits are less at middle of bays
  #2b) This is reduced in areas of pollinator treatments
#3a) In general, Leafcutter ~ -Honeybees
  #3b) Leafcutter visits increase with distance from beehives
#4) Leafcutter visits higher in male bays - FALSE
#5) Leafcutter visits differ among treatments

#lbeeVis~1 - No real effect of shelter distance, either with rank order (p=0.95) or mixed models (p=0.76). Perhaps I need to get closer? No effect of bay (p=0.80). Adding ZI terms doesn't help.
temp=filter(lbees,EdgeCent=='Edge') %>% select(-EdgeCent)

pr=list(R = list(V = 1, nu = 0.002),
        G = list(G1 = list(V = 1, nu = 0.002, alpha.mu=0, alpha.V=1000)))

#lbeeShelter1=MCMCglmm(lbeeVis~minShelter*factor(Distance),family='poisson',data=temp,random=~Field,verbose=F,nitt=100000+5000,burnin=5000,thin=100,prior=pr)

DIC(lbeeShelter1,lbeeShelter2,lbeeShelter3,lbeeShelter4,lbeeShelter5)

pred=with(temp,data.frame(minShelter,Distance=factor(Distance),
                          predict(lbeeShelter1,interval='confidence'))) %>%
  distinct() %>%
  arrange(Distance,minShelter)

p6=ggplot(pred,aes(minShelter,fit*12))+geom_line(aes(col=Distance),size=1)+
  geom_point(data=temp,aes(minShelter,lbeeVis*12,col=factor(Distance)),size=2)+
  geom_ribbon(aes(ymax=upr*12,ymin=lwr*12,fill=factor(Distance)),alpha=0.3)+
  scale_colour_manual(values=c('blue','steelblue','salmon','red'))+
  scale_fill_manual(values=c('blue','steelblue','salmon','red'))+
  labs(x='Distance to nearest shelter (m)',y='Leafcutter Visits/hr',col='Honeybee \n Distance',
       fill='Honeybee \n Distance')
ggsave(paste(figs,'lbeeShelter.png',sep='\\'),p6,width=8,height=6)

ggplot(temp,aes(minShelter,lbeeVis*12,col=factor(Distance)))+geom_point()+geom_smooth(method='glm.nb',se=F)+
  labs(x='Distance to nearest shelter (m)',y='Leafcutter Visits/hr',col='Distance')+
  scale_colour_manual(values=c('blue','steelblue','salmon','red'))

#lbeeVis~1 - No real effect of distance from honeybees(p=0.71), or bay (p=0.9)
temp=filter(lbees,EdgeCent=='Edge') %>% select(-EdgeCent)

lbeeDist2=MCMCglmm(lbeeVis~Distance+Bay,family='poisson',random=~Field,data=temp,nitt=23000,thin=20,verbose=F,pl=T)

pred=with(temp,data.frame(Distance,Bay,predict(lbeeDist2,interval='confidence'))) %>%
  distinct()%>%
  arrange(Bay,Distance)

p7=ggplot(pred,aes(Distance,fit*12))+geom_line(aes(col=Bay),size=1)+
  geom_point(data=lbees,aes(Distance,lbeeVis*12,col=Bay),position=position_jitter(width=5,height=3),size=2)+
  geom_ribbon(aes(ymax=upr*12,ymin=lwr*12,fill=Bay),alpha=0.3)+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  labs(x='Distance to beehives(m)',y='Leafcutter Visits/hr')+ylim(0,600)
ggsave(paste(figs,'lbeeDist.png',sep='\\'),p7,width=8,height=6)

#Number of visits, broken down by bay position across distance
#lbeeVis~Distance*EdgeCenter
temp=filter(lbees,Bay=='F',Distance!=20,Distance!=100) %>%
  select(-Bay)%>%
  mutate(Distance=factor(Distance,labels=c('Near','Far')))

ggplot(temp,aes(Distance,lbeeVis,fill=EdgeCent))+geom_boxplot()+ylim(0,30)+scale_fill_manual(values=c('purple','forestgreen'))+facet_wrap(~Treatment)

#No difference using regular poisson. Data is slightly zero-inflated, but ZA models give the same answer as regular poisson
# pr=list(R = list(V = diag(1), nu = 0.002),
#        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior
# lbeeVisFbay1=MCMCglmm(lbeeVis~trait*(Distance*EdgeCent*Treatment),
#                        rcov=~trait:units,
#                        random=~trait:Field,
#                        family='zapoisson',prior=pr,
#                        nitt=150000+5000,thin=150,burnin=5000,verbose=F,
#                        data=temp)

lbeeVisFbay1=MCMCglmm(lbeeVis~Distance*EdgeCent,random=~Field,family='poisson',
                      verbose=F,data=temp,nitt=55000,burnin=5000,thin=50)

pred=with(temp,data.frame(Distance,EdgeCent,Treatment,
                          predict(lbeeVisFbay1,interval='confidence'))) %>%
  distinct()%>%
  arrange(EdgeCent,Distance,Treatment)

p5=ggplot(temp,aes(Distance,lbeeVis*12,col=EdgeCent))+
  geom_point(position=position_jitterdodge(jitter.width=0.3,dodge.width=0.5),size=2)+
  scale_colour_manual(values=c('purple','forestgreen'))+
  geom_pointrange(data=pred,aes(x=Distance,y=fit*12,ymax=upr*12,ymin=lwr*12,col=NULL,group=EdgeCent),position=position_dodge(width=.5),size=1)+
  labs(col='Position',fill='Position',y='Leafcutter Visits/hr')+ylim(0,500)
ggsave(paste(figs,'lbeeVisFbay.png',sep='\\'),p5,width=8,height=6)

#lbeeVis~Distance*EdgeCenter+Treatment
temp=filter(lbees,Bay=='F',Distance!=20,Distance!=100) %>%
  select(-Bay)%>%
  mutate(Distance=factor(Distance,labels=c('Near','Far')))

ggplot(temp,aes(Distance,lbeeVis*12,fill=EdgeCent))+
  geom_boxplot()+
  scale_fill_manual(values=c('purple','forestgreen'))+facet_wrap(~Treatment,ncol=3)+
  labs(col='Position',y='Leafcutter Visits/hr')+ylim(0,300)

#lbeeVis~-Honeybees - No real effect
ggplot(lbees,aes(hbeeVis,lbeeVis))+geom_point()+geom_smooth(method='lm')
lbeeHbee1=MCMCglmm(lbeeVis~hbeeVis,family='poisson',random=~Field,data=lbees,verbose=F)
pred=with(lbees,data.frame(hbeeVis,predict(lbeeHbee1,interval='confidence'))) %>%
  distinct()%>%
  arrange(hbeeVis)
ggplot(pred,aes(hbeeVis,y=fit,ymax=upr,ymin=lwr))+geom_ribbon(alpha=0.3)+geom_line()+geom_point(data=lbees,aes(x=hbeeVis,y=lbeeVis,ymax=NULL,ymin=NULL))

#lbeeVis~HoneybeeDistance - No real effect
ggplot(lbees,aes(Distance,lbeeVis))+geom_point()+geom_smooth(method='lm')
lbeeHbeeDist1=MCMCglmm(lbeeVis~Distance,family='poisson',random=~Field,data=lbees)
pred=with(lbees,data.frame(Distance,predict(lbeeHbeeDist1,interval='confidence'))) %>%
  distinct()%>%
  arrange(Distance)
ggplot(pred,aes(Distance,y=fit,ymax=upr,ymin=lwr))+geom_ribbon(alpha=0.3)+geom_line()+geom_point(data=lbees,aes(x=Distance,y=lbeeVis,ymax=NULL,ymin=NULL))

#lbeeVis~Treatment - slightly more visits in DoubleTent+DoubleBees
temp=filter(lbees,EdgeCent=='Edge') %>% select(-EdgeCent) %>%
  mutate(Treatment=factor(Treatment,labels=c('Control','Double Tent','Double Tent \n Double Bees')))

#ZA version
pr=list(R = list(V = 1, nu = 0.002, fix=2), #Fixes residual variance for ZI process at 1
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process
lbeeTreatment1=MCMCglmm(lbeeVis~trait*(Treatment),family='zapoisson',data=temp,
                      rcov=~trait:units,
                      random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                      verbose=F,nitt=100000+5000,
                      burnin=5000,thin=100,prior=pr)

pred=with(temp,data.frame(Treatment,predict(lbeeTreatment1,interval='confidence'))) %>%
  distinct()

p=ggplot(pred,aes(Treatment,fit*12))+
  geom_point(data=temp,aes(Treatment,lbeeVis*12),col='darkgrey',position=position_jitter(width=0.1))+
  geom_pointrange(aes(ymax=upr*12,ymin=lwr*12),size=1)+
  labs(y='Leafcutter Visits/hr')+theme(axis.title.x=element_blank())
ggsave(paste(figs,'lbeeVisTreat.png',sep='\\'),p,width=8,height=6)

# IFD-related stuff (for Riley and Ralph) ---------------------------------

#Honeybee visits and IFD (i.e.% of total honeybees vs % of total fls for each plot)
temp=select(survey,Field,Distance,Bay,EdgeCent,minShelter,FlDens,Honeybee) %>%
  filter(EdgeCent=='Edge') %>%
  select(-EdgeCent) %>%
  unite(FlDens_Honeybee,FlDens,Honeybee) %>%
  spread(Bay,FlDens_Honeybee) %>%
  separate(F,c('FlDensF','HoneybeeF'),'_',convert=T) %>%
  separate(M,c('FlDensM','HoneybeeM'),'_',convert=T) %>%
  # mutate(totHbee=HoneybeeF*6+HoneybeeM,totFls=FlDensF*6+FlDensM,
  #        percHbeeF=HoneybeeF*6/totHbee,percHbeeM=HoneybeeM/totHbee,
  #        percFlsF=FlDensF*6/totFls,percFlsM=FlDensM/totFls)
  mutate(totHbee=HoneybeeF+HoneybeeM,totFls=FlDensF+FlDensM,
         percHbeeF=HoneybeeF/totHbee,percHbeeM=HoneybeeM/totHbee,
         percFlsF=FlDensF/totFls,percFlsM=FlDensM/totFls) %>%
  mutate(lDist=cut(minShelter,4)) %>%
  filter(complete.cases(.))

#% visits dependent on honeybee distance
ggplot(temp,aes(percFlsF,percHbeeF))+geom_point()+labs(x='% Female Flowers',y='% Visits to Females')+geom_smooth(method='lm')+geom_abline(intercept=0,slope=1,linetype='dashed')+lims(x=c(0,1),y=c(0,1))+facet_wrap(~Distance)+geom_text(data=data.frame(x=rep(0.1,4),y=rep(0.9,4),label=paste('n=',unname(summary(as.factor(temp$Distance))),sep=''),Distance=unique(temp$Distance)),aes(x=x,y=y,label=label))

#% visits dependent on leafcutter distance
ggplot(temp,aes(percFlsF,percHbeeF))+geom_point()+labs(x='% Female Flowers',y='% Visits to Females')+geom_smooth(method='lm')+geom_abline(intercept=0,slope=1,linetype='dashed')+lims(x=c(0,1),y=c(0,1))+facet_wrap(~lDist)+geom_text(data=data.frame(x=rep(0.1,4),y=rep(0.9,4),label=paste('n=',unname(summary(temp$lDist)),sep=''),lDist=unique(temp$lDist)),aes(x=x,y=y,label=label)) #THIS ISN"T WORKING


a=lm(percHbeeF~percFlsF*minShelter,data=temp)
b=lm(percHbeeF~percFlsF+minShelter,data=temp)
c=lm(percHbeeF~percFlsF,data=temp)
d=lm(percHbeeF~minShelter,data=temp)
e=lm(percHbeeF~1,data=temp)
AIC(a,b,c,d,e) #Distance*PercFlsF interaction, and minShelter*PercFlsF interaction

#Leafcutter visits and IFD (i.e.% of total leafcutters vs % of total fls for each plot)
temp=select(survey,Field,Distance,Bay,EdgeCent,minShelter,FlDens,Leafbee) %>%
  filter(EdgeCent=='Edge') %>%
  select(-EdgeCent) %>%
  unite(FlDens_Leafbee,FlDens,Leafbee) %>%
  spread(Bay,FlDens_Leafbee) %>%
  separate(F,c('FlDensF','LeafbeeF'),'_',convert=T) %>%
  separate(M,c('FlDensM','LeafbeeM'),'_',convert=T) %>%
  mutate(totLbee=LeafbeeF+LeafbeeM,totFls=FlDensF+FlDensM,
         percLbeeF=LeafbeeF/totLbee,percLbeeM=LeafbeeM/totLbee,
         percFlsF=FlDensF/totFls,percFlsM=FlDensM/totFls) %>%
  filter(complete.cases(.))

ggplot(temp,aes(percFlsF,percLbeeF))+geom_point()+labs(x='% Female Flowers',y='% Visits to Females')+geom_smooth(method='lm')+geom_abline(intercept=0,slope=1,linetype='dashed')+lims(x=c(0,1),y=c(0,1))

a=lm(percLbeeF~percFlsF*minShelter,data=temp)
b=lm(percLbeeF~percFlsF+minShelter,data=temp)
c=lm(percLbeeF~percFlsF,data=temp)
d=lm(percLbeeF~percFlsF,data=temp)
AIC(a,b,c,d) #No real effect of hbee distance or leafcutter distance

# Pollination deposition -----------------------------------------------------

#PREDICTIONS:
#1) Pollen deposition decreases with distance from bee colony (TRUE)
#2) Pollen is lower at center of F bay (TRUE)
#3) Pollen is lower away from leafcutter shelters
#4) Pollen changes with pollination treatment (FALSE)

#Pollen~Distance
pr=list(R = list(V = diag(1), nu = 0.002),
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for plot and field random effects
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

polDist1=MCMCglmm(Pollen~Distance,random=~Field+FieldPlot,family='poisson',
                  verbose=F,data=pollen)

pred6=unique(with(pollen,data.frame(Distance,predict(polDist1,interval='confidence'))))
p6=ggplot(pred6,aes(Distance,fit))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(ymax=NULL,ymin=NULL),size=1)+
  geom_point(data=pollen,aes(Distance,Pollen),position=position_jitter(width=5))+
  ylim(0,150)+labs(y='Pollen count',x='Distance to beehives(m)')
ggsave(paste(figs,'polDist.png',sep='\\'),p6,width=8,height=6)

#Pollen~1 - no effect of treatments
temp=mutate(pollen,Treatment=factor(Treatment,labels=c('Control','Double Tent','Double Tent \n Double Bees')))

polTreat1=MCMCglmm(Pollen~Treatment,random=~Field+FieldPlot,family='poisson',verbose=F,data=temp)
pred=unique(with(temp,data.frame(Treatment,predict(polTreat1,interval='confidence'))))

p=ggplot(pred,aes(Treatment,fit))+
  geom_point(data=temp,aes(Treatment,Pollen),col='darkgrey',position=position_jitter(width=0.1))+
  geom_pointrange(aes(ymax=upr,ymin=lwr),size=1)+
  labs(y='Pollen count')+theme(axis.title.x=element_blank())+ylim(0,150)
ggsave(paste(figs,'polTreat.png',sep='\\'),p,width=8,height=6)


#Pollen~Distance + EdgeCent? - weak effect of EdgeCent(p=0.1,DIC=0.5), no effect of EdgeCent:Distance (p=0.52) or Treatment (DIC=0.17)

temp=filter(pollen,Distance!=20,Distance!=100) %>% #Only edge and middle of field
  mutate(Distance=factor(Distance,labels=c('Near','Far'))) #Convert Distance to factor

ggplot(temp,aes(Distance,Pollen,fill=EdgeCent))+geom_boxplot()+ylim(0,150)+
  labs(fill='Position',y='Pollen count')+scale_fill_manual(values=c('purple','green'))+
  facet_wrap(~Treatment)

pr=list(R = list(V = diag(1), nu = 0.02),
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for plot and field random effects
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

polPos1=MCMCglmm(Pollen~EdgeCent+Distance,
                 random=~Field+FieldPlot,family='poisson',
                 prior=pr,nitt=50000+5000,burnin=5000,thin=50,
                 verbose=F,data=temp)

polPos2=MCMCglmm(Pollen~Distance,
                 random=~Field+FieldPlot,family='poisson',
                 prior=pr,nitt=50000+5000,burnin=5000,thin=50,
                 verbose=F,data=temp)
DIC(polPos1,polPos2)

pred7=unique(with(temp,data.frame(Distance,EdgeCent,predict(polPos1,interval='confidence'))))
p7=ggplot(pred7,aes(Distance,fit))+
  scale_colour_manual(values=c('purple','forestgreen'))+
  labs(y='Pollen count',col='Position',fill='Position')+ylim(0,150)+
  geom_point(data=temp,aes(x=Distance,y=Pollen,col=EdgeCent,fill=EdgeCent),position=position_jitterdodge(jitter.width=0.2,dodge.width=0.5),size=2)+
  geom_pointrange(aes(ymax=upr,ymin=lwr,col=NULL,fill=EdgeCent),position=position_dodge(width=.5),size=1)
ggsave(paste(figs,'polDistFbay.png',sep='\\'),p7,width=8,height=6)

#Pollen~lbeeVisits+hbeeVisits

ggplot(pollen,aes(Leafbee,Pollen))+geom_point()



#Pod/seed effects

#PREDICTIONS:
#1) Floral success declines with distance from beehives/leafcutter tents
#2) Floral success is lower at center of middle bay
#3) Seed weight:plant weight declines with distance from beehives/leafcutter tents
#4) Seed weight:plant weight is lower at center of middle bay
#5) #seeds/pod declines with distance from beehives/leafcutter tents
#6) #seeds/pod is lower at center of middle bay
#7) Seed size declines with distance from beehives/leafcutter tents
#8) Seed size is lower at center of middle bay


# Flower success -------------------------------------------------------------

#Pod success~Distance+VegMass -- no effect of Treatment (-1.21) or minShelter (-0.84, p=0.98)
temp=select(plants,Distance,Field,Plant:Missing,minShelter,EdgeCent,Treatment) %>%
  filter(EdgeCent=='Edge',!is.na(VegMass)) %>% #Removes center plots
  unite(FieldPlot,Field,Distance,EdgeCent,remove=F)

pr=list(R=list(V=1,nu=0.2),
       G=list(G1=list(V = diag(1), nu = 1),G2=list(V = diag(1), nu = 1))) #Priors
       #G=list(G1=list(V = diag(1), nu = 1))) #Priors
#flWinsDist1=MCMCglmm(cbind(Pods,Missing)~Distance,random=~Field+FieldPlot,

flWinsDist1=MCMCglmm(cbind(Pods,Missing)~VegMass+Distance,random=~Field+FieldPlot,
                     family='multinomial2',data=temp,verbose=F,prior=pr)
flWinsDist3=MCMCglmm(cbind(Pods,Missing)~poly(VegMass,2,raw=TRUE)+Distance,random=~Field+FieldPlot,
                     family='multinomial2',data=temp,verbose=F,prior=pr,nitt=50000+5000,burnin=5000,thin=50)

varComp(flWinsDist1,units=T) #Field=42% (25-61), Plot=28% (19-47)

#Can't figure out how to do partial plots, so doing this instead:
flWinsDist1a=MCMCglmm(cbind(Pods,Missing)~log(VegMass),random=~Field+FieldPlot,
                     family='multinomial2',data=temp,verbose=F,prior=pr)
flWinsDist1b=MCMCglmm(cbind(Pods,Missing)~Distance,random=~Field+FieldPlot,
                     family='multinomial2',data=temp,verbose=F,prior=pr)

pred=predict(flWinsDist1b,interval='confidence',marginal=~Field+FieldPlot) %>% #Distance effects
  data.frame(.) %>%
  mutate(Distance=temp$Distance,total=temp$Pods+temp$Missing) %>%
  mutate(fitprob=fit/total,lwrprob=lwr/total,uprprob=upr/total)%>%
  select(Distance,fitprob,lwrprob,uprprob) %>%
  round(.,4) %>%
  distinct()
p8b=ggplot(temp,aes(Distance,Pods/(Pods+Missing)))+
  geom_point()+
  geom_ribbon(data=pred,aes(x=Distance,y=NULL,ymax=uprprob,ymin=lwrprob),alpha=0.3)+
  geom_line(data=pred,aes(x=Distance,y=fitprob))+
  labs(x='Distance to Honeybees (m)',y='Proportion flower success')
ggsave(paste(figs,'flWinsDist.png',sep='\\'),p8b,width=8,height=6)

pred=predict(flWinsDist1a,interval='confidence',marginal=~Field+FieldPlot) %>% #VegMass effects
  data.frame(.) %>%
  mutate(VegMass=temp$VegMass,total=temp$Pods+temp$Missing) %>%
  mutate(fitprob=fit/total,lwrprob=lwr/total,uprprob=upr/total)%>%
  select(VegMass,fitprob,lwrprob,uprprob) %>%
  round(.,4) %>%
  distinct()
p8a=ggplot(temp,aes(VegMass,Pods/(Pods+Missing)))+
  geom_point()+
  geom_ribbon(data=pred,aes(x=VegMass,y=NULL,ymax=uprprob,ymin=lwrprob),alpha=0.3)+
  geom_line(data=pred,aes(x=VegMass,y=fitprob))+scale_x_log10()+
  labs(x='Plant Size (g)',y='Proportion flower success')
ggsave(paste(figs,'flWinsVeg.png',sep='\\'),p8a,width=8,height=6)




#Pod Success~1 - No effect of leafcutter distance...

temp=select(plants,Distance,Field,Plant:Missing,minShelter,EdgeCent,Treatment) %>%
  filter(EdgeCent=='Edge',!is.na(VegMass)) %>% #Removes center plots
  unite(FieldPlot,Field,Distance,EdgeCent,remove=F)

p9=ggplot(temp,aes(minShelter,Pods/(Pods+Missing)))+geom_point()+geom_smooth(method='glm',method.args=list(family='binomial'),se=T)+
  labs(x='Leafcutter Distance(m)',y='Proportion flower success')
ggsave(paste(figs,'flWinsShelter.png',sep='\\'),p9,width=8,height=6)

pr=list(R=list(V=1,fix=1),
        G=list(G1=list(V = diag(1), nu = 1),G2=list(V = diag(1), nu = 1))) #Priors
flWinsShelter1=MCMCglmm(cbind(Pods,Missing)~minShelter,random=~Field+FieldPlot,
                     family='multinomial2',data=temp,verbose=F,prior=pr,pr=T)
varComp(flWinsShelter1,units=F) #Field = 20-59%, Plot = 21-53%, Plant = 15-39%

#Pod success~Distance+EdgeCenter - Edge is weakly significant (p=0.1), but Distance is very significant (p<0.001)
temp=filter(plants,Distance!=20,Distance!=100) %>%
  select(Distance,Field,FieldPlot,Plant:Missing,minShelter,EdgeCent,Treatment) %>%
  mutate(Distance=factor(Distance,labels=c('Near','Far'))) #Convert distance to factor

avg=mutate(temp,prop=Pods/(Pods+Missing)) %>%
  group_by(Distance,EdgeCent)%>%
  summarize(mean=mean(prop),sd=sd(prop))

p10=ggplot(temp,aes(Distance,Pods/(Pods+Missing),col=EdgeCent))+
  geom_point(position=position_jitterdodge(jitter.width=0.3,dodge.width=0.5),size=2)+
  labs(x='Distance (m)',y='Proportion flower success',col='Bay position')+
  scale_colour_manual(values=c('purple','forestgreen'))+
  geom_pointrange(data=avg,aes(y=mean,ymin=mean-(sd*1.96),ymax=mean+(sd*1.96),group=EdgeCent),position=position_dodge(width=0.5),col='black',size=1)
ggsave(paste(figs,'flWinsFbay.png',sep='\\'),p10,width=8,height=6)


pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1))) #Priors
flWinsFbay1=MCMCglmm(cbind(Pods,Missing)~Distance*EdgeCent,random=~Field,
                        family='multinomial2',data=temp,verbose=F,prior=pr,pr=T)
varComp(flWinsFbay1) #Fieldplot is singular?? Can't estimate variance properly. Using Field for now.

pred=predict(flWinsFbay1)


# Seed weight : veg weight ------------------------------------------------

#Seed weight~Plant weight * Distance?
temp=select(plants,Distance,Field,FieldPlot,Plant:Missing,minShelter,EdgeCent) %>%
  mutate(seedVeg=SeedMass/VegMass) %>%
  filter(EdgeCent=='Edge') #Removes center plots

#Seed weight~plant weight *distance - strong distance:plant weight interaction
pr=list(R=list(V=diag(1),nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1), G2=list(V = diag(1), nu = 1))) #Priors
seedVegDist1=MCMCglmm(SeedMass~VegMass*as.factor(Distance),family='gaussian',
                      random=~Field+FieldPlot,data=temp,prior=pr,verbose=F,
                      nitt=30000+3000,thin=30,burnin=3000)
varComp(seedVegDist1)

pred=unique(with(temp,data.frame(Distance=as.factor(Distance),VegMass,predict(seedVegDist1,interval='confidence'))))

p11=ggplot(pred,aes(VegMass,fit,col=Distance))+geom_line()+
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=Distance,col=NULL),alpha=0.3)+
  geom_point(data=temp,aes(y=SeedMass,col=as.factor(Distance)))+
  labs(x='Veg Mass(g)',y='Seed Mass(g)',col='Distance')+geom_smooth(method='lm',se=F)+
  scale_fill_manual(values=c('blue','steelblue','salmon','red'))+
  scale_colour_manual(values=c('blue','steelblue','salmon','red'))
ggsave(paste(figs,'seedVegDist.png',sep='\\'),p11,width=8,height=6)

#SeedMass~VegMass*Distance+minShelter - weak minShelter:VegMass interaction (p=0.07, i.e. vegmass is lower towards shelters), but DIC suggest simpler model is better (DIC=-0.88).

temp=select(plants,Distance,Field,FieldPlot,Plant:Missing,minShelter,EdgeCent) %>%
  mutate(seedVeg=SeedMass/VegMass) %>%
  filter(EdgeCent=='Edge') #Removes center plots

pr=list(R=list(V=diag(1),nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1), G2=list(V = diag(1), nu = 1))) #Priors
seedVegShelter2=MCMCglmm(SeedMass~VegMass+VegMass:as.factor(Distance)+minShelter-1,family='gaussian',
                      random=~Field+FieldPlot,data=temp,prior=pr,verbose=F,
                      nitt=30000+3000,thin=30,burnin=3000) #Actual model

pred=unique(with(temp,data.frame(VegMass,Distance=as.factor(Distance),minShelter=cut(minShelter,4),predict(seedVegShelter2,interval='confidence'))))

#Seed:veg relationship with bee distance and lbee distance
ggplot(pred,aes(VegMass,fit,col=Distance))+
  facet_wrap(~minShelter)+
  geom_line()+
  geom_point(data=mutate(temp,minShelter=cut(minShelter,4)),
             aes(x=VegMass,y=SeedMass,col=as.factor(Distance)))+
  geom_line(size=1)+
  labs(x='Veg mass(g)',y='Seed mass(g)',col='Distance')+
  scale_colour_manual(values=c('blue','steelblue','salmon','red'))

seedVegShelter3=MCMCglmm(seedVeg~as.factor(Distance)+minShelter-1,family='gaussian',
                         random=~Field+FieldPlot,data=temp,prior=pr,verbose=F,
                         nitt=30000+3000,thin=30,burnin=3000) #Used to display relationship

pred=unique(with(temp,data.frame(Distance=as.factor(Distance),minShelter,predict(seedVegShelter3,interval='confidence'))))

p12=ggplot(pred,aes(minShelter,fit,col=Distance))+
  #geom_ribbon(aes(ymin=lwr,ymax=upr,col=NULL,fill=Distance),alpha=0.3)+
  geom_point(data=temp,aes(x=minShelter,y=SeedMass/VegMass,col=as.factor(Distance)))+
  geom_line(size=1)+
  labs(x='Leafcutter Distance(m)',y='Seed:Veg',col='Distance')+
  #scale_fill_manual(values=c('blue','steelblue','salmon','red'))+
  scale_colour_manual(values=c('blue','steelblue','salmon','red'))
ggsave(paste(figs,'seedVegShelter.png',sep='\\'),p12,width=8,height=6)


#Seed weight~Plant weight * Distance + EdgeCenter + PlantWeight:EdgeCenter
temp=select(plants,Distance,Field,FieldPlot,Plant:Missing,minShelter,EdgeCent) %>%
  filter(Distance!=20,Distance!=100) %>% #Only edge and middle of field
  mutate(Distance=factor(Distance,labels=c('Near','Far'))) #Convert Distance to factor

pr=list(R=list(V=diag(1),nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1), G2=list(V = diag(1), nu = 1))) #Priors
seedVegFbay1=MCMCglmm(SeedMass~VegMass*Distance+EdgeCent+VegMass:EdgeCent,
                        random=~Field+FieldPlot,family='gaussian',data=temp,
                        nitt=30000+3000,burnin=3000,thin=30,verbose=F,prior=pr)

pred=unique(with(temp,data.frame(Distance,EdgeCent,VegMass,predict(seedPlantFbay1,interval='confidence')))) %>%unite(DistEdge,Distance,EdgeCent,remove=F)

p13=ggplot(pred,aes(x=VegMass,y=fit,ymax=upr,ymin=lwr,group=DistEdge))+
  geom_line(aes(linetype=Distance,col=EdgeCent),size=1)+
  geom_ribbon(aes(fill=EdgeCent),alpha=0.3)+
  geom_point(data=temp,aes(x=VegMass,y=SeedMass,col=EdgeCent,shape=Distance,ymax=NULL,ymin=NULL,group=NULL))+
  scale_colour_manual(values=c('purple','forestgreen'))+
  scale_fill_manual(values=c('purple','forestgreen'))+
  labs(x='Veg Mass(g)',y='Seed Mass(g)')+scale_shape_manual(values=c(19,1))
ggsave(paste(figs,'seedVegFbay.png',sep='\\'),p13,width=8,height=6)


#Less obvious for leafcutters
p14=ggplot(temp,aes(minShelter,SeedMass/VegMass,col=EdgeCent,linetype=Distance))+
  geom_point(aes(shape=Distance))+geom_smooth(method='lm',se=F)+
  labs(x='Leafcutter Distance (m)',y='Seed:Veg',col='Bay position')+
  scale_colour_manual(values=c('purple','forestgreen'))+scale_shape_manual(values=c(19,1))
ggsave(paste(figs,'seedPlantShelterFbay.png',sep='\\'),p14,width=8,height=6)


# Seeds/pod ---------------------------------------------------------------

#PodCount~Distance + Treatment -- No effect of Shelter distance (p=0.85), Shelter Dist:Honeybee Dist (p=0.96), or Treatment:Dist (p=0.15 & 0.56)
temp=filter(seeds,!is.na(PodCount),!is.na(PodMass),EdgeCent=='Edge')
pr=list(R=list(V=diag(1),nu=0.02),
        G=list(G1=list(V = diag(1),nu=0.02),
               G2=list(V = diag(1), nu =0.02, alpha.mu=0, alpha.V=1000),
               G3=list(V = diag(1), nu =0.02, alpha.mu=0, alpha.V=1000))) #Priors
seedsPod1=MCMCglmm(PodCount~Treatment+Distance,family='poisson',random=~Field+FieldPlot+Plant,data=temp,
                   prior=pr,nitt=50000+5000,burnin=5000,thin=50,verbose=F)

pred=unique(with(temp,data.frame(Distance,Treatment,predict(seedsPod1,interval='confidence'))))

p15=ggplot(pred,aes(x=Distance,y=fit,colour=Treatment,fill=Treatment,ymax=upr,ymin=lwr))+
  geom_line(size=1)+geom_ribbon(aes(colour=NULL),alpha=0.3)+
  scale_colour_manual(values=c('blue','darkorange','darkgreen'))+scale_fill_manual(values=c('blue','orange','green'))+
  geom_point(data=temp,aes(x=Distance,y=PodCount,ymax=NULL,ymin=NULL),position=position_jitter(),size=1)+
  labs(x='Distance from Honeybees (m)',y='Seeds per pod')
ggsave(paste(figs,'seedsPodDistanceTreatment.png',sep='\\'),p15,width=8,height=6)

#Podcount~Edge*Distance -- no effect of Distance:Treatment (0.89), Treatment:EdgeCent (0.51) or Distance:EdgeCent (0.72). Treatment doesn't improve DIC (0.56), but Intercept for Double Tent+Bees is marginally different (p=0.10).
temp=filter(seeds,Distance!=20,Distance!=100,!is.na(PodCount),PodCount<40) %>%
  select(Distance,Field,FieldPlot,Plant,PodCount,minShelter,EdgeCent,Treatment) %>%
  mutate(Distance=factor(Distance,labels=c('Near','Far'))) #Convert distance to factor

pr=list(R=list(V=diag(1),nu=0.02),
        G=list(G1=list(V = diag(1),nu=0.02),
               G2=list(V = diag(1), nu =0.02, alpha.mu=0, alpha.V=1000),
               G3=list(V = diag(1), nu =0.02, alpha.mu=0, alpha.V=1000))) #Priors

seedsPodFbay1=MCMCglmm(PodCount~Distance+Treatment+EdgeCent,family='poisson',
                       random=~Field+FieldPlot+Plant,data=temp,
                       prior=pr,nitt=50000+5000,burnin=5000,thin=50,verbose=F)

pred=unique(with(temp,data.frame(Distance,EdgeCent,Treatment,
                                 predict(seedsPodFbay1,interval='confidence'))))

p16=ggplot(pred,aes(x=Distance,y=fit,colour=EdgeCent))+
  geom_pointrange(aes(ymax=upr,ymin=lwr,group=EdgeCent),position=position_dodge(width=.5),size=1)+
  scale_colour_manual(values=c('purple','forestgreen'))+
  labs(col='Position',colour='Position',y='Seeds per pod')+
  facet_wrap(~Treatment)+
  geom_point(data=temp,aes(x=Distance,y=PodCount,ymax=NULL,ymin=NULL,group=EdgeCent),
             size=1,position=position_jitterdodge(jitter.width=0.2,dodge.width=0.5))

ggsave(paste(figs,'seedsPodFbayTreatment.png',sep='\\'),p16,width=8,height=6)

# Seed size ---------------------------------------------------------------



# Seed weight -------------------------------------------------------------


