#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN COMMODITY CANOLA FIELDS (2014+2015)

# Libraries and ggplot theme ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

#Big-text theme, no grid lines (used for Bayer 2016 presentation)
prestheme=theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15), 
                axis.title=element_text(size=20),
                title=element_text(size=20),
                panel.grid.major=element_line(size=0.5,colour='black',linetype='dotted'),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(size=1,colour='black'),
                strip.text=element_text(size=15))
theme_set(theme_bw()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

load("~/Projects/UofC/canola_yield_project/Commodity field analysis/commodityfieldDataAll.RData")
rm(AICc,brix2mgul,deltaAIC,DIC,plotFixedTerms,predFixef,se,varComp,zeroCheck,conversion)
setwd('~/Projects/UofC/canola_yield_project')

# Pod-level model ------------------------------------------

temp=filter(seedsAll)

# ggplot(temp,aes(PodCount))+geom_histogram(bins=45)+labs(x='Seeds per pod',y='Count')
# ggplot(temp,aes(PodMass))+geom_histogram(bins=45)+labs(x='Weight of pod',y='Count')+xlim(0,0.2)
# ggplot(temp,aes(PodCount,PodMass))+geom_point()+labs(x='Seeds per pod',y='Weight of pod',y='Count')
# ggplot(temp,aes(PodCount,PodMass/PodCount))+geom_point()+labs(x='Seeds per pod',y='Weight of single seed',y='Count')+geom_smooth(method='lm')+ylim(0,0.007)

#Pod count model
# datalistPod=with(temp[!is.na(temp$PodCount)&temp$PodCount>0,],list(
#   Npod=length(PodCount), #Length of unique seed counts
#   Nunique=length(unique(PodCount)), #Unique seed counts
#   uniquePodCount=sort(unique(PodCount)),
#   SeedCount=PodCount
#   ))

# #test
# invlogit <- function(x) exp(x)/(1+exp(x))
# logit <- function(x) log(x/(1-x))
# 
# Novules <- 100
# Nseeds <- 20
# Nfert <- c(Nseeds:Novules)
# fertProp <- seq(0.1,0.9,0.1)
# seedProp <- seq(0.1,0.9,0.1)
# prob1 <- array(rep(0,length(fertProp)*length(seedProp)),dim=c(length(fertProp),length(seedProp)))
# 
# for(m in 1:length(fertProp)){
#   for(n in 1:length(seedProp)){
#     for(i in 1:length(Nfert)){
#       prob1[m,n] <- dbinom(Nfert[i],Novules,fertProp[m],log=T)+dbinom(Nseeds,Nfert[i],seedProp[n],log=T) + prob1[m,n]
#     }
#   }
# }
# 
# contour(x=fertProp,y=seedProp,z=prob1,xlab='FertProp',ylab='SeedProp')

# mod1=jags(data=datalistPod,inits=startlist,
#           parameters.to.save=c('int.Pol','slope.Fert','fit','fit.new'),
#           model.file='pod_level.txt',
#           n.chains=3,n.adapt=1000,n.iter=11000,n.burnin=1000,n.thin=10,parallel=T)
# summary(mod1)
# xyplot(mod1)
# traceplot(mod1)
# pp.check(mod1,'fit','fit.new') #Fit is good.

#Pod weight-count modelling
library(jagsUI)

#Pod count model
datalistPod=with(temp[!is.na(temp$PodCount)&temp$PodMass>0,],list(
  Npod=length(PodCount),
  SeedCount=PodCount,
  PodMass=PodMass
))

startlist=function() list(
  int.Pol= rnorm(1,-3.5,0.3),
  slope.Fert = rnorm(1,0.1,0.3),
  int.PodMass = rnorm(1,0,1), 
  slope.PodMass = rnorm(1,0,1), 
  prec.PodMass = rgamma(1,0.01,0.01) 
)

mod2=jags(data=datalistPod,inits=startlist,
          parameters.to.save=c(
            'int.Pol','slope.Fert', #Seed count variables
            'int.PodMass','slope.PodMass','prec.PodMass', #Seed weight-count variables
            'fit.SeedCount','fit.SeedCount.new', #Post. pred checks
            'fit.PodMass','fit.PodMass.new'),
          model.file='pod_count_weight.txt',
          n.chains=3,n.adapt=1000,n.iter=2000,n.thin=2,parallel=T)

summary(mod2)
print(mod2)
xyplot(mod2) #Good convergence
# traceplot(mod2)

#Test
par(mfrow=c(3,1))
simFertCount <- rpois(datalistPod$Npod,120*exp(mod2$q50$int.Pol)) #Simulate fert ov count
simSeedCount <- rpois(datalistPod$Npod,(exp(mod2$q50$slope.Fert*simFertCount)*simFertCount)+0.01) #Simulate seed count
hist(simFertCount,main=NULL,xlab='Fert Ov Count',xlim=c(0,60),breaks=30) 
hist(simSeedCount,main=NULL,xlab='Sim Seed Count',xlim=c(0,60),breaks=30)
hist(datalistPod$SeedCount,main=NULL,xlab='Actual Seed Count',xlim=c(0,60),breaks=30)

pp.check(mod2,'fit.SeedCount','fit.SeedCount.new') #Good
pp.check(mod2,'fit.PodMass','fit.PodMass.new') #Weird. Very close to 1:1, but offset by constant value. Does the variance change with podcount too?

par(mfrow=c(2,1))
hist(temp$PodMass,breaks=50,xlim=c(0,0.7),main='Actual Pod Mass',xlab='g')
hist(rlnorm(1000,-3.04,1/sqrt(1.964)),breaks=100,xlim=c(0,0.7),main='Simulated Pod Mass',xlab='g')

detach("package:jagsUI", unload=TRUE)

#Same thing, but trying in Stan for speed purposes

library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

#Rows without missing seed counts/weights
noMiss <- with(temp,rowSums(is.na(cbind(PodCount,PodMass)))==0 & PodCount!=0)

#List structure for Stan

datalist <- with(temp[!is.na(temp$PodCount)&temp$PodCount>0,],list( #Strips NAs and 0 counts
  Npod=length(PodCount), #Number of observed seed counts
  Nunique=length(unique(PodCount)), #Number of unique seed counts
  uniquePodCount=sort(unique(PodCount)), 
  uniqueMatches=match(PodCount,sort(unique(PodCount))), #Matching index
  SeedCount=PodCount
))
str(datalist)

modPodcount = stan(file='pod_level.stan',data=datalist,iter=2000,chains=3,
             control=list(adapt_delta=0.8))

#Check output
pars <- c('intPol','intFert','slopeFert')
print(modPodcount,pars=pars)
traceplot(modPodcount,pars=pars)
pairs(modPodcount,pars=pars)

#Simulate seed count
simcounts <- lapply(extract(modPodcount,pars=pars),median)
invLogit <- function(x) exp(x)/(1+exp(x))
par(mfrow=c(3,1))
simFertCount <- rbinom(datalist$Npod,120,invLogit(simcounts[[1]])) #Simulate fert ov count
simSeedCount <- rbinom(datalist$Npod,simFertCount,invLogit(simcounts[[2]]+simcounts[[3]]*simFertCount)) 
hist(simFertCount,main=NULL,xlab='Fert Ov Count',xlim=c(0,100)) 
hist(simSeedCount,main=NULL,xlab='Sim Seed Count',xlim=c(0,100),breaks=length(unique(simSeedCount)))
hist(datalist$SeedCount,main=NULL,xlab='Actual Seed Count',xlim=c(0,100),breaks=length(unique(datalist$SeedCount)))

#Doesn't look very good. Perhaps try an arrival model of pollen?

#Fit distributions
poisFun <- function(x,dat) -sum(dpois(dat,x,log=T)) #Poisson
logLikPois <- optimize(f=poisFun,lower=0,upper=1000,dat=flowersAll$Pollen[!is.na(flowersAll$Pollen)])
geomFun <- function(x,dat) -sum(dgeom(dat,x,log=T)) #Geometric
logLikGeom <- optimize(f=geomFun,c(0,1),dat=flowersAll$Pollen[!is.na(flowersAll$Pollen)])
negbinFun <- function(x,dat) -sum(dnbinom(dat,mu=x[1],size=x[2],log=T)) #Neg Bin
logLikNB <- optim(c(1,1),negbinFun,dat=flowersAll$Pollen[!is.na(flowersAll$Pollen)])

2*2+(2*logLikNB$value) #AIC for NB - best
2*1+(2*logLikGeom$objective) #AIC for geom - worse
2*1+(2*logLikPois$objective) #AIC for poisson - much worse

par(mfrow=c(4,1)) #Plot results
hist(rnbinom(datalist$Npod,mu=logLikNB[[1]][1],size=logLikNB[[1]][2]),main=NULL,xlab='Neg Bin',xlim=c(0,4000),breaks=seq(0,4000,20))
hist(rgeom(datalist$Npod,logLikGeom[[1]][1]),main=NULL,xlab='Geom',xlim=c(0,4000),breaks=seq(0,4000,20))
hist(rpois(datalist$Npod,logLikPois[[1]][1]),main=NULL,xlab='Pois',xlim=c(0,4000),breaks=seq(0,4000,20))
hist(flowersAll$Pollen[!is.na(flowersAll$Pollen)],xlim=c(0,4000),breaks=seq(0,4000,20),main=NULL,xlab='Actual pollen counts')

# Test model to simulate seed counts --------------------------------------
invLogit <- function(x){ 
  i <- exp(x)/(1+exp(x)) 
  i[is.nan(i)] <- 1 #If x is large enough, becomes Inf/Inf = NaN, but Lim(invLogit) x->Inf = 1
  return(i)
}
#Linear breakpoint function - two lines with intersection "b"
bpoint <- function(x,int1,slope1,b,slope2) ifelse(x<b,int1+slope1*x,b*slope1+(x-b)*slope2+int1)

#Likelihood method: (non-vectorized version)

#p(Seed)*p(Ovule)*p(Pollen)
# prob <- array(NA,c(maxPolNum,maxOvNum,length(seedCount)))
# # prob <- matrix(NA,nrow=maxPolNum,ncol=maxOvNum)
# for(seed in 1:length(seedCount)){ #For each seed count
#   for(polCount in seedCount[seed]:maxPolNum){ #For each potential pollen count
#     for(ovCount in seedCount[seed]:min(maxOvNum,polCount)){ #For each potential ovule count
#       prob[polCount,ovCount,seed] <- dnbinom(polCount,mu=292.9317,size=0.6136358,log=T)+
#         dbinom(ovCount,polCount,0.05,log=T)+
#         dbinom(seedCount[seed],ovCount,0.01,log=T)
#     }
#   }
# }

#Values for pollination process:
#mu=292.9316584
#size=0.6136358

margProb <- function(coefs,dat,maxOvNum,maxPolNum,mu,size){
  
  seedCount <- c(1:52) #Number of seeds to try
  mu <- 292.9316584 #Pollen dist
  size <- 0.6136358
  lambda <- 32 #Ovule dist
  
  maxOvNum <- max(datalist$SeedCount) #Number of potential ovules
  maxPolNum <- qnbinom(.95,mu=mu,size=0.6136358) #Upper number of potential pollen grains (95% of density)
  
  coefs <- c(0.1,-1,0.1) #10% pollen survival, invLogit(-6+0.5*seedNum) pod abortion rate
  
  #Replicate all levels of seed/ovule/pollen counts
  prob <- expand.grid(seedCount=seedCount,ovCount=1:maxOvNum,polCount=1:maxPolNum)
  #Remove impossible categories
  prob <- prob[(prob$ovCount>=prob$seedCount) & (prob$polCount>=prob$seedCount),]
  #Pollen survial prob
  prob$probPol <- coefs[1]
  #Pod survival prob
  prob$probSeed <- invLogit(coefs[2]+coefs[3]*prob$seedCount)
  
  #Calculate log-prob for each
  prob$lp <- dnbinom(prob$polCount,mu=mu,size=size,log=T)+ #Pollen prob
    dpois(prob$ovCount,lambda=lambda,log=T)+ #Ovule prob
    ifelse(prob$seedCount==prob$ovCount, #Fert prob
           #if seeds==ovules, successful pollen count could have been anything between #seeds and #pollen
           pbinom(prob$seedCount-1,prob$polCount,prob$probPol,log=T,lower.tail=F),
           #if seeds<ovules, only one way to get seed count
           dbinom(prob$seedCount,prob$polCount,prob$probPol,log=T))+
    dbinom(1,1,prob$probSeed,log=T) #Pod abort prob
  #Sum by seed counts
  seedLp <- with(prob,tapply(lp,seedCount,function(x) log(sum(exp(x))))) 
  barplot(exp(seedLp),pch=19,ylab='p(x)',xlab='Seed count')
  #p(seeds=1:52)
  sum(exp(seedLp)) #Problem: sum(seedprob)>1 : 0.55
  #p(seeds!=1:52)
  
  dnbinom(0,mu=mu,size=size) #Prob of 0 pollen - 0.02
  
  with(expand.grid(pol=1:maxPolNum,ov=1:maxOvNum),{ #Prob of >1 pollen, but no fertilization - 0.07
    sum(exp(dnbinom(pol,mu=mu,size=size,log=T)+
        dpois(ov,lambda=lambda,log=T)+
        dbinom(0,pol,coefs[1],log=T)))})
  
  
  
  prob$lp_ab <- dnbinom(prob$polCount,mu=mu,size=size,log=T)+ #Pollen prob
    dpois(prob$ovCount,lambda=lambda,log=T)+ #Ovule prob
    ifelse(prob$seedCount==prob$ovCount, #Fert prob
           #if seeds==ovules, successful pollen count could have been anything between #seeds and #pollen
           pbinom(prob$seedCount-1,prob$polCount,prob$probPol,log=T,lower.tail=F),
           #if seeds<ovules, only one way to get seed count
           dbinom(prob$seedCount,prob$polCount,prob$probPol,log=T))+ 
    dbinom(0,1,prob$probSeed,log=T) #Pod abort prob
  
  abLp <- with(prob,tapply(lp_ab,seedCount,function(x) log(sum(exp(x))))) 
  sum(exp(abLp))
  
  
  
  
  
  
  
  
  #This appears to work. Double-check the math and try to figure out p(0 seeds)
  
  #Match -lp to actual counts, sum, return
  return(-sum(seedLp[match(dat,names(seedLp))])) #Match 
}

simSeeds <- function(maxNovules,coefs,plotResults=T,actualSeeds=NA,propDiff=F){
  #Idea: 
  #Pollination is a Negative Binomial process with mu=292.9,size=0.61
  mu <- 292.9316584 
  size <- 0.6136358
  #Ovule number per flower is a Poisson process with lambda =~32. Wang et al 2011 found that there are cross-season difference, but they aren't that big (~30-35 ovules).
  lambda <- 32
  #Fertilization is a Binomial process, where Fert Ovules ~ Bin(Pollen,p)
  #If more ovules are fertilized than available, fert ovules = total ovules (integrate all outcomes < Novules)
  
  #Pod abortion is a Bernoulli process where Abortion ~ Bern(invLogit(int1+slope1*Nfert))
  p <- coefs[1]
  int1 <- coefs[2]
  slope1 <- coefs[3]
  
  set.seed(1)
  #Generate pollen counts for 1000 flowers
  simPol <- rnbinom(1000,mu=mu,size=size)
  #Generate ovule counts for 1000 flowers
  simOv <- rpois(1000,32)
  simOv[simOv>maxNovules] <- maxNovules #Set all ovules > maxNovules to maxNovules
  
  #Ovules are fertilized by pollen
  simFert <- rbinom(1000,simPol,p) #Fixed perc of pollen makes it to ovules
  simFert[simFert>simOv] <- simOv[simFert>simOv] #If Npollen>Novules, reverts to ov number "first come first served"
  
  #Fertilized ovules become seeds
  simSeed <- simFert[simFert>0] #Fertilized (suriving pollen>1) ovules become seeds
  probSeed <- invLogit(int1+slope1*simSeed) #Prob of pod abortion ~ # of fertilized ovules
  simSeed <- simSeed[rbinom(1000,1,probSeed)==1] #Plant aborts pods depending on number of ovules
  
  if(plotResults==T & (length(actualSeeds)!=1 & !is.na(actualSeeds[1]))){
    #Plots
    par(mfrow=c(5,1),mar=c(5,5,2,5))
    #Pollen
    hist(simPol,main=NULL,xlab='Pollen grains',breaks=50,xlim=c(0,max(simPol)))
    par(new=T) #Pollen survival
    curve(invLogit(p+0*x),0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Pollen survival ',col='red')
    #Ovules
    hist(simOv,xlab='Ovules',xlim=c(0,max(simOv)),breaks=seq(0,max(simOv),1),main=NULL)
    #Fertilized ovules
    hist(simFert,main=NULL,xlab='Fertilized Ovules',breaks=seq(0,max(simFert),1),xlim=c(0,max(simOv)))
    par(new=T) #Survival prob
    curve(invLogit(int1+slope1*x),0,max(simFert),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Pod survival',col='red')
    #Seed counts
    hist(simSeed,main=NULL,xlab='Seeds/pod',xlim=c(0,max(simOv)),breaks=seq(0,max(simOv),1))
    #Actual seeds
    hist(actualSeeds,main=NULL,xlab='Actual Seed Count',xlim=c(0,max(simOv)),breaks=length(unique(datalist$SeedCount)))
    par(mfrow=c(1,1),mar= c(5, 4, 4, 2) + 0.1)
  } else if(propDiff==T) {
    uniqueCounts <- sort(unique(c(actualSeeds,simSeed))) #Unique counts (from sim or actual)
    props <- matrix(0,nrow=length(uniqueCounts),ncol=2) #Matrix to store proportions
    propsActual <- table(actualSeeds)/length(actualSeeds)
    propsSim <- table(simSeed)/length(simSeed)
    props[match(names(propsActual),uniqueCounts),1] <- propsActual
    props[match(names(propsSim),uniqueCounts),2] <- propsSim
    return(sum(abs(props[,1]-props[,2])))
  } else {
    return(list(simPol=simPol,simFert=simFert,simSeed=simSeed))
  }
}

simSeeds(maxNovules=50,coefs=c(0.2,-6,0.5),plotResults=T,actualSeeds=datalist$SeedCount)


#Proportion missing seeds. Should be related to 
plantsAll %>% mutate(PropMis=Missing/(Pods+Missing)) %>% filter(!is.na(PropMis),PropMis>0) %>% 
  ggplot(aes(PropMis))+geom_histogram()


# Field-level visitation, nectar, and pollen deposition model -------------

datalistField <- with(fieldsAll,list(
  NFields=nrow(fieldsAll),
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

datalistPlot <- with(filter(surveyAll,Distance!=1),list(
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

datalistFlw <- with(filter(flowersAll,Distance!=1),list(
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
mean(fieldsAll$Temp) #24 C is mean temperature

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
  geom_jitter(data=filter(surveyAll,Distance!=1),aes(Distance,Honeybee),col='black',size=1)+
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
  geom_point(data=flowersAll,aes(x=Honeybee/(FlDens*Time/60),y=Vol,col=NULL))+
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
  geom_point(data=flowersAll,aes(x=Honeybee/(FlDens*Time/60),y=Vol,col=NULL))+
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
  geom_point(data=flowersAll,aes(x=Honeybee,y=Pollen))+
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
ggplot(surveyAll,aes(Distance,Honeybee,col=BeeYard))+geom_jitter()+
  geom_smooth(method='glm.nb',formula=y~log(x),se=F)+
  scale_colour_manual(values=c('orange','green'))+
  ylim(0,20)

filter(surveyAll,Distance!=1,BeeYard=='Stocked'|Honeybee<20) %>% 
  ggplot(aes(Distance,Honeybee))+geom_jitter()+
  facet_wrap(~BeeYard,ncol=1)+
  geom_smooth(method='glm.nb',formula=y~log(x),se=F)+
  scale_colour_manual(values=c('orange','green'))#+
  #ylim(0,20)

mod1ML<-glm.nb(Honeybee~log(Distance)*NumHives,data=surveyAll)
summary(mod1ML)



# Unused code -------------------------------------------------------------

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
# temp=filter(seedsAll,Year==2014,PodMass>0)
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
