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

#Double-linear survival, with number of ovules included
simSeeds <- function(mu,size,Novules,coefs,plotResults=T,actualSeeds=NA,propDiff=F){
  
  intFert <- coefs[1]
  slopeFert <- coefs[2]
  intSeed <- coefs[3]
  slopeSeed <- coefs[4]
  
  set.seed(1)
  #Generate pollen counts for 1000 flowers
  simPol <- rnbinom(5000,mu=mu,size=size)
  #"Abort" flowers with 0 pollen grains
  simPol <- simPol[simPol>0]
  
  #Proportion of pollen germinates and makes it to ovules
  propFert <- invLogit(intFert+slopeFert*simPol) #Proportion surviving
  simFert <- rbinom(length(simPol),simPol,propFert)
  #If more than 120 pollen grains make it, cut off fertilization at max ovules
  simFert[simFert>Novules] <- Novules 
  #If there are zero fertilized ovules, "abort" flowers
  simFert <- simFert[simFert>0]
  
  #Proportion of fertilized ovules become seeds
  propSeed <- invLogit(intSeed+slopeSeed*simFert) #Proportion surviving
  simSeed <- rbinom(length(simFert),simFert,propSeed)
  #If there are zero seeds, "abort" flowers
  simSeed <- simSeed[simSeed>0]
  
  # #Plant aborts flowers according to seed production ('rolls the dice' again)
  # propSurv <- invLogit()
  # simSeed <- simSeed

  if(plotResults==T & (length(actualSeeds)!=1 & !is.na(actualSeeds[1]))){
    #Plots
    par(mfrow=c(4,1),mar=c(5,5,2,5))
    #Pollen
    hist(simPol,main=NULL,xlab='Pollen grains',breaks=50,xlim=c(0,max(simPol)))
    par(new=T) #Survival prob
    curve(invLogit(intFert+slopeFert*x),0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Prob survival',col='red')
    #Ovules
    hist(simFert,main=NULL,xlab='Fertilized ovules',breaks=50,xlim=c(0,Novules))
    par(new=T) #Survival prob
    curve(invLogit(intSeed+slopeSeed*x),
          0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Prob survival',col='red')
    #Seeds
    hist(simSeed,main=NULL,xlab='Seeds',breaks=length(unique(simSeed)),xlim=c(0,Novules))
    hist(actualSeeds,main=NULL,xlab='Actual Seed Count',xlim=c(0,Novules),breaks=length(unique(datalist$SeedCount)))
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

# testcoef <- c(-10,0.5,5,-0.3)
testcoef <- c(0,0.05,40,-.1)
simSeeds(coefs=testcoef,
                 mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=120,plotResults=T,
                 actualSeeds=datalist$SeedCount,propDiff=F) #Works
testcoef <- c(0,0.05,40,-.1)
par(mfrow=c(2,1))
curve(bpoint(x,int1=testcoef[1],slope1=testcoef[2],b=testcoef[3],slope2=testcoef[4]),0,120,ylab='f(x)')
curve(invLogit(bpoint(x,int1=testcoef[1],slope1=testcoef[2],b=testcoef[3],slope2=testcoef[4])),0,120,ylab='p(x)')
par(mfrow=c(1,1))

#This set of parameters works reasonably well:
# c(1,0.2,1,-.02,0)

#Optimize
bestCoefs <- optim(par=c(0,0.05,40,-.1),simSeeds,
      mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=120,plotResults=F,
      actualSeeds=datalist$SeedCount,propDiff=T)

#plot estimate
simSeeds(coefs=bestCoefs$par,
         mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=120,plotResults=T,
         actualSeeds=datalist$SeedCount)
#Still doesn't look very good. 
par(mfrow=c(2,1))
curve(bpoint(x,bestCoefs$par[1],bestCoefs$par[2],bestCoefs$par[3],bestCoefs$par[4]),0,120,ylab='f(x)')
curve(invLogit(bpoint(x,bestCoefs$par[1],bestCoefs$par[2],bestCoefs$par[3],bestCoefs$par[4])),0,120,ylab='p(x)')
par(mfrow=c(1,1))


#Likelihood method:

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

maxPolNum <- 500 #Maximum number of pollen grains (for now)
maxOvNum <- 100 #Maximum number of ovules
seedCount <- sort(unique(datalist$SeedCount)) #Number of seeds

prob <- expand.grid(seedCount=seedCount,ovCount=1:maxOvNum,polCount=1:maxPolNum,lp=NA)
prob <- prob[(prob$ovCount>=prob$seedCount) & (prob$polCount>=prob$ovCount),]

margProb <- function(coefs,dat,maxOvNum,maxPolNum,mu,size){
  int1 <- coefs[1]
  slope1 <- coefs[2]
  int2 <- coefs[3]
  slope2 <- coefs[4]
  
  #Replicate all levels of seed/ovule/pollen counts
  prob <- expand.grid(seedCount=seedCount,ovCount=1:maxOvNum,polCount=1:maxPolNum,lp=NA)
  #Remove impossible categories
  prob <- prob[(prob$ovCount>=prob$seedCount) & (prob$polCount>=prob$ovCount),]
  
  prob$probOv <- invLogit(int1+polCount*slope1)
  prob$probSeed <- invLogit(int2+ovCount*slope2)
  #Calculate log-prob for each
  prob$lp <- dnbinom(prob$polCount,mu=mu,size=size,log=T)+
    dbinom(prob$ovCount,prob$polCount,prob$probOv,log=T)+
    dbinom(prob$seedCount,prob$ovCount,prob$probSeed,log=T)
  #Sum by seed counts
  seedLp <- with(prob,tapply(lp,seedCount,function(x) log(sum(exp(x))))) 
  #Match -lp to actual counts, sum, return
  return(-sum(seedLp[match(dat,names(seedLp))])) #Match 
}

#Breakpoint model with no intermediate fertilized ovules
margProbBP <- function(coefs,dat,maxPolNum,mu,size){
  int1 <- coefs[1]
  slope1 <- coefs[2]
  bp <- coefs[3]
  slope2 <- coefs[4]
  
  #Replicate all levels of seed/ovule/pollen counts
  prob <- expand.grid(seedCount=seedCount,polCount=1:maxPolNum,lp=NA)
  #Remove impossible categories
  prob <- prob[(prob$polCount>=prob$seedCount),]
  
  prob$probSeed <- invLogit(bpoint(prob$polCount,int1,slope1,bp,slope2))
  #Calculate log-prob for each
  prob$lp <- dnbinom(prob$polCount,mu=mu,size=size,log=T)+
    dbinom(prob$seedCount,prob$polCount,prob$probSeed,log=T)
  #Sum by seed counts
  seedLp <- with(prob,tapply(lp,seedCount,function(x) log(sum(exp(x))))) 
  #Match -lp to actual counts, sum, return
  return(-sum(seedLp[match(dat,names(seedLp))])) #Match 
}

simSeedsBp <- function(mu,size,Novules,coefs,plotResults=T,actualSeeds=NA,propDiff=F){
  
  int1 <- coefs[1]
  slope1 <- coefs[2]
  bp <- coefs[3]
  slope2 <- coefs[4]
  
  set.seed(1)
  #Generate pollen counts for 1000 flowers
  simPol <- rnbinom(5000,mu=mu,size=size)
  #"Abort" flowers with 0 pollen grains
  simPol <- simPol[simPol>0]
  
  #Proportion of fertilized ovules become seeds
  propSeed <- invLogit(bpoint(simPol,int1,slope1,bp,slope2)) #Proportion surviving
  simSeed <- rbinom(length(simPol),simPol,propSeed)
  #If there are zero seeds, "abort" flowers
  simSeed <- simSeed[simSeed>0]
  
  if(plotResults==T & (length(actualSeeds)!=1 & !is.na(actualSeeds[1]))){
    #Plots
    par(mfrow=c(3,1),mar=c(5,5,2,5))
    #Pollen
    hist(simPol,main=NULL,xlab='Pollen grains',breaks=50,xlim=c(0,max(simPol)))
    par(new=T) #Survival prob
    curve(invLogit(bpoint(x,int1,slope1,bp,slope2)),0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Prob survival',col='red')
    #Ovules
    # hist(simFert,main=NULL,xlab='Fertilized ovules',breaks=50,xlim=c(0,Novules))
    # par(new=T) #Survival prob
    # curve(invLogit(intSeed+slopeSeed*x),
    #       0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    # axis(side=4,cex.axis=1,col='red',col.axis='red') 
    # mtext(side=4,line=3,'Prob survival',col='red')
    #Seeds
    hist(simSeed,main=NULL,xlab='Seeds',breaks=length(unique(simSeed)),xlim=c(0,Novules))
    hist(actualSeeds,main=NULL,xlab='Actual Seed Count',xlim=c(0,Novules),breaks=length(unique(datalist$SeedCount)))
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

#Breakpoint model with intermediate fertilized ovules - linear model for pollen survival, breakpoint for ovule survival
margProbBP2 <- function(coefs,dat,maxPolNum,mu,size){
  polInt <- coefs[1]
  polSlope <- coefs[2]
  ovInt1 <- coefs[3]
  ovSlope1 <- coefs[4]
  ovBP <- coefs[5]
  ovSlope2 <- coefs[6]
  
  #Replicate all levels of seed/ovule/pollen counts
  prob <- expand.grid(seedCount=seedCount,ovCount=1:maxOvNum,polCount=1:maxPolNum,lp=NA)
  #Remove impossible categories
  prob <- prob[(prob$ovCount>=prob$seedCount) & (prob$polCount>=prob$ovCount),]
  
  prob$probOv <- invLogit(polInt+polSlope*prob$polCount)
  prob$probSeed <- invLogit(bpoint(prob$ovCount,ovInt1,ovSlope1,ovBP,ovSlope2))
  #Calculate log-prob for each
  prob$lp <- dnbinom(prob$polCount,mu=mu,size=size,log=T)+
    dbinom(prob$ovCount,prob$polCount,prob$probOv,log=T)+
    dbinom(prob$seedCount,prob$ovCount,prob$probSeed,log=T)
  #Sum by seed counts
  seedLp <- with(prob,tapply(lp,seedCount,function(x) log(sum(exp(x))))) 
  #Match -lp to actual counts, sum, return
  return(-sum(seedLp[match(dat,names(seedLp))])) #Match 
}

simSeedsBP2 <- function(mu,size,Novules,coefs,plotResults=T,actualSeeds=NA,propDiff=F){
  
  polInt <- coefs[1]
  polSlope <- coefs[2]
  ovInt1 <- coefs[3]
  ovSlope1 <- coefs[4]
  ovBP <- coefs[5]
  ovSlope2 <- coefs[6]
  
  set.seed(1)
  #Generate pollen counts for 1000 flowers
  simPol <- rnbinom(5000,mu=mu,size=size)
  #"Abort" flowers with 0 pollen grains
  simPol <- simPol[simPol>0]
  
  #Proportion of pollen germinates and makes it to ovules
  propFert <- invLogit(polInt+polSlope*simPol) #Proportion surviving
  simFert <- rbinom(length(simPol),simPol,propFert)
  #If more than 120 pollen grains make it, cut off fertilization at max ovules
  simFert[simFert>Novules] <- Novules 
  #If there are zero fertilized ovules, "abort" flowers
  simFert <- simFert[simFert>0]
  
  #Proportion of fertilized ovules become seeds
  propSeed <- invLogit(bpoint(simFert,ovInt1,ovSlope1,ovBP,ovSlope2)) #Proportion surviving
  simSeed <- rbinom(length(simFert),simFert,propSeed)
  #If there are zero seeds, "abort" flowers
  simSeed <- simSeed[simSeed>0]
  
  # #Plant aborts flowers according to seed production ('rolls the dice' again)
  # propSurv <- invLogit()
  # simSeed <- simSeed
  
  if(plotResults==T & (length(actualSeeds)!=1 & !is.na(actualSeeds[1]))){
    #Plots
    par(mfrow=c(4,1),mar=c(5,5,2,5))
    #Pollen
    hist(simPol,main=NULL,xlab='Pollen grains',breaks=50,xlim=c(0,max(simPol)))
    par(new=T) #Survival prob
    curve(invLogit(polInt+polSlope*x),0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Prob survival',col='red')
    #Ovules
    hist(simFert,main=NULL,xlab='Fertilized ovules',breaks=50,xlim=c(0,Novules))
    par(new=T) #Survival prob
    curve(invLogit(bpoint(x,ovInt1,ovSlope1,ovBP,ovSlope2)),
          0,max(simFert),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Prob survival',col='red')
    #Seeds
    hist(simSeed,main=NULL,xlab='Seeds',breaks=length(unique(simSeed)),xlim=c(0,Novules))
    hist(actualSeeds,main=NULL,xlab='Actual Seed Count',xlim=c(0,Novules),breaks=length(unique(datalist$SeedCount)))
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


margProb(coefs=c(-0.475,-0.002,0.004,0.000),dat=datalist$SeedCount,maxOvNum=maxOvNum,maxPolNum=maxPolNum,mu=292.9317,size=0.6136358)
margProbBP(coefs=c(0,0.1,20,-0.1),dat=datalist$SeedCount,maxPolNum=maxPolNum,mu=292.9317,size=0.6136358)
margProbBP2(coefs=c(-0.475,-0.002,-20,1,30,-.3),dat=datalist$SeedCount,maxPolNum=maxPolNum,mu=292.9317,size=0.6136358)

simSeedsBP2(coefs=c(-0.475,-0.002,-20,1,30,-.3),
         mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=maxOvNum,plotResults=T,
         actualSeeds=datalist$SeedCount,propDiff=F) #Works, but not well (looks like seed field distribution)

curve(invLogit(bpoint(x,0,0.01,30,-1)),0,100)



#best coefs for simple model -> 0.155, 0.579, -LL=17563.75
#double-linear model -> -0.475 -0.002  0.004  0.000, -LL=17563.94
bestCoefs <- optim(par=c(0,0,0,0),fn=margProb,
      dat=datalist$SeedCount,maxOvNum=maxOvNum,maxPolNum=maxPolNum,mu=292.9317,size=0.6136358,
      method='BFGS')
      # method="L-BFGS-B",lower=c(-50,0),upper=c(1,1))

simSeeds(coefs=bestCoefs$par,
         mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=maxOvNum,plotResults=T,
         actualSeeds=datalist$SeedCount,propDiff=F) #Works, but not well (looks like seed field distribution)

#Breakpoint model skipping ovules
bestCoefs2 <- optim(par=c(-3,0.1,10,-0.1),fn=margProbBP,
                    dat=datalist$SeedCount,maxPolNum=1500,mu=292.9317,size=0.6136358)

margProbBP(coefs=bestCoefs2$par,dat=datalist$SeedCount,maxPolNum=1500,mu=292.9317,size=0.6136358)

simSeedsBp(coefs=bestCoefs2$par,
         mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=maxOvNum,plotResults=T,
         actualSeeds=datalist$SeedCount,propDiff=F) #Doesn't work properly - same as first model

#Breakpoint model with ovules
# -1.008786435  -0.003902827 -18.149592217   6.460526408  26.304270993  -1.327666477
bestCoefs3 <- optim(par=c(-0.475,-0.002,-20,1,30,-.3),fn=margProbBP2,
                    dat=datalist$SeedCount,maxPolNum=1000,mu=292.9317,size=0.6136358)

simSeedsBP2(#coefs=bestCoefs3$par,
            coefs=c(-1.008786435,-0.003902827,-18.149592217,6.460526408,26.304270993,-1.327666477),
            mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=maxOvNum,plotResults=T,
            actualSeeds=datalist$SeedCount,propDiff=F) #Works, but not well (looks like seed field distribution)


#Linear model with intermediate fertilized ovules - linear for pollen survival, linear for ovule survival, linear for flower survival
margProb3 <- function(coefs,dat,maxOvNum,maxPolNum,mu,size){
  polInt <- coefs[1] #Pollen survival
  polSlope <- coefs[2]
  ovInt <- coefs[3] #Ovule survival
  ovSlope <- coefs[4]
  flwInt <- coefs[5] #Flower survival
  flwSlope <- coefs[6]
  
  #Replicate all levels of seed/ovule/pollen counts
  prob <- data.matrix(expand.grid(seedCount=sort(unique(dat)),ovCount=1:maxOvNum,polCount=1:maxPolNum,
                                  lp=NA,probOv=NA,probSeed=NA,probPod=NA))
  #Remove impossible categories
  prob <- prob[(prob[,'ovCount']>=prob[,'seedCount']) & (prob[,'polCount']>=prob[,'ovCount']),]
  
  prob[,'probOv'] <- invLogit(polInt+polSlope*prob[,'polCount']) #Prob of pollen -> fert ovule
  prob[,'probSeed'] <- invLogit(ovInt+ovSlope*prob[,'ovCount']) #Prob of fert ovule -> seed
  prob[,'probPod'] <- invLogit(flwInt+flwSlope*prob[,'ovCount']) #Prob of pod survival
  
  #Calculate log-prob for each
  prob[,'lp'] <- dnbinom(prob[,'polCount'],mu=mu,size=size,log=T)+ # Pollen count ~ NegBin(mu,size)
    dbinom(prob[,'ovCount'],prob[,'polCount'],prob[,'probOv'],log=T)+ #Ovule count ~ Binomial(Pollen count,probOv)
    dbinom(1,1,prob[,'probPod']) + #Pod survival | Fertilization ~ Bernoulli(propPod)
    dbinom(prob[,'seedCount'],prob[,'ovCount'],prob[,'probSeed'],log=T) #Seed count
    
  #Sum by seed counts
  seedLp <- tapply(prob[,'lp'],prob[,'seedCount'],function(x) log(sum(exp(x)))) 
  #Match -lp to actual counts, sum, return
  return(-sum(seedLp[match(dat,sort(unique(dat)))])) #Match 
}

simSeeds3 <- function(mu,size,Novules,coefs,plotResults=T,actualSeeds=NA,propDiff=F){
  
  polInt <- coefs[1] #Pollen survival
  polSlope <- coefs[2]
  ovInt <- coefs[3] #Ovule survival
  ovSlope <- coefs[4]
  flwInt <- coefs[5] #Flower survival
  flwSlope <- coefs[6]
  
  set.seed(1)
  #Generate pollen counts for 15000 flowers
  simPol <- rnbinom(15000,mu=mu,size=size)
  #"Abort" flowers with 0 pollen grains
  simPol <- simPol[simPol>0]
  
  #Proportion of pollen germinates and makes it to ovules
  propFert <- invLogit(polInt+polSlope*simPol) #Proportion surviving
  simFert <- rbinom(length(simPol),simPol,propFert)
  #If more than 120 pollen grains make it, cut off fertilization at max ovules
  simFert[simFert>Novules] <- Novules 
  #If there are zero fertilized ovules, "abort" flowers
  simFert <- simFert[simFert>0]
  
  #Proportion of fertilized ovules become seeds
  propSeed <- invLogit(ovInt+ovSlope*simFert) #Proportion surviving
  simSeed <- rbinom(length(simFert),simFert,propSeed)
  #If there are zero seeds, "abort" flowers
  simSeed <- simSeed[simSeed>0]
  
  #Plant aborts flowers according to seed production ('rolls the dice' again)
  propAbort <- invLogit(flwInt+flwSlope*propSeed) #Prob of flower suriving
  simAbort <- rbinom(length(propSeed),1,propAbort) #Bernoulli dist.
  simSeed <- simSeed[simAbort==1]
  
  if(plotResults==T & (length(actualSeeds)!=1 & !is.na(actualSeeds[1]))){
    #Plots
    par(mfrow=c(4,1),mar=c(5,5,2,5))
    #Pollen
    hist(simPol,main=NULL,xlab='Pollen grains',breaks=50,xlim=c(0,max(simPol)))
    par(new=T) #Survival prob
    curve(invLogit(polInt+polSlope*x),0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Prob survival',col='red')
    #Ovules
    hist(simFert,main=NULL,xlab='Fertilized ovules',breaks=50,xlim=c(0,Novules))
    par(new=T) #Survival prob
    curve(invLogit(ovInt+ovSlope*x), #Ovule survival
          0,max(simFert),xlab=NA,ylim=c(range(c(propAbort,propSeed))),ylab=NA,col='red',lty=2,axes=F)
    curve(invLogit(flwInt+flwSlope*x), #Flower survival
          0,max(simFert),xlab=NA,ylab=NA,col='darkred',lty=1,add=T)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Prob survival',col='red')
    #Seeds
    hist(simSeed,main=NULL,xlab='Seeds',breaks=length(unique(simSeed)),xlim=c(0,Novules))
    hist(actualSeeds,main=NULL,xlab='Actual Seed Count',xlim=c(0,Novules),breaks=length(unique(datalist$SeedCount)))
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

simSeeds3(coefs=c(0.36650299,-0.00370369,2.01436660,-0.03364391,-5.37489532,1.80187306),
            mu=logLikNB[[1]][1],size=logLikNB[[1]][2],Novules=120,plotResults=T,
            actualSeeds=datalist$SeedCount,propDiff=F)

margProb3(coefs=c(-1,-0.003,4,-0.1,-5,0.1),dat=datalist$SeedCount,
            maxOvNum=100,maxPolNum=1000,
            mu=logLikNB[[1]][1],size=logLikNB[[1]][2]) 

#LL = 10014.17 - This appears to work OK.
#0.36650299 -0.00370369  2.01436660 -0.03364391 -5.37489532  1.80187306
bestCoefs4 <- optim(par=c(-1,-0.003,4,-0.1,-5,0.1),fn=margProb3,maxOvNum=120,maxPolNum=1000,
                    dat=datalist$SeedCount,mu=292.9317,size=0.6136358)

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
