#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN COMMODITY CANOLA FIELDS (2014+2015)

# Libraries and ggplot theme ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(beepr)

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

#Commodity field data
setwd('~/Projects/UofC/canola_yield_project')
load("./Commodity field analysis/commodityfieldDataAll.RData") 
rm(AICc,brix2mgul,deltaAIC,DIC,plotFixedTerms,predFixef,se,varComp,zeroCheck,conversion,visitorsAll,visitors2015)
fieldsAllComm <- fieldsAll; flowersAllComm <- flowersAll; plantsAllComm <- plantsAll;
seedsAllComm <- seedsAll; surveyAllComm <- surveyAll;
rm(fieldsAll,flowersAll,plantsAll,seedsAll,surveyAll)
#Seed field data
load("./Seed field analysis/seedFieldDataAll.RData")
fieldsAllSeed <- allFields; plantsAllSeed <- allPlants; pollenAllSeed <- allPollen; seedsAllSeed <- allSeeds;
surveyAllSeed <- allSurvey; 
plantsAllSeed$Field <- gsub('Unrah','Unruh',plantsAllSeed$Field) #Fixes spelling error
seedsAllSeed$Field <- gsub('Unrah','Unruh',seedsAllSeed$Field)
seedsAllSeed$EdgeCent <- ifelse(seedsAllSeed$EdgeCent=='Cent','Center',seedsAllSeed$EdgeCent) #Fixes cent/center
rm(allFields,allPlants,allPollen,allSeeds,allSurvey,behav2015,visitors2016,nectar2016,folder)
#Remove Warnock 4 - plants/seeds observed but not field
plantsAllSeed <- filter(plantsAllSeed,Field!='Warnock 4')
seedsAllSeed <- filter(seedsAllSeed,Field!='Warnock 4')
#Set 'negative' missing pods (mistake in counting) to NA.
plantsAllComm <- mutate(plantsAllComm,Missing=ifelse(Missing<0,NA,Missing))
plantsAllSeed <- mutate(plantsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
seedsAllComm <- mutate(seedsAllComm,Missing=ifelse(Missing<0,NA,Missing))
seedsAllSeed <- mutate(seedsAllSeed,Missing=ifelse(Missing<0,NA,Missing))
#Extra seed field data (from Riley W.)
rileyExtra <- read.csv('./Seed field analysis/rileyExtra.csv')
setwd('~/Projects/UofC/canola_yield_project')

#Convenience functions
logit <- function(x){
  return(log(x/(1-x)))
}
invLogit <- function(x){ 
  i <- exp(x)/(1+exp(x)) 
  i[is.nan(i)] <- 1 #If x is large enough, becomes Inf/Inf = NaN, but Lim(invLogit) x->Inf = 1
  return(i)
}
coefs <- function(L){ #Shorter summary from list of L coefficients from a Stan list
  upr=sapply(L,function(x) round(quantile(x,0.975),3))
  lwr=sapply(L,function(x) round(quantile(x,0.025),3))
  meanCoef=sapply(L, mean); stdev=sapply(L, sd)
  data.frame(median=sapply(L,function(x) round(median(x),3)),
             lwr=lwr,upr=upr,mean=round(meanCoef,3),sd=round(stdev,3),z=round(meanCoef/stdev,3),
             overlap=sapply(L,function(x) overlap(x)),
             pval=round(2*(1-pnorm(abs(meanCoef/stdev),0,1)),4)) }
#Linear breakpoint function - two lines with intersection "b"
bpoint <- function(x,int1,slope1,b,slope2) ifelse(x<b,int1+slope1*x,b*slope1+(x-b)*slope2+int1)
#Effect size for Posterior samples
effSize <- function(x) unname(median(x)/diff(quantile(x,c(0.025,0.975))))
#Does 95% of posterior overlap zero?
overlap <- function(x) {r <- quantile(x,c(0.025,0.975))>=0; xor(r[1],r[2]);}
temp <- filter(seedsAllComm)

#Everything above this line run to start

# Data from Wang et al 2011 - ovule counts --------------------------------
ovDat <- read.csv('wang2011ovData.csv') %>% 
  mutate(Freq=ifelse(Freq<0,0,Freq))

ovDatNew <- expand.grid(Freq=NA,Count=20:50,Day=unique(ovDat$Day))
for(i in unique(ovDat$Day)){
  ovDatNew[ovDatNew$Day==i,'Freq'] <- with(ovDat,approx(x=Count[Day==i],y=Freq[Day==i],xout=c(20:50)))$y #Interpolate
  ovDatNew$Freq[is.na(ovDatNew$Freq[ovDatNew$Day==i])] <- 0 #Turn NAs to zeros
  ovDatNew[ovDatNew$Day==i,'Freq'] <- ovDatNew[ovDatNew$Day==i,'Freq']/sum(ovDatNew[ovDatNew$Day==i,'Freq']) #Normalize
}
ggplot(ovDatNew,aes(Count,Freq,col=Day))+geom_line()

#Negative binomial/poisson from data has too much variance. Trying normal dist
par(mfrow=c(5,1))
for(i in unique(ovDatNew$Day)){
  pars <- data.frame(mean=rep(NA,1000),sd=rep(NA,1000))
  for(rep in 1:1000){
    #Simulated count data
    simCount <- with(ovDatNew[ovDatNew$Day==i,],sample(ovDatNew$Count,100,prob=ovDatNew$Freq,replace=T))
    simPar <- optim(c(0.5,1),function(x,count) -sum(dnorm(x=count,mean=x[1],sd=x[2],log=T)),count=simCount)$par
    pars$mean[rep] <- simPar[1]
    pars$sd[rep] <- simPar[2]
  }
  plot(c(20:45),dnorm(c(20:45),mean=median(pars$mean),sd=median(pars$sd)),ylab='Density',type='l',ylim=c(0,max(ovDatNew$Freq)),
       main=i)
  with(ovDatNew[ovDatNew$Day==i,],lines(Count,Freq,col='red'))
}

#overall
par(mfrow=c(1,1))
pars <- data.frame(mean=rep(NA,1000),sd=rep(NA,1000))
for(rep in 1:1000){
  #Simulated count data
  simCount <- with(ovDatNew,sample(ovDatNew$Count,100,prob=ovDatNew$Freq,replace=T))
  simPar <- optim(c(0.5,1),function(x,count) -sum(dnorm(x=count,mean=x[1],sd=x[2],log=T)),count=simCount)$par
  pars$mean[rep] <- simPar[1]
  pars$sd[rep] <- simPar[2]
}
plot(c(20:45),dnorm(c(20:45),mean=median(pars$mean),sd=median(pars$sd)),ylab='Density',type='l',ylim=c(0,0.15))
ovDatNew %>% group_by(Count) %>% summarize(Freq=sum(Freq)) %>% ungroup() %>% mutate(Freq=Freq/sum(Freq)) %>% 
  with(.,lines(Count,Freq,col='red'))
points(c(20:45),dpois(c(20:45),32.15),col='blue') #Poisson is too wide

#Mean=32.15, SD=2.66 for overall ovule distribution. Discrete version:
plot(10:55,log(pnorm(11:56,32.15,2.66)-pnorm(10:55,32.15,2.66)),ylab='LogLik')



#Pod weight-count modelling (JAGS) -----------------------------
library(jagsUI)
setwd('./Models')
#Pod count model
datalist <- with(temp[!is.na(temp$PodCount)&temp$PodCount>0,],list( #Strips NAs and 0 counts
  Npod=length(PodCount), #Number of observed seed counts
  Nunique=length(unique(PodCount)), #Number of unique seed counts
  uniquePodCount=sort(unique(PodCount)),
  NuniquePodCount=as.vector(table(PodCount)),
  uniqueMatches=match(PodCount,sort(unique(PodCount))), #Matching index
  SeedCount=PodCount
))

datalist <- c(datalist,with(flowersAllComm[!is.na(flowersAllComm$Pollen),],{list(
  Npollen=length(Pollen),
  PollenCount=Pollen
)})) #Append pollen data to list

datalist <- c(datalist,with(plantsAllComm[!is.na(plantsAllComm$Pods)&!is.na(plantsAllComm$Missing)&plantsAllComm$Missing>0,],{list(
  Nplants=length(Pods),
  Pods=Pods,PodsMissing=Missing
)})) #Append flower success data to list

datalist$SeedCount2 <- datalist$SeedCount

startlist <- function() list(
  #lambda=5.6, r=0.6,
  intPol=-1.15,
  slopePol1=-3.1,
  intPod=1.16,
  slopePod=-0.001,
  simPollenCount=datalist$SeedCount+10,
  simOvCount=datalist$SeedCount+10
)

params <- c(#'lambda','r', #Params for negbin pollination process
  # 'simPollenCount',
  # 'simOvCount',
  'simPollenZero','simFertFailPod','simPodFailPod','totalPodFail',
  #'LLpod','LLseed','lik',
  'simSeedCount','fitSeedCount','fitSimSeedCount',
  'intPol','slopePol'
  #'intPod','slopePod'
  )

tick <- Sys.time()
mod2 <- jags(data=datalist,
          inits=startlist,
          parameters.to.save=params,
          model.file='pod_count_weight.txt',
          n.chains=3,n.adapt=200,n.iter=300,n.burnin=100,n.thin=1,parallel=T)
tock <- Sys.time()
tock-tick 
beep(1)
# save(mod2,file='./podCountResults.RData')
# load('podCountResults.Rdata') 

print(mod2)
traceplot(mod2,parameters=c('intPol','slopePol'))
# traceplot(mod2,parameters='totalPodFail')
pp.check(mod2,actual='fitSeedCount',new='fitSimSeedCount') #Variance = higher in actual dataset. Poisson performs OK, but is weirdly above 1:1 line

mod2samp <- mod2$sims.list #Samples from mod2
str(mod2samp)

plot(datalist$SeedCount,apply(mod2samp$LLpod,2,median),pch=19) #LL decreases with SeedCount. This should be in the opposite direction
points(datalist$SeedCount,apply(mod2samp$LLpod,2,max),pch=3) 
points(datalist$SeedCount,apply(mod2samp$LLpod,2,min),pch=3) 

hist(apply(mod2samp$LLseed,2,median))



#Predicted (median) seed counts vs actual
data.frame(med=apply(mod2samp$simSeedCount,2,median),
           upr=apply(mod2samp$simSeedCount,2,function(x) quantile(x,0.95)),
           lwr=apply(mod2samp$simSeedCount,2,function(x) quantile(x,0.05)),
           actual=datalist$SeedCount) %>%
  ggplot(aes(actual,med,ymax=upr,ymin=lwr))+geom_pointrange()+
  geom_abline(intercept=0,slope=1)+labs(x='Actual',y='Predicted')+
  coord_cartesian(xlim=c(0,55),ylim=c(0,70))

#Simulated seed count
hist(mod2samp$simSeedCount[10,])
table(mod2samp$simSeedCount[10,])

with(mod2samp,{ #Shape of likelihood functions
  lambda <- 5.685
  r <- 0.615
  p <- r/(r+exp(lambda))
  meanNegBin <- (r*(1-p))/p
  
  par(mfrow=c(3,1))
  # curve(invLogit(median(intPol)+(median(slopePol1)*((x-meanNegBin)/1000))+(median(slopePol2)*((x-meanNegBin)/1000)^2)),1,2500,
  #       ylab='p(pollenSurv)')
  # plot(invLogit(median(intPol)+(median(slopePol1)*((c(1:2500)-meanNegBin)/1000))+(median(slopePol2)*((c(1:2500)-meanNegBin)/1000)^2))*
  #         dnbinom(c(1:2500),prob=p,size=r),ylab='p(pollenSurv)*p(pollen)',pch=19)
  curve(invLogit(bpoint(x/1000,median(intPol),median(slopePol1),median(polBp),median(slopePol2))),1,2500,ylab='p(pollenSurv)')
  plot(invLogit(bpoint(1:2500/1000,median(intPol),median(slopePol1),median(polBp),median(slopePol2)))*dnbinom(c(1:2500),
        prob=p,size=r),ylab='p(pollenSurv)*p(pollen)',pch=19)
  
  plot(1:2500,((c(1:2500)-meanNegBin)/1000),pch=19,type='l')
  print(((c(1)-meanNegBin)/1000))
  par(mfrow=c(1,1))
})

#Pairplot of mean seeds vs other params
with(mod2samp,data.frame(naSeeds=apply(simSeedCount,1,function(x) sum(is.na(x))),
                         meanSeeds=apply(simSeedCount,1,function(x) mean(x,na.rm=T)),
                         #meanOv=apply(simOvCount,1,mean),
                         deviance,
                         intPol,slopePol1,slopePol2,polBp)) %>% pairs()

# #Ovule count compared to seed count
# data.frame(SeedCount=datalist$SeedCount,simOvMed=apply(mod2samp$simOvCount,2,median),
#            upr=apply(mod2samp$simOvCount,2,max),
#            lwr=apply(mod2samp$simOvCount,2,min)) %>% 
#   ggplot(aes(SeedCount,simOvMed))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
#   geom_abline(intercept=0,slope=1)
# 
# #Pollen count compared to seed count
# data.frame(SeedCount=datalist$SeedCount,simPolMed=apply(mod2samp$simPollenCount,2,median),
#            upr=apply(mod2samp$simPollenCount,2,max),
#            lwr=apply(mod2samp$simPollenCount,2,min)) %>%
#   ggplot(aes(SeedCount,simPolMed))+geom_pointrange(aes(ymax=upr,ymin=lwr))+
#   geom_abline(intercept=0,slope=1)
# 
# data.frame(isna=is.na(mod2samp$simSeedCount[5,]),simSeed=mod2samp$simSeedCount[5,],
#            simPol=mod2samp$simPollenCount[5,],#simOv=mod2samp$simOvCount[5,],
#            seedCount=datalist$SeedCount) %>% 
#   ggplot(aes(simOv,simSeed))+geom_point()+
#   geom_point(aes(simOv,seedCount),col='red')


hist(mod2samp$simSeedCount[100,])

with(mod2samp,{
  par(mfrow=c(2,1))
  iter <- 51
  hist(simPollenCount[iter,],main='Simulated pollen',xlab='Count') #Simulated pollen (50th iteration)
  hist(simOvCount[iter,],main='Simulated ovules',xlab='Count') #Simulated fertilized ovules
  # hist(simFertSeed[iter,],main='Simulated fertilized seeds',xlab='Count',breaks=seq(0,52,1)) #Simulated fertilized ovules 
  hist(datalist$SeedCount,main='Acutual seed count',xlab='Count',breaks=seq(1,52,1))
  print(c(intPol[iter],slopePol1[iter],slopePol2[iter],polBp[iter]))
  par(mfrow=c(1,1))
})



#Test
par(mfrow=c(2,1))
simPolCount <- 
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


#Pod weight-count modeling - custom likelihood (Stan) ------------------------------------

library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

#Rows without missing seed counts/weights
# noMiss <- with(temp,rowSums(is.na(cbind(PodCount,PodMass)))==0 & PodCount!=0)

#List structure for Stan

datalist <- with(temp[!is.na(temp$PodCount)&temp$PodCount>0,],list( #Strips NAs and 0 counts
  Npod=length(PodCount), #Number of observed seed counts
  Nunique=length(unique(PodCount)), #Number of unique seed counts
  uniquePodCount=sort(unique(PodCount)),
  NuniquePodCount=as.vector(table(PodCount)),
  uniqueMatches=match(PodCount,sort(unique(PodCount))), #Matching index
  SeedCount=PodCount
))

datalist <- c(datalist,with(flowersAllComm[!is.na(flowersAllComm$Pollen),],{list(
  Npollen=length(Pollen),
  PollenCount=Pollen
)})) #Append pollen data to list

datalist <- c(datalist,with(plantsAllComm[!is.na(plantsAllComm$Pods)&!is.na(plantsAllComm$Missing)&plantsAllComm$Missing>0,],{list(
  Nplants=length(Pods),
  Pods=Pods,PodsMissing=Missing
)})) #Append flower success data to list

modPodcount = stan(file='pod_level.stan',data=datalist,iter=100,chains=1,control=list(adapt_delta=0.8),
                   init=function() list(intPolSurv=0,slopePolSurv=-3,intPodSurv=-2,slopePodSurv=3))
print(modPodcount)
traceplot(modPodcount) #Very poor traces, esp for pollen mu and phi. Likely what is happening is the mu and phi terms are being "influenced" by seed counts. Not sure how to separate these two, aside from heirarchical terms.
#Try running model without pollen included, and use results from below for pollen input.
#Same results occur even when using pollen. Try to run this in JAGS instead?

modPolcount <- stan(file='pollenMod.stan',data=datalist[c('Npollen','PollenCount')],iter=1000,chains=3)
print(modPolcount) #mu=293.75, phi=0.61
traceplot(modPolcount)

# podMod <- stan_model(file='pod_level.stan')
# optFit <- optimizing(podMod,data=datalist,init=list(
#   intPolSurv=-1,slopePolSurv=-3.5,intPodSurv=-2.9,slopePodSurv=3.5))
# round(optFit$par,4) #Weird. Chooses super-low survival probs.

#Check output
pars <- c('intPol','intFert','slopeFert')
print(modPodcount,pars=pars)
traceplot(modPodcount,pars=pars)
pairs(modPodcount,pars=pars)

# #Simulate seed count
# simcounts <- lapply(extract(modPodcount,pars=pars),median)
# invLogit <- function(x) exp(x)/(1+exp(x))
# par(mfrow=c(3,1))
# simFertCount <- rbinom(datalist$Npod,120,invLogit(simcounts[[1]])) #Simulate fert ov count
# simSeedCount <- rbinom(datalist$Npod,simFertCount,invLogit(simcounts[[2]]+simcounts[[3]]*simFertCount)) 
# hist(simFertCount,main=NULL,xlab='Fert Ov Count',xlim=c(0,100)) 
# hist(simSeedCount,main=NULL,xlab='Sim Seed Count',xlim=c(0,100),breaks=length(unique(simSeedCount)))
# hist(datalist$SeedCount,main=NULL,xlab='Actual Seed Count',xlim=c(0,100),breaks=length(unique(datalist$SeedCount)))

#Doesn't look very good. Perhaps try an arrival model of pollen?
# #Fit distributions
# poisFun <- function(x,dat) -sum(dpois(dat,x,log=T)) #Poisson
# logLikPois <- optimize(f=poisFun,lower=0,upper=1000,dat=flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)])
# geomFun <- function(x,dat) -sum(dgeom(dat,x,log=T)) #Geometric
# logLikGeom <- optimize(f=geomFun,c(0,1),dat=flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)])
# negbinFun <- function(x,dat) -sum(dnbinom(dat,mu=x[1],size=x[2],log=T)) #Neg Bin
# logLikNB <- optim(c(1,1),negbinFun,dat=flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)])
# 
# 2*2+(2*logLikNB$value) #AIC for NB - best
# 2*1+(2*logLikGeom$objective) #AIC for geom - worse
# 2*1+(2*logLikPois$objective) #AIC for poisson - much worse
# 
# par(mfrow=c(4,1)) #Plot results
# hist(rnbinom(datalist$Npod,mu=logLikNB[[1]][1],size=logLikNB[[1]][2]),main=NULL,xlab='Neg Bin',xlim=c(0,4000),breaks=seq(0,4000,20))
# hist(rgeom(datalist$Npod,logLikGeom[[1]][1]),main=NULL,xlab='Geom',xlim=c(0,4000),breaks=seq(0,4000,20))
# hist(rpois(datalist$Npod,logLikPois[[1]][1]),main=NULL,xlab='Pois',xlim=c(0,4000),breaks=seq(0,4000,20))
# hist(flowersAllComm$Pollen[!is.na(flowersAllComm$Pollen)],xlim=c(0,4000),breaks=seq(0,4000,20),main=NULL,xlab='Actual pollen counts')


# Test model to simulate seed counts --------------------------------------

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

margProb <- function(coefs,dat,maxOvNum,maxPolNum,mu,size,lambda){
  
  seedCount <- c(1:maxOvNum) #Number of seeds to try
  
  #Replicate all levels of seed/ovule/pollen counts
  prob <- expand.grid(seedCount=seedCount,ovCount=1:maxOvNum,polCount=1:maxPolNum)
  #Remove impossible categories
  prob <- prob[(prob$ovCount>=prob$seedCount) & (prob$polCount>=prob$seedCount),]
  #Pollen survial prob
  prob$probPol <- invLogit(coefs[1]+coefs[2]*(prob$polCount/1000))
  #Pod survival prob
  # prob$probSeed <- invLogit(coefs[3]+coefs[4]*prob$seedCount)
  
  lpPollen <- dnbinom(1:maxPolNum,mu=mu,size=size,log=T) #Pollen LP
  lpOv <- dpois(1:maxOvNum,lambda=lambda,log=T) #Ovule LP
  # lpPolSurv <- matrix(NA,nrow=maxPolNum,ncol=maxPolNum) #Pollen survival LP
  # for(i in 1:maxPolNum){ #For each potential pollen receipt
  #   lpPolSurv[i,1:i] <- dbinom(1:i,i,invLogit(coefs[1]+coefs[2]*i)) #Pollen survival lp
  # }
  lpPodSurv <- dbinom(1,1,invLogit(coefs[3]+coefs[4]*(seedCount/10)),log=T) #Pod survival prob
  
  #Calculate log-prob for each
  prob$lp <- lpPollen[prob$polCount] +
    lpOv[prob$ovCount]+
    ifelse(prob$seedCount==prob$ovCount, #Fert prob
           #if seeds==ovules, successful pollen count could have been anything between #seeds and #pollen
           pbinom(prob$seedCount-1,prob$polCount,prob$probPol,log=T,lower.tail=F),
           #if seeds<ovules, only one way to get seed count
           # lpPolSurv[prob$seedCount,prob$polCount]
           dbinom(prob$seedCount,prob$polCount,prob$probPol,log=T)
           )+
    lpPodSurv[prob$seedCount]
  #Sum by seed counts
  seedLp <- with(prob,tapply(lp,seedCount,function(x) log(sum(exp(x))))) 
  # barplot(c(1-sum(exp(seedLp)),exp(seedLp)),pch=19,ylab='p(x)',xlab='Seed count')
  
  
  #Match -lp to actual counts, sum
  lp <- sum(seedLp[match(dat$SeedCount,names(seedLp))])+ #Seed count lp
    sum(dbinom(dat$FailFlw,dat$AllFlw,1-sum(exp(seedLp)),log=T)) #Flw abortion lp
  
  return(-lp) 
}

margProb(c(-2,-1,-3.0506205,4.580386),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32)

library(doSNOW)
cl <- makeCluster(3,type='SOCK')
clusterExport(cl,list('margProb','invLogit','datalist'))

#Try across various pollen survival probs (1st arg)
ll3 <- parSapply(cl,seq(-3,3,length.out=15),function(x) 
  margProb(c(x,-3.5249114,-2.8723591,0.3460556),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))

#Pollen survival slopes (2nd arg)
ll4 <- parSapply(cl,seq(-6,-2,length.out=15),function(x) 
  margProb(c(-1.0314562,x,-2.8723591,0.3460556),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))

#Try across various pod abortion intercepts (3rd arg)
ll2 <- parSapply(cl,seq(-6,0,length.out=15),function(x) 
  margProb(c(-1.0314562,-3.5249114,x,0.3460556),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))

#Try across various pod abortion slopes (4th arg)
ll <- parSapply(cl,seq(-3,3,length.out=15),function(x)
  margProb(c(-1.0314562,-3.5249114,-2.8723591,x),dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),52,1800,292.9316584,0.6136358,32))

par(mfrow=c(2,2))
plot(seq(-3,3,length.out=15),ll3,xlab='Pollen surv int | Pod surv Intercept,Slope',ylab='-Log Lik',pch=19,main='Arg1: pol surv int')
plot(seq(-6,-2,length.out=15),ll4,xlab='Pollen surv int | Pod surv Intercept,Slope',ylab='-Log Lik',pch=19,main='Arg2: pol surv slope') 
plot(seq(-6,0,length.out=15),ll2,xlab='Intercept | Slope values, Pollen surv',ylab='-Log Lik',pch=19,main='Arg3: pod surv int')
plot(seq(-3,3,length.out=15),ll,xlab='Pod surv slope | Pod surv int, Pol surv int, slope',ylab='-Log Lik',pch=19,main='Arg4: pod surv slope')
par(mfrow=c(1,1))
beep(1)

stopCluster(cl)

library(optimParallel)
cl <- makeCluster(3,type='SOCK')
clusterExport(cl,list('margProb','invLogit','datalist'))
setDefaultCluster(cl=cl)

# -1.031621,-3.524373,-2.871942,3.460029; LL:23946.96
optFit <- optimParallel(c(-1.031,-3.524,-2.872,3.46),margProb,
                dat=list(SeedCount=datalist$SeedCount,FailFlw=datalist$PodsMissing,AllFlw=datalist$PodsMissing+datalist$Pods),
                maxOvNum=52,maxPolNum=1500,mu=292.9316584,size=0.6136358,lambda=32,
                method='L-BFGS-B',lower=c(-10,-10,-20,-1),upper=c(10,10,10,13),parallel=cl)
# optFit <- c(-1.031621,-3.524373,-2.871942,3.460029)

stopCluster(cl)
detach('package:optimParallel')

simSeeds <- function(maxNovules,coefs,plotResults=T,actualSeeds=NA,propDiff=F){
  #Idea: 
  #Pollination is a Negative Binomial process with mu=292.9,size=0.61
  mu <- coefs[1]
  size <- coefs[2]
  #Ovule number per flower is a Poisson process with lambda =~32. Wang et al 2011 found that there are cross-season difference, but they aren't that big (~30-35 ovules).
  lambda <- 32
  #Fertilization is a Binomial process, where Fert Ovules ~ Bin(Pollen,p)
  #If more ovules are fertilized than available, fert ovules = total ovules (integrate all outcomes < Novules)
  
  #Pod abortion is a Bernoulli process where Abortion ~ Bern(invLogit(int1+slope1*Nfert))
  int1 <- coefs[3]
  slope1 <- coefs[4]
  
  int2 <- coefs[5]
  slope2 <- coefs[6]
  
  set.seed(1)
  #Generate pollen counts for 1000 flowers
  simPol <- rnbinom(5000,mu=mu,size=size)
  #Generate ovule counts for flowers
  simOv <- rpois(5000,32)
  simOv[simOv>maxNovules] <- maxNovules #Set all ovules > maxNovules to maxNovules
  
  #Ovules are fertilized by pollen
  probFert <- invLogit(int1+slope1*simPol/1000)
  simFert <- rbinom(length(simPol),simPol,probFert) #Fixed perc of pollen makes it to ovules
  simFert[simFert>simOv] <- simOv[simFert>simOv] #If Npollen>Novules, reverts to ov number "first come first served"
  
  #Fertilized ovules become seeds
  simSeed <- simFert[simFert>0] #Fertilized (suriving pollen>1) ovules become seeds
  probSeed <- invLogit(int2+slope2*simSeed/10) #Prob of pod abortion ~ # of fertilized ovules
  simSeed <- simSeed[rbinom(length(simSeed),1,probSeed)==1] #Plant aborts pods depending on number of ovules
  
  if(plotResults==T & (length(actualSeeds)!=1 & !is.na(actualSeeds[1]))){
    #Plots
    par(mfrow=c(5,1),mar=c(5,5,2,5))
    #Pollen
    hist(simPol,main=NULL,xlab='Pollen grains',breaks=50,xlim=c(0,max(simPol)))
    par(new=T) #Pollen survival
    curve(invLogit(int1+slope1*x/1000),0,max(simPol),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
    axis(side=4,cex.axis=1,col='red',col.axis='red') 
    mtext(side=4,line=3,'Pollen survival ',col='red')
    #Ovules
    hist(simOv,xlab='Ovules',xlim=c(0,max(simOv)),breaks=seq(0,max(simOv),1),main=NULL)
    #Fertilized ovules
    hist(simFert,main=NULL,xlab='Fertilized Ovules',breaks=seq(0,max(simFert),1),xlim=c(0,max(simOv)))
    par(new=T) #Survival prob
    curve(invLogit(int2+slope2*x/10),0,max(simFert),xlab=NA,ylab=NA,col='red',lty=2,axes=F)
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

simSeeds(maxNovules=52,coefs=c(292.9316584,0.6136358,optFit$par),plotResults=T,actualSeeds=datalist$SeedCount)
#Looks pretty good, actually

#Proportion missing seeds. Should be related to 
plantsAllComm %>% mutate(PropMis=Missing/(Pods+Missing)) %>% filter(!is.na(PropMis),PropMis>0) %>% 
  ggplot(aes(PropMis))+geom_histogram()



# Commodity field visitation and pollen deposition (Stan) ---------------------------------
library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
library(shinystan)

#List structure for Stan
datalistField <- with(arrange(fieldsAllComm,as.factor(paste(Year,Field))),list( #Field-level measurements
  Nfield=length(Year), #Number of fields
  #fieldIndex = as.numeric(as.factor(paste(Year,Field))), #Index for each field
  numHives=NumHives, #Number of hives/field
  is2015=Year==2015, #Is year 2015?
  isGP=Area!='Lethbridge', #Is area Grand Prairie?
  isIrrigated=Irrigated=='Irrigated' #Is field irrigated?
))

datalistPlot <- with(arrange(surveyAllComm,factor(paste(Year,Field,Distance))),list( #Plot-level measurements
  Nplot=length(Distance), #Number of plots
  plotIndex=as.numeric(as.factor(paste(Year,Field))), #Index for field (which field?)
  dist=Distance, #Distance from edge
  hbeeVis=Honeybee, #Visits by honeybees
  flyVis=Fly, #Visits by flies
  totalTime=TotalTime/10, #Total time (mins/10)
  #Planting density
  Nplot_densObs=sum(!is.na(PlDens)), #Number of plots with observed flower density
  Nplot_densMiss=sum(is.na(PlDens)), #Plots with unobserved flower density
  plDens_obs=log(PlDens[!is.na(PlDens)])-mean(log(PlDens[!is.na(PlDens)])), #Scaled (log) Planting density
  obsPlDens_ind=which(!is.na(PlDens)),
  missPlDens_ind=which(is.na(PlDens)),
  flDens=sqrt(FlDens)-21 #Flower density - square root and centered
))

datalistFlw <- flowersAllComm %>% 
  mutate(flowerIndex=as.numeric(factor(paste(Year,Field,Distance)))) %>% 
  arrange(flowerIndex) %>% 
  filter(!is.na(Pollen)) %>% 
  with(list(Nflw=sum(!is.na(Pollen)), #Number of pollen samples
  flowerIndex=flowerIndex, #Index for flower (which plot?)
  pollenCount=Pollen
))

datalistPlant <- plantsAllComm %>% ungroup() %>% 
  filter(SeedMass!=0) %>%
  filter(!is.na(Pods),!is.na(Missing)) %>% 
  mutate(plantIndex=as.numeric(factor(paste(Year,Field,Distance)))) %>% #Index for plant (which plot?)
  arrange(plantIndex) %>% 
  with(list(Nplant=length(VegMass), #Number of plant samples (some missing)
  # Nplant_obs=sum(!is.na(VegMass)), #Observed plants
  # Nplant_miss=sum(is.na(VegMass)), #Missing plants
  podCount=Pods, #Successful pods
  flwCount=Pods+Missing, #Pods + Missing (total flw production)
  plantIndex=plantIndex, #Index for plant (which plot?)
  plantSize=log(VegMass[!is.na(VegMass)])-mean(log(VegMass[!is.na(VegMass)])), #log weight of veg mass (g), centered
  yield=SeedMass[!is.na(SeedMass)] #observed weight of all seeds (g)
  # obsPl_ind=which(!is.na(VegMass)),missPl_ind=which(is.na(VegMass)), #Indices for missing plant size
  # obsYield_ind=which(!is.na(SeedMass)), missYield_ind=which(is.na(SeedMass)) #Indices for missing yield
))

#Problem: plots exist at the pod level which do not exist at plant level
a <- plantsAllComm %>% ungroup() %>% filter(!is.na(Pods),!is.na(Missing)) %>%  filter(SeedMass!=0) %>%
  transmute(index=factor(paste(Year,Field,Distance,Plant))) %>% distinct() #Index for plant (which plot?)
b <- seedsAllComm %>% ungroup() %>% filter(!is.na(Plant)&!is.na(PodCount)&PodCount>0&!is.na(PodMass)) %>%
  filter(!is.na(Pods),!is.na(Missing)) %>% #Remove plants from plant level
  transmute(index=factor(paste(Year,Field,Distance,Plant))) %>% distinct()
keep <- which(a$index %in% b$index) #Plants to keep from seeds plants dataset

datalistPod <- seedsAllComm %>% ungroup() %>% 
  mutate(podIndex=as.numeric(factor(paste(Year,Field,Distance,Plant)))) %>%   
  arrange(podIndex) %>% 
  filter(podIndex %in% keep) %>%
  filter(!is.na(Plant)&!is.na(PodCount)&PodCount>0&!is.na(PodMass)) %>% 
  filter(!is.na(Pods),!is.na(Missing)) %>% #Remove plants from plant level
  with(list(Npod=length(Distance), #Number of seeds measured
  seedCount=PodCount, #Number of seeds per pod
  seedMass=(PodMass/PodCount)*1000, #Weight per seed (mg)
  podIndex=podIndex) #Index for pod (which plant?)
)
datalistPod$seedMass[datalistPod$seedMass>8] <- with(datalistPod,seedMass[seedMass>8]/10) #Fixes weird outliers
datalist <- c(datalistField,datalistPlot,datalistFlw,datalistPlant,datalistPod)
rm(datalistField,datalistPlot,datalistFlw,datalistPlant,datalistPod,a,b,keep) #Cleanup
str(datalist)
 
inits <- function() { with(datalist,
  list(plDens_miss=rep(0,Nplot_densMiss),intPlDens=1, slopeGPPlDens=0, #Plant density
    slopeDistPlDens=0,sigmaPlDens=0.5,sigmaPlDens_field=0.5,
    intPlSize=0,slopePlDensPlSize=0, #Plant size
    slopeDistPlSize=0,slopeGpPlSize=0,slopeIrrigPlSize=0,slope2015PlSize=0,
    sigmaPlSize_field=0.5,sigmaPlSize_plot=0.5,sigmaPlSize=0.5,
    intPlSize_field=rep(0,Nfield),intPlSize_plot=rep(0,Nplot),
    intFlDens=1,slopePlSizeFlDens=0, slopeHbeeDistFlDens=0, #Flower density
    sigmaFlDens=0.5,sigmaFlDens_field=0.5,intFlDens_field=rep(0,Nfield),
    intVisit=-1,slopeYearVis=0,slopeGpVis=0,slopeYearGpVis=1.5,slopeDistVis=0, #Visitation
    slopeHiveVis=0.5,slopeFlDens=0,sigmaVisField=2,visitHbeePhi=0.7,intVisit_field=rep(0.5,Nfield),lambdaVisField=1,
    intPollen=5.5,slopeVisitPol=0,slopeHbeeDistPollen=0,slopeStockingPollen=0, #Pollen deposition
    slopeStockingHbeeDistPollen=0,sigmaPolField=0.5,sigmaPolPlot=0.4,pollenPhi=0.7,
    intPollen_field=rep(0,datalist$Nfield),intPollen_plot=rep(0,datalist$Nplot),
    intFlwCount=1,slopePlSizeFlwCount=0,sigmaFlwCount_field=0.5,sigmaFlwSurv_plot=0.5, #Flower count per plant
    intFlwSurv_plot=rep(0,Nplot),intFlwCount_field=rep(0,Nfield),flwCountPhi=1,
    #Flower survival
    intFlwSurv=1,slopeVisitSurv=0,slopePolSurv=0,slopePlSizeSurv=0.02, 
    slopeIrrigSurv=0,slope2015Surv=0,slopeIrrig2015Surv=0,slopePlSizeIrrigSurv=0,
    slopeFlwCountSurv=0,slopePlDensSurv=0,sigmaFlwSurv_plot=0.3,sigmaFlwSurv_field=0.3, 
    intFlwSurv_field=rep(0,Nfield),intFlwSurv_plot=rep(0,Nplot),flwSurvPhi=1,
    #Seed count
    intSeedCount=3.15,slopeVisitSeedCount=0.05,slopePolSeedCount=0,seedCountPhi=21, 
    slopeIrrigSeedCount=0,slopeIrrig2015SeedCount=0,slopePlSizeIrrigSeedCount=0,
    slopePlSize2015SeedCount=0,slopePlSizeIrrig2015SeedCount=0,
    slope2015SeedCount=0,slopePlSizeCount=0, sigmaSeedCount_plant=0.14,
    sigmaSeedCount_plot=0.07,sigmaSeedCount_field=0.1,
    intSeedCount_field=rep(0,datalist$Nfield),intSeedCount_plot=rep(0,datalist$Nplot),
    intSeedCount_plant=rep(0,datalist$Nplant),
    #Seed weight
    intSeedWeight=0,slopeVisitSeedWeight=0,slopePolSeedWeight=0,slopeSeedCount=0, 
    slopePlSizeWeight=0,slopeIrrigSeedWeight=0,
    slope2015SeedWeight=0,slope2015IrrigSeedWeight=0,slopeSeedCountVegMassSeedWeight=0,
    slopePlSizeIrrigSeedWeight=0,slopeSeedCount2015SeedWeight=0,slopePlSize2015SeedWeight=0,
    slopePlSizeIrrig2015SeedWeight=0,
    sigmaSeedWeight=0.5,sigmaSeedWeight_plant=0.5,sigmaSeedWeight_plot=0.5,sigmaSeedWeight_field=0.5,
    intSeedWeight_field=rep(0,datalist$Nfield),slopePlSizeSeedWeight=0,
    intSeedWeight_plot=rep(0,datalist$Nplot),intSeedWeight_plant=rep(0,datalist$Nplant),
    lambdaSeedWeight=1.5,
    #Yield
    intYield=0,slopeYield=1,sigmaYield=0.5 
   ))
}

# #Claims list
# claims6 <- stan(file='./Commodity model claims 1/commodity_claims50.stan',data=datalist,iter=800,chains=4,
#                 control=list(adapt_delta=0.8),init=inits)
# beep(1)
# newpar='slopeFlwCountFlDens'
# stan_hist(claims6,pars=c(pars,newpar))+geom_vline(xintercept=0,linetype='dashed')
# traceplot(claims6,pars=c(pars,newpar),inc_warmup=F)+geom_hline(yintercept=0,linetype='dashed')
# mod1 <- extract(claims6)
# 2*(1-pnorm(abs(mean(mod1[[1]])/sd(mod1[[1]])),0,1)) #p-val for claim (2-tailed)
# print(claims1,pars=c(pars,newpar))

#Full model - 1.7 hrs for 1000 iter
modPodcount3 <- stan(file='visitation_pollen_model.stan',data=datalist,iter=1000,chains=3,
                   control=list(adapt_delta=0.8),init=inits)
library(loo)
loo_modPodcount1 <- loo(modPodcount,pars='log_lik_plant') #Model without plDens:stocking
loo_modPodcount2 <- loo(modPodcount2,pars='log_lik_plant') #model with plDens:stocking
loo_modPodcount3 <- loo(modPodcount3,pars='log_lik_plant') #model without plDens:stocking or stocking

compare(loo_modPodcount1,loo_modPodcount2,loo_modPodcount3) #Very little difference between model with
compare(loo_modPodcount1,loo_modPodcount2) #Very little difference between model with
compare(loo_modPodcount1,loo_modPodcount3) #Very little difference between model with 


# load('modPodcount.Rdata') #Seed weight, size, and yield
# load('modPodcount2.Rdata') #All other coefficients (plot/plant level) - 22 mins for 1000 iter

# print(modPodcount)
pars=c('intPlDens','slope2015PlDens','slopeIrrigPlDens','slope2015IrrigPlDens',
       'slopeDistPlDens','slopeGPPlDens','sigmaPlDens','sigmaPlDens_field') #Planting density
pars=c('intPlSize','slopePlDensPlSize','slopeDistPlSize','slopeGpPlSize', #Plant size
       'slopeIrrigPlSize','slope2015PlSize','slopeStockingPlSize',
       'slopePlDensStockingPlSize',
       'sigmaPlSize_field','sigmaPlSize_plot','sigmaPlSize')
pars=c('intFlDens','slopePlSizeFlDens','slopeHbeeDistFlDens','sigmaFlDens','sigmaFlDens_field') #Flower density
pars=c('intVisit','slopeYearVis','slopeGpVis','slopeYearGpVis','slopeIrrigVis',
       'slopeDistVis','slopeHiveVis','slopeFlDens', #Visitation
       'sigmaVisField','lambdaVisField','visitHbeePhi')
pars=c('intPollen','slopeVisitPol','slopeHbeeDistPollen',#'slopeFlyVisPol',
       'sigmaPolField','pollenPhi') #Pollen deposition
pars=c('intFlwCount','slopePlSizeFlwCount', #Flower count per plant
       'sigmaFlwCount_field','sigmaFlwCount_plot','flwCountPhi')
pars=c('intFlwSurv','slopeVisitSurv','slopePolSurv','slopePlSizeSurv',
       'slopePlDensSurv','slopeIrrigSurv','slope2015Surv',
       'slopeIrrig2015Surv','slopePlSizeIrrigSurv',
       'sigmaFlwSurv_field','flwSurvPhi') #Flower survival
pars=c('intSeedCount','slopeVisitSeedCount','slopePolSeedCount','slopePlSizeCount',
       'slope2015SeedCount','slopeIrrigSeedCount','slopeIrrig2015SeedCount','slopePlSizeIrrigSeedCount',
       'slopePlSize2015SeedCount','slopePlSizeIrrig2015SeedCount',
       'seedCountPhi','sigmaSeedCount_plant','sigmaSeedCount_field') #Seed count
pars=c('intSeedWeight','slopeVisitSeedWeight','slopePolSeedWeight',#Seed weight
       'slopeSeedCount','slopePlSizeWeight','slopeIrrigSeedWeight',
       'slope2015SeedWeight','slope2015IrrigSeedWeight','slopeSeedCountPlSizeSeedWeight',
       'slopePlSizeIrrigSeedWeight','slopeSeedCount2015SeedWeight','slopePlSize2015SeedWeight',
       'slopePlSizeIrrig2015SeedWeight',
       'sigmaSeedWeight','sigmaSeedWeight_plant','sigmaSeedWeight_field','lambdaSeedWeight')
pars=c('intYield','slopeYield','sigmaYield',
       'sigmaYield_field','sigmaYield_plot','L_field','L_plot')
stan_hist(modPodcount3,pars=pars)+geom_vline(xintercept=0,linetype='dashed')
traceplot(modPodcount,pars=c(pars),inc_warmup=F)
# launch_shinystan(modPodcount)
print(modPodcount2,pars=pars)
# pairs(modPodcount,pars=pars) #Takes way too long

mod3 <- extract(modPodcount3)
coefs(mod3[pars])
library(xtable)
print(xtable(coefs(mod3[pars]),digits=c(0,3,3,3,3,3,3,0,4)))

#Faster Pairplots
pairs(mod3[c(pars,'lp__')],lower.panel=function(x,y){
  par(usr=c(0,1,0,1))
  text(0.5, 0.5, round(cor(x,y),2), cex = 1 * exp(abs(cor(x,y))))})

#Plot of random intercepts
t(apply(mod3$ranEffPlSize_plot[,2,],2,function(x) quantile(x,c(0.5,0.975,0.025)))) %>%
  as.data.frame() %>% rename(median='50%',upr='97.5%',lwr='2.5%') %>% arrange(median) %>%
  mutate(row=1:nrow(.)) %>% 
  ggplot(aes(row,median))+geom_pointrange(aes(ymax=upr,ymin=lwr))+geom_hline(yintercept=0,col='red')
  # ggplot(aes(median))+geom_density()

qqnorm(apply(mod3$intPollen_plot,2,median));qqline(apply(mod3$intPollen_plot,2,median));
mean(apply(mod3$intVisit_field,2,median))

#Plots of random intercepts/slopes for yield
par(mfrow=c(2,1))
plot(apply(mod3$ranEffYield_field[,1,],2,mean),apply(mod3$ranEffYield_field[,2,],2,mean),
     xlab='Intercept',ylab='Slope',main='Field',pch=19) #Weakly correlated (r=-0.55)
abline(h=0,lty='dashed');abline(v=0,lty='dashed')
plot(apply(mod3$ranEffYield_plot[,1,],2,mean),apply(mod3$ranEffYield_plot[,2,],2,mean),
     xlab='Intercept',ylab='Slope',main='Plot',pch=19) #Highly correlated (r=-0.98)
abline(h=0,lty='dashed');abline(v=0,lty='dashed')
par(mfrow=c(1,1))

#Check model fit:
par(mfrow=c(2,1))
#planting density - good
with(mod3,plot(apply(plDens_resid,1,function(x) sum(abs(x))), #PP plot - good
               apply(predPlDens_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Planting density'))
abline(0,1,col='red') 
plot(datalist$plDens_obs,apply(mod3$predPlDens,2,median)[datalist$obsPlDens_ind], #Predicted vs Actual
     ylab='Predicted plant density',xlab='Actual plant density'); abline(0,1,col='red');

#plant size - good
with(mod3,plot(apply(plSize_resid,1,function(x) sum(abs(x))),
               apply(predPlSize_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Plant size')) 
abline(0,1,col='red'); #PP plot - good
plot(datalist$plantSize,apply(mod3$predPlSize,2,median), #Predicted vs Actual - good
     ylab='Predicted plant size',xlab='Actual plant size'); abline(0,1,col='red')

#flower density per plot
with(mod3,plot(apply(flDens_resid,1,function(x) sum(abs(x))),
               apply(predFlDens_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Flower density per plot')) 
abline(0,1,col='red'); #PP plot - good
plot(datalist$flDens,apply(mod3$predFlDens,2,median), #Predicted vs Actual - good
     ylab='Predicted flower density',xlab='Actual flower density'); abline(0,1,col='red')

#hbee visits - good
with(mod3,plot(apply(hbeeVis_resid,1,function(x) sum(abs(x))),
               apply(predHbeeVis_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Hbee visits')) 
abline(0,1,col='red')
legend('topleft',paste('p=',round(mean(with(mod3,apply(hbeeVis_resid,1,function(x) sum(abs(x)))>apply(predHbeeVis_resid,1,function(x) sum(abs(x))))),3)))
#PP plot - OK
plot(with(datalist,hbeeVis/totalTime),jitter(apply(mod3$predHbeeVis,2,median)), #Predicted vs Actual - good
     ylab='Predicted visits',xlab='Actual visits')
abline(0,1,col='red')

#pollen - OK, but plot level random effects throw off PP checks
with(mod3,plot(apply(pollen_resid,1,function(x) sum(abs(x))),
               apply(predPollen_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Pollen counts')) 
abline(0,1,col='red') #PP plot
with(mod3,{x <- sum(apply(pollen_resid,1,function(x) sum(abs(x)))<
             apply(predPollen_resid,1,function(x) sum(abs(x))))/nrow(predPollen_resid)
  legend('topleft',paste('p =',round(min(x,1-x),3)))})
plot(datalist$pollenCount,apply(mod3$predPollenCount,2,median), #Predicted vs Actual 
     ylab='Predicted pollen',xlab='Actual pollen')
abline(0,1,col='red') #PP plot

#seeds per pod - bad
with(mod3,plot(apply(seedCount_resid,1,function(x) sum(abs(x))),
               apply(predSeedCount_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Seeds per pod')) 
abline(0,1,col='red') #PP plot - not good
plot(datalist$seedCount,apply(mod3$predSeedCount,2,median), #Predicted vs Actual - good
     ylab='Predicted seed count',xlab='Actual seed count'); abline(0,1,col='red')

#weight per seed - good
with(mod3,plot(apply(seedMass_resid,1,function(x) sum(abs(x))),
               apply(predSeedMass_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Weight per seed'))
abline(0,1,col='red') #PP plot 
plot(datalist$seedMass,apply(mod3$predSeedMass,2,median), #Predicted vs Actual - OK, but weird outliers: check data
     ylab='Predicted seed weight',xlab='Actual seed weight'); abline(0,1,col='red'); 

#Partial plots
temp <- with(datalist,data.frame(int=rep(1,Npod),seedCount, #Dataframe to hold info
  logHbeeVis=with(datalist,log((hbeeVis/totalTime)+0.5))[plantIndex[podIndex]],
  pollen=apply(mod3$pollenPlot,2,mean)[plantIndex[podIndex]], #Mean pollen at plot
  lambdaInt=rep(1,Npod),
  resid=apply(mod3$seedMass_resid,2,mean))) %>% as.matrix()

#Model matrix
modMat <- as.matrix(with(mod3,data.frame(intSeedWeight,slopeSeedCount,
                               slopeVisitSeedWeight,slopePolSeedWeight,
                               lambdaSeedWeight=1/lambdaSeedWeight)))
#Seedcount
temp2 <- as.data.frame(temp) %>% #Marginalize across all vars except seedCount
  mutate_at(vars(logHbeeVis,pollen),mean) %>% as.matrix()
temp2 <- cbind(temp2,t(apply(temp2[,c(1:5)] %*% t(modMat),1,function(x) quantile(x,c(0.5,0.975,0.025))))) %>% 
  as.data.frame() %>% rename(median='50%',upr='97.5%',lwr='2.5%') %>% 
  mutate(predUpr=apply(mod3$predSeedMass,2,quantile,0.975)) %>% 
  mutate(predLwr=apply(mod3$predSeedMass,2,quantile,0.025)) %>% 
  arrange(seedCount) %>% filter(resid<6,seedCount<40)

ggplot(temp2,aes(seedCount,median))+geom_ribbon(aes(ymax=upr,ymin=lwr),fill='red',alpha=0.3)+
  geom_point(aes(seedCount,(median+resid)),alpha=0.3)+
  # geom_point(aes(seedCount,predUpr),alpha=0.3,col='red')+
  # geom_point(aes(seedCount,predLwr),alpha=0.3,col='red')+
  geom_line(size=1,col='red')+ylab('Weight per seed (mg)')+xlab('Seeds per pod')
  
  
#Flower count - good predictions, but slightly overdispersed
with(mod3,plot(apply(flwCount_resid,1,function(x) sum(abs(x))),
               apply(predFlwCount_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Flower count per plant'))
abline(0,1,col='red');  #PP plot 
plot(datalist$flwCount,apply(mod3$predFlwCount,2,median), #Predicted vs Actual - good
     ylab='Predicted flower count',xlab='Actual flower count')  
abline(0,1,col='red')

#flower survival (pod count) - better with a beta-binomial, but still not the best
with(mod3,plot(apply(podCount_resid,1,function(x) sum(abs(x))),
               apply(predPodCount_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Pods per plant'))
abline(0,1,col='red') #PP plot - good
plot(datalist$podCount,apply(mod3$predPodCount,2,median), #Predicted vs Actual - good
     ylab='Predicted number of pods',xlab='Actual number of pods'); abline(0,1,col='red')

#Partial plots
temp <- with(datalist,data.frame(int=rep(1,Nplant), #Dataframe to hold info
             logHbeeVis=with(datalist,log((hbeeVis/totalTime)+0.5))[plantIndex],
             pollen=apply(mod3$pollenPlot,2,mean)[plantIndex], #Mean pollen at plot
             plDens=apply(mod3$plDens,2,mean)[plantIndex], #Planting density
             irrig=isIrrigated[plotIndex[plantIndex]],is2015=is2015[plotIndex[plantIndex]],
             irrig2015=isIrrigated[plotIndex[plantIndex]]*is2015[plotIndex[plantIndex]],
             plSize=plantSize,plSizeIrrig=plantSize*isIrrigated[plotIndex[plantIndex]],
             logflwCount=log(flwCount),
             resid=apply(mod3$podCount_resid,2,mean),podCount,flwCount=flwCount)) %>% as.matrix()

#Model matrix
modMat <- as.matrix(with(mod3,data.frame(intFlwSurv,slopeVisitSurv,
            slopePolSurv,slopePlDensSurv,
            slopeIrrigSurv,slope2015Surv,slopeIrrig2015Surv,
            slopePlSizeSurv,slopePlSizeIrrigSurv,slopeFlwCountSurv)))
#FlowerCount
temp2 <- as.data.frame(temp) %>% #Marginalize across all vars except (log)FlowerCount
  mutate_at(vars(int:plSizeIrrig),mean) %>% as.matrix()
temp2 <- cbind(temp2,t(apply(temp2[,c(1:10)] %*% t(modMat),1,function(x) quantile(x,c(0.5,0.975,0.025))))) %>% 
  as.data.frame() %>% rename(pred='50%',upr='97.5%',lwr='2.5%') %>% 
  mutate(predUpr=apply(mod3$predPodCount,2,quantile,0.975)) %>%
  mutate(predLwr=apply(mod3$predPodCount,2,quantile,0.025)) %>%
  mutate_at(vars(pred,upr,lwr),invLogit) %>% #Transform to proportion
  mutate(pred=pred*flwCount,upr=upr*flwCount,lwr=lwr*flwCount) %>% #Tranform to predicted count
  arrange(flwCount)
ggplot(data=temp2)+
  geom_ribbon(aes(x=flwCount,ymax=upr/flwCount,ymin=lwr/flwCount),alpha=0.2)+
  geom_point(aes(flwCount,(pred-resid)/flwCount))+
  geom_line(aes(flwCount,pred/flwCount),size=1,col='red')+
  # scale_x_log10()+
  labs(x='Flowers per plant',y='Flower survival')+ylim(NA,1)

#Plant size
temp2 <- as.data.frame(temp) %>% #Marginalize across all vars except Plant size
  mutate(flwCountMeas=flwCount) %>% 
  mutate_at(vars(int:irrig2015,plSizeIrrig:flwCount),mean) %>% as.matrix()
temp2 <- cbind(temp2,t(apply(temp2[,c(1:10)] %*% t(modMat),1,function(x) quantile(x,c(0.5,0.975,0.025))))) %>% 
  as.data.frame() %>% rename(pred='50%',upr='97.5%',lwr='2.5%') %>% 
  mutate(predUpr=apply(mod3$predPodCount,2,quantile,0.975)) %>%
  mutate(predLwr=apply(mod3$predPodCount,2,quantile,0.025)) %>%
  mutate_at(vars(pred,upr,lwr),invLogit) %>% #Transform to proportion
  mutate(pred=pred*flwCountMeas,upr=upr*flwCountMeas,lwr=lwr*flwCountMeas) %>% #Tranform to predicted count
  mutate(plSize=exp(plSize+2.62)) %>% #Transform to original scale
  arrange(plSize)
ggplot(data=temp2)+
  geom_ribbon(aes(x=plSize,ymax=upr/flwCountMeas,ymin=lwr/flwCountMeas),alpha=0.2)+
  geom_point(aes(plSize,(pred-resid)/flwCountMeas))+
  geom_line(aes(plSize,pred/flwCountMeas),size=1,col='red')+
  # scale_x_log10()+
  labs(x='Plant size',y='Flower survival')+ylim(NA,1)


#Yield per plant
with(mod3,plot(apply(yield_resid,1,function(x) sum(abs(x))),
               apply(predYield_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Yield per plant'))
abline(0,1,col='red') #PP plot - good
plot(log(datalist$yield),apply(mod3$predYield,2,median), #Predicted vs Actual - good
     ylab='log Predicted yield',xlab='log Actual yield'); abline(0,1,col='red')
par(mfrow=c(1,1))

# Partial effects plots for commodity fields -----------------------------
mod3 <- extract(modPodcount) #Get coefficients from stan model

#Plant-level data (792 rows)
plantDat <- with(datalist,data.frame(
  field=plotIndex[plantIndex],
  plot=plantIndex,
  visit=(hbeeVis/totalTime)[plotIndex[plantIndex]],
  logvisit=log((hbeeVis/totalTime)+0.5)[plotIndex[plantIndex]],
  # pollen=apply(mod3$pollenPlot,2,median)[plotIndex[plantIndex]],
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
  flwCount=flwCount
  # flwCountResid=log(flwCount)-apply(mod3$flwCountMu,2,median), #Log-resid for flower count
  # podCountResid=logit(podCount/flwCount)-apply(mod3$flwSurv,2,median) #Logit-resid for flower survival
  )) 

## Model averaging
# library(MuMIn)
# library(lme4)
# temp<- with(plantDat,data.frame(plSize,plDens,logdist,GP,year2015,stocking,irrigation,field=factor(field),plot=factor(plot)))
# tempmod <- lmer(plSize~plDens*stocking*logdist-plDens:stocking:logdist+GP*year2015+stocking*year2015+
#                   year2015*irrigation+(plDens|field)+(1|plot),data=temp,na.action=na.fail,REML=F)
# summary(tempmod)
# tempmod_d <- dredge(tempmod)
# summary(model.avg(tempmod_d))

#Model matrix for plant size
MM_plSize <- with(plantDat,data.frame(int=1,plDens,logdist,GP,year2015,
                    stocking,irrigation,plDensStocking)) %>% as.matrix()
#Coefficient matrix for plant size
coef_plSize <- with(mod3,data.frame(intPlSize,slopePlDensPlSize,slopeDistPlSize,
            slopeGpPlSize,slope2015PlSize,slopeStockingPlSize,slopeIrrigPlSize,
            slopePlDensStockingPlSize))

#Partial effect of plant density and stocking
MM_temp <- MM_plSize %>% as.data.frame() %>% mutate_at(vars(-plDens,-stocking,-plDensStocking),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_plSize),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plantDat$plSizeResid+pred) %>%
  filter(stocking==0|stocking==40) %>%
  mutate(stocking=factor(stocking,labels=c('No hives','40 hives'))) %>%
  mutate(plDens=exp(plDens+mean(log(surveyAllComm$PlDens),na.rm=T))) %>% #Untransform plant density
  mutate_at(vars(upr,lwr,pred,resid),function(x,y) exp(x+mean(log(plantsAllComm$VegMass),na.rm=T))) %>% 
  ggplot(aes(x=plDens))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=stocking),alpha=0.3)+
  geom_jitter(aes(y=resid,col=stocking))+
  geom_line(aes(y=pred,col=stocking),size=1)+
  labs(x=expression(paste('Plants per ',m^2)),y='Plant size (g)')+xlim(10,100)+ylim(0,40)+
  scale_colour_manual(values=c('forestgreen','darkorange'))+
  scale_fill_manual(values=c('forestgreen','darkorange'))+
  theme(legend.position=c(0.8,0.9),legend.background=element_rect(fill='white',colour='black',linetype='solid'),
        legend.title=element_blank())

#Model matrix for flower production (per plant)
MM_flwCount <- with(plantDat,data.frame(int=1,plSize)) %>% as.matrix()

#Coefficient matrix for flower production
coef_flwCount <- with(mod3,data.frame(intFlwCount,slopePlSizeFlwCount))

#Effect of plant size on flower production
MM_temp <- MM_flwCount 

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
  mutate(resid=plantDat$flwCountResid+pred) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=plantDat$plSize))+
  # geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid),alpha=0.5)+
  geom_line(aes(y=pred),size=1,col='red')+
  labs(x='Plant size',y='Flowers per plant')

#Model matrix for flower survival (pod count)
MM_flwSurv <- with(plantDat,data.frame(int=1,logvisit,pollen,plSize,plDens,year2015,irrigation,
                                         irrig2015,plSizeIrrig)) %>% as.matrix()
#Coefficient matrix for flower survival
coef_flwSurv <- with(mod3,data.frame(intFlwSurv,slopeVisitSurv,slopePolSurv,slopePlSizeSurv,slopePlDensSurv,
                                     slopeIrrigSurv,slope2015Surv,slopeIrrig2015Surv,slopePlSizeIrrigSurv))

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

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwSurv),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$podCountResid+pred) %>% mutate_at(vars(pred:resid),invLogit) %>%
  ggplot(aes(x=plSize))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid),alpha=0.5)+ylim(0.5,0.9)+xlim(-2,2)+
  geom_line(aes(y=pred),size=1,col='red')+
  labs(x='Plant size',y='Pod survival')

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

#Pod-level data (3872 rows)
podDat <- with(datalist,data.frame(
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
  seedCountResid=
    apply(log(matrix(rep(datalist$seedCount,nrow(mod3$seedCountMu)),ncol=nrow(mod3$seedCountMu)))-t(mod3$seedCountMu),1,mean)
  )) 

#Model matrix for seedCount
MM_seedCount <- with(podDat,data.frame(int=1,visit,pollen,plSize,year2015,irrigation,irrig2015,
                                       plSizeIrrig,plSize2015,plSizeIrrig2015)) %>% as.matrix()
#Coefficient matrix for seedCount
coef_seedCount <- with(mod3,data.frame(intSeedCount,slopeVisitSeedCount,slopePolSeedCount,slopePlSizeCount,
            slope2015SeedCount,slopeIrrigSeedCount,slopeIrrig2015SeedCount,slopePlSizeIrrigSeedCount,
            slopePlSize2015SeedCount,slopePlSizeIrrig2015SeedCount))

#Partial effects plot for year
MM_temp <- MM_seedCount %>% as.data.frame() %>% mutate_at(vars(-year2015),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedCountResid+pred) %>% mutate(year2015=factor(year2015,labels=c('2014','2015'))) %>% 
  ggplot()+geom_jitter(aes(x=year2015,y=exp(resid)),width=0.25,alpha=0.3)+
  geom_pointrange(aes(x=year2015,y=exp(pred),ymax=exp(upr),ymin=exp(lwr)),col='red',size=1)+
  labs(x='Year',y='Seeds per pod')+ylim(0,40)

#Model matrix for seedMass 
MM_seedMass <- with(podDat,data.frame(int=1,visit,pollen,seedCount,plSize,irrigation,year2015,irrig2015,seedCountPlSize,
                                      plSizeIrrig,seedCount2015,plSize2015,plSizeIrrig2015)) %>% as.matrix()

#Coefficent matrix for seedMass
coef_seedMass <- with(mod3,data.frame(intSeedWeight,slopeVisitSeedWeight,slopePolSeedWeight,slopeSeedCount,
                      slopePlSizeWeight,slopeIrrigSeedWeight,slope2015SeedWeight,slope2015IrrigSeedWeight,
                      slopeSeedCountPlSizeSeedWeight,slopePlSizeIrrigSeedWeight,slopeSeedCount2015SeedWeight,
                      slopePlSize2015SeedWeight,slopePlSizeIrrig2015SeedWeight)) %>% as.matrix()

#Partial effect of seed count/plant size
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-seedCount,-plSize),mean) %>% 
  mutate(plSize=rep(c(-1,0,1),length.out=n())) %>% #10,50,90 percentile of plant size
  mutate(seedCount=rep(c(1:40),length.out=n())) %>% #1-40 seeds
  mutate(seedCountPlSize=seedCount*plSize,plSizeIrrig=plSize*irrigation,
         plSize2015=plSize*year2015,plSizeIrrig2015=plSize*irrigation*year2015) %>% 
  distinct() %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(plSize=factor(plSize,labels=c('Small','Medium','Large'))) %>% 
  ggplot()+
  geom_ribbon(aes(seedCount,ymax=upr,ymin=lwr,fill=plSize),alpha=0.3)+
  geom_line(aes(seedCount,pred,col=plSize),size=1)+
  labs(x='Seeds per pod',y='Seed Weight')

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


# Seed field visitation and pollen deposition (Stan) -----------------------------
library(rstan)
setwd('./Models')
rstan_options(auto_write = TRUE)
options(mc.cores = 4)
library(shinystan)

#List structure for Stan

#Organize names of extra sites from Riley
rileyFields <- levels(rileyExtra$site)[!(levels(rileyExtra$site) %in% levels(fieldsAllSeed$Field))]
samFields <- levels(fieldsAllSeed$Field)
fieldsAllSeed$Field <- factor(fieldsAllSeed$Field,levels=c(samFields,rileyFields)) #Apply new ordering of fields to dataframes
rileyExtra$site <- factor(rileyExtra$site,levels=c(samFields,rileyFields))  
surveyAllSeed$Field <- factor(surveyAllSeed$Field,levels=c(samFields,rileyFields))
pollenAllSeed$Field <- factor(pollenAllSeed$Field,levels=c(samFields,rileyFields))
plantsAllSeed$Field <- factor(plantsAllSeed$Field,levels=c(samFields,rileyFields))
seedsAllSeed$Field <- factor(seedsAllSeed$Field,levels=c(samFields,rileyFields))
plantsAllSeed <- plantsAllSeed %>% group_by(Field,Distance,EdgeCent) %>% mutate(Plant=1:n()) %>% ungroup()
#Fix date structure in Riley's fields
rileyExtra$date <- as.character(rileyExtra$date)
rileyExtra$date[grepl('-',rileyExtra$date)] <- as.character(as.Date(rileyExtra$date,format='%d-%b-%y')[grepl('-',rileyExtra$date)])
rileyExtra$date[grepl('/',rileyExtra$date)] <- as.character(as.Date(rileyExtra$date,format='%m/%d/%Y')[grepl('/',rileyExtra$date)])
rileyExtra$date <- as.Date(rileyExtra$date,format='%F')

#Joins plant ID to seedsAllSeed
temp <- select(plantsAllSeed,Year,Field,Distance,EdgeCent,Branch,Pods,Missing,Plant) %>% 
  unite(ID,Year:Missing)
seedsAllSeed <- seedsAllSeed %>%
  unite(ID,Year,Field,Distance,EdgeCent,Branch,Pods,Missing,remove=F) %>% 
  left_join(temp,by='ID') %>%  select(-ID)
rm(temp)

#Add average pollen per plot
surveyAllSeed <- surveyAllSeed %>% 
  unite(plot,Field,Distance,EdgeCent,remove=F) %>%
  left_join(summarize(group_by(unite(pollenAllSeed,plot,Field,Distance,EdgeCent),plot),
                      polCountPlot=log(mean(Pollen))),by='plot') 

datalistField <- list( #Field-level measurements
  Nfield=length(fieldsAllSeed$Year), #Number of fields
  Nfield_extra=length(rileyFields) #Number of extra fields from Riley
)

datalistPlot <- with(surveyAllSeed,list( #Plot-level measurements
  Nplot=length(Distance), #Number of plots
  plotIndex=as.numeric(Field), #Index for field (which field is plot from?)
  lbeeStocking=Treatment=='Double tent', #Half leafcutter stocking?
  is2016=Year==2016, #Is field from 2016?
  hbee_dist=Distance, #Distance from honeybees
  hbeeVis=hbee, #Visits by honeybees
  lbee_dist=minDist, #Centered distance from leafcutters
  lbeeVis=lbee,
  isCent=EdgeCent=='Center', #Is plot from center?
  isMBay=Bay=='M', #Is plot from M bay?
  totalTime=TotalTime/10, #Total time (mins/10)
  plotList=paste(Field,Distance,Bay,EdgeCent),
  flDens=sqrt(FlDens*4)-22, #Flower density/m2 - sqrt transform and centered
  plDens=log(PlDens)-mean(log(PlDens),na.rm=T), #Plant density - sqrt transformed and centered
  #Day of year centered around July 9
  surveyDay=as.numeric(format(fieldsAllSeed$Surveyed[match(surveyAllSeed$Field,fieldsAllSeed$Field)],format='%j'))-190,
  # plSizePlot=avgPlSize, #Average plant size (plot-level)
  polCountPlot=polCountPlot, #Average pollen count (plot-level)
  Nplot_F=sum(Bay=='F'), #Number of female plots
  plotIndex_F=match(1:length(Distance),which(Bay=='F')), #Index for female-only plots (which female plot j does plot i belong to?)
  plotIndex_F2=match(which(Bay=='F'),1:length(Distance)) #Reverse index (which plot i does female plot j belong to?)
))
datalistPlot$plotIndex_F[is.na(datalistPlot$plotIndex_F)] <- 0 #Set male plot indices to zero
datalistPlot$totalTime[is.na(datalistPlot$totalTime)] <- 0.5 #Fix one missing time point
#Join in extra data from Riley
datalistPlot <- c(datalistPlot,with(rileyExtra,list(
  Nplot_extra=length(ldist), #Number of extra plots
  plotIndex_extra=as.numeric(site), #Index for field (which field?)
  is2016_extra=Year==2016,
  lbeeStocking_extra=treatment=='Double tent',
  hbee_dist_extra=hdist,
  hbeeVis_extra=hbee_vis,
  lbee_dist_extra=ldist,
  lbeeVis_extra=lbee_vis,
  isCent_extra=rep(FALSE,length(ldist)), #Riley's plots were at edge of bay
  isMBay_extra=Bay=='Male', #Is plot from M bay?
  totalTime_extra=rep(10,length(ldist))/10, #Riley used 10 mins for everything
  flDens_extra=sqrt(flDens)-22, #Flower density - sqrt transformed and centered
  surveyDay_extra=as.numeric(format(rileyExtra$date,format='%j'))-190 #Day of year centered around July 9
)))
datalistFlw <- with(pollenAllSeed,list( #Pollen samples
  Nflw=length(Pollen), #Number of pollen samples
  flowerIndex=match(paste(Field,Distance,'F',EdgeCent),datalistPlot$plotList),#Index for flower (which plot?)
  pollenCount=Pollen
))

datalistPlant <- plantsAllSeed %>% filter(!is.na(Pods),!is.na(Missing)) %>% #Filter out plants with missing pods/flw counts
  with(.,list(
  Nplant=length(Distance), #Number of plant samples
  podCount=Pods, #Successful pods
  flwCount=Pods+Missing, #Pods + Missing (total flw production)
  totalSeedMass=SeedMass, #Weight of all seeds (g)
  # avgSeedCount=AvPodCount,avgSeedMass=AvPodMass/AvPodCount, #Average seeds per pod and weight per seed
  plantIndex=match(paste(Field,Distance,'F',EdgeCent),datalistPlot$plotList), #Index for plant (which plot?)
  plantSize=log(VegMass)-mean(log(VegMass),na.rm=T), #log weight of veg mass (g), centered
  plantList=paste(Field,Distance,'F',EdgeCent,Plant)
))

datalistPod <- seedsAllSeed %>% ungroup() %>% 
  mutate(plantList=paste(Field,Distance,'F',EdgeCent,Plant)) %>% select(-Field,-Distance,-EdgeCent,-Plant) %>% 
  mutate(inPlantList=plantList %in% datalistPlant$plantList) %>% 
  filter(inPlantList) %>% #Filter out pods where plants had missing flower counts
  filter(!is.na(PodCount),!is.na(PodMass)) %>% #Filter out pods with missing seed counts
  with(.,list(
  Npod=length(PodCount), #Number of seeds measured
  seedCount=PodCount, #Number of seeds per pod
  seedMass=1000*PodMass/PodCount, #Weight per seed (mg)
  #Index for pod (which plant?)
  podIndex=match(plantList,datalistPlant$plantList) 
))
datalistPlant$plantList <- datalistPlot$plotList <-  NULL

#No NAs in pollen

#If any plot-level extra measurements are missing
naPlot <- with(datalistPlot,is.na(hbee_dist_extra)|is.na(hbeeVis_extra)|is.na(lbee_dist_extra)) #Strip NAs from extra plot data
datalistPlot <- within(datalistPlot,{
  plotIndex_extra <- plotIndex_extra[!naPlot]
  lbeeStocking_extra <- lbeeStocking_extra[!naPlot]
  is2016_extra <- is2016_extra[!naPlot]
  hbee_dist_extra <- hbee_dist_extra[!naPlot]
  hbeeVis_extra <- hbeeVis_extra[!naPlot]
  lbee_dist_extra <- lbee_dist_extra[!naPlot]
  lbeeVis_extra <- lbeeVis_extra[!naPlot]
  isCent_extra <- isCent_extra[!naPlot]
  isMBay_extra <- isMBay_extra[!naPlot]
  totalTime_extra <- totalTime_extra[!naPlot]
  flDens_extra <- flDens_extra[!naPlot]
  surveyDay <- surveyDay[!naPlot]
  #Parameters for data imputation
  #Flower density
  Nplot_flsObs <- sum(!is.na(flDens)) #Number of plots observed
  Nplot_flsMiss <- sum(is.na(flDens)) #Number of plots missing
  obsFls_ind <- which(!is.na(flDens)) #Index for observed
  missFls_ind <- which(is.na(flDens)) #Index for missing
  flDens_obs <- flDens[!is.na(flDens)] #Observed flower density
  #Plant density
  Nplot_densObs <- sum(!is.na(plDens)) #Number of plots observed
  Nplot_densMiss <- sum(is.na(plDens)) #Number of plots missing
  obsPlDens_ind <- which(!is.na(plDens)) #Index for observed 
  missPlDens_ind <- which(is.na(plDens)) #Index for missing
  plDens_obs <- plDens[!is.na(plDens)] #Observed number of plants
  #Extra flower density
  Nplot_flsObs_extra <- sum(!is.na(flDens_extra))
  Nplot_flsMiss_extra <- sum(is.na(flDens_extra))
  obsFls_ind_extra <- which(!is.na(flDens_extra))
  missFls_ind_extra <- which(is.na(flDens_extra))
  flDens_obs_extra <- flDens_extra[!is.na(flDens_extra)] 
  # #Plot-level average plant size
  # Nplot_avgPlSizeObs <- sum(!is.na())
  # Nplot_avgPlSizeMiss
  # obsAvgPlSizeMiss_ind
  # missAvgPlSizeMiss_ind
  # avgPlSizeObs_obs
  Nplot_extra <- sum(!naPlot) #Revises number of plots
})

# #Impute missing plant data
# naPlant <- with(datalistPlant,is.na(podCount)|is.na(flwCount)|
#                   podCount>flwCount|is.na(plantSize_obs)|is.na(totalSeedMass))
# datalistPlant <-(within(datalistPlant,{
#   podCount <- podCount[!naPlant]
#   flwCount <- flwCount[!naPlant]
#   totalSeedMass <- totalSeedMass[!naPlant]
#   plantSize_obs <- plantSize_obs[!naPlant]
#   Nplant_obs <- sum(!naPlant) #Number of plants with complete measurements
#   plantSurvIndex <- plantIndex[!naPlant]
#   Nplant_miss <- sum(naPlant) #Number of plants with missing measurements
#   missPlant_ind <- which(naPlant) #Index of missing plants
#   obsPlant_ind <- which(!naPlant) #Index of non-missing plants
# }))

# #If any pod measurements are NA, or missing plant above |(datalistPod$podIndex %in% which(naPlant))
# naPod <- with(datalistPod,is.na(seedCount)|is.na(seedMass))
# datalistPod <- within(datalistPod,{
#   seedCount <- seedCount[!naPod]
#   seedMass <- seedMass[!naPod]
#   podIndex <- podIndex[!naPod]
#   # podIndex <- as.numeric(factor(podIndex)) #Re-index
#   Npod <- sum(!naPod)
# })
# datalistPlot$flDens <- datalistPlot$plDens <- datalistPlot$flDens_extra <- NULL #Remove NA fields
# datalistPlant$avgSeedCount <- datalistPlant$avgSeedMass <- NULL

datalist <- c(datalistField,datalistPlot,datalistFlw,datalistPlant,datalistPod)
rm(datalistField,datalistPlot,datalistFlw,datalistPlant,datalistPod,naPlant,naPod,rileyFields,samFields,naPlot) #Cleanup
str(datalist)

# print(modPodcount_claims21,pars=pars)
# traceplot(modPodcount_claims28,pars=c(pars,'slopeStockingFlwSurv'))
# claims <- extract(modPodcount_claims28)
# 2*(1-pnorm(abs(mean(claims[[1]])/sd(claims[[1]])),0,1)) #p-val
# with(claims,plot(flwSurv_resid,predFlwSurv_resid));abline(0,1) #Posterior predictive checks. Beta-binomial much better than binomial. Predicted slightly higher, but looks OK for now.

inits <- function() {with(datalist,list(
  #Planting density
  plDens_miss=rep(0,Nplot_densMiss),plDens_miss_extra=rep(0,Nplot_extra),
  intPlDens=0,slopeHbeeDistPlDens=0.06,
  sigmaPlDens=0.3,sigmaPlDens_field=0.4,intPlDens_field=rep(0,Nfield+Nfield_extra),
  #Plant size,
  intPlSize=0.04,slopePlDensPlSize=-0.75,
  slopeDistPlSize=0.08,sigmaPlSize=0.6,
  #Flower density per plot
  flDens_miss=rep(0,Nplot_flsMiss),flDens_extra_miss=rep(0,Nplot_flsMiss_extra),
  intFlDens=0.5,slopePlSizeFlDens=1,slopeBayFlDens=1,slope2016FlDens=-3,slopeDistFlDens=1.2,
  intFlDens_field=rep(0,Nfield+Nfield_extra),sigmaFlDens=5,sigmaFlDens_field=3.5,
  #Hbees
  intVisitHbee=1.9,slopeHbeeDistHbee=-0.25,slopeLbeeDistHbee=0.4,slopeLbeeHbeeDistHbee=0,
  slopeLbeeVisHbee=0,slopeCentHbee=0.35,slopeFlDensHbee=0,slopeFBayHbee=0.1,
  visitHbeePhi=0.7,zeroVisHbeeTheta=0.3,
  #Lbees
  intVisitLbee=4.5,slopeHbeeDistLbee=-0.3,slopeLbeeDistLbee=-0.8,slopeCentLbee=-0.6,
  slopeFBayLbee=0,slopeStocking=0,slopeCentHbeeDistLbee=-0.2,slopeStockingHbeeDistLbee=0.3,
  slopeFlDensLbee=0.03,sigmaLbeeVisField=1,visitLbeePhi=0.5,intVisitLbee_field=rep(0,Nfield+Nfield_extra),
  #Pollen
  intPol=2.5,slopeHbeePol=0.03,slopeLbeePol=0.2,slopeCentPol=-0.4,slopeHbeeDistPol=-0.12, 
  slopeFlDensPol=-0.03,sigmaPolField=0.8,sigmaPolPlot=0.5,pollenPhi=0.7,
  intPol_field=rep(0,Nfield),intPol_plot=rep(0,Nplot_F),
  #Flower count per plant
  intFlwCount=6,slopePlSizeFlwCount=0.56,slopeCentFlwCount=0,
  slopePolFlwCount=0,slopeLbeeVisFlwCount=0,slopeFlDensFlwCount=0,
  sigmaFlwCount_field=0.26,intFlwCount_field=rep(0,Nfield),flwCountPhi=5.1,
  #Flower survival
  intFlwSurv=-5,slopePolSurv=0.5,slopePlSizeSurv=0,slopeEdgeCentSurv=0.4,
  slopeSeedCountSurv=2.3,slopeHbeeDistSurv=0,slopeLbeeDistSurv=0,
  slopeFlwCountSurv=0,slopeFlwDensSurv=0,sigmaFlwSurv_field=0.3,
  sigmaFlwSurv_plot=0.5,intFlwSurv_field=rep(0,Nfield),intFlwSurv_plot=rep(0,Nplot_F),
  #Seed count
  intSeedCount=3,slopePolSeedCount=0.06,slopePlSizeCount=0.2,slopeEdgeCentSeedCount=-0.15,
  slopeHbeeDistSeedCount=-0.01,slopeFlDensSeedCount=-0.005,slopeFlwCountSeedCount=-3e-4,slopeSurvSeedCount=0.15,
  seedCountPhi=3.5,sigmaSeedCount_field=0.1,sigmaSeedCount_plot=0.1,sigmaSeedCount_plant=0.1,
  intSeedCount_field=rep(0,Nfield),intSeedCount_plot=rep(0,Nplot_F),intSeedCount_plant=rep(0,Nplant),
  #Seed weight
  intSeedWeight=3.5,slopePolSeedWeight=-0.05,slopeSeedCount=-0.035,slopePlSizeSeedWeight=0.25,
  slope2016SeedWeight=0.5,slopeLbeeDistSeedWeight=0.1,slopePlDensSeedWeight=0.3,
  slopeStockingSeedWeight=0.2,
  sigmaSeedWeight=1,sigmaSeedWeight_field=0.3,sigmaSeedWeight_plot=0.3,
  sigmaSeedWeight_plant=0.5,intSeedWeight_field=rep(0,Nfield),intSeedWeight_plot=rep(0,Nplot_F),
  intSeedWeight_plant=rep(0,Nplant)),lambdaSeedWeight=4
)}

#Full model
modPodcount_seed <- stan(file='visitation_pollen_model_seed.stan',data=datalist,
                         iter=150,chains=3,control=list(adapt_delta=0.8),init='inits')
# save(modPodcount_seed,file='modPodcount_seed.Rdata') #Seed count/weight model
modPodcount_seed2 <- stan(file='visitation_pollen_model_seed.stan',data=datalist,
                         iter=1000,chains=4,control=list(adapt_delta=0.8),init='inits')
# save(modPodcount_seed2,file='modPodcount_seed2.Rdata') #Flower count/pod success
# Setting max_treedepth=15 takes about 2-3x as long to run model. Use with care.

# #ML version, no variance estimates
# modPodcount_seed2 <- optimizing(stan_model(file='visitation_pollen_model_seed.stan'),data=datalist,init=inits)
# str(modPodcount_seed2)
# modPodcount_seed2$par[names(modPodcount_seed2$par) %in% pars] #Get optim results


pars=c('intPlDens','slopeHbeeDistPlDens',#'slopeHbeeDistSqPlDens', #Planting density
       'sigmaPlDens','sigmaPlDens_field') 
pars=c('intPlSize','slopePlDensPlSize','slopeDistPlSize','sigmaPlSize') #Plant size
pars=c('intFlDens','slopePlSizeFlDens','slopeBayFlDens',
       'slope2016FlDens','slopeDistFlDens',
       'sigmaFlDens','sigmaFlDens_field') #Flower density
pars=c('intVisitLbee','slopeHbeeDistLbee','slopeLbeeDistLbee','slopeCentLbee','slopeFBayLbee', #Lbee vis
       'slopeStockingLbee','slope2016Lbee','slopeCentHbeeDistLbee','slopeStockingHbeeDistLbee',
       'slopeFlDensLbee','sigmaLbeeVisField','visitLbeePhi')
pars=c('intVisitHbee','slopeHbeeDistHbee','slopeLbeeDistHbee','slopeLbeeHbeeDistHbee','slopeLbeeVisHbee', #Hbee vis 
       'slopeCentHbee','slopeFlDensHbee','slopeFBayHbee', 
       'visitHbeePhi','zeroVisHbeeTheta') 
pars=c('intPol','slopeHbeePol','slopeLbeePol','slopeCentPol','slopeHbeeDistPol','slopeFlDensPol', #Pollen
       'pollenPhi','sigmaPolField')
pars=c('intFlwCount','slopePlSizeFlwCount', #Flower count per plant
       'slopeCentFlwCount','slopePolFlwCount','slopeLbeeVisFlwCount',
       'slopeFlDensFlwCount','sigmaFlwCount_field','flwCountPhi')
pars=c('intFlwSurv','slopePolSurv','slopePlSizeSurv', #Flower survival
       'slopeEdgeCentSurv','slopeHbeeDistSurv','slopeLbeeDistSurv','slopeFlwCountSurv',
       'slopeFlwDensSurv','sigmaFlwSurv_field','sigmaFlwSurv_plot')
pars=c('intSeedCount','slopePolSeedCount','slopePlSizeCount', #Seeds per pod
       'slopeEdgeCentSeedCount','slopeHbeeDistSeedCount','slopeFlDensSeedCount',
       'slopeSurvSeedCount','seedCountPhi','sigmaSeedCount_field','sigmaSeedCount_plot','sigmaSeedCount_plant')
pars=c('intSeedWeight','slopePolSeedWeight','slopeSeedCount', #Weight per seed
       'slopePlSizeSeedWeight',
       'slope2016SeedWeight','slopeLbeeDistSeedWeight','slopePlDensSeedWeight',
       'slopeStockingSeedWeight','lambdaSeedWeight',
       'sigmaSeedWeight','sigmaSeedWeight_field', #Random effects for plot don't converge well
       'sigmaSeedWeight_plot')#'sigmaSeedWeight_plant' 
stan_hist(modPodcount_seed,pars=pars)+geom_vline(xintercept=0,linetype='dashed')
traceplot(modPodcount_seed,pars=pars,inc_warmup=F)+geom_hline(yintercept=0,linetype='dashed')
print(modPodcount_seed,pars=pars)

#Check model fit:
mod3 <- extract(modPodcount_seed)
coefs(mod3[pars])

#Distribution of random effects intercepts
t(apply(mod3$intSeedCount_plot,2,function(x) quantile(x,c(0.5,0.025,0.975)))) %>% 
  as.data.frame() %>% rename(med='50%',lwr='2.5%',upr='97.5%') %>% 
  arrange(med) %>% mutate(plot=1:n()) %>% 
  ggplot()+geom_pointrange(aes(x=plot,y=med,ymax=upr,ymin=lwr))+
  geom_hline(yintercept=0,col='red')

#Faster pair plots
pairs(mod3[c(pars,'lp__')],lower.panel=function(x,y){
  par(usr=c(0,1,0,1))
  text(0.5, 0.5, round(cor(x,y),2), cex = 1 * exp(abs(cor(x,y))))})

par(mfrow=c(2,1))
#planting density - good
with(mod3,plot(apply(plDens_resid,1,function(x) sum(abs(x))), #PP plot - good
               apply(predPlDens_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Planting density'))
abline(0,1,col='red') 
plot(datalist$plDens_obs,apply(mod3$predPlDens,2,median)[datalist$obsPlDens_ind], #Predicted vs Actual
     ylab='Predicted plant density',xlab='Actual plant density'); abline(0,1,col='red');

#plant size - OK
with(mod3,plot(apply(plSize_resid,1,function(x) sum(abs(x))),
               apply(predPlSize_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Plant size')) 
abline(0,1,col='red'); #PP plot - good
plot(datalist$plantSize_obs,apply(mod3$predPlSize,2,median)[datalist$obsPlant_ind], #Predicted vs Actual - good
     ylab='Predicted plant size',xlab='Actual plant size'); abline(0,1,col='red')

#hbee visits - not the best, but OK
with(mod3,plot(apply(hbeeVis_resid,1,function(x) sum(abs(x))),
               apply(predHbeeVis_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Hbee visits')) 
abline(0,1,col='red') #PP plot - OK
plot(with(datalist,c(hbeeVis,hbeeVis_extra)),apply(mod3$predHbeeVis_all,2,median), #Predicted vs Actual - good
     ylab='Predicted visits',xlab='Actual visits')
abline(0,1,col='red')

#lbee visits - not the best, but OK
#PP plots - OK, but not the best
# Negbin (ZI or regular) is not good, but Poisson (ZI or regular) is far worse, and traces for intercept are bad. Sticking with regular NB for now.
with(mod3,plot(apply(lbeeVis_resid,1,function(x) sum(abs(x))),
               apply(predLbeeVis_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Lbee visits')); 
abline(0,1,col='red') #PP plot - not so good
plot(with(datalist,c(lbeeVis,lbeeVis_extra)),apply(mod3$predLbeeVis_all,2,median), #Predicted vs Actual - OK
     ylab='Predicted visits',xlab='Actual visits')
abline(0,1,col='red') #PP plot - not so good

#pollen - good
with(mod3,plot(apply(pollen_resid,1,function(x) sum(abs(x))),
              apply(predPollen_resid,1,function(x) sum(abs(x))),
              xlab='Sum residuals',ylab='Sum simulated residuals',main='Pollen counts')) 
abline(0,1,col='red') #PP plot
plot(datalist$pollenCount,apply(mod3$predPollenCount,2,median), #Predicted vs Actual 
     ylab='Predicted pollen',xlab='Actual pollen')
abline(0,1,col='red') #PP plot

#seeds per pod - not good, but close
with(mod3,plot(apply(seedCount_resid,1,function(x) sum(abs(x))),
               apply(predSeedCount_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Seeds per pod')) 
abline(0,1,col='red') #PP plot - not good
plot(datalist$seedCount,apply(mod3$predSeedCount,2,median), #Predicted vs Actual - good
     ylab='Predicted seed count',xlab='Actual seed count'); abline(0,1,col='red')

#weight per seed - good
with(mod3,plot(apply(seedMass_resid,1,function(x) sum(abs(x))),
               apply(predSeedMass_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Weight per seed'))
abline(0,1,col='red') #PP plot - bad
plot(datalist$seedMass,apply(mod3$predSeedMass,2,median), #Predicted vs Actual - OK, but weird outliers: check data
     ylab='Predicted seed weight',xlab='Actual seed weight'); abline(0,1,col='red'); 





#Flower count - good, but close
with(mod3,plot(apply(flwCount_resid,1,function(x) sum(abs(x))),
               apply(predFlwCount_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Flower count per plant'))
abline(0,1,col='red');  #PP plot - OK
plot(datalist$flwCount,apply(mod3$predFlwCount,2,median), #Predicted vs Actual - good, but variance increases
     ylab='Predicted flower count',xlab='Actual flower count',ylim=c(0,3000))  
abline(0,1,col='red')
#flower survival (pod count) - good
with(mod3,plot(apply(podCount_resid,1,function(x) sum(abs(x))),
               apply(predPodCount_resid,1,function(x) sum(abs(x))),
               xlab='Sum residuals',ylab='Sum simulated residuals',main='Pods per plant'))
abline(0,1,col='red') #PP plot - good
plot(datalist$podCount,apply(mod3$predPodCount,2,median), #Predicted vs Actual - good
     ylab='Predicted number of pods',xlab='Actual number of pods'); abline(0,1,col='red')

par(mfrow=c(1,1))

# Partial effects plots for seed fields -----------------------------------

#Residual plots (normalQQ,equal var)
with(mod3,plot(apply(seedWeightMu,2,mean),apply(seedMass_resid,2,mean)))
with(mod3,qqnorm(apply(seedMass_resid,2,mean)));with(mod3,qqline(apply(seedMass_resid,2,mean)))

#Pod-level data (2885 rows)
podDat <- with(datalist,data.frame(
  pollen=tapply(apply(mod3$pollenMu-matrix(rep(mod3$intPol,1050),ncol=1050),2,median),flowerIndex,mean)[
    plotIndex_F[plantIndex[podIndex]]],
  plSize=plantSize[podIndex], #Plant size
  edge=isCent[plantIndex[podIndex]], 
  lbeeDist=(log(c(lbee_dist,lbee_dist_extra))-mean(log(c(lbee_dist,lbee_dist_extra))))[plantIndex[podIndex]],
  flDens=apply(mod3$flDens,2,mean)[plantIndex[podIndex]], plDens=apply(mod3$plDens,2,mean)[plantIndex[podIndex]],
  stocking=lbeeStocking[plantIndex[podIndex]],surv=logit(podCount/flwCount)[plantIndex[podIndex]],
  year2016=is2016[plantIndex[podIndex]],seedCount=seedCount,seedWeight=seedMass,
  seedMassResid=apply(mod3$seedMass_resid,2,median) #Median resid for seedMass
))

#Model matrix for seedMass 
MM_seedMass <- with(podDat,data.frame(int=1,pollen,seedCount,plSize,year2016,lbeeDist,plDens,stocking)) %>% 
  as.matrix()
#Coefficent matrix for seedMass
coef_seedMass <- with(mod3,data.frame(intSeedWeight,slopePolSeedWeight,slopeSeedCount,slopePlSizeSeedWeight,slope2016SeedWeight,
                      slopeLbeeDistSeedWeight,slopePlDensSeedWeight,slopeStockingSeedWeight)) %>% as.matrix()

#Partial effect of seed count
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-seedCount),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred) %>% arrange(seedCount) %>% 
  ggplot()+
  geom_ribbon(aes(x=seedCount,ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(seedCount,pred),size=1)+
  geom_point(aes(seedCount,resid),alpha=0.5)+
  ylim(0,7.5)+labs(x='Seeds per pod',y='Weight per seed')

#Partial effect of plant size
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-plSize),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred) %>% arrange(plSize) %>% 
  ggplot()+
  geom_ribbon(aes(x=plSize,ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(plSize,pred),size=1)+
  geom_point(aes(plSize,resid),alpha=0.5)+
  ylim(0,7.5)+labs(x='Plant size',y='Weight per seed')

#Partial effect of year
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-year2016),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% mutate(year2016=factor(year2016,labels=c('2015','2016'))) %>% 
  mutate(resid=podDat$seedMassResid+pred) %>% 
  ggplot()+
  # geom_jitter(aes(year2016,resid),alpha=0.5,width=0.25)+
  geom_violin(aes(year2016,resid),fill='gray50')+
  geom_pointrange(aes(x=year2016,y=pred,ymax=upr,ymin=lwr),size=1)+
  ylim(0,7.5)+labs(x='Year',y='Weight per seed')

#Partial effect of plant density
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-plDens),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=podDat$seedMassResid+pred) %>% arrange(plSize) %>% 
  ggplot(aes(x=plDens))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  geom_point(aes(y=resid),alpha=0.5)+
  ylim(0,7.5)+labs(x='Planting Density',y='Weight per seed')

#Partial plot of both plant density and plant size (I think they're related)
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-plSize,-plDens),mean) %>% 
  as.matrix()
data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
  mutate(resid=podDat$seedMassResid+pred) %>%
  arrange(plSize,plDens) %>%
  # ggplot()+geom_point(aes(x=plSize,y=plDens,col=pred))
  ggplot()+geom_point(aes(x=plSize,y=plDens,col=pred))+xlim(-2,2)+
  labs(x='Plant size',y='Plant density',col='Weight per seed')

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
