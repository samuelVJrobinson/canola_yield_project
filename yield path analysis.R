#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN COMMODITY CANOLA FIELDS (2014+2015)

# Libraries and ggplot theme ---------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
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
#Faster pair plot
fastPairs <- function(l){ #List l
  pairs(l,lower.panel=function(x,y){
  par(usr=c(0,1,0,1))
  text(0.5, 0.5, round(cor(x,y),2), cex = 1 * exp(abs(cor(x,y))))})
}


#Linear breakpoint function - two lines with intersection "b"
bpoint <- function(x,int1,slope1,b,slope2) ifelse(x<b,int1+slope1*x,b*slope1+(x-b)*slope2+int1)
#Effect size for Posterior samples
effSize <- function(x) unname(median(x)/diff(quantile(x,c(0.025,0.975))))
#Does 95% of posterior overlap zero?
overlap <- function(x) {r <- quantile(x,c(0.025,0.975))>=0; xor(r[1],r[2]);}
#Posterior predictive check plots
PPplots <- function(resid,predResid,actual,pred,main=NULL){
  par(mfrow=c(2,1))
  plot(resid,predResid,xlab='Sum residuals',ylab='Sum simulated residuals',main=main)
  x <- sum(resid<predResid)/length(resid)
  legend('topleft',paste('p =',round(min(x,1-x),3)))
  abline(0,1,col='red') #PP plot
  plot(actual,pred, #Predicted vs Actual - good
       ylab=paste('Predicted',main),xlab=paste('Actual',main)) 
  abline(0,1,col='red')
  par(mfrow=c(1,1))
}
#Convert g/m2 to bushels/acre
g2bushels <- function(x){ 
  x*4046.86/22679.6 
  # 453.592*50/bushel - using 50 lbs/bushel estimate
  # 44.0920 #bushels per tonne
  # 2204.62 #lbs per tonne
}

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

datalistPlant <- 
  plantsAllComm %>% ungroup() %>% 
  filter(SeedMass!=0) %>%
  filter(!is.na(Pods),!is.na(Missing),!is.na(AvPodCount),!is.na(AvPodMass)) %>%
  # filter(!is.na(Pods),!is.na(Missing)) %>%
  mutate(plantIndex=as.numeric(factor(paste(Year,Field,Distance)))) %>% #Index for plant (which plot?)
  arrange(plantIndex) %>% 
  with(list(Nplant=length(VegMass), #Number of plant samples (some missing)
  # Nplant_obs=sum(!is.na(VegMass)), #Observed plants
  # Nplant_miss=sum(is.na(VegMass)), #Missing plants
  podCount=Pods, #Successful pods
  flwCount=Pods+Missing, #Pods + Missing (total flw production)
  plantIndex=plantIndex, #Index for plant (which plot?)
  plantSize=log(VegMass[!is.na(VegMass)])-mean(log(VegMass[!is.na(VegMass)])), #log weight of veg mass (g), centered
  #Averaged seeds per pod and weight per seed
  avgSeedCount=AvPodCount,
  avgSeedMass=(AvPodMass/AvPodCount)*1000, #Weight per seed (mg)
  yield=SeedMass[!is.na(SeedMass)] #observed weight of all seeds (g)
  # obsPl_ind=which(!is.na(VegMass)),missPl_ind=which(is.na(VegMass)), #Indices for missing plant size
  # obsYield_ind=which(!is.na(SeedMass)), missYield_ind=which(is.na(SeedMass)) #Indices for missing yield
)) 

#Problem: plots exist at the pod level which do not exist at plant level
a <- plantsAllComm %>% ungroup() %>% filter(!is.na(Pods),!is.na(Missing),!is.na(AvPodCount),!is.na(AvPodMass)) %>%  
  filter(SeedMass!=0) %>%  transmute(index=factor(paste(Year,Field,Distance,Plant))) %>% 
  distinct() #Index for plant (which plot?)
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
#Adds average flower count per plant, and average plant weight (plot-level), used in claim 6, 12, 15
datalistPlot <- c(datalistPlot,with(datalistPlant,
  list( #Both terms are missing from plots 248-271, so use the same N and indices
    Nplot_flwCountObs = length(unique(plantIndex)),
    Nplot_flwCountMiss = datalistPlot$Nplot-length(unique(plantIndex)),
    flwCountPlot_obs = log(unname(tapply(flwCount,plantIndex,mean))), #Log-transformed average flower count
    plSizePlot_obs = unname(tapply(plantSize,plantIndex,mean)), #Plant size
    obsFlwCount_ind = unique(plantIndex),
    missFlwCount_ind = (1:datalistPlot$Nplot)[!(1:datalistPlot$Nplot %in% unique(plantIndex))]
  )
))

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
    #Pollen deposition
    intPollen=5.5,slopeVisitPol=0,slopeHbeeDistPollen=0,slopeStockingPollen=0, 
    slopeStockingHbeeDistPollen=0,sigmaPolField=0.5,sigmaPolPlot=0.4,pollenPhi=0.7,
    intPollen_field=rep(0,datalist$Nfield),intPollen_plot=rep(0,datalist$Nplot),
    #Flower count per plant
    intFlwCount=1,slopePlSizeFlwCount=0,slopeSurvFlwCount=0,slope2015FlwCount=0,
    sigmaFlwCount_field=0.5,sigmaFlwSurv_plot=0.5, 
    phiFlwCount_field=0.1,intPhiFlwCount=0,slopePlSizePhiFlwCount=0,sigmaPhiFlwCount_field=0.1,
    intPhiFlwCount_field=rep(0,Nfield),
    intFlwSurv_plot=rep(0,Nplot),intFlwCount_field=rep(0,Nfield),flwCountPhi=1,
    #Flower survival
    intFlwSurv=1,slopeVisitSurv=0,slopePolSurv=0,slopePlSizeSurv=0.02, 
    slopeIrrigSurv=0,slope2015Surv=0,slopeIrrig2015Surv=0,slopePlSizeIrrigSurv=0,
    slopeFlwCountSurv=0,slopePlDensSurv=0,sigmaFlwSurv_plot=0.3,sigmaFlwSurv_field=0.3, 
    intFlwSurv_field=rep(0,Nfield),intFlwSurv_plot=rep(0,Nplot),flwSurvPhi=1,
    slopePlSizePlDensSurv=0,intPhiFlwSurv=0,slopePlSizePhiFlwSurv=0,sigmaPhiFlwSurv_field=0.1,
    intPhiFlwSurv_field=rep(0,Nfield),
    #Seed count
    intSeedCount=3.15,slopeVisitSeedCount=0,slopePolSeedCount=0,
    slopePlSizeCount=0.01,slope2015SeedCount=0.15,
    seedCountPhi=21, sigmaSeedCount_plant=0.14,sigmaSeedCount_plot=0.08,sigmaSeedCount_field=0.06,
    intSeedCount_field=rep(0,datalist$Nfield),intSeedCount_plot=rep(0,datalist$Nplot),
    intSeedCount_plant=rep(0,datalist$Nplant),
    #Seed weight
    intSeedWeight=2,slopeVisitSeedWeight=0,slopePolSeedWeight=0,slopeSeedCount=0.02, 
    slopePlSizeWeight=0,slopeIrrigSeedWeight=0,slope2015SeedWeight=0.5,
    slope2015IrrigSeedWeight=0,slopeSeedCountVegMassSeedWeight=0,
    slopePlSizeIrrigSeedWeight=0,slopeSeedCount2015SeedWeight=0,slopePlSize2015SeedWeight=0,
    slopePlSizeIrrig2015SeedWeight=0,
    sigmaSeedWeight=0.5,sigmaSeedWeight_plant=0.2,sigmaSeedWeight_plot=0.5,sigmaSeedWeight_field=0.4,
    intSeedWeight_field=rep(0,datalist$Nfield),slopePlSizeSeedWeight=0,
    intSeedWeight_plot=rep(0,datalist$Nplot),intSeedWeight_plant=rep(0,datalist$Nplant),
    lambdaSeedWeight=1.5,
    #Yield
    intYield=-0.3,slopeYield=1,sigmaYield=0.25,
    sigmaYield_field=c(0.1,0.05),sigmaYield_plot=c(0.4,0.17),
    L_field=t(chol(matrix(c(1,-0.5,-0.5,1),ncol=2))),
    L_plot=t(chol(matrix(c(1,-0.5,-0.5,1),ncol=2))),
    zYield_field=matrix(rep(0,datalist$Nfield*2),nrow=2),
    zYield_plot=matrix(rep(0,datalist$Nplot*2),nrow=2)
   ))
}

# #Feed datalist into stan_rdump for use in CmdStan
# with(datalist,stan_rdump(names(datalist),'tempDat.data.R'))
# #Feed inits into stan_rdump
# temp <- inits()
# with(temp,stan_rdump(names(temp),'inits.data.R'))
# 
# #Claims list
# claims1 <- stan(file='./Commodity model claims 2/commodity_claims39.stan',
#                 data=datalist,iter=10,chains=1,control=list(adapt_delta=0.8),init=inits)
# beep(1)
# pars <- c('intPollen','slopeVisitPol','slopeHbeeDistPollen',
#           'sigmaPolField','pollenPhi','sigmaPolPlot')	
# newpar <- 'slopeStockingPollen'

stan_trace(claims1,pars=c(pars,newpar))
# pairs(claims1,pars=c(pars,newpar)) #Takes a long time
mod1 <- extract(claims1)
coefs(mod1[c(pars,newpar)])

#Full model - 1.7 hrs for 1000 iter
modPodcount <- stan(file='visitation_pollen_model.stan',data=datalist,iter=1,chains=1,
                   control=list(adapt_delta=0.8),init=inits)
beep(1)
# save(modPodcount,file='modPodcount.Rdata')
load('modPodcount.Rdata') #Load all parameters

# print(modPodcount)
pars=c('intPlDens','slope2015PlDens','slopeIrrigPlDens','slope2015IrrigPlDens',
       'slopeDistPlDens','slopeGPPlDens','sigmaPlDens','sigmaPlDens_field') #Planting density
pars=c('intPlSize','slopePlDensPlSize','slopeDistPlSize','slopeGpPlSize', #Plant size
       'slopeIrrigPlSize','slope2015PlSize',
       'sigmaPlSize_field','sigmaPlSize_plot','sigmaPlSize')
pars=c('intFlDens','slopePlSizeFlDens','slopeHbeeDistFlDens','sigmaFlDens','sigmaFlDens_field') #Flower density
pars=c('intVisit','slopeYearVis','slopeGpVis','slopeYearGpVis','slopeIrrigVis',
       'slopeDistVis','slopeHiveVis','slopeFlDens', #Visitation
       'sigmaVisField','lambdaVisField','visitHbeePhi')
pars=c('intPollen','slopeVisitPol','slopeHbeeDistPollen',#'slopeFlyVisPol',
       'sigmaPolField','sigmaPolPlot','pollenPhi') #Pollen deposition
pars=c('intFlwCount','slopePlSizeFlwCount','slopeSurvFlwCount','slope2015FlwCount', #Flower count per plant
       'phiFlwCount_field','intPhiFlwCount','slopePlSizePhiFlwCount','sigmaPhiFlwCount_field')
pars=c('intFlwSurv','slopeVisitSurv','slopePolSurv','slopePlSizeSurv',
       'slopePlDensSurv','slopeIrrigSurv','slope2015Surv','sigmaFlwSurv_field',
       # 'flwSurvPhi') #Flower survival
       'intPhiFlwSurv','slopePlSizePhiFlwSurv','sigmaPhiFlwSurv_field')
pars=c('intSeedCount','slopeVisitSeedCount','slopePolSeedCount','slopePlSizeCount',
       'slope2015SeedCount','seedCountPhi','sigmaSeedCount_plant','sigmaSeedCount_field') #Seed count
pars=c('intSeedWeight','slopeVisitSeedWeight','slopePolSeedWeight',#Seed weight
       'slopeSeedCount','slopePlSizeWeight','slopeIrrigSeedWeight',
       'slope2015SeedWeight','slope2015IrrigSeedWeight','sigmaSeedWeight',
       'sigmaSeedWeight_plant','sigmaSeedWeight_field','lambdaSeedWeight')
pars=c('intYield','slopeYield','sigmaYield',
       'sigmaYield_field[1]','sigmaYield_field[2]','sigmaYield_plot[1]','sigmaYield_plot[2]',
       'L_field','L_plot')
stan_hist(modPodcount,pars=pars)+geom_vline(xintercept=0,linetype='dashed')
traceplot(modPodcount,pars=c(pars),inc_warmup=F)
# launch_shinystan(modPodcount)
# print(modPodcount,pars=pars) #Takes way too long
# pairs(modPodcount,pars=pars) #Takes way too long

mod3 <- extract(modPodcount) #Extract values
#Replace Cholesky matrices with covariance term
storage <- rep(NA,nrow(mod3$L_field))
for(i in 1:nrow(mod3$L_field)){
  storage[i] <- (mod3$L_field[i,,] %*% t(mod3$L_field[i,,]))[1,2]
}
mod3$L_field <- storage #Field level term
for(i in 1:nrow(mod3$L_plot)){
  storage[i] <- (mod3$L_plot[i,,] %*% t(mod3$L_plot[i,,]))[1,2]
}
mod3$L_plot <- storage


#Coefficients in table form for LaTeX
(mod3coefs <- data.frame(par='Plant size',parname=rownames(coefs(mod3[pars])),coefs(mod3[pars]),row.names=NULL))
library(xtable)
print(xtable(mod3coefs,digits=c(0,0,0,3,3,3,3,3,3,0,4)),include.rownames=F)

#Faster Pairplots
pairs(mod3[c(pars,'lp__')],lower.panel=function(x,y){
  par(usr=c(0,1,0,1))
  text(0.5, 0.5, round(cor(x,y),2), cex = 1 * exp(abs(cor(x,y))))})

#Plot of random intercepts
t(apply(mod3$intSigmaFlwCount_field,2,function(x) quantile(x,c(0.5,0.975,0.025)))) %>%
  as.data.frame() %>% rename(median='50%',upr='97.5%',lwr='2.5%') %>% arrange(median) %>%
  mutate(row=1:nrow(.)) %>% 
  ggplot(aes(row,median))+geom_pointrange(aes(ymax=upr,ymin=lwr))+geom_hline(yintercept=0,col='red')
  # ggplot(aes(median))+geom_density()

plot(apply(mod3$intFlwCount_field,2,median),apply(mod3$intSigmaFlwCount_field,2,median))
qqnorm(apply(mod3$intPollen_plot,2,median));qqline(apply(mod3$intPollen_plot,2,median));
mean(apply(mod3$intVisit_field,2,median))

#Plots of random intercepts/slopes for yield; not sure if the negative correlation actually means much
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
with(mod3,PPplots(apply(plDens_resid,1,function(x) sum(abs(x))),
                  apply(predPlDens_resid,1,function(x) sum(abs(x))),
                  apply(plDens,2,mean),apply(predPlDens,2,median),main='Plant density'))

#plant size - good
with(mod3,PPplots(apply(plSize_resid,1,function(x) sum(abs(x))),
                  apply(predPlSize_resid,1,function(x) sum(abs(x))),
                  datalist$plantSize,apply(predPlSize,2,median),main='Plant size'))

#flower density per plot
with(mod3,PPplots(apply(flDens_resid,1,function(x) sum(abs(x))),
                  apply(predFlDens_resid,1,function(x) sum(abs(x))),
                  datalist$flDens,apply(predFlDens,2,median),main='Flower density'))

#hbee visits - good
with(mod3,PPplots(apply(hbeeVis_resid,1,function(x) sum(abs(x))),
                  apply(predHbeeVis_resid,1,function(x) sum(abs(x))),
                  with(datalist,hbeeVis/totalTime),apply(predHbeeVis,2,median),main='Visits per 10mins'))

#pollen - OK, but plot level random effects throw off PP checks
with(mod3,PPplots(apply(pollen_resid,1,function(x) sum(abs(x))),
                  apply(predPollen_resid,1,function(x) sum(abs(x))),
                  datalist$pollenCount,apply(predPollenCount,2,median),main='Pollen per stigma'))

#Flower count - good
with(mod3,PPplots(apply(flwCount_resid,1,function(x) sum(abs(x))),
                  apply(predFlwCount_resid,1,function(x) sum(abs(x))),
                  datalist$flwCount,apply(predFlwCount,2,median),main='Flowers per plant'))

#flower survival (pod count) - good
with(mod3,PPplots(apply(podCount_resid,1,function(x) sum(abs(x))),
                  apply(predPodCount_resid,1,function(x) sum(abs(x))),
                  datalist$podCount,apply(predPodCount,2,median),main='Pods per plant'))

#seeds per pod - bad
with(mod3,PPplots(apply(seedCount_resid,1,function(x) sum(abs(x))),
                  apply(predSeedCount_resid,1,function(x) sum(abs(x))),
                  datalist$seedCount,apply(predSeedCount,2,median),main='Seeds per pod'))

#weight per seed - good
with(mod3,PPplots(apply(seedMass_resid,1,function(x) sum(abs(x))),
                  apply(predSeedMass_resid,1,function(x) sum(abs(x))),
                  datalist$seedMass,apply(mod3$predSeedMass,2,median),main='Weight per seed (TKW)'))

#Yield per plant - good
with(mod3,PPplots(apply(yield_resid,1,function(x) sum(abs(x))),
                  apply(predYield_resid,1,function(x) sum(abs(x))),
                  datalist$yield,exp(apply(mod3$predYield,2,median)),main='Total yield'))

# Partial effects plots for commodity fields -----------------------------
mod3 <- extract(modPodcount) #Get coefficients from stan model
# load('modPodcount3.Rdata') #All extracted coefficients in list form (mod3)

#Plot-level data (271 rows)
plotDat <- with(datalist,data.frame( 
  field=plotIndex,year=is2015[plotIndex],irrigation=isIrrigated[plotIndex],
  dist=log(dist)-mean(log(dist)), #Centered on 3.412705
  GP=isGP[plotIndex],hives=log(numHives[plotIndex]+1), #Log-number of hives
  plDens=apply(mod3$plDens,2,mean),flDens=flDens,
  hbeeResid=log(hbeeVis+1/totalTime)-apply(mod3$visitHbeeMu,2,median) #Log-residuals
))

#Model matrix for hbee visits
MM_hbee <- with(plotDat,data.frame(int=1,year,GP,yearGP=year*GP,dist,hives,flDens,irrigation)) %>% as.matrix()
#Coefficient matrix for plant size
coef_hbee <- with(mod3,data.frame(intVisit,slopeYearVis,slopeGpVis,slopeYearGpVis,
                                  slopeDistVis,slopeHiveVis,slopeFlDens,slopeIrrigVis))

#Partial effect of year/area, distance
MM_temp <- MM_hbee %>% as.data.frame() %>% mutate_at(vars(-year,-GP,-yearGP,-dist),mean) %>% 
   as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_hbee),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=plotDat$hbeeResid+pred) %>% mutate_at(vars(pred:resid),exp) %>% 
  mutate(year=factor(year,labels=c('2014','2015')),GP=factor(GP,labels=c('Lethbridge','Grand Prairie'))) %>% 
  mutate(yearGP=factor(paste(GP,year))) %>%
  ggplot(aes(x=exp(dist+3.41)))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_jitter(aes(y=resid))+
  geom_line(aes(y=pred),size=1)+facet_grid(year~GP)+ylim(0,10)+
  labs(x='Distance from edge(m)',y='Visits/10 mins')
ggsave('../Figures/Commodity/slopeYearGPVis.png',p1)


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
  seedWeight=seedMass,seedWeightResid=apply(mod3$seedMass_resid,2,median) #Median resid for seedMass
  # seedCountResid=apply(mod3$seedCount_resid,2,median) #Residuals on real scale
  # seedCountResid=apply(log(matrix(rep(datalist$seedCount,nrow(mod3$seedCountMu)),ncol=nrow(mod3$seedCountMu)))-t(mod3$seedCountMu),1,mean)
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
MM_seedMass <- with(podDat,data.frame(int=1,visit,pollen,seedCount=seedCount,plSize,irrigation,year2015,irrig2015)) %>% as.matrix()

#Coefficent matrix for seedMass
coef_seedMass <- with(mod3,data.frame(intSeedWeight,slopeVisitSeedWeight,slopePolSeedWeight,slopeSeedCount,
                      slopePlSizeWeight,slopeIrrigSeedWeight,slope2015SeedWeight,slope2015IrrigSeedWeight)) %>% as.matrix()

#Partial effect of seed count
MM_temp <- MM_seedMass %>% as.data.frame() %>% mutate_at(vars(-seedCount),mean) %>% 
  as.matrix()

temp <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_seedMass),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>%
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%
  mutate(resid=podDat$seedWeightResid+pred) 
temp2 <- data.frame(x=rep(temp$seedCount,each=nrow(mod3$predSeedMass_resid)),y=rep(temp$pred,each=nrow(mod3$predSeedMass_resid))+as.vector(mod3$predSeedMass_resid)) %>% 
  group_by(x) %>% summarize(upr=quantile(y,0.975),lwr=quantile(y,0.025)) %>% ungroup()
  
ggplot(temp,aes(x=seedCount))+
  geom_ribbon(data=temp2,aes(x=x,ymax=upr,ymin=lwr),alpha=0.3)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  geom_jitter(aes(y=resid),width=0.2,alpha=0.3)+xlim(0,40)+
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
  plDens=log(PlDens)-mean(log(PlDens),na.rm=T), #Plant density - log transformed and centered
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
  yield=SeedMass, #Weight of all seeds (g)
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
datalist <- datalist[!sapply(datalist,function(x) any(is.na(x)))] #Get rid of indices with NAs in them
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
  intFlDens=0.5,slopePlSizeFlDens=1,slopeMBayFlDens=1,slope2016FlDens=-3,slopeDistFlDens=1.2,
  intFlDens_field=rep(0,Nfield+Nfield_extra),sigmaFlDens=5,sigmaFlDens_field=3.5,
  #Hbees
  intVisitHbee=1.9,slopeHbeeDistHbee=-0.25,slopeLbeeDistHbee=0.4,slopeLbeeHbeeDistHbee=0,
  slopeLbeeVisHbee=0,slopeCentHbee=0.35,slopeFlDensHbee=0,slopeMBayHbee=0.1,
  visitHbeePhi=0.7,zeroVisHbeeTheta=0.3,
  #Lbees
  intVisitLbee=4.5,slopeHbeeDistLbee=-0.3,slopeLbeeDistLbee=-0.8,slopeCentLbee=-0.6,
  slopeMBayLbee=0,slopeStocking=0,slopeCentHbeeDistLbee=-0.2,slopeStockingHbeeDistLbee=0.3,
  slopeFlDensLbee=0.03,sigmaLbeeVisField=1,visitLbeePhi=0.5,intVisitLbee_field=rep(0,Nfield+Nfield_extra),
  #Pollen
  intPol=2.4,slopeHbeePol=0.03,slopeLbeePol=0.2,slopeCentPol=-0.6,slopeHbeeDistPol=-0.15, 
  slopeFlDensPol=-0.023,sigmaPolField=0.8,sigmaPolPlot=0.6,pollenPhi=0.82,
  intPol_field=rep(0,Nfield),intPol_plot=rep(0,Nplot_F),
  #Flower count per plant
  intFlwCount=5.9,slopePlSizeFlwCount=0.9,slopeCentFlwCount=0.1,
  # slopePolFlwCount=-0.035,slopeFlDensFlwCount=0.005,
  slopeFlwSurvFlwCount=0,
  sigmaFlwCount_field=0.09,sigmaFlwCount_plot=0.1,
  intFlwCount_field=rep(0,Nfield),intFlwCount_plot=rep(0,Nplot_F),flwCountPhi=35,
  #Flower survival
  intFlwSurv=0.75,slopePolSurv=0.13,slopePlSizeSurv=0.2,slopeEdgeCentSurv=-0.24,
  slopeHbeeDistSurv=-0.1,slopeLbeeDistSurv=-0.16,slopeFlwDensSurv=-0.01,sigmaFlwSurv_field=0.35,
  sigmaFlwSurv_plot=0.26,intFlwSurv_field=rep(0,Nfield),intFlwSurv_plot=rep(0,Nplot_F),
  #Seed count
  intSeedCount=2.7,slopePolSeedCount=0.05,slopePlSizeCount=0.07,slopeEdgeCentSeedCount=-0.15,
  slopeHbeeDistSeedCount=-0.01,slopeFlDensSeedCount=-0.006,
  slopeSurvSeedCount=0.15,seedCountPhi=3.5,sigmaSeedCount_field=0.1,
  sigmaSeedCount_plot=0.1,sigmaSeedCount_plant=0.1,
  intSeedCount_field=rep(0,Nfield),intSeedCount_plot=rep(0,Nplot_F),intSeedCount_plant=rep(0,Nplant),
  #Seed weight
  intSeedWeight=3.6,slopePolSeedWeight=0,slopeSeedCount=-0.035,slopePlSizeSeedWeight=0.25,
  slope2016SeedWeight=0.5,slopeLbeeDistSeedWeight=0.1,slopePlDensSeedWeight=0.3,
  slopeStockingSeedWeight=0.2,slopePlDensPlSizeSeedWeight=0,
  sigmaSeedWeight=1.03,sigmaSeedWeight_field=0.3,sigmaSeedWeight_plot=0.3,
  sigmaSeedWeight_plant=0.6,intSeedWeight_field=rep(0,Nfield),intSeedWeight_plot=rep(0,Nplot_F),
  intSeedWeight_plant=rep(0,Nplant)),lambdaSeedWeight=5,
  #Yield
  intYield=-0.3,slopeYield=1,sigmaYield=0.32,sigmaYield_field=c(0.2,0.04),sigmaYield_plot=c(0.7,0.24) 
)}

#Feed datalist into stan_rdump for use in CmdStan
with(datalist,stan_rdump(names(datalist),'tempDatSeed.data.R'))
#Feed inits into stan_rdump
temp <- inits()
with(temp,stan_rdump(names(temp),'initsSeed.data.R'))

#Full model
modPodcount_seed <- stan(file='visitation_pollen_model_seed.stan',data=datalist,
                         iter=1,chains=1,control=list(adapt_delta=0.8),init=inits)
beep(1)
# Setting max_treedepth=15 takes about 2-3x as long to run model. Use with care.

# save(modPodcount_seed,file='modPodcount_seed.Rdata') #Flw count, seed, yield model - 4hrs for 2000 iter
# save(modPodcount_seed,file='modPodcount_seed2.Rdata') #All other params (visits,plSize,etc)
# load('modPodcount_seed.Rdata')

# #ML version, no variance estimates
# modPodcount_seed2 <- optimizing(stan_model(file='visitation_pollen_model_seed.stan'),data=datalist,init=inits)
# str(modPodcount_seed2)
# modPodcount_seed2$par[names(modPodcount_seed2$par) %in% pars] #Get optim results

pars=c('intPlDens','slopeHbeeDistPlDens',#'slopeHbeeDistSqPlDens', #Planting density
       'sigmaPlDens','sigmaPlDens_field') 
pars=c('intPlSize','slopePlDensPlSize','slopeDistPlSize','sigmaPlSize') #Plant size
pars=c('intFlDens','slopePlSizeFlDens',
       'slope2016FlDens','slopeDistFlDens',
       'sigmaFlDens','sigmaFlDens_field') #Flower density
pars=c('intVisitLbee','slopeHbeeDistLbee','slopeLbeeDistLbee','slopeCentLbee','slopeMBayLbee', #Lbee vis
       'slopeStockingLbee','slope2016Lbee','slopeCentHbeeDistLbee','slopeStockingHbeeDistLbee',
       'slopeFlDensLbee','sigmaLbeeVisField','visitLbeePhi')
pars=c('intVisitHbee','slopeHbeeDistHbee','slopeLbeeDistHbee','slopeLbeeHbeeDistHbee',
       'slopeLbeeVisHbee', #Hbee vis 
       'slopeCentHbee','slopeFlDensHbee','slopeMBayHbee', 
       'visitHbeePhi','zeroVisHbeeTheta') 
pars=c('intPol','slopeHbeePol','slopeLbeePol','slopeCentPol','slopeHbeeDistPol','slopeFlDensPol', #Pollen
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
       'slopeSurvSeedCount','seedCountPhi','sigmaSeedCount_field','sigmaSeedCount_plot',
       'sigmaSeedCount_plant')
pars <- c('intSeedWeight','slopePolSeedWeight','slopeSeedCount', #Weight per seed
       'slopePlSizeSeedWeight',
       'slope2016SeedWeight','slopeLbeeDistSeedWeight','slopePlDensSeedWeight',
       'slopeStockingSeedWeight',
       'lambdaSeedWeight',
       'sigmaSeedWeight','sigmaSeedWeight_field', #Random effects for plot don't converge well
       'sigmaSeedWeight_plant')
pars <- c('intYield','slopeYield','sigmaYield','sigmaYield_field','sigmaYield_plot',
          'L_field','L_plot')
stan_hist(modPodcount_seed,pars=pars)+geom_vline(xintercept=0,linetype='dashed')
traceplot(modPodcount_seed,pars=pars,inc_warmup=F)
print(modPodcount_seed,pars=pars)

#Check model fit:
mod3 <- extract(modPodcount_seed)

fastPairs(mod3[pars])

#Coefficients
(mod3coefs <- data.frame(par='NA',parname=rownames(coefs(mod3[pars])),coefs(mod3[pars]),row.names=NULL))
print(xtable(mod3coefs,digits=c(0,0,0,3,3,3,3,3,3,0,4)),include.rownames=F)

#Distribution of random effects intercepts
t(apply(mod3$intPlDens_field,2,function(x) quantile(x,c(0.5,0.025,0.975)))) %>% 
  as.data.frame() %>% rename(med='50%',lwr='2.5%',upr='97.5%') %>% 
  arrange(med) %>% mutate(plot=1:n()) %>% 
  ggplot()+geom_pointrange(aes(x=plot,y=med,ymax=upr,ymin=lwr),alpha=0.5)+
  geom_hline(yintercept=0,col='red')

par(mfrow=c(2,1))
#planting density - good
with(mod3,PPplots(apply(plDens_resid,1,function(x) sum(abs(x))),
                  apply(predPlDens_resid,1,function(x) sum(abs(x))),
                  datalist$plDens_obs,apply(mod3$predPlDens,2,median)[datalist$obsPlDens_ind],
                  main='Plant density'))

#plant size - OK, but weird curvilinear shape
with(mod3,PPplots(apply(plSize_resid,1,function(x) sum(abs(x))),
               apply(predPlSize_resid,1,function(x) sum(abs(x))),
               datalist$plantSize,apply(mod3$predPlSize,2,median),
               main='Plant size')) 

#flower density - OK, but underpredicted at high densities
with(mod3,PPplots(apply(flDens_resid,1,function(x) sum(abs(x))),
                  apply(predFlDens_resid,1,function(x) sum(abs(x))),
                  apply(flDens,2,median),
                  apply(predFlDens,2,median),
                  main='Flower density')) 

#hbee visits - not the best, but OK
with(mod3,PPplots(apply(hbeeVis_resid,1,function(x) sum(abs(x))),
               apply(predHbeeVis_resid,1,function(x) sum(abs(x))),
               with(datalist,c(hbeeVis,hbeeVis_extra)),apply(mod3$predHbeeVis_all,2,median),
               main='Hbee visits')) 

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
with(mod3,PPplots(apply(pollen_resid,1,function(x) sum(abs(x))),
              apply(predPollen_resid,1,function(x) sum(abs(x))),
              datalist$pollenCount,apply(mod3$predPollenCount,2,median),
              main='Pollen counts')) 

#Flower count - good
with(mod3,PPplots(apply(flwCount_resid,1,function(x) sum(abs(x))),
               apply(predFlwCount_resid,1,function(x) sum(abs(x))),
               datalist$flwCount,apply(mod3$predFlwCount,2,median),
               main='Flower count per plant'))

#flower survival - good
with(mod3,PPplots(apply(podCount_resid,1,function(x) sum(abs(x))),
               apply(predPodCount_resid,1,function(x) sum(abs(x))),
               datalist$podCount,apply(mod3$predPodCount,2,median),main='Pods per plant'))

#seeds per pod - not good
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

mod3 <- extract(modPodcount_seed)

#Plot-level data (647 rows) - visitation data
plotDat <- with(datalist,data.frame(
  plot=1:(Nplot+Nplot_extra),
  hbeeDist=c(hbee_dist,hbee_dist_extra),lbeeDist=c(lbee_dist,lbee_dist_extra), #Centered on 4.641845
  lbeeVis=c(lbeeVis,lbeeVis_extra),isMBay=c(isMBay,isMBay_extra), #Centered on 3.113714
  flDens=apply(mod3$flDens,2,mean),isCent=c(isCent,isCent_extra), #Sqrt tranformed - 22
  halfStock=c(lbeeStocking,lbeeStocking_extra), #half-stocking
  year=c(is2016,is2016_extra),
  hbeeResid=log((c(hbeeVis,hbeeVis_extra)+1)/c(totalTime,totalTime_extra))-apply(mod3$visitMu_hbee,2,median), #Log-residuals
  lbeeResid=log((c(lbeeVis,lbeeVis_extra)+1)/c(totalTime,totalTime_extra))-apply(mod3$visitMu_lbee,2,median)
  )) %>% 
  mutate(logHbeeDist=log(hbeeDist)-mean(log(hbeeDist)),logLbeeDist=log(lbeeDist)-mean(log(lbeeDist)))

#Model matrix for hbee distance/stocking
MM_lbeeVis <- with(plotDat,data.frame(int=1,logLbeeDist,logHbeeDist,isMBay,isCent,halfStock,year,flDens,
                         centHbeeDist=isCent*logHbeeDist,stockHbeeDist=halfStock*logHbeeDist)) %>% as.matrix()

#Coefficent matrix for lbee visits
coef_lbeeVis <- with(mod3,data.frame(intVisitLbee,slopeLbeeDistLbee,slopeHbeeDistLbee,slopeFBayLbee,slopeCentLbee,
                    slopeStockingLbee,slope2016Lbee,slopeFlDensLbee,slopeCentHbeeDistLbee,
                    slopeStockingHbeeDistLbee)) %>% as.matrix()

#Partial effect of hbee distance/stocking
MM_temp <- MM_lbeeVis %>% as.data.frame() %>% mutate_at(vars(-logHbeeDist,-halfStock,-stockHbeeDist),mean) %>% as.matrix()

p1 <- data.frame(MM_temp,t(apply(MM_temp %*% t(coef_lbeeVis),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>%  
  mutate(resid=plotDat$hbeeResid+pred) %>% 
  mutate(halfStock=factor(halfStock,labels=c('Full','Half'))) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=exp(logHbeeDist+4.641)))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=halfStock),alpha=0.3,show.legend=F)+
  geom_line(aes(y=pred,col=halfStock),size=1)+
  geom_jitter(aes(y=resid,col=halfStock),alpha=0.5,size=1,width=5)+
  xlim(0,400)+ylim(0,50)+
  labs(x='Distance from edge',y='Visits/10 mins',col='Stocking')+
  scale_colour_manual(values=c('red','blue'))+scale_fill_manual(values=c('red','blue'))+
  theme(legend.position=c(0.8,0.85),legend.background=element_rect(fill='white',colour='black',linetype='solid'),
        legend.title=element_text(size=18))
ggsave('../Figures/slopeStockingHbeeDistLbee.png',p1)


#Plant-level data (582 rows)
plantDat <- with(datalist,data.frame(
  plant=1:Nplant,
  plDens=apply(mod3$plDens,2,median)[plotIndex_F[plantIndex]], #Centered around 3.573633
  plSize=plantSize, #Centered around 3.197572
  hbeeDist=(log(hbee_dist)-mean(log(c(hbee_dist,hbee_dist_extra))))[plotIndex_F[plantIndex[1:Nplant]]], #Centered around 4.641845
  lbeeDist=(log(lbee_dist)-mean(log(c(lbee_dist,lbee_dist_extra))))[plotIndex_F[plantIndex[1:Nplant]]], #Centered around 3.113714
  edgeCent=isCent[plotIndex_F[plantIndex[1:Nplant]]],
  pol=apply(mod3$pollenMu_plot,2,median)[plotIndex_F[plantIndex]],
  lbeeVis=log(1+hbeeVis/totalTime)[plotIndex_F[plantIndex[1:Nplant]]],
  hbeeVis=log(1+lbeeVis/totalTime)[plotIndex_F[plantIndex[1:Nplant]]],
  flDens=apply(mod3$flDens,2,median)[plotIndex_F[plantIndex[1:Nplant]]],
  flwCount=flwCount,podCount=podCount,
  flwCount_resid=apply(mod3$flwCountMu,2,median)-log(flwCount) #Log-residuals
))

#Model matrix for flower count 
MM_flwCount <- with(plantDat,data.frame(int=1,edgeCent,pol,lbeeVis,flDens,plSize)) %>% 
  as.matrix()
#Coefficent matrix for seedMass
coef_flwCount <- with(mod3,data.frame(intFlwCount,slopeCentFlwCount,slopePolFlwCount,slopeLbeeVisFlwCount,
                                      slopeFlDensFlwCount,slopePlSizeFlwCount)) %>% as.matrix()

#Partial effect of plant size
MM_temp <- MM_flwCount %>% as.data.frame() %>% mutate_at(vars(-plSize),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$flwCount_resid+pred) %>% arrange(flDens) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=exp(plSize+3.204388)))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  geom_point(aes(y=resid),alpha=0.5)+#ylim(250,650)+
  labs(x='Plant size',y='Flowers per plant')


#Partial effect of flower density
MM_temp <- MM_flwCount %>% as.data.frame() %>% mutate_at(vars(-flDens),mean) %>% as.matrix()

data.frame(MM_temp,t(apply(MM_temp %*% t(coef_flwCount),1,function(x) quantile(x,c(0.5,0.025,0.975))))) %>% 
  rename(pred=X50.,upr=X97.5.,lwr=X2.5.) %>% 
  mutate(resid=plantDat$flwCount_resid+pred) %>% arrange(flDens) %>% 
  mutate_at(vars(pred:resid),exp) %>% 
  ggplot(aes(x=(flDens+22)^2))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  geom_point(aes(y=resid),alpha=0.5)+ylim(250,650)+
  labs(x='Flower density',y='Flowers per plant')

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


#Pod-level data (2885 rows)
podDat <- with(datalist,data.frame(
  field=plotIndex[plantIndex[podIndex[1:Npod]]],plot=plantIndex[podIndex[1:Npod]],plant=podIndex[1:Npod],
  pollen=tapply(apply(mod3$pollenMu-matrix(rep(mod3$intPol,1050),ncol=1050),2,median),flowerIndex,mean)[
    plotIndex_F[plantIndex[podIndex]]],
  plSize=plantSize[podIndex], #Plant size - centered on 3.197572
  edge=isCent[plantIndex[podIndex]], 
  lbeeDist=(log(c(lbee_dist,lbee_dist_extra))-mean(log(c(lbee_dist,lbee_dist_extra))))[plantIndex[podIndex]],
  hbeeDist=(log(c(hbee_dist,hbee_dist_extra))-mean(log(c(hbee_dist,hbee_dist_extra))))[plantIndex[podIndex]],
  flDens=apply(mod3$flDens,2,mean)[plantIndex[podIndex]], 
  plDens=apply(mod3$plDens,2,mean)[plantIndex[podIndex]], #centered on 3.573633
  stocking=lbeeStocking[plantIndex[podIndex]],surv=logit(podCount/flwCount)[plantIndex[podIndex]],
  year2016=is2016[plantIndex[podIndex]],seedCount=seedCount,seedWeight=seedMass,
  seedMassResid=seedMass-apply(mod3$seedWeightMu,2,median) #Median resid for seedMass
))

temp <- podDat %>% mutate_at(vars(field,plot,plant),factor) %>% 
  lmer(seedWeight~pollen+year2016+lbeeDist+log(seedCount)+plDens*plSize+stocking+(1|field/plot/plant),data=.)

summary(temp)
plot(temp) #Residuals
plot(podDat$seedWeight,predict(temp),xlab='Actual',ylab='Predicted');abline(0,1,col='red')



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
  ggplot(aes(x=exp(plSize+3.197572)))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(y=pred),size=1)+
  geom_point(aes(y=resid),alpha=0.5)+
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
  ggplot(aes(x=exp(plDens+3.573633)))+
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
  ggplot()+geom_point(aes(x=exp(plSize+3.197572),y=exp(plDens+3.573633),size=pred))+
  labs(x='Plant size',y='Plant density',col='Weight per seed')

# Total effects for seed fields -------------------------------------------

# #Load all coefficients from models
setwd('./Models')
# load('modPodcount_seed.Rdata') #Seed count/weight model  + flw count + yield model
# mod1 <- extract(modPodcount_seed)
# rm(modPodcount_seed); gc()
# load('modPodcount_seed2.Rdata') #Everything else
# mod2 <- extract(modPodcount_seed)
# rm(modPodcount_seed)
# # load('modPodcount_seed3.Rdata') #All other params (visits,plSize,etc)
# # mod3 <- extract(modPodcount_seed)
# # rm(modPodcount_seed)
# mod3 <- c(mod1,mod2)
# mod3 <- mod3[unique(c(names(mod1),names(mod2)))] #Removes doubles
# mod3$lp__ <- NULL #Removes LP
# rm(mod1,mod2); gc();
# save(mod3,file='modPodcount_seed3.Rdata')
load('modPodcount_seed3.Rdata')
#Fixes Fbay-Mbay names
names(mod3)[grepl('FBay',names(mod3))] <- sub('FBay','MBay',names(mod3)[grepl('FBay',names(mod3))])

# transdens <- log(plDens)-3.573633 #Plant density
# transPlSize <- log(plSize)-3.197572 #Plant size
# transHbeeDist <- log(hbeeDist)-4.641845 #Hbee distance
# transLbeeDist <- log(lbeeDist)-3.113714 #Lbee distance
# transFlDens <- sqrt(flDens)-22 #Flw dens
# transLbeeVis <- log(1+lbeeVis)

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

#Generate flower counts from plant sizes
genFlwCount <- function(plSize,isCent,avgPol,flDens,
                        int,slopePlSize,slopeCent,slopePol,slopeFlDens,
                        intPhi,slopePlSizePhi,addVar=T){
  N <- getN(list(plSize,isCent,avgPol,flDens))
  
  #plSize=plant size, intercept,slope of plant size, phi
  transPlSize <- log(plSize)-3.197572 #Plant size
  transFlDens <- sqrt(flDens)-22 #Flw dens
  mu <- int+slopePlSize*transPlSize+slopeCent*isCent+
              slopePol*avgPol+slopeFlDens*transFlDens #Expected value
  if(!addVar) return(exp(mu)) #No variance
  phi <- exp(intPhi+slopePlSizePhi*transPlSize) #Dispersion
  a <- round(exp(rnorm(N,mu,phi)))
  if(any(is.na(a))) stop('NAs from rnorm')
  return(a)
}

#Generate flower survival (pod count) from flower counts
genPodCount <- function(flwCount,pol,plSize,isCent,hbeeDist,lbeeDist,flwDens,
                        int,slopePol,slopePlSize,slopeCent,slopeHbeeDist,slopeLbeeDist,slopeFlwDens,
                        intPhi,slopePlSizePhi,addVar=T){
  #noBeta removes beta-binomial step
  N <- getN(list(flwCount,pol,plSize,isCent,hbeeDist,lbeeDist,flwDens))
  transPlSize <- log(plSize)-3.197572 #Plant size
  transHbeeDist <- log(hbeeDist)-4.641845 #Hbee distance
  transLbeeDist <- log(lbeeDist)-3.113714 #Lbee distance
  transFlDens <- sqrt(flwDens)-22 #Flw dens
  #Expected value
  mu <- invLogit(int+slopePol*pol+slopePlSize*transPlSize+slopeCent*isCent+
                   slopeHbeeDist*transHbeeDist+slopeLbeeDist*transLbeeDist+
                   slopeFlwDens*transFlDens)
  if(!addVar) return(mu*flwCount) #Without dispersion
  phi <- exp(intPhi+slopePlSizePhi*transPlSize) #Dispersion
  theta <- rbeta(N,mu*phi,(1-mu)*phi) #Adds beta
  a <- rbinom(N,flwCount,theta)
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
  
  # #Test values
  # hbeeDist=seq(1,401,50)
  # lbeeDist=10
  # isCent=0
  # is2016=0
  # isHalfStock=0
  # dat=mod3
  # returnAll=F
  # useMean=F
  # plotVar=F
  # rm(hbeeDist,lbeeDist,isCent,is2016,returnAll,useMean,isHalfStock) #Cleanup
  # rm(simPlDens,simPlSize,simFlwDens,simLbeeVis,simHbeeVis,simAvgPol,simFlwCount,
  #    simPodCount,simSurv,simSeedCount,simSeedWeight,simYield)
  
  if(useMean) dat <- lapply(dat,mean)
  
  #Simulated plant density
  simPlDens <- round(with(dat,genPlDens(hbeeDist=hbeeDist,int=sample(intPlDens,1),
                      slopeDist=sample(slopeHbeeDistPlDens,1),
                      stDev=ifelse(plotVar,sample(sigmaPlDens,1),0),transform=T)))
  
  #Compare to actual - looks OK, but a bit too sparse at the high end. Mainly caused by random effect
  # par(mfrow=c(2,1))
  # hist(surveyAllSeed$PlDens,breaks=seq(0,150,5),xlab='Actual Plant Density',main=NULL)
  # hist(simPlDens,breaks=seq(0,150,5),xlab='Sim Plant Density',main=NULL)
  
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
  
  # #Looks good
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
  # hist(with(dat,genAvgPol(hbeeVis=simHbeeVis,lbeeVis=simLbeeVis,isCent=isCent,hbeeDist=hbeeDist,
  #           flDens=simFlwDens,int=sample(intPol,1),slopeHbeeVis=sample(slopeHbeePol,1),
  #           slopeLbeeVis=sample(slopeLbeePol,1),
  #           slopeCent=sample(slopeCentPol,1),slopeHbeeDist=sample(slopeHbeeDistPol,1),
  #           slopeFlDens=sample(slopeFlDensPol,1),
  #           stDev=sample(sigmaPolPlot,1),center=F,transform=T)),xlab='Sim Avg pollen/plot',
  #      main=NULL,breaks=seq(0,200,5))
  # hist(with(pollenAllSeed,tapply(Pollen,paste(Field,Year,Distance),mean)),
  #      xlab='Actual Avg Pollen/plot',main=NULL,breaks=seq(0,200,5))
  
  #Simulate flower count per plant
  simFlwCount <- with(dat,mapply(genFlwCount,plSize=simPlSize,isCent=isCent,
                     avgPol=simAvgPol,flDens=simFlwDens,int=sample(intFlwCount,1),
                     slopePlSize=sample(slopePlSizeFlwCount,1),slopeCent=sample(slopeCentFlwCount,1),
                     slopePol=sample(slopePolFlwCount,1),slopeFlDens=sample(slopeFlDensFlwCount,1),
                     intPhi=sample(intPhiFlwCount,1),slopePlSizePhi=sample(slopePlSizePhiFlwCount,1),
                     addVar=plotVar,SIMPLIFY=F))
  
  # #Compare to actual - looks good
  # par(mfrow=c(2,1))
  # upr <- max(unlist(simFlwCount),with(plantsAllSeed,Pods+Missing),na.rm=T)
  # hist(unlist(simFlwCount),xlab='Sim Flowers/plant',main=NULL,breaks=seq(0,upr+50,50))
  # hist(with(plantsAllSeed,Pods+Missing),xlab='Actual Flowers/plant',main=NULL,breaks=seq(0,upr+50,50))
  
  #Simulate flower survival to pod
  simPodCount <- with(dat,mapply(genPodCount,flwCount=simFlwCount,pol=simAvgPol,plSize=simPlSize,isCent=isCent,
                           hbeeDist=hbeeDist,lbeeDist=lbeeDist,flwDens=simFlwDens,int=sample(intFlwSurv,1),
                           slopePol=sample(slopePolSurv,1),slopePlSize=sample(slopePlSizeSurv,1),
                           slopeCent=sample(slopeEdgeCentSurv,1),slopeHbeeDist=sample(slopeHbeeDistSurv,1),
                           slopeLbeeDist=sample(slopeLbeeDistSurv,1),slopeFlwDens=sample(slopeFlwDensSurv,1),
                           intPhi=sample(intPhiFlwSurv,1),slopePlSizePhi=sample(slopePlSizePhiFlwSurv,1),
                           addVar=plotVar,SIMPLIFY=F))
  
  # #Flower survivorship (proportion) - looks OK
  # par(mfrow=c(2,1))
  # hist(unlist(simPodCount)/unlist(simFlwCount),xlab='Sim flower survival',main=NULL,breaks=seq(0,1,0.05))
  # hist(with(plantsAllSeed,Pods/(Pods+Missing)),xlab='Actual flower survival',main=NULL,breaks=seq(0,1,0.05))
  
  #Proportion survival
  simSurv <- mapply(function(x,y) x/y,x=simPodCount,y=simFlwCount,SIMPLIFY=F) 
  
  #Simulate avg seeds per pod
  simSeedCount <- with(dat,mapply(genSeedCount,pol=simAvgPol,plSize=simPlSize,isCent=isCent,hbeeDist=hbeeDist,
                    flDens=simFlwDens,flwSurv=simSurv,int=sample(intSeedCount,1),
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
  # abline(v=mean(unlist(simSeedWeight)),col='red')
  # abline(v=median(unlist(simSeedWeight)),col='red',lty='dashed')
  # hist(with(plantsAllComm,AvPodMass*1000/AvPodCount),xlab='Actual weight per seed (mg)',main=NULL,breaks=seq(min,upr,0.2))
  # abline(v=mean(with(plantsAllComm,AvPodMass*1000/AvPodCount),na.rm=T),col='red')
  # abline(v=median(with(plantsAllComm,AvPodMass*1000/AvPodCount),na.rm=T),col='red',lty='dashed')
  
  #Simulate total yield for each plot
  simYield <- with(dat,mapply(genYield,pods=simPodCount,seedCount=simSeedCount,seedWeight=simSeedWeight,
                              int=sample(intYield,1),slopeYield=sample(slopeYield,1)))
  
  # #Looks OK. Tails aren't as long, but this is probably due to tails from sim seed count
  # par(mfrow=c(2,1))
  # hist(simYield,xlab='Sim yield/m2',main=NULL,breaks=seq(0,1000,20))
  # hist(with(datalist,tapply(yield*exp(plDens[plantIndex]+3.573633),plantIndex,mean)),
  #      xlab='Actual yield/m2',main=NULL,breaks=seq(0,1000,20))
  
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

#Simulate pollination effects at generic fields
scenario <- expand.grid(hDist=c(20,400),lDist=seq(1,51,5),cent=c(0,0.5,1),halfStock=c(0,1))
results <- replicate(500,with(scenario,simSeed(hbeeDist=hDist,lbeeDist=lDist,isCent=cent,
                                              is2016=0,isHalfStock=halfStock,dat=mod3,plotVar=F)))
beep(2)
#Simulated yield distribution
results2 <- data.frame(scenario,t(apply(results,1,function(x) quantile(x,c(0.5,0.05,0.95),na.rm=T)))) %>%
  rename('pred'='X50.','lwr'='X5.','upr'='X95.')

#Results in bu/acre
results2 %>%
  mutate_at(vars(pred,upr,lwr),g2bushels) %>%
  mutate(cent=factor(cent,labels=c('Edge','Halfway','Center')),hDist=factor(hDist,labels=c('Near','Far'))) %>%
  mutate(halfStock=factor(halfStock,labels=c('Regular','Half-stocking'))) %>% 
  ggplot(aes(x=lDist,y=pred,col=cent))+
  # geom_ribbon(aes(ymax=upr,ymin=lwr,col=NULL,fill=cent),alpha=0.2,show.legend=F)+
  geom_line(aes(linetype=halfStock),size=1)+
  facet_grid(~hDist)+
  labs(y='Predicted yield (bu/acre)',x='Shelter Distance',col='Bay\nposition')+
  scale_colour_manual(values=c('darkorange','red','blue'))+scale_fill_manual(values=c('darkorange','red','blue'))






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
