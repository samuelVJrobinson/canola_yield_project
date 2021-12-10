#HELPER FUNCTIONS FOR PATH ANALYSIS

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
#Function to calculate R2 for GLMMs. fixVar, ranVar and resVar are either variable names (vector of var names for random effects), or NA, indicating no effect specified
r2calc <- function(mod,fixVar,ranVar,resVar){
  fixed <- ifelse(!is.na(fixVar),sd(apply(mod[[fixVar]],2,median)),0) #fixed effect sd
  if(!is.na(ranVar[1])){
    ran <- 0
    for(i in 1:length(ranVar)){ #Adds variances from all levels of random effects
      ran <- ran + median(mod[[ranVar[i]]])
    }
  } else ran <- 0
  res <- ifelse(!is.na(resVar),median(mod[[resVar]]),0) #residual variance
  print(paste('Fixed effect var =',round(fixed,4)))
  print(paste('Random effect var =',round(ran,4)))
  print(paste('Residual var =',round(res,4)))
  print(paste('Conditional R2 =',round((fixed+ran)/(fixed+ran+res),4)))
  print(paste('Marginal R2 =',round((fixed)/(fixed+ran+res),4)))
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
  resRange <- range(resid)
  pResRange <- range(predResid)
  xl <- yl <- c(pmin(resRange[1],pResRange[1]),pmax(resRange[2],pResRange[2]))
  x <- sum(resid<predResid)/length(resid)
  plot(resid~predResid,ylab='Sum actual residuals',xlab='Sum simulated residuals',xlim=xl,ylim=yl,
       main=paste(main,' (p =',round(min(x,1-x),3),')'))
  # legend('topleft',)
  abline(0,1,col='red') #PP plot
  plot(actual~pred, #Predicted vs Actual
       xlab=paste('Predicted',main),ylab=paste('Actual',main)) 
  abline(0,1,col='red')
  abline(lm(actual~pred),col='red',lty=2)
  par(mfrow=c(1,1))
}

#Plot of random intercepts
compareRE <- function(mod,parSet){
  require(ggpubr)
  p1 <- extract(mod,pars=parSet)[[1]] %>% 
    apply(.,2,function(x) quantile(x,c(0.1,0.5,0.9))) %>% 
    t() %>% data.frame() %>% setNames(c('lwr','med','upr')) %>% 
    ggplot(aes(sample=med))+ ##q-q plot
    geom_qq()+
    geom_qq_line()  
  p2 <- extract(mod,pars=parSet)[[1]] %>% 
    apply(.,2,function(x) quantile(x,c(0.1,0.5,0.9))) %>% 
    t() %>% data.frame() %>% setNames(c('lwr','med','upr')) %>% 
    tibble::rownames_to_column() %>%
    mutate(rowname=as.numeric(rowname)) %>% 
    ggplot(aes(x=rowname,y=med))+
    geom_pointrange(aes(ymax=upr,ymin=lwr))+
    geom_hline(yintercept = 0,col='red',linetype='dashed')+
    labs(x='Intercept',y='Posterior Dist.')
  ggarrange(p1,p2,ncol=1,nrow=2)
}

#Convert g/m2 to bushels/acre
g2bushels <- function(x){ 
  x*4046.86/22679.6 
  # 453.592*50/bushel - using 50 lbs/bushel estimate
  # 44.0920 #bushels per tonne
  # 2204.62 #lbs per tonne
}