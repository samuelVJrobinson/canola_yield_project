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
#Get table of parameters to go into kable
parTable <- function(mod,pars){
  require(dplyr); require(tidyr)
  summary(mod,pars=pars)$summary %>% 
    data.frame() %>% 
    transmute(mean,sd,Z=mean/sd,
              median=X50.,lwr=X2.5.,upr=X97.5.,
              pval=pnorm(0,abs(mean),sd)*2,
              n_eff,Rhat) %>% 
    mutate(across(c(mean:upr,Rhat),~round(.x,3))) %>%
    mutate(pval=round(pval,4),n_eff=round(n_eff)) %>% 
    tibble::rownames_to_column('param')
}
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

# #Posterior predictive check plots - older version
# PPplots <- function(resid,predResid,actual,pred,main=NULL){
#   par(mfrow=c(2,1))
#   resRange <- range(resid)
#   pResRange <- range(predResid)
#   xl <- yl <- c(pmin(resRange[1],pResRange[1]),pmax(resRange[2],pResRange[2]))
#   x <- sum(resid<predResid)/length(resid)
#   plot(resid~predResid,ylab='Sum actual residuals',xlab='Sum simulated residuals',xlim=xl,ylim=yl,
#        main=paste(main,' (p =',round(min(x,1-x),3),')'))
#   # legend('topleft',)
#   abline(0,1,col='red') #PP plot
#   plot(actual~pred, #Predicted vs Actual
#        xlab=paste('Predicted',main),ylab=paste('Actual',main)) 
#   abline(0,1,col='red')
#   abline(lm(actual~pred),col='red',lty=2)
#   par(mfrow=c(1,1))
# }

#Posterior predictive check plots
PPplots <- function(mod,actual=NULL,pars=c('pred','resid','predResid'),main='',index=NA,jitterX=NA){
  require(ggpubr)
  oldtheme <- theme_get() #Get theme
  theme_set(theme_classic())

  modVals <- extract(mod)
  nam <- names(modVals)
  
  if(main=='') main <- deparse(substitute(actual))
  
  pred <- apply(modVals[[which(nam==pars[1])]],2,median) #Median predicted value
  resid <- apply(modVals[[which(nam==pars[2])]],1,function(x) sum(abs(x))) #Sum of actual residuals
  predResid <- apply(modVals[[which(nam==pars[3])]],1,function(x) sum(abs(x))) #Sum of simulated residuals
  
  if(length(index)>1 || !is.na(index)) pred <- pred[index] #If some predictions are imputed, show only observed
  
  d1 <- data.frame(resid,predResid) #Data frame for p1
  
  resRange <- range(resid) #Range limits
  pResRange <- range(predResid)
  xl <- yl <- c(pmin(resRange[1],pResRange[1]),pmax(resRange[2],pResRange[2]))
  x <- sum(resid<predResid)/length(resid) #Posterior p-val
  
  p1 <- d1 %>% ggplot(aes(x=predResid,y=resid)) + 
    geom_point() +
    geom_abline(intercept = 0,slope=1,col='blue',linetype='solid')+ #1:1 line
    labs(y='Sum actual residuals',x='Sum simulated residuals', 
         title=paste0(main,' (p = ',round(min(x,1-x),3),')'))+
    xlim(xl)+ylim(yl)
  
  d2 <- data.frame(actual,pred) #Data for predicted vs actual plot
  
  p2 <- d2 %>% ggplot(aes(x=pred,y=actual))
  
  if(is.na(jitterX)){
    p2 <- p2 + geom_point()
  } else {
    p2 <- p2 + geom_point(position=position_jitter(width = jitterX,height=0))
  }
  
  p2 <- p2 +
    geom_smooth(method='lm',se=FALSE,col='blue',linetype='dashed',formula = y~x)+ #Regression line
    geom_abline(intercept = 0,slope=1,col='blue',linetype='solid')+ #1:1 line
    labs(x=paste('Predicted',main),y=paste('Actual',main)) #Axis labels
  p <- ggarrange(p1,p2,ncol=1,nrow=2)
  
  theme_set(oldtheme)
  
  return(p)
  
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