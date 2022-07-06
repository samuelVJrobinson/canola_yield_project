#Test estimation of ZINB models
N <- 300
x <- runif(N,0,5)
Beta <- c(0.5,0.2)
X <- model.matrix(~x)
mu <- exp(X%*%Beta)
phi <- 0.7 #NB phi
theta <- 0.3 #Bernoulli theta 
ziProc <- rbinom(N,1,theta) #(extra zeros)
nbProc <- rnbinom(N,mu=mu,size=phi)
y <- ifelse(ziProc==1,0,nbProc)
# cbind(ziProc,nbProc,y)

datalist <- list(N=N,Nbeta=length(Beta),X=X,y=y)

testMod <- stan(file='./zinbTest.stan',data=datalist,iter=2000,chains=4,init=0) #Run model

p <- c('Beta','phi','theta')
p2 <- c('Beta[1]','Beta[2]','phi','theta')
stan_trace(testMod,pars=p)
stan_dens(testMod,pars=p)
fastPairs(testMod,p2)
PPplots(testMod,datalist$y,c('generated','resid','generated_resid'),'y',ZIpar='theta') #The actual vs generated plot sucks. Can it be improved?

#Something like this:
t(apply(exp(X %*% t(extract(testMod,'Beta')[[1]]))*(1-outer(rep(1,datalist$N),extract(testMod,'theta')[[1]])),
        1,function(x) quantile(x,c(0.5,0.1,0.9)))) %>% 
  data.frame() %>% setNames(c('med','lwr','upr')) %>% 
  mutate(x=x,y=datalist$y) %>% 
  ggplot(aes(x=med,y=y))+
  geom_point() + 
  geom_errorbarh(aes(xmax=upr,xmin=lwr),height=0,alpha=0.1)+
  geom_smooth(method='lm',se=FALSE,col='blue',linetype='dashed',formula = y~x)+ #Regression line
  geom_abline(intercept = 0,slope=1,col='blue',linetype='solid')+ #1:1 line
  labs(x='Predicted',y='Actual') #Axis labels

#Actual plot looks OK
t(apply(exp(X %*% t(extract(testMod,'Beta')[[1]]))*(1-outer(rep(1,datalist$N),extract(testMod,'theta')[[1]])),
        1,function(x) quantile(x,c(0.5,0.1,0.9)))) %>% 
  data.frame() %>% setNames(c('med','lwr','upr')) %>% 
  mutate(x=x,y=datalist$y) %>% 
  arrange(x) %>% 
  ggplot(aes(x=x,y=y))+
  geom_line(aes(x=x,y=med))+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_point()
