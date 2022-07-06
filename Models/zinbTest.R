#Test estimation of ZINB models
N <- 300
x <- runif(N,0,5)
Beta <- c(0.5,0.1)
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
PPplots(testMod,datalist$y,c('generated','resid','generated_resid'),'y') 

plot(apply(extract(testMod,'generated')[[1]],2,max),
      datalist$y,xlab='Generated',ylab='Actual')
