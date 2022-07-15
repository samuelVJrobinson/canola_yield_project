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
parTable <- function(x){
  require(dplyr); require(tidyr)
  if(is.null(x)) return(NA)
  n <- names(x) #Model parameters
  n <- n[(!grepl('(\\[[0-9]+,*[0-9]*\\]|lp)',n))|grepl('([sS]igma|slope)',n)] #Drops parameter vectors, unless name contains "sigma" (variance term) or "slope"
  n <- n[!grepl('_miss',n)] #Gets rid of imputed values
  
  parTable <- summary(x,pars=n)$summary %>% 
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
fastPairs <- function(l,pars=NULL,N=500){ 
  library(rstan)
  if(class(l)=='stanfit'){
    l <- rstan::extract(l,pars=pars)
    if(any(sapply(l,function(x) length(dim(x)))>1)){
      stop('Parameters have more than 1 dimension: ',paste(pars[sapply(l,function(x) length(dim(x))>1)],collapse=', '))
    }
    l <- do.call('cbind',l)
    if(nrow(l)>N) l <- l[round(seq(1,nrow(l),length.out=N)),,drop=FALSE]
  }
  if(ncol(l)==1){
    warning('Only 1 parameter found. Displaying histogram instead of biplot')
    hist(l,main='',xlab=pars)
  } else {
    if(length(ncol(l))>10) warning('More than 10 parameters. Pairplot may take some time.')
    pairs(l,lower.panel=function(x,y){
      par(usr=c(0,1,0,1))
      text(0.5, 0.5, round(cor(x,y),2), cex = 1 * exp(abs(cor(x,y))))},
      panel=function(x,y){ points(x,y,pch=19)}
    )
  }
  
  
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
PPplots <- function(mod,actual=NULL,pars=c('pred','resid','predResid'),main='',index=NA,jitterX=NA,scale='',ZIpar=NA){
  if(is.null(mod)) stop('Model not found')
  require(ggpubr)
  oldtheme <- theme_get() #Get theme
  theme_set(theme_classic())

  modVals <- extract(mod)
  nam <- names(modVals)
  
  if(any(!pars %in% nam)){
    stop(paste0('Parameters not found: ',paste(pars[!pars %in% nam],collapse=', ')))
  }
  
  if(main=='') main <- deparse(substitute(actual))
  
  genVals <- modVals[[which(nam==pars[1])]] #Generated values
  pred <- apply(genVals,2,median) #Median predicted value
  predUpr <- apply(genVals,2,function(x) quantile(x,0.9)) #Upper predicted
  predLwr <- apply(genVals,2,function(x) quantile(x,0.1)) #Lower predicted
  
  resid <- apply(modVals[[which(nam==pars[2])]],1,function(x) sum(abs(x))) #Sum of actual residuals
  predResid <- apply(modVals[[which(nam==pars[3])]],1,function(x) sum(abs(x))) #Sum of simulated residuals
  
  if(length(index)>1 || !is.na(index)){
    pred <- pred[index] #If some predictions are imputed, show only observed
    predUpr <- predUpr[index]; predLwr <- predLwr[index]
  } 
  
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
  
  #Actual vs Generated plot
  
  d2 <- data.frame(actual,pred,predUpr,predLwr) #Data for generated vs actual plot
  
  if(scale=='log'){
    d2 <- d2 %>% mutate(across(everything(),~.x+0.1)) #If log-log, add 0.1 to generated and actual
  }
  
  p2 <- d2 %>% ggplot(aes(x=pred,y=actual))
  
  if(is.na(jitterX)){
    p2 <- p2 + geom_point() + geom_errorbarh(aes(xmax=predUpr,xmin=predLwr),height=0,alpha=0.1)
  } else {
    p2 <- p2 + geom_point(position=position_jitter(width = jitterX,height=0)) +
      geom_errorbarh(aes(xmax=predUpr,xmin=predLwr),height=0,alpha=0.1,
                     position=position_jitter(width = jitterX,height=0))
  }
  
  p2 <- p2 +
    geom_smooth(method='lm',se=FALSE,col='blue',linetype='dashed',formula = y~x)+ #Regression line
    geom_abline(intercept = 0,slope=1,col='blue',linetype='solid')+ #1:1 line
    labs(x=paste('Generated',main),y=paste('Actual',main)) #Axis labels
  
  if(scale=='log'){
    p2 <- p2 + scale_x_log10()+scale_y_log10()
  } else if(scale=='sqrt') {
    p2 <- p2 + scale_x_sqrt()+scale_y_sqrt()
  }
  
  # Density plots of predicted datasets
  
  xMax <- max(actual) #Upper and lower lims
  xMin <- min(actual)
  
  yDens <- density(actual,from=xMin,to=xMax) #Kernel density of actual data
  predDens <- t(apply(modVals[[which(nam==pars[1])]],1,function(i) density(i,from=xMin,to=xMax)$y))
  predDens <- t(apply(predDens,2,function(i) quantile(i,c(0.05,0.25,0.75,0.95))))
  predDens <- ifelse(predDens<.Machine$double.eps,0,predDens)
  colnames(predDens) <- c('lwr2','lwr1','upr1','upr2')
  
  p3 <- data.frame(x=yDens$x,y=yDens$y,predDens) %>% 
    ggplot(aes(x=x))+
    geom_ribbon(aes(ymin=lwr2,ymax=upr2),alpha=0.2)+
    geom_ribbon(aes(ymin=lwr1,ymax=upr1),alpha=0.2)+
    geom_line(aes(y=y))+
    labs(x='Value',y='Density',title='Actual vs Simulated Density')
  
  if(is.na(ZIpar)){
    
    p <- ggarrange(p1,p2,p3,ncol=1,nrow=3)  
  } else {
    #Zero-inflation plots
    propActualZero <- sum(actual==0)/length(actual)
    propSimZero <- apply(modVals[[which(nam==pars[1])]],1,function(x) sum(x==0)/length(x))
    p4 <- data.frame(type=rep(c('Simulated Zeros','ZI parameter'),each=length(modVals[[which(nam==ZIpar)]])),
               val=c(propSimZero,modVals[[which(nam==ZIpar)]])) %>% 
      ggplot()+
      geom_density(aes(x=val,fill=type),alpha=0.5)+
      geom_vline(xintercept = propActualZero,linetype='dashed')+
      labs(x='Proportion',y='Density',fill=NULL)+
      scale_fill_manual(values=c('black','red'))+
      theme(legend.position = c(0.25,0.9))
    p <- ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)  
  }
  
  theme_set(oldtheme)
  return(p)
  
}

#Plot of random intercepts
compareRE <- function(mod,parSet,colNum=1,alpha=1){
  if(is.null(mod)) stop('Model not found')
  pars <- extract(mod,pars=parSet)[[1]]
  
  if(length(dim(pars))>2){
    pars <- pars[,colNum,]
    print(paste0('Display column: ',colNum))
  } 
  
  p <- pars %>% 
    apply(.,2,function(x) quantile(x,c(0.1,0.5,0.9))) %>% 
    t() %>% data.frame() %>% setNames(c('lwr','med','upr'))
  
  require(ggpubr)
  p1 <- p %>% 
    ggplot(aes(sample=med))+ ##q-q plot
    geom_qq(alpha=alpha)+
    geom_qq_line(alpha=alpha)  
  p2 <- p %>% 
    tibble::rownames_to_column() %>%
    mutate(rowname=as.numeric(rowname)) %>% 
    ggplot(aes(x=rowname,y=med))+
    geom_pointrange(aes(ymax=upr,ymin=lwr),alpha=alpha)+
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

#Extract fits for marginal plots
# mod = stanfit model or dataframe with parameters 
# parList = named list with vectors of parameters to predict at (passed through expand.grid)
# otherPars = character vector of other parameters to add to the dataframe
# q = quantiles to predict at
# nrep = number of replicates of parameters to generate (for data.frame parameters)

getPreds <- function(mod,parList=NULL,offset=NULL,ZIpar=NULL,otherPars=NULL,q=c(0.5,0.025,0.975),nrep=4000,trans=NULL){
  if(is.null(parList)) stop('Specify parameters')
  dfNames <- c("param","mean","sd","Z","median",
               "lwr","upr","pval","n_eff","Rhat")
  
  if(class(mod)=='stanfit'){
    m <- extract(mod,par=names(parList))
    m <- do.call('cbind',m) #All draws
    m_mean <- apply(m,1,mean) #Mean values
    if(!is.null(otherPars)){
      extras <- extract(mod,par=otherPars) %>% sapply(.,median)
    }
  } else if(class(mod)=='data.frame' & !any(!names(mod)==dfNames)){ #If using a dataframe
    if(any(!names(parList) %in% mod$param)){
      stop(paste0('Parameters not found: ',paste(names(parList)[!names(parList) %in% mod$param],collapse=', ')))
    } 
    chooseRows <- mod$param %in% names(parList) 
    m <- mod[chooseRows,] 
    m_mean <- setNames(m$mean,names(parList)) #Mean value
    m <- do.call('cbind',lapply(1:nrow(m),function(i) runif(nrep,m$lwr[i],m$upr[i]))) #Generate draws b/w upr and lower
    colnames(m) <- names(parList)
    
    if(!is.null(otherPars)){
      extras <- mod$median[mod$param %in% otherPars] 
      names(extras) <- otherPars
    }
  } else {
    warning('mod is not stanfit object or dataframe with correct column names')
    return(NA)
  }
  
  mm <- expand.grid(parList) #Model matrix
  mmMat <- as.matrix(mm)
  
  if(!is.null(offset)){ #Add offset term if needed
    m <- cbind(m,offset=rep(1,nrow(m))) #Slope fixed at 1
    m_mean <- c(m_mean,setNames(1,'offset'))
    mmMat <- cbind(mmMat,offset=rep(log(offset),nrow(mmMat)))
  } 
  
  if(!is.null(ZIpar)){ #Add ZI term if needed (on log scale)
    i <- mod$param==ZIpar
    # m <- cbind(m,ZI=log(1-rnorm(nrep,mod$mean[i],mod$sd[i])))
    # m <- cbind(m,ZI=log(1-c(mod$lwr[i],mod$med[i],mod$upr[i])))
    m <- cbind(m,ZI=log(1-runif(nrep,mod$lwr[i],mod$upr[i]))) #Add log(1-ZI) param
    m_mean <- c(m_mean,setNames(log(1-mod$mean[i]),'ZI'))
    mmMat <- cbind(mmMat,ZI=rep(1,nrow(mmMat)))
  }
  
  pred <- mmMat %*% t(m) #Multiply parameters by model matrix
  pred <- data.frame(mmMat %*% m_mean,t(apply(pred,1,function(x) quantile(x,q))))
  names(pred) <- c('mean','med','lwr','upr')
  if(!is.null(trans)){ #If transformation needed
    pred <- eval(call(trans,pred))
  }
  
  if(!is.null(otherPars)){ #If other parameters needed
    pred <- cbind(pred,sapply(extras,function(x) rep(x,nrow(pred))))
  }
  return(cbind(mm,pred))
}

# d-separation claims list from Shipley 2009; shows original terms in brackets with claim & conditioning set outside
# g = dag object, or list of formulae to be passed to ggdag
# form = Convert to R-style formulas from dagitty structure?
# 

shipley.test <- function(g,form=FALSE){ 
  require(dagitty)
  if(class(g)=='list'){
    require(ggdag); g <- do.call(dagify,g)
  }
  
  claim2Form <- function(x){ #Converts claim to formula, if needed
    # dep <- x$X
    # ind <- paste(x$Y,paste(x$Z,collapse=' + '),sep = ' + ')
    # paste0(dep,' ~ ',ind)
    dep <- x$X #Dependent variable
    cl <- x$Y #Claim
    c1 <- x$Z[!x$Z %in% parents(g,x$X)] #Parent(s) of claim
    c1 <- paste0(c(cl,c1),collapse=' + ') #Claim + parent of claim
    c2 <- x$Z[x$Z %in% parents(g,x$X)] #Original parents of X
    c2 <- paste0('(',paste0(c2,collapse=' + '),')') #Add brackets
    ind <- paste0(c(c1,c2),collapse=' + ') #Add claim + conditioning set
    paste0(dep,' ~ ',ind) #Return string
  }
  
  # From Shipley 2009
  #1. Express causal relationship as DAG
  #2. List each of the k pairs of variables that do not have an arrow b/w them
  claims <- impliedConditionalIndependencies(g) #All variables without arrows between them
  for(i in 1:length(claims)) claims[[i]]$Z <- list() #Remove conditioning sets
  claims <- claims[!duplicated(claims)] #Remove duplicate claims
  #Get exogenous variables
  exoVars <- names(g)[sapply(names(g),function(x,dag) length(parents(dag,x))==0,dag=g)] #Exogenous variables only
  isExo <- logical(length(claims)) #Is claim b/w only exogenous variables?
  for(i in 1:length(claims)) isExo[i] <- claims[[i]]$X %in% exoVars & claims[[i]]$Y %in% exoVars 
  claims <- claims[!isExo] #Remove claims b/w only exogenous variables
  
  for(i in 1:length(claims)) { #This isn't part of Shipley's test, but makes modeling more "realistic", and doesn't change results (independence claims _should_ be the same in either direction)
    if(length(ancestors(g,claims[[i]]$X))<length(ancestors(g,claims[[i]]$Y))){ #If x-value has fewer causal ancestors
      X <- claims[[i]]$Y # Switch position of x and y
      claims[[i]]$Y <- claims[[i]]$X
      claims[[i]]$X <- X
    }
  }
  if(exists('X')) rm(X)
  
  #3. For each of the k pairs of variables (Xi, Xj), list the set of other
  #variables, {Z} in the graph that are direct causes of either Xi or Xj. 
  for(i in 1:length(claims)){
    var1 <- claims[[i]]$X
    var2 <- claims[[i]]$Y
    parents1 <- parents(g,var1)
    parents2 <- parents(g,var2)
    claims[[i]]$Z <- unique(c(parents1,parents2))
  }
  
  # Sort by causal rank (more ancestors = lower on the list)
  
  #Names of nodes in order of ancestry length
  nodeOrder <- names(sort(sapply(names(g),function(v) length(ancestors(g,v))-1))) 
  
  depOrder <- match(sapply(claims,function(x) x$X),nodeOrder) #Order of dependent vars
  claimOrder <- match(sapply(claims,function(x) x$Y),nodeOrder) #Order of claims
  claims <- claims[order(depOrder,claimOrder)] #Sort by order of dep vars, then order of claims
  
  for(i in 1:length(claims)) claims[[i]]$Z <- sort(claims[[i]]$Z) #Sort conditioning set
  claims <- unname(claims) #Get rid of names
  if(form) claims <- lapply(claims,claim2Form) #Convert to formula if needed
  return(claims)
}

#Shipley's d-sep test, using a dataframe d with a column of p-values p and a column of labels lab

shipley.dSep <- function(d,p,lab){
  require(ggplot2); require(dplyr)

    pvals <- d %>% pull({{p}}) #p-values
  
  kval <- 2*length(pvals) #k-val (df)
  if(any(pvals==0)){
    warning('Zero p-value: C-stat will be infinity')
    ret <- c(Inf,kval,0)
  } else {
    cstat <- -2*sum(log(pvals)) #C-statistic
    cp <- 1-pchisq(cstat,kval) 
    ret <- c(cstat,kval,cp)
  }
  
  tCol <- rev(ifelse(pvals<0.05,'red','black'))
  ttl <- paste0('C-stat: ',round(ret[1],2),', k: ',ret[2],', p: ',formatC(ret[3],digits=3))
  
  #(Approximately) how many terms need to be added to make p-value<0.05?
  removeP <- function(p){ #Which p-values should be removed?
    cstatPval <- function(x) 1-pchisq(-2*sum(log(x)),2*length(x)) 
    remove <- rep(FALSE,length(p))
    while(cstatPval(p[!remove])<0.05) remove[which.min(ifelse(remove,Inf,p))] <- TRUE 
    return(remove)  
  }
  
  d %>% mutate(param=factor({{lab}},levels=rev({{lab}}))) %>% 
    mutate(tooLow={{p}}<=0.05) %>% 
    mutate(rem=removeP({{p}})) %>%
    ggplot(aes(x={{p}},y={{lab}},col=tooLow))+
    geom_point(aes(size=rem),show.legend = FALSE)+
    geom_point(show.legend = FALSE)+
    labs(x='p-value',y=NULL,title=ttl) +
    scale_x_log10(breaks=c(0.001,0.01,0.05,0.5))+
    scale_colour_manual(values=c('black','red'))+
    theme(axis.text.y = element_text(colour=tCol,size=10),
          axis.text.x = element_text(size=10),
          title = element_text(size=15))
}

capFirst <- function(x,decap=FALSE){ #Capitalize/decapitalize first letter of a string
  if(!is.character(x)) stop('Not a character or character vector')
  if(decap){
    substr(x, 1, 1) <- tolower(substr(x, 1, 1))
  } else {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  }
  return(x)
}

