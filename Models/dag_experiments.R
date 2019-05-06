#R script to look at how basis sets are constructed for conditional independence claims
library(dagitty)
# Tests -------------------------------------------------------------------
#Model 1 from Shipley 2009
g1 <- dagitty("dag{x1 -> x2 -> x3 -> x4 -> x5}")
#Stated claims list from Shipley 2009
# x1 _||_ x3 | x2 
# x1 _||_ x4 | x3 
# x1 _||_ x5 | x4 
# x2 _||_ x4 | x1, x3 #Includes parents of x2 as well as x4
# x2 _||_ x5 | x1, x4
# x3 _||_ x5 | x2, x4 

print(impliedConditionalIndependencies(g1,type='missing.edge')) #This doesn't look right
# "missing.edge"
# x1 _||_ x3 | x2 Same as paper
# x1 _||_ x4 | x3 Same as paper
# x1 _||_ x4 | x2 Why is this here? x2 not a parent of either term
# x1 _||_ x5 | x4 Same as paper
# x1 _||_ x5 | x3 ??
# x1 _||_ x5 | x2 ??
# x2 _||_ x4 | x3 Same as paper, but missing x1 (parent of x2)
# x2 _||_ x5 | x4 Same as paper, but missing x1
# x2 _||_ x5 | x3 ??
# x3 _||_ x5 | x4 Same as paper, but missing x2

print(impliedConditionalIndependencies(g1,type='basis.set')) 
#This looks better, but only uses parents of the main variable in the conditioning set
# "basis.set"
# x3 _||_ x1 | x2 
# x4 _||_ x1, x2 | x3
# x5 _||_ x1, x2, x3 | x4
# "basis.set" revised
# x3 _||_ x1 | x2 Same as paper 
# x4 _||_ x1 | x3 Same as paper 
# x5 _||_ x1 | x4 Same as paper 
# x4 _||_ x2 | x3 Same as paper, but missing x1
# x5 _||_ x2 | x4 Same as paper, but missing x1
# x5 _||_ x3 | x4 Same as paper, but missing x2

# Trying the same thing in pSEM -------------------------------------------
library(piecewiseSEM)
library(nlme)
library(lme4)
data(shipley)

shipley.psem <- psem(
  lme(DD ~ lat, random = ~ 1 | site / tree, na.action = na.omit, data = shipley),
  lme(Date ~ DD, random = ~ 1 | site / tree, na.action = na.omit, data = shipley),
  lme(Growth ~ Date, random = ~ 1 | site / tree, na.action = na.omit, data = shipley),
  glmer(Live ~ Growth + (1 | site) + (1 | tree), family = binomial(link = "logit"), data = shipley) 
)
(new.summary <- summary(shipley.psem, .progressBar = F,conditional=T))
new.summary$IC$AIC

summary(shipley.psem)

# Model (correct)
# x1  -> x2 -> x3   -> x4     -> x5
# lat -> DD -> Date -> Growth -> Live

# lat _||_ date | DD 
# lat _||_ Growth | Date 
# lat _||_ Live | Growth 
# DD _||_ Growth | lat, Date #Includes parents of DD as well as Growth
# DD _||_ Live | lat, Growth
# Date _||_ Live | DD, Growth 

# Independ.Claim 
# Date  ~  lat  + DD - correct
# Growth  ~  lat  + Date - correct
# Live  ~  lat  + Growth - correct
# Growth  ~  DD  + lat  + Date - correct
# Live  ~  DD  + lat  + Growth - correct
# Live  ~  Date  + DD  + Growth - correct

# Writing my own d-sep claims list. Appears that dagitty doesn't have this. --------------
shipley.test <- function(g){ 
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
  for(i in 1:length(claims)) {
    if(length(ancestors(g,claims[[i]]$X))<length(ancestors(g,claims[[i]]$Y))){ #If x-value has fewer causal ancestors
      X <- claims[[i]]$Y # Switch position of x and y
      claims[[i]]$Y <- claims[[i]]$X 
      claims[[i]]$X <- X
    }
  } 
  rm(X) #Cleanup
  
  #3. For each of the k pairs of variables (Xi, Xj), list the set of other
  #variables, {Z} in the graph that are direct causes of either Xi or Xj. 
  for(i in 1:length(claims)){
    var1 <- claims[[i]]$X
    var2 <- claims[[i]]$Y
    parents1 <- parents(g,var1)
    parents2 <- parents(g,var2)
    claims[[i]]$Z <- unique(c(parents1,parents2))
  }
  
  # Sort by causal rank
  claims <- claims[order(sapply(sapply(claims,function(x) x$X),function(x,dag) length(ancestors(dag,x)), dag=g),
                         sapply(claims,function(x) x$X),
                         sapply(sapply(claims,function(x) x$Y),function(x,dag) length(ancestors(dag,x)), dag=g))]
  
  for(i in 1:length(claims)) claims[[i]]$Z <- sort(claims[[i]]$Z)
  return(claims)
}

g1 <- dagitty("dag{x1 -> x2 -> x3 -> x4 -> x5}")
shipley.test(g1)

#Stated claims list from Shipley 2009
# x1 _||_ x3 | x2 
# x1 _||_ x4 | x3 
# x1 _||_ x5 | x4 
# x2 _||_ x4 | x1, x3 #Includes parents of x2 as well as x4
# x2 _||_ x5 | x1, x4
# x3 _||_ x5 | x2, x4 

#Results from my function. Looks OK.
# x3 _||_ x1 | x2
# x4 _||_ x1 | x3
# x4 _||_ x2 | x1, x3
# x5 _||_ x1 | x4
# x5 _||_ x2 | x1, x4
# x5 _||_ x3 | x2, x4


# Dagitty claims list for commodity fields --------------------------------
commDAG <- dagitty(" dag {
        plDens <- is2015; plDens <- isIrrig; plDens <- hbeeDist; plDens <- isGP;
        plSize <- plDens; plSize <- hbeeDist; plSize <- isGP; plSize <- is2015; plSize <- isIrrig;
        flDens <- plSize; flDens <- hbeeDist;
        hbeeVis <- is2015; hbeeVis <- isGP; hbeeVis <- hbeeDist; hbeeVis <- numHives; 
        hbeeVis <- flDens; hbeeVis <- isIrrig;
        pollen <- hbeeVis; pollen <- hbeeDist;
        flwSurv <- hbeeVis; flwSurv <- pollen; flwSurv <- plDens; flwSurv <- isIrrig; flwSurv <- is2015;
        flwSurv <- plSize; 
        flwCount <- is2015; flwCount <- plSize; flwCount <- flwSurv;
        seedCount <- hbeeVis; seedCount <- pollen;  seedCount <- is2015; seedCount <- plSize;
        seedWeight <- hbeeVis; seedWeight <- pollen; seedWeight <- isIrrig; seedWeight <- is2015
        seedWeight <- plSize; seedWeight <- seedCount
}")

plot(graphLayout(commDAG))

shipley.test(commDAG)


# Dagitty claims list for seed fields -------------------------------------

seedDAG <- dagitty( "dag {
              plDens <- hbeeDist 
              plSize <- hbeeDist; plSize <- plDens
              flDens <- plSize; flDens <- isMbay; flDens <- is2016; flDens <- hbeeDist;
              hbeeVis <- flDens; hbeeVis <- hbeeDist; hbeeVis <- lbeeDist;
              hbeeVis <- lbeeVis; hbeeVis <- isMbay; hbeeVis <- isCent
              lbeeVis <- lbeeDist; lbeeVis <- hbeeDist; lbeeVis <- isMbay
              lbeeVis <- isCent; lbeeVis <- lbeeStocking; lbeeVis <- is2016; lbeeVis <- flDens;
              pol <- hbeeVis; pol <- lbeeVis; pol <- isCent; pol <- hbeeDist; pol <- flDens
              flwCount <- plSize; flwCount <- isCent; flwCount <- surv;
              surv <- pol; surv <- plSize; surv <- isCent; surv <- hbeeDist; surv <- lbeeDist; surv <- flDens
              seedCount <- pol; seedCount <- plSize; seedCount <- isCent; seedCount <- hbeeDist;
              seedCount <- flDens; seedCount <- surv; 
              seedWeight <- pol; seedWeight <- seedCount; seedWeight <- plSize; seedWeight <- plDens;
              seedWeight <- is2016; seedWeight <- lbeeDist; seedWeight <- lbeeStocking
              }"
)
plot(graphLayout(seedDAG))

shipley.test(seedDAG)













