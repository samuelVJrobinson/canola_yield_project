#R script to look at how basis sets are constructed for conditional independence claims
library(dagitty)
source('../helperFunctions.R')
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

g1 <- dagitty("dag{x1 -> x2 -> x3 -> x4 -> x5}")
# debugonce(shipley.test)
shipley.test(g1)
g2 <- list(x2 ~ x1, x3 ~ x2, x4 ~ x3, x5 ~ x4)
shipley.test(g2,form=TRUE)

#Stated claims list from Shipley 2009
# x1 _||_ x3 | x2 
# x1 _||_ x4 | x3 
# x1 _||_ x5 | x4 
# x2 _||_ x4 | x1, x3 #Includes parents of x2 as well as x4
# x2 _||_ x5 | x1, x4
# x3 _||_ x5 | x2, x4 

#Results from my function. Looks OK. 
#Reversed dependencies doesn't matter for dSep claims (i.e. x3 _||_ x1 = x1 _||_ x3)
# x3 _||_ x1 | x2
# x4 _||_ x1 | x3
# x4 _||_ x2 | x1, x3
# x5 _||_ x1 | x4
# x5 _||_ x2 | x1, x4
# x5 _||_ x3 | x2, x4
















