#R script to look at how basis sets are constructed for conditional independence claims

# Tests -------------------------------------------------------------------
library(dagitty)
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
(new.summary <- summary(shipley.psem, .progressBar = F))
new.summary$IC$AIC

summary(shipley.psem)


# Dagitty claims list for seed fields -------------------------------------


g1 <- dagitty( "dag {
               flDens <- plSize <- plDens <- hbeeDist
               plSize <- hbeeDist
               flDens <- isMbay; flDens <- is2016; flDens <- hbeeDist
               hbeeVis <- flDens; hbeeVis <- hbeeDist; hbeeVis <- lbeeDist;
               hbeeVis <- lbeeVis; hbeeVis <- isMbay; hbeeVis <- isCent
               lbeeVis <- lbeeDist; lbeeVis <- hbeeDist; lbeeVis <- isMbay
               lbeeVis <- isCent; lbeeVis <- lbeeStocking; lbeeVis <- is2016; lbeeVis <- flDens;
               pol <- hbeeVis; pol <- lbeeVis; pol <- isCent; pol <- hbeeDist; pol <- flDens
               flwCount <- plSize; flwCount <- isCent; flwCount <- pol; flwCount <- flDens
               surv <- pol; surv <- plSize; surv <- isCent; surv <- hbeeDist; surv <- lbeeDist; surv <- flDens
               seedCount <- pol; seedCount <- plSize; seedCount <- isCent; seedCount <- hbeeDist;
               seedCount <- flDens; seedCount <- surv; 
               seedWeight <- pol; seedWeight <- seedCount; seedWeight <- plSize; seedWeight <- is2016
               seedWeight <- lbeeDist; seedWeight <- lbeeStocking
               }"
)
plot(graphLayout(g1))

print(impliedConditionalIndependencies(g1,type='basis.set'))

