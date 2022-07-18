#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN COMMODITY CANOLA FIELDS (2014+2015)

# Libraries and ggplot theme ---------------------------------------------------------
library(ggplot2)
library(MCMCglmm)
library(nlme)
library(mgcv)
library(dplyr)
library(tidyr)

#Big-text theme, no grid lines (used for Bayer 2016 presentation)
prestheme=theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15),
                axis.title=element_text(size=20),
                title=element_text(size=20),
                axis.line.x=element_line(colour='black'),
                axis.line.y=element_line(colour='black'),
                #panel.grid.major=element_line(size=0.5,colour='black',linetype='dotted'),
                #panel.grid.minor=element_blank,
                panel.border=element_blank(),
                strip.text=element_text(size=15))
theme_set(theme_classic()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

#Load and organize everything -----------

#Load from Rdata file
load("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\canola_yield_project\\Commodity field analysis\\commodityfieldDataAll.RData")

# fields2014=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\field info 2014.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# survey2014=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\survey data 2014.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# nectar2014=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\nectar calibration.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# plants2014=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\plant info 2014 commodity.csv",stringsAsFactors=T,strip.white=T,na.strings= c("",'NA'))
# 
# fields2015=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\field info 2015.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# survey2015=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\survey data 2015.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# plants2015=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\plant info 2015 commodity.csv",stringsAsFactors=T,strip.white=T,na.strings= c("",'NA'))
# 
# AICc=function(model,...){ #Second-order bias corrected AIC
#   #require(plyr)
#   corAIC=function(x){
#     k=length(coef(x))
#     n=length(residuals(x))
#     return(AIC(x)+(2*k*k-1)/(n-k-1))
#   }
#   if(missing(...)){
#     return(corAIC(model))
#   }
#   else{
#     mlist=list(model,...)
#     len=plyr::laply(mlist,function(x) length(residuals(x)))
#     if(length(len)!=sum(len %in% len[1]))
#       warning('Models use different # data points. Comparing AICc values is invalid')
#     forms=sapply(mlist, function(x) paste(formula(x)))
#     data.frame(row.names=paste(forms[2,],forms[1,],forms[3,]),AICc=plyr::laply(mlist,corAIC))
#   }
# }
# 
# deltaAIC=function(m1,m2,...,correct=TRUE){
#   if(correct) a=AICc(m1,m2,...) else a=AIC(m1,m2,...)
#   a$delta_AICc=a$AICc-min(a$AICc)
#   a$w_i=exp(-0.5*a$delta_AICc)/sum(exp(-0.5*a$delta_AICc)) #Model Probabilities (Akaike weights)
#   a
# }
# 
# DIC=function(m1,...){ #Convenience function to extract DIC values from MCMCglmm models
#   modlist=list(m1,...)
#   dicvals=unlist(lapply(modlist,'[[',c('DIC')))
#   fixeffs=as.character(lapply((lapply(modlist,'[[','Fixed')),'[[','formula'))
#   raneffs=as.character(lapply((lapply(modlist,'[[','Random')),'[[','formula'))
#   data.frame(Fixed=fixeffs,Random=raneffs,DIC=dicvals,deltaDIC=dicvals-min(dicvals))
# }
# 
# plotFixedTerms=function(m1){ #Convenience function to plot 95% distributions for MCMCglmm fixed effects
#   dat=as.data.frame(summary(m1)$solutions)
#   dat$params=as.factor(rownames(dat))
#   colnames(dat)[c(2,3)]=c('lwr','upr')
#   ggplot(dat,aes(x=params,y=post.mean,ymax=upr,ymin=lwr))+geom_pointrange()+labs(x='Fixed Effects',y='Posterior Mean')+geom_hline(yintercept=0,colour='red')+theme(axis.text.x=element_text(angle=90))
# }
# 
# se=function(x) sd(x)/sqrt(length(x)) #Convenience SE function
# 
# zeroCheck=function(mod,dat) { #Posterior predictive checks for zero-inflation (MCMCglmm course notes pg. 104)
#   require(MCMCglmm)
#   require(ggplot2)
#   try(if(is.null(mod$X)) stop("Needs fixed effect design matrix (use saveX=T in MCMCglmm arguments)"))
#   nz=1:nrow(mod$Sol) #Number of samples to hold
#   oz=sum(dat[as.character(mod$Fixed$formula[2])]==0) #Actual number of 0s in the dataset
#   for (i in 1:nrow(mod$Sol)) {
#     pred.l <- rnorm(nrow(dat), (mod$X %*% mod$Sol[i,])@x, sqrt(mod$VCV[i])) #Predicted mean
#     nz[i] <- sum(rpois(nrow(dat), exp(pred.l)) == 0)} #Expected number under poisson
#   qplot(nz,xlab='Expected number of zeros')+geom_vline(aes(xintercept=oz))
#   #Expected zeros, given the overdispersed poisson model
# }
# 
# varComp=function(mod,units=T) { #Tabulates variance components of random effects for MCMCglmm models
#   VCV=mod$VCV
#   if(units==F) VCV=VCV[,c(1:ncol(VCV)-1)] #Strips out "units" column (for binomial models)
#   round(data.frame(mode=posterior.mode(VCV/rowSums(VCV)),HPDinterval(VCV/rowSums(VCV))),3)
# }
# 
# #Predictions on the data scale for MCMCglmm
# predFixef=function(mod,X=NA,link=NA) { #Model, X coefficients, and scale ('gaussian','poisson','binomial')
#   if(length(X)!=ncol(mod$Sol)||is.na(X)) stop('All coefficients must have an X-value')
#   XB=matrix(mod$Sol,ncol=ncol(mod$Sol))%*%X
#   sigma2=matrix(mod$VCV,ncol=ncol(mod$VCV))%*%rep(1,ncol(mod$VCV))
#   #For eqns, see pg 46/47 of Hadfield's course notes. Also, see R-sig thread for discussion of back-transformation, and how to deal with variance from R and G-structures:
#   #https://stat.ethz.ch/pipermail/r-sig-mixed-models/2016q1/024543.html
#   #In Hadfield's eqn: E[y]=exp(XB+Zu+e), E_u,e[y]=exp(XB+0.5*sigma^2) (Poisson link)
#   if(link=='gaussian') {
#     return(mean(XB+sigma2))
#   } else if(link=='poisson') {
#     return(exp(mean(XB+0.5*sigma2)))
#   } else if(link=='binomial'){
#     invlogit=function(a) exp(a)/(1+exp(a))
#     return(invlogit(mean(XB-0.5*sigma2*tanh(XB*(1+2*exp(-0.5*sigma2))/6)))) #McCulloch and Searle 2001 formula
#     #invlogit(XB/sqrt(1+sigma2*(16*sqrt(3)/15*pi)^2)) #Diggle et al 2004 formula
#   } else stop('Link must be gaussian, poisson, or binomial')
# } #Appears to work
# 
# #General field data
# fields2014=
#   fields2014 %>%
#   transmute(Year=2014,Field=paste(Grower,X.), #Concatenates grower name and field number
#        Area=as.factor(Area), #Converts Area to factor
#        Lat=Deg.N+Min.N/60,Lon=-(Deg.W+Min.W/60), #Lat/Lon
#        FieldSize_ha=Field.Size..ha., #Field size (ha)
#        Variety=as.factor(Variety), #Converts Area to factor
#        Irrigated=factor(Irrigated,labels=c('Unirrigated','Irrigated')),
#        Irrigated=factor(Irrigated,levels=c('Irrigated','Unirrigated')),
#        #Converts Irrigated from logical to factor
#        BeeYard=factor(Number.of.Hives>0,labels=c('Unstocked','Stocked')),
#        BeeYard=factor(BeeYard,levels=c('Stocked','Unstocked')),
#        NumHives=Number.of.Hives, #Converts BeeYard from logical to factor
#        Surveyed=as.POSIXct(Surveyed,format='%b %d, %Y'),
#        StartTime=as.POSIXct(paste(Surveyed,Start.Time),format='%Y-%m-%d %H:%M'),
#        EndTime=as.POSIXct(paste(Surveyed,End.Time),format='%Y-%m-%d %H:%M'),
#        Temp=as.numeric(Temperature),WindSp=Wind.Speed.m.s.*3.6, #Weather (wind speed in km/hr)
#        RH=as.numeric(Humidity),Cloud=Cloud.Cover,
#        Cloud=ifelse(grepl('overcast',Cloud),'8/8',Cloud),
#        Cloud=as.numeric(substr(Cloud,0,1))/8) %>%
#   filter(!is.na(Surveyed)) #Strips unsurveyed fields
# 
# fields2015=filter(fields2015,Type=='Commodity') %>%
#   transmute(Year=2015,Field,
#             Area=as.factor(Area), #Converts Area to factor
#             Lat=Deg.N,Lon=-Deg.W, #Lat/Lon
#             FieldSize_ha=Hectares, #Field size (ha)
#             Variety=as.factor(Variety), #Converts Area to factor
#             Irrigated=factor(Irrigated,labels=c('Unirrigated','Irrigated')), #Converts Irrigated from logical to factor
#             BeeYard=factor(Number.of.Hives>0,labels=c('Unstocked','Stocked')), #Converts BeeYard from logical to factor
#             BeeYard=factor(BeeYard,levels=c('Stocked','Unstocked')),
#             NumHives=Number.of.Hives,
#             Surveyed=as.POSIXct(Surveyed,format='%b %d, %Y'),
#             StartTime=as.POSIXct(paste(Surveyed,Start.Time),format='%Y-%m-%d %H:%M'),
#             EndTime=as.POSIXct(paste(Surveyed,End.Time),format='%Y-%m-%d %H:%M'),
#             BowlStart=as.POSIXct(paste(Surveyed,BowlStart),format='%Y-%m-%d %H:%M'),
#             BowlEnd=as.POSIXct(paste(Surveyed,BowlEnd),format='%Y-%m-%d %H:%M'),
#             Temp=Temperature,WindSp=Wind.Speed.km.hr.,RH=Humidity,Cloud=Cloud.Cover) %>% #Wind speed in km/hr
#   filter(!is.na(Surveyed)) #Strips unsurveyed fields
# 
# matches2014=match(survey2014$Field,fields2014$Field)
# 
# survey2014=transmute(survey2014,Field,Year=2014,Distance,Lat=DegN+MinN/60,Lon=-(DegW+MinW/60),
#                  Irrigated=fields2014$Irrigated[matches2014], #Matches Irrigation/Stocking from fields to survey
#                  BeeYard=fields2014$BeeYard[matches2014],
#                  NumHives=fields2014$NumHives[matches2014],
#                  Area=fields2014$Area[matches2014],
#                  Variety=fields2014$Variety[matches2014], #Matches crop varieties
#                  Temp=as.numeric(fields2014$Temp[matches2014]), #Air Temperature
#                  WindSp=fields2014$WindSp[matches2014], #Wind Speed (km/hr)
#                  RH=as.numeric(fields2014$RH[matches2014]), #Relative Humidity
#                  FlDens=Fls.50cm.2*4,  #Flower Density, converting fls/50cm2 to fls/m2
#                  Date=fields2014$Surveyed[matches2014],
#                  StartTime=as.POSIXct(paste(Date,StartTime),format='%Y-%m-%d %H:%M'),
#                  EndTime=as.POSIXct(paste(Date,EndTime),format='%Y-%m-%d %H:%M'),
#                  TotalTime=as.numeric(EndTime-StartTime,units='mins'), #Total time spent on survey in minutes
#                  Honeybee,Leafcutterbee,Bumblebee,Otherbee,Hoverfly,Butterfly,
#                  Fly=LargerFly+SmallerFly+CalliphoridFly, #Combines all "fly" categories (-hoverfly))
#                  Total=Honeybee+Leafcutterbee+Bumblebee+Otherbee+Hoverfly+Butterfly+Fly, #Total visits category
#                  Len1,Len2,Len3,Len4,Len5,Perc1,Perc2,Perc3,Perc4,Perc5,Microcap.volume,
#                  Pollen1,Pollen2,Pollen3,Pollen4,Pollen5,MountingOK,PlDens=Plants.m.2)
# 
# survey2015=filter(survey2015,Type=='Commodity') %>%
#   select(-Type:-EdgeDir,-LeafShelter1:-LeafShelter4,-BagNumbers) #Strips seed field data
# matches2015=match(survey2015$Field,fields2015$Field)
# 
# survey2015=transmute(survey2015,Field,Year=2015,Distance,Lat=DegN,Lon=-DegW,
#                      Irrigated=fields2015$Irrigated[matches2015], #Matches Irrigation/Stocking from fields to survey
#                      BeeYard=fields2015$BeeYard[matches2015],
#                      NumHives=fields2015$NumHives[matches2015],
#                      Area=fields2015$Area[matches2015],
#                      Variety=fields2015$Variety[matches2015], #Matches crop varieties
#                      Temp=fields2015$Temp[matches2015], #Air Temperature
#                      WindSp=fields2015$WindSp[matches2015], #Wind Speed (km/hr)
#                      RH=fields2015$RH[matches2015], #Relative Humidity
#                      FlDens=Fls.50cm.2*4,  #Flower Density, converting fls/50cm2 to fls/m2
#                      Date=fields2015$Surveyed[matches2015],
#                      StartTime=as.POSIXct(paste(Date,StartTime),format='%Y-%m-%d %H:%M'),
#                      EndTime=as.POSIXct(paste(Date,EndTime),format='%Y-%m-%d %H:%M'),
#                      TotalTime=10, #Total time spent on survey in minutes (all plots in 2015 used 10 mins)
#                      HoneybeeTopPollen,HoneybeeTopNectar,HoneybeeSidePollen,HoneybeeSideNectar,
#                      NumHTP,NumHTN,NumHSP,NumHSN,Leafcutterbee,NumLeafcutterbee,Bumblebee,NumBumblebee,
#                      Otherbee,NumOtherbee,Hoverfly,NumHoverfly,Fly,NumFly,
#                      Honeybee=HoneybeeTopPollen+HoneybeeTopNectar+HoneybeeSidePollen+HoneybeeSideNectar,#Total honeybees
#                      NumHoneybee=NumHTP+NumHTN+NumHSP+NumHSN,
#                      Total=HoneybeeTopPollen+HoneybeeTopNectar+HoneybeeSidePollen+HoneybeeSideNectar+
#                        Leafcutterbee+Bumblebee+Otherbee+Hoverfly+Fly, #Total visits
#                      NumTotal=NumHTP+NumHTN+NumHSP+NumHSN+NumLeafcutterbee+NumBumblebee+NumOtherbee+
#                        NumHoverfly+NumFly, #Total number of visits
#                      Len1,Len2,Len3,Len4,Len5,Perc,Microcap.volume,
#                      Pollen1,Pollen2,Pollen3,Pollen4,Pollen5,MountingOK,PlDens=Plants.m.2,Soilwater)
# 
# #Insect data only
# visitors2014=rowwise(survey2014) %>%
#   mutate(AvgNec=mean(c(Len1,Len2,Len3,Len4,Len5),na.rm=T)/Microcap.volume,
#          AvgPerc=mean(c(Perc1,Perc2,Perc3,Perc4,Perc5),na.rm=T),
#          AvgPol=mean(c(Pollen1,Pollen2,Pollen3,Pollen4,Pollen5),na.rm=T))%>%
#   select(Field,Year,Area,Distance,Irrigated,BeeYard,NumHives,Temp,WindSp,RH, #Selects survey2014 data
#          TotalTime,Honeybee:Butterfly,Fly,Total, #Flower visitors (minus "other")
#          FlDens,AvgNec,AvgPerc,AvgPol) %>%
#   filter(!is.na(TotalTime)) %>% #Removes NA rows (from extra nectar measurements)
#   gather('Clade','Visits',Honeybee:Total) %>% #Reshapes dataframe into long form, with pollinator class as a factor.
#   mutate(Field=as.factor(Field),Clade=as.factor(Clade)) #Field and Clade as factors
# 
# visitors2015=rowwise(survey2015) %>%
#   mutate(AvgNec=mean(c(Len1,Len2,Len3,Len4,Len5),na.rm=T)/Microcap.volume,
#          AvgPol=mean(c(Pollen1,Pollen2,Pollen3,Pollen4,Pollen5),na.rm=T))%>%
#   select(Field,Year,Area,Distance,Irrigated,BeeYard,NumHives,Temp,WindSp,RH, #Selects survey2015 data
#          StartTime,TotalTime,HoneybeeTopPollen:NumTotal, #Flower visitors (minus "other")
#          FlDens,AvgNec,Perc,AvgPol,Soilwater) %>%
#   filter(!is.na(TotalTime)) %>% #Removes NA rows (from extra nectar measurements)
#   unite(HoneybeeTopPollen,HoneybeeTopPollen,NumHTP) %>% unite(HoneybeeTopNectar,HoneybeeTopNectar,NumHTN) %>%
#   unite(HoneybeeSidePollen,HoneybeeSidePollen,NumHSP) %>% unite(HoneybeeSideNectar,HoneybeeSideNectar,NumHSN) %>%
#   unite(Leafcutterbee,Leafcutterbee,NumLeafcutterbee) %>% unite(Bumblebee,Bumblebee,NumBumblebee) %>%
#   unite(Otherbee,Otherbee,NumOtherbee) %>% unite(Hoverfly,Hoverfly,NumHoverfly) %>%
#   unite(Fly,Fly,NumFly) %>% unite(Honeybee,Honeybee,NumHoneybee) %>% unite(Total,Total,NumTotal) %>%
#   gather('Clade','Count',HoneybeeTopPollen:Total) %>% #Reshapes dataframe into long form, with pollinator class as a factor.
#   separate(Count,c('Visits','Count'),convert=T,sep='_') %>%
#   mutate(Field=as.factor(Field),Clade=as.factor(Clade)) #Field and Clade as factors
# 
# #Reshapes nectar length, %, and stigma pollen count
# flowers2014=select(survey2014,Field,Year,Area,Variety,Distance,Irrigated,BeeYard,NumHives,StartTime,Len1:MountingOK,FlDens) %>%
#   unite(Samp1,Len1,Perc1,Pollen1) %>% unite(Samp2,Len2,Perc2,Pollen2) %>% unite(Samp3,Len3,Perc3,Pollen3) %>%
#   unite(Samp4,Len4,Perc4,Pollen4) %>% unite(Samp5,Len5,Perc5,Pollen5) %>% #Pairs observations
#   gather(sample,values,Samp1:Samp5) %>% #Gathers samples into 1 column
#   separate(values,c('Len','Perc','Pollen'),sep='_',convert=T) %>% #Separates out observations
#   separate(sample,c('sample','Flower'),sep=4,convert=T) %>% #Separates observation number
#   mutate(Vol=Len*Microcap.volume/32) %>%
#   select(-sample,-Len) #Cleanup
# 
# flowers2015=select(survey2015,Field,Year,Area,Variety,Distance,Irrigated,BeeYard,NumHives,StartTime,Len1:MountingOK,FlDens,Soilwater) %>%
#   unite(Samp1,Len1,Pollen1) %>% unite(Samp2,Len2,Pollen2) %>% unite(Samp3,Len3,Pollen3) %>%
#   unite(Samp4,Len4,Pollen4) %>% unite(Samp5,Len5,Pollen5) %>% #Pairs observations
#   gather(sample,values,Samp1:Samp5) %>% #Gathers samples into 1 column
#   separate(values,c('Len','Pollen'),sep='_',convert=T) %>% #Separates out observations
#   separate(sample,c('sample','Flower'),sep=4,convert=T) %>% #Separates observation number
#   mutate(Vol=Len*Microcap.volume/32) %>%
#   select(-sample,-Len) #Cleanup
# 
# # Microcap calibration
# nectar2014=gather(nectar2014,'Microcap','Len',2:3,convert=T) %>% #Gathers values
#   mutate(Microcap=factor(Microcap,labels=c('2','5')),Microcap=as.numeric(as.character(Microcap)),
#          Measured=(Len/32)*Microcap,Microcap=as.factor(Microcap)) %>%
#   #Converts to numeric, calculates measured volume, then converts back to factor
#   filter(!is.na(Len)) #Strips NAs
# 
# #Conversion formula for %brix (g/100g solution) to mg/uL sugar
# brix2mgul=function(B){ #B is corrected %brix read from refractometer
#   1.177*((0.9988*B) +(0.389*(B^2)/100) +(0.113*(B^3)/10000) +(0.091*(B^4)/1000000) -(0.039*(B^5)/100000000))/100} #Pyke's formula (from Ralph)
# 
# a=lm(Vol~Measured+Measured:Microcap,nectar2014) #Model actual volume based on measured volume
# flowers2014=mutate(flowers2014,Vol=predict(a,data.frame(Measured=Vol,Microcap=as.factor(Microcap.volume))),
#                #Correct the nectar volume for field measurements ("flowers" dataframe)
#                Field=as.factor(Field), #Convert field to factor
#                Sugar=Vol*brix2mgul(Perc), #Calculates mg sugar per nectar reading
#                FullPollen=Pollen>159) %>% #Is flower fully pollinated? (Mequida & Renard 1984 estimate)
#   filter(!(is.na(Pollen)&is.na(Perc)&is.na(Sugar))) #Removes extra nectar measurements
# 
# flowers2015=mutate(flowers2015,Vol=predict(a,data.frame(Measured=Vol,Microcap=as.factor(Microcap.volume))),
#                    #Correct the nectar volume for field measurements ("flowers" dataframe)
#                    Field=as.factor(Field), #Convert field to factor
#                    Sugar=Vol*brix2mgul(Perc), #Calculates mg sugar per nectar reading
#                    FullPollen=Pollen>159) %>% #Is flower fully pollinated? (Mequida & Renard 1984 estimate)
#   filter(!(is.na(Pollen)&is.na(Perc)&is.na(Sugar)))
# 
# #Matches visitation rates to distance
# matches2014=match(paste(flowers2014$Field,flowers2014$Distance),paste(visitors2014$Field, visitors2014$Distance))
# 
# flowers2014=mutate(flowers2014,Time=as.numeric(visitors2014$TotalTime)[matches2014],
#                Honeybee=filter(visitors2014,Clade=='Honeybee')$Visits[matches2014],
#                Fly=filter(visitors2014,Clade=='Fly')$Visits[matches2014],
#                Hoverfly=filter(visitors2014,Clade=='Hoverfly')$Visits[matches2014],
#                Total=filter(visitors2014,Clade=='Total')$Visits[matches2014]) %>%
#   select(-c(Microcap.volume)) #Cleanup field measurements
# 
# matches2015=match(paste(flowers2015$Field,flowers2015$Distance),paste(visitors2015$Field, visitors2015$Distance))
# 
# flowers2015=mutate(flowers2015,Time=as.numeric(visitors2015$TotalTime)[matches2015],
#          Honeybee=filter(visitors2015,Clade=='Honeybee')$Visits[matches2015],
#          Fly=filter(visitors2015,Clade=='Fly')$Visits[matches2015],
#          Hoverfly=filter(visitors2015,Clade=='Hoverfly')$Visits[matches2015],
#          Total=filter(visitors2015,Clade=='Total')$Visits[matches2015]) %>%
#   select(-c(Microcap.volume)) #Cleanup field measurements
# 
# survey2014=filter(survey2014,!is.na(StartTime)) #Removes extra plots (2uL measurements)
# rm(a,matches2014,matches2015,nectar2014)
# 
# #Plant data 2014
# plants2014=arrange(plants2014,FieldName,Plot,Plant) #Re-orders plant measurements
# 
# seeds2014 <- transmute(plants2014,Field=FieldName,Plot,Plant,
#                 VegMass=TotalMass-SeedMass,SeedMass,Branch,
#                 Pods=Pods+Bag+Bag_Tscar,Missing=Missing-Bag,
#                 Pod1,Pod2,Pod3,Pod4,Pod5,
#                 Weigh1,Weigh2,Weigh3,Weigh4,Weigh5,SeedCount) %>% #Selects correct columns
#   gather('Pod','Value',9:18) %>% #Melts pod count and pod weight
#   separate(Pod,into=c('Param','PodNum'),sep=-1) %>%
#   mutate(Param=as.factor(Param),PodNum=as.numeric(PodNum)) %>%
#   spread(Param,Value) %>% #Casts pod count and pod weight back
#   rename(PodCount=Pod,PodMass=Weigh,Pod=PodNum) %>%
#   mutate(FieldPlot=paste(Field,Plot,sep='_')) %>%
#   filter(paste(Field,Plot,Plant,sep='.')!='McKee 5.1.5')   #Removes McKee5.1.5 (DISPUTE IN LABELLING)
# 
# 
# plants2014 <- group_by(seeds2014,Field,Plot,Plant) %>%
#   summarise_at(vars(-FieldPlot),funs(mean(.,na.rm=T),se)) %>% #Mean and SE of metrics for each plant (all mean values except for pod measurements will equal plant-level metrics)
#   select(Field,Plot,Plant,VegMass=VegMass_mean,SeedMass=SeedMass_mean,
#          Branch=Branch_mean,Pods=Pods_mean,Missing=Missing_mean,SeedCount=SeedCount_mean,
#          AvPodCount=PodCount_mean,SEPodCount=PodCount_se,AvPodMass=PodMass_mean,SEPodMass=PodMass_se) %>%
#   mutate(FieldPlot=as.factor(paste(Field,Plot,sep='_'))) #Creates FieldPlot factor
# 
# #Merge Plot to Distance, merging from some other dataframe. Merge other info (variety, stocking, etc.) as well
# dist <- select(survey2014,Field,Distance,Area,Variety,Irrigated,BeeYard,NumHives,StartTime,EndTime,Temp,WindSp,RH,PlDens) %>%
#    unite(FieldPlot,Field,Distance)
# 
# #Merges Distances and other data into seed and plant dataframes
# seeds2014 <- left_join(seeds2014,dist,by='FieldPlot') %>%
#   select(-FieldPlot) %>%
#   rename(Distance=Plot) %>%
#   arrange(Field,Distance,Plant,Pod)
# plants2014 <- left_join(plants2014,dist,by='FieldPlot') %>%
#   select(-FieldPlot) %>%
#   rename(Distance=Plot) %>%
#   mutate(PropMissing=Missing/(Pods+Missing)) %>%
#   arrange(Field,Distance,Plant)
# rm(dist) #Remove distance
# 
# #Plant data 2015
# plants2015 <- arrange(plants2015,FieldName,Plot,Plant) #Re-orders plant measurements
# seeds2015 <- transmute(plants2015,Field=FieldName,Plot,Plant,
#                     VegMass=TotalMass-SeedMass,SeedMass,Branch,
#                     Pods=Pods+Bag+Bag_Tscar,Missing=Missing-Bag,
#                     Pod1,Pod2,Pod3,Pod4,Pod5,
#                     Weigh1,Weigh2,Weigh3,Weigh4,Weigh5,SeedCount) %>% #Selects correct columns
#   gather('Pod','Value',9:18) %>% #Melts pod count and pod weight
#   separate(Pod,into=c('Param','PodNum'),sep=-1) %>%
#   mutate(Param=as.factor(Param),PodNum=as.numeric(PodNum)) %>%
#   spread(Param,Value) %>% #Casts pod count and pod weight back
#   rename(PodCount=Pod,PodMass=Weigh,Pod=PodNum) %>%
#   mutate(FieldPlot=paste(Field,Plot,sep='_'))
# 
# plants2015 <- group_by(seeds2015,Field,Plot,Plant) %>%
#   summarise_at(vars(-FieldPlot),funs(mean(.,na.rm=T),se))%>% #Mean and SE of metrics for each plant (all mean values except for pod measurements will equal plant-level metrics)
#   select(Field,Plot,Plant,VegMass=VegMass_mean,SeedMass=SeedMass_mean,
#          Branch=Branch_mean,Pods=Pods_mean,Missing=Missing_mean,SeedCount=SeedCount_mean,
#          AvPodCount=PodCount_mean,SEPodCount=PodCount_se,AvPodMass=PodMass_mean,SEPodMass=PodMass_se) %>%
#   mutate(FieldPlot=as.factor(paste(Field,Plot,sep='_'))) #Creates FieldPlot factor
# 
# #Merge Plot to Distance, merging from some other dataframe. Merge other info (variety, stocking, etc.) as well
# dist <- select(survey2015,Field,Distance,Area,Variety,Irrigated,BeeYard,NumHives,StartTime,EndTime,Temp,WindSp,RH,PlDens) %>%
#   unite(FieldPlot,Field,Distance)
# 
# #Merges Distances and other data into seed and plant dataframes
# seeds2015 <- left_join(seeds2015,dist,by='FieldPlot') %>%
#   select(-FieldPlot) %>%
#   rename(Distance=Plot) %>%
#   arrange(Field,Distance,Plant,Pod)
# plants2015 <- left_join(plants2015,dist,by='FieldPlot') %>%
#   select(-FieldPlot) %>%
#   rename(Distance=Plot) %>%
#   mutate(PropMissing=Missing/(Pods+Missing)) %>%
#   arrange(Field,Distance,Plant)
# rm(dist) #Remove distance
# 
# #Merge 2014/2015 dataframes
# fieldsAll <- bind_rows(fields2014,fields2015) %>%
#   mutate(Field=factor(Field),Variety=factor(Variety))
# surveyAll <- bind_rows(survey2014,survey2015)%>%
#   mutate(Field=factor(Field))
# flowersAll <- bind_rows(flowers2014,flowers2015)%>%
#   mutate(Field=factor(Field),Variety=factor(Variety))
# visitorsAll <- bind_rows(visitors2014,visitors2015)
# 
# seedsAll <- bind_rows(mutate(seeds2014,Year=2014),mutate(seeds2015,Year=2015))
# plantsAll <- bind_rows(mutate(plants2014,Year=2014),mutate(plants2015,Year=2015))
# 
# #Filters out honeybee visit classifications (side,pollen,nectar) from 2015
# visitorsAll <- visitorsAll %>%
#   filter(Clade!='Total',!grepl('Honeybee.',Clade))
# 
# #Cleanup
# rm(fields2014,fields2015,flowers2014,flowers2015,survey2014,survey2015,visitors2014,seeds2014,seeds2015,plants2014,plants2015)
# 
# #Calculate plot yield in bu/ac (only for plantsAll)
# #50 lbs/bushel * 453.492 g/lb = 22679.6 g/bushel = 4.409249e-05 bushels/g
# #1 acre = 4046.86m2 = 0.0002471052 acres/m2
# #Grams/m2 * (4.409249e-05/0.0002471052) = g/m2 * 0.1784361 = bushels/acre
# conversion <- (1/(50*453.492))/(1/4046.86)
# plantsAll <- mutate(plantsAll,Yield=SeedMass*PlDens*conversion)
# 
# # Save workspace
# save.image("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\canola_yield_project\\Commodity field analysis\\commodityfieldDataAll.RData")
# 
# #Data frame descriptions:
# #FIELDS: general field information (not really used)
# #SURVEY: plot-level info taken from data sheets (not really used)
# #FLOWERS: flower-level nectar and pollen counts
# #VISITORS: plot-level visitation, with mean nectar, %sugar and pollen counts taken from flowers. Broken down by Clade.
# #PLANTS: plant-level yield metrics, matched with Honeybee, Fly, and Total visit
# #SEEDS: flower-level seed mass & seed count, matched with field- and plot-level measurements
# 
# #Folder to save images
# folder="C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Commodity field analysis\\Figures"

# General visitation plots --------------------------------------------------------

#OVERALL VISITORS

#Type of honeybee visits (side-working, etc. from 2015 data)
temp=filter(visitors2015,grepl('Honeybee.',Clade),Distance!=400) %>%
  mutate(Behav=ifelse(grepl('Side',Clade),'Side','Top'),Type=ifelse(grepl('Nectar',Clade),'Nectar','Pollen')) %>%
  mutate(Type=factor(Type),Behav=factor(Behav)) %>%
  group_by(Area,Type,Behav) %>% #,Distance
  summarize(sumVisits=sum(Visits),avgVisits=mean(Visits),sumCount=sum(Count),avgCount=mean(Count))

#Total Visits
p=ggplot(temp,aes(x=Behav,y=sumVisits,fill=Type))+geom_bar(stat='identity',position=position_dodge())+
  facet_grid(~Area)+
  labs(y='Total Visits',fill='Forager\nType',x='Behaviour')+
  scale_fill_manual(values=c('forestgreen','orange'))
ggsave(paste(folder,'hbeeVisType.png',sep='\\'),p,width=8,height=6)

# #Avg visits/hr
# ggplot(temp,aes(x=Behav,y=avgVisits*6,fill=Type))+geom_bar(stat='identity',position=position_dodge())+
#   facet_grid(Distance~Area)+
#   labs(y='Average Visits/hr',fill='Type',x=NULL)+
#   scale_fill_manual(values=c('forestgreen','orange'))
# 
# #Total Count
# ggplot(temp,aes(x=Behav,y=sumCount,fill=Type))+geom_bar(stat='identity',position=position_dodge())+
#   facet_grid(Distance~Area)+
#   labs(y='Total Counts',fill='Type',x=NULL)+
#   scale_fill_manual(values=c('forestgreen','orange'))
# 
# #Avg Count/hr
# ggplot(temp,aes(x=Behav,y=avgCount*6,fill=Type))+geom_bar(stat='identity',position=position_dodge())+
#   facet_grid(Distance~Area)+
#   labs(y='Average Count/hr',fill='Type',x=NULL)+
#   scale_fill_manual(values=c('forestgreen','orange'))

#Total counts
cladeOrder=group_by(visitorsAll,Clade) %>%
  summarize(sumVisits=sum(Visits)) %>%
  arrange(sumVisits)

p=group_by(visitorsAll,Area,BeeYard,Year,Clade) %>%
  summarize(Visits=sum(Visits)) %>%
  mutate(Clade=factor(Clade,levels=cladeOrder$Clade))%>%
  ggplot(aes(Clade,Visits,fill=BeeYard))+
  geom_bar(stat='identity',position=position_dodge())+
  facet_grid(Year~Area)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_fill_manual(values=c('orange','blue'))+labs(fill='Stocking',y='Total Visits')
ggsave(paste(folder,'totalVis.png',sep='\\'),p,width=8,height=6)

#% of total for that Area/Year/Stocking - kind of deceptive, IMO
group_by(temp,Area,BeeYard,Year,Clade) %>%
  summarize(Visits=sum(Visits)) %>%
  spread(Clade,Visits,fill=0) %>%
  rowwise() %>%
  mutate(Total=sum(Bumblebee,Butterfly,Fly,Honeybee,Hoverfly,Leafcutterbee,Otherbee)) %>%
  gather(Clade,Visits,Bumblebee:Otherbee) %>%
  mutate(Visits=Visits*100/Total) %>%
  ggplot(aes(Clade,Visits,fill=Area))+
  geom_bar(stat='identity',position=position_dodge())+
  facet_grid(BeeYard~Year)+
  theme(axis.text.x=element_text(angle=90,hjust=0.5))+labs(y='% Visits')+
  scale_fill_manual(values=c('orange','blue'))

#Visitation 2014,Stocked - heavy zero-inflation, and these models should be run for longer
pr=list(R = list(V = diag(2), nu = 0.02, fix=2), #Fixes residual variance for ZI process at 1
        #G = list(G1 = list(V = diag(2), nu = 0.002))) #Old prior (unused)
        G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

a=MCMCglmm(Visits~trait*(Area*Clade),
           rcov=~idh(trait):units,
           random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
           family='zapoisson',prior=pr,
           data=filter(visitorsAll,BeeYard=='Stocked',Year==2014,
                       Clade!='Bumblebee',Clade!='Leafcutterbee',Clade!='Otherbee'))

b=MCMCglmm(Visits~trait*(Area*Clade),
           rcov=~idh(trait):units,
           random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
           family='zapoisson',prior=pr,
           data=filter(visitorsAll,BeeYard=='Unstocked',Year==2014,
                       Clade!='Bumblebee',Clade!='Leafcutterbee'))

c=MCMCglmm(Visits~trait*(Area*Clade),
           rcov=~idh(trait):units,
           random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
           family='zapoisson',prior=pr,
           data=filter(visitorsAll,BeeYard=='Stocked',Year==2015,
                       Clade!='Bumblebee',Clade!='Butterfly',Clade!='Leafcutterbee'))

d=MCMCglmm(Visits~trait*(Area*Clade),
           rcov=~idh(trait):units,
           random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
           family='zapoisson',prior=pr,
           data=filter(visitorsAll,BeeYard=='Unstocked',Year==2015,
                       Clade!='Butterfly',Clade!='Leafcutterbee'))

predA=with(filter(visitorsAll,BeeYard=='Stocked',Year==2014,
       Clade!='Bumblebee',Clade!='Leafcutterbee',Clade!='Otherbee'),
     data.frame(predict(a,interval='confidence')*6,Area,Clade,BeeYard='Stocked',Year=2014)) %>%unique()
predB=with(filter(visitorsAll,BeeYard=='Unstocked',Year==2014,
                  Clade!='Leafcutterbee',Clade!='Bumblebee'),
       data.frame(predict(b,interval='confidence')*6,Area,Clade,BeeYard='Unstocked',Year=2014))%>%
  unique()
predC=with(filter(visitorsAll,BeeYard=='Stocked',Year==2015,
                  Clade!='Bumblebee',Clade!='Butterfly',Clade!='Leafcutterbee'),
       data.frame(predict(c,interval='confidence')*6,Area,Clade,BeeYard='Stocked',Year=2015))%>%
  unique()
predD=with(filter(visitorsAll,BeeYard=='Unstocked',Year==2015,
                  Clade!='Butterfly',Clade!='Leafcutterbee'),
       data.frame(predict(d,interval='confidence')*6,Area,Clade,BeeYard='Unstocked',Year=2015))%>%
  unique()

pred=bind_rows(predA,predB,predC,predD) %>%
  mutate(fit=round(fit,2),upr=round(upr,2),lwr=round(lwr,2)) %>%
  complete(Area,Clade,BeeYard,Year,fill=list(fit=0,lwr=0,upr=0)) %>% #Fills in non-estimated
  mutate(nonEst=ifelse(fit==0,T,F)) #Non estimated vector

cladeOrder=group_by(pred,Clade) %>%
  summarize(avgfit=mean(fit)) %>%
  arrange(avgfit)

pred=mutate(pred,Clade=factor(Clade,levels=cladeOrder$Clade))

p=ggplot(filter(pred,!nonEst),aes(Clade,fit,col=BeeYard))+
  geom_pointrange(aes(ymax=upr,ymin=lwr),position=position_dodge(width=0.5))+
  facet_grid(Year~Area)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  scale_colour_manual(values=c('orange','blue'))+
  geom_point(data=filter(pred,nonEst),aes(x=Clade,y=fit),col='black')+
  labs(y='Visits/hr',x=NULL,col='Stocking')
ggsave(paste(folder,'visAreaClade.png',sep='\\'),p,width=8,height=6) 

#VISITATION VS DISTANCE

#Honeybee visitation rate with distance
ggplot(filter(visitors,Clade=='Honeybee'&!is.na(Distance)&!is.na(Rate)),
       aes(Distance,round(Rate),colour=Area))+
  geom_smooth(method='glm.nb',se=T,size=1.5)+
  labs(y='Visits/hr')+
  #theme(legend.position='none')+ 
  scale_colour_manual(values=c('darkcyan','darkorange'))+facet_wrap(~BeeYard,ncol=1)


# Visitation with distance -------------------------------------------------------
 
# #OVERALL VISITORS - greater number of flies in Lethbridge, greater number of other wild pollinators in GP. DIC supports Clade*Area model.
# temp=filter(visitors,Clade!='Total')
# 
# vis_a=MCMCglmm(Count~as.numeric(Total.Time,units='hours')+Clade*Area,
#                random=~Field,
#                data=temp,family='poisson',prior=pr)
# 
# pred=mutate(cbind(temp[,c('Clade','Area','Total.Time')],
#                   predict(vis_a,marginal=~Field,interval='confidence'))
#             ,Total.Time=as.numeric(Total.Time))%>%
#   filter(Total.Time==10) %>%
#   unique() %>%
#   select(-Total.Time) %>%
#   mutate(fit=fit*6,lwr=lwr*6,upr=upr*6) #Visits per hour
# 
# p=ggplot(pred,aes(Clade,fit,colour=Area))+ #Way more flies in Lethbridge
#   geom_pointrange(aes(ymin=lwr,ymax=upr),size=1.1,position=position_dodge(width=0.3))+
#   theme(axis.text.x = element_text(size=20,angle = 90,hjust = 1,vjust=0.5))+
#   labs(y='Visits/hr')+ 
#   scale_colour_manual(values=c('darkcyan','darkorange'))
# ggsave(paste(folder,'visitsClade.png',sep='\\'),p,width=8,height=6)
# 
# #DOMINANT VISITORS - only Hoverflies, Flies, & Honeybees
# temp=filter(visitors,Clade=='Hoverfly'|Clade=='Flies'|Clade=='Honeybee')
# 
# pr_var=diag(13)*1e+6
# pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
# pr=list(B=list(mu=matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0),13),V=pr_var),
#         R=list(V=1,nu=0.02),
#         G=list(G1=list(V = 1, nu = 1,alpha.mu=0,alpha.V=2000)))
# 
# vis_b=MCMCglmm(Count~as.numeric(Total.Time,units='hours')+Clade*Area*BeeYard,
#                random=~Field,data=temp,family='poisson',prior=pr,nitt=30000,thin=30)
# 
# pred=mutate(cbind(temp[,c('Clade','Area','BeeYard','Total.Time')],
#                   predict(vis_b,marginal=~Field,interval='confidence'))
#             ,Total.Time=as.numeric(Total.Time))%>%
#   filter(Total.Time==10) %>%
#   unique() %>%
#   select(-Total.Time) %>%
#   mutate(fit=fit*6,lwr=lwr*6,upr=upr*6) #Visits per hour
# 
# p=ggplot(pred,aes(Clade,fit,colour=Area))+ 
#   geom_pointrange(aes(ymin=lwr,ymax=upr),size=1.1,position=position_dodge(width=0.3))+
#   theme(axis.text.x = element_text(size=20,angle = 90,hjust = 1,vjust=0.5),strip.text=element_text(size=15))+
#   labs(y='Visits/hr')+facet_wrap(~BeeYard,ncol=2)+
#   scale_colour_manual(values=c('darkcyan','darkorange'))
# ggsave(paste(folder,'visitsCladeStocking.png',sep='\\'),p,width=8,height=6)
# 
# #VISITATION VS DISTANCE - DIC supports Distance*Area*Clade model, also supports idh(1+Distance):Field term
# 
# temp=droplevels(filter(visitors,Clade!='Total',Clade!='Leafcutterbee',Clade!='Bumblebee',Clade!='Otherbee')) %>%
#   select(Area,Field,Distance,Clade,Count,Total.Time) %>%
#   mutate(Time=as.numeric(Total.Time))
#   
# pr_var=diag(17)*1e+6
# pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
# pr_mu=matrix(rep(0,17))
# pr_mu[2]=1
# pr=list(B=list(mu=pr_mu,V=pr_var),
#         R=list(V=1,nu=0.02),
#         G=list(G1=list(V = diag(2), nu = 1,alpha.mu=c(0,0),alpha.V=diag(2)*2000)))
# vis_c=MCMCglmm(Count~as.numeric(Total.Time,units='hours')+Distance*Area*Clade,
#                random=~idh(1+Distance):Field,
#                data=temp,family='poisson',prior=pr,nitt=50000,thin=100)
# 
# pred=mutate(cbind(temp[,c('Clade','Distance','Area','Total.Time')],
#                   predict(vis_c,interval='confidence')) #Marginalizing across random effects
#             ,Total.Time=as.numeric(Total.Time))%>%
#   filter(Total.Time==10,Distance!=325) %>%
#   unique() %>%
#   select(-Total.Time) %>%
#   mutate(fit=fit*6,lwr=lwr*6,upr=upr*6) #Per hour estimates
# #write.csv(pred,"C:\\Users\\Samuel\\Desktop\\visits_dist_clade.csv") #Summary for Shelley
# 
# p=ggplot(pred,aes(Distance,y=fit,ymax=upr,ymin=lwr,colour=Area))+ #
#   geom_pointrange(position=position_dodge(width=10),size=1)+
#   geom_smooth(method='glm',formula=round(y)~x,family='poisson',se=F,size=1.1)+
#   facet_wrap(~Clade)+
#   theme(strip.text=element_text(size=15))+
#   scale_colour_manual(values=c('darkcyan','darkorange'))+
#   ylim(0,45)+ylab('Visits/hr') 
# ggsave(paste(folder,'visitsDistClade.png',sep='\\'),p,width=8,height=6)

#VISITATION VS DISTANCE (HONEYBEES ONLY) - not enough 20-hive fields to accurately assess pattern
temp=filter(visitorsAll,Clade=='Honeybee',Distance!=1) %>% #Removes 1m plot
  #Lumps stocking rate into categories
  mutate(StockingRate=cut(NumHives,c(-1,19,31,40),labels=c('None(0)','Low(20-30)','High(36-40)')))

tempA=filter(temp,Year==2014)
pr=list(R = list(V = diag(2), nu = 0.2, fix=2), #Fixes residual variance for ZI process at 1
        #G = list(G1 = list(V = diag(2), nu = 0.002))) #Old prior (unused)
        G=list(G1=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000), #Prior for poisson process
               G2=list(V=1, nu=0.002, alpha.mu=0, alpha.V=1000))) #Prior for ZI process

vis_a=MCMCglmm(Visits~trait*(Distance),
                     rcov=~idh(trait):units,
                     random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
                     family='zapoisson',prior=pr,
                     #nitt=50000+5000,thin=50,burnin=5000,verbose=T,
                     data=tempA)

predA=with(tempA,data.frame(round(predict(vis_a,interval='confidence')*6,2),Year,Distance)) %>%
  distinct() %>% arrange(Distance)

tempB=filter(temp,Year==2015)

vis_b=MCMCglmm(Visits~trait*(Distance),
               rcov=~idh(trait):units,
               random=~idh(at.level(trait,1)):Field+idh(at.level(trait,2)):Field,
               family='zapoisson',prior=pr,
               #nitt=50000+5000,thin=50,burnin=5000,verbose=T,
               data=tempB)

predB=with(tempB,data.frame(round(predict(vis_b,interval='confidence')*6,2),Year,Distance)) %>%
  distinct() %>% arrange(Distance)

pred=bind_rows(predA,predB)

p=ggplot(pred,aes(Distance,fit))+geom_line()+facet_wrap(~Year,ncol=1)+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_point(data=temp,aes(x=Distance,y=Visits*6))+labs(y='Visits/hr',x='Distance from honeybees(m)')
ggsave(paste(folder,'hbeeVisDistAll.png',sep='\\'),p,width=8,height=6)

# Nectar volume plots -----------------------------------------------------
temp=filter(flowersAll,!is.na(Vol),Distance!=1,BeeYard=='Stocked',Area=='Lethbridge')

#NECTAR VOLUME VS DISTANCE VS STOCKING VS REGION
pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(2),nu =0.002,alpha.mu=c(0,0),alpha.V=diag(2)*2000)))

vol_a=MCMCglmm(Vol~Distance*factor(Year),random=~idh(1+Distance):Field,
               #prior=pr,nitt=50000,thin=50,
               data=temp,verbose=F)
summary(vol_a)

pred=with(temp,data.frame(round(predict(vol_a,interval='confidence'),2),Year,Distance,BeeYard,Area)) %>%
  distinct() %>% arrange(Year,Area,BeeYard,Distance)

p=ggplot(pred,aes(Distance,fit))+geom_ribbon(aes(ymax=upr,ymin=lwr,fill=BeeYard),alpha=0.3)+
  geom_line(aes(col=BeeYard),size=1)+facet_grid(Year~Area)+
  scale_fill_manual(values=c('darkorange','blue'))+
  scale_colour_manual(values=c('darkorange','blue'))+
  labs(x='Distance from honeybees(m)',fill='Stocking',col='Stocking',y=expression(paste('Nectar/flower(',mu,'L)',sep='')))+
  geom_point(data=temp,aes(x=Distance,y=Vol,col=BeeYard),position=position_jitter(width=8,height=0.02))+ylim(0,3)
ggsave(paste(folder,'nectarVolDistBeeyard.png',sep='\\'),p,width=8,height=6)

#NECTAR VOLUME VS SOIL WATER CONTENT (2015 ONLY)
temp=filter(flowersAll,Year==2015,!is.na(Vol),!is.na(Soilwater)) %>% 
  select(Field,Area,Distance,BeeYard,Vol,Soilwater,Irrigated) %>%
  arrange(Area,Irrigated,Field,Distance)

ggplot(temp,aes(Soilwater,Vol,col=BeeYard))+geom_point()+facet_wrap(~Area,ncol=1)

#In 2014, Unirrigated (p=0.05) and Distance (p=0.004) are significant, but no interaction. 
#in 2015, Distance (p<0.001), but not Irrigation (p=0.912) and no interaction
vol_soilwater=MCMCglmm(Vol~Year*Irrigated+Distance,
                       random=~idh(1+Distance):Field,
                       data=filter(flowersAll,Area=='Lethbridge'))

pred=with(filter(flowersAll,Area=='Lethbridge'),
     data.frame(predict(vol_soilwater,interval='confidence')^2,Year,Irrigated,Distance)) %>%
  distinct() 

p=ggplot(pred,aes(Distance,fit,fill=Irrigated))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line(aes(col=Irrigated),size=1)+facet_wrap(~Year,ncol=1)+
  labs(col='Irrigation',fill='Irrigation',y=expression(paste('Nectar/flower(',mu,'L)',sep='')),x='Distance from honeybees(m)')+
  scale_fill_manual(values=c('blue','red'))+
  scale_colour_manual(values=c('blue','red'))+
  geom_point(data=filter(flowersAll,Area=='Lethbridge'),aes(x=Distance,y=Vol,col=Irrigated),
             position=position_jitter(width=8,height=0.1))+
  ylim(0,3)
ggsave(paste(folder,'nectarVolDistIrrigation.png',sep='\\'),p,width=8,height=6)
# Nectar concentration plots ----------------------------------------------
ylab='Nectar concentration (%)'

ggplot(flowersAll,aes(Field,Perc))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+labs(y=ylab) #Nectar concentration split by field. Doesn't appear to be any large differences

ggplot(flowersAll,aes(Distance,Perc,colour=Irrigated))+geom_point()+labs(y=ylab)+geom_smooth(method='lm',se=F)+scale_colour_manual(values=c('blue','red'))+ylim(c(50,70)) #Nectar concentration vs Distance. 

ggplot(flowersAll,aes(Vol,Perc))+geom_point()+labs(y=ylab,x='Nectar Volume(uL)')+geom_smooth(method='loess')+facet_wrap(~Irrigated) #Nectar concentration vs Volume split by Irrigation. Positive relationship? Possibly some kind of piecewise model? Weird...

flowersAll %>% 
  summarize(Perc=mean(Perc,na.rm=T)) %>% 
  mutate(mgul=brix2mgul(Perc))

# Nectar concentration models ---------------------------------------------
library(Kendall)

#Perc vs Distance
Perc_a=lm(Perc~log(Distance),data=flowersAll) 
qqnorm(resid(Perc_a))#Residuals non-normal

#Mixed effects
Perc_b=lme(Perc~Distance,random=~Distance|Field,data=flowersAll,na.action=na.omit) #SHOULD CONTAIN SOME KIND OF TEST OF AIR TEMP/HEAT INDEX
qqnorm(resid(Perc_b))
qqline(resid(Perc_b)) #Not perfect residuals, but better than lm()

plot(lm(Perc~Distance+Vol+Irrigated,data=flowersAll,na.action=na.omit),which=2) 

#Look into this more... appears to be a significant positive relationship, but residuals are non-normal. Apply centering, play around with mixed models, GAM, etc.

Perc_c=MCMCglmm(Perc~Irrigated,random=~Field,data=flowersAll)


# Nectar production -------------------------------------------------------
# ylab=expression(paste('Nectar(',mu,'L)',sep='')) #Nectar volume
# ylabhr=expression(paste('Nectar(',mu,'L)/hr',sep='')) #Nectar volume/hr
# temp= flowersAll %>%
#   filter(Distance==500) %>%
#   select(Field,Year,Area,Irrigated,Variety,Distance,StartTime,Honeybee,Total,TotalTime=Time,Flower,Vol) %>%
#   left_join(select(fieldsAll,Field,Temp:Cloud),by='Field') %>%
#   mutate(Time=format(StartTime,format='%H:%M'),
#          Time=as.POSIXct(Time,format='%H:%M')-as.POSIXct("8:00",format='%H:%M'),
#          Time=as.numeric(Time,units='hours')) %>% #Time in hours after 7AM (sunrise ~5:45AM during mid-July)
#   mutate(Cloud=ifelse(grepl('overcast',Cloud),'8/8',Cloud),
#          Cloud=as.numeric(substr(Cloud,0,1))/8) 
# 
# ggplot(filter(temp,Field!='Hofer 3',Irrigated=='Unirrigated'),aes(Time,Vol))+
#   geom_point(position=position_jitter(width=.1,height=0.1))+facet_wrap(Year~Area)+
#   labs(y=ylab,x='Hours since 7:00AM')+geom_smooth(method='lm',se=F) #Vol vs Time for low-visited plots
# 
# ggplot(filter(temp,Total<2,Field!='Hofer 3'),aes(Time,Vol,colour=Irrigated))+geom_point()+facet_wrap(Year~Area)+
#   scale_colour_manual(values=c('blue','red'))+
#   labs(y=ylab,x='Hours since 7:00AM')+geom_smooth(method='lm',se=F) #Vol vs Time with Irrigation
# 
# ggplot(filter(temp,Field!='Hofer 3',Irrigated=='Unirrigated'),aes(Time,Vol,colour=factor(Honeybee)))+geom_point()+facet_wrap(~Area)+
#   #scale_colour_manual(values=c('black','darkgrey'))+
#   labs(y=ylab,x='Hours since 7:00AM',colour='Honeybees \n /10mins')+geom_smooth(method='lm',se=F) #Vol vs Time showing hbee visits
# 
# a=lme4::lmer(Vol~Time*Irrigated+(1|Field),data=filter(temp,Field!='Hofer 3',Year=='2014',Area=='Lethbridge'),REML=T)
# fixef(a)[2]+fixef(a)[4] #For Lethbridge 2014 unirrigated, seems like about 0.25 uL/hr is reasonable. Not sure what happens in the irrigated fields (negative trend)
# 
# a=lm(Vol~Time,data=filter(temp,Year=='2015',Area=='Grand Prairie'))
# a=lme4::lmer(Vol~Time+(1|Field),data=filter(temp,Year=='2015',Area=='Grand Prairie'),REML=F)
# fixef(a)[2] #For Grand Prairie 2015, nectar production has a slope of 0.1 uL/hr (p=0.07). 
# 
# ggplot(temp,aes(Temp,Vol,Irrigated=='Unirrigated'))+geom_point()+facet_wrap(~Area) #Vol vs Temperature
# ggplot(temp,aes(RH,Vol,Irrigated=='Unirrigated'))+geom_point()+facet_wrap(~Area) #Vol vs RH
# ggplot(temp,aes(Cloud,Vol,Irrigated=='Unirrigated'))+geom_point()+facet_wrap(~Area) #Vol vs Cloud cover
# 
# ggplot(filter(temp,Field!='Hofer 3',Irrigated=='Unirrigated'),aes(Total,Vol))+geom_point()+facet_wrap(~Area)+
#   labs(y=ylab,x='Total visits/10mins')+geom_smooth(method='lm',se=F) #Vol vs Visits
# 
# 
# pr=list(R=list(V=1,nu=0.2),
#         G=list(G1=list(V = diag(1), nu = 1,alpha.mu=c(0),alpha.V=diag(1)*2000)))
# 
# a=MCMCglmm(Vol~Time,random=~Field,data=filter(temp,Area=='Lethbridge',Field!='Hofer 3',Irrigated=='Unirrigated'),
#             prior=pr,thin=100,burnin=5000,nitt=5000+100000,verbose=F)
# #Time is significant (p=0.02, p=0.10 if Hofer3 included), but Total Visits (p=0.21) doesn't matter (no honeybees were present) 
# 
# b=MCMCglmm(Vol~Total,random=~Field,data=filter(temp,Area=='Grand Prairie',Irrigated=='Unirrigated'),
#             prior=pr,thin=100,burnin=5000,nitt=5000+100000,verbose=F) 
# #Neither Time (p=0.45), Honeybee visits (p=0.15), nor Total visits (p=0.84) are important for Vol in the  G.P.500m plots 

#Pyke's formula (from Ralph) to convert %brix to mg/uL of sugar
brix2mgul=function(B){ #B is corrected %brix read from refractometer
  1.177*((0.9988*B) +(0.389*(B^2)/100) +(0.113*(B^3)/10000) +(0.091*(B^4)/1000000) -(0.039*(B^5)/100000000))/100} 
#Data from Mohr and Jay(1990), using B. napus var. Regent
#Nectar production from the "continuous removal" method - nectar removed from individual flowers throughout the day

Vol=c(0.69,1.32,2.01,2.73,3.40,4.01,4.58)
Vol[2:7]=Vol[2:7]-Vol[1:6]
Perc=c(36.3,49.1,54.1,50.5,49.6,49.9,43.8)
Sugar=Vol*brix2mgul(Perc)

mohrData=data.frame(Time=as.POSIXlt(paste(seq(8,20,2),rep(':00',7),sep=''),format='%H:%M'),Vol=Vol,Perc=Perc,Sugar=Sugar) %>%
  mutate(Vol=cumsum(Vol),Sugar=cumsum(Sugar)) %>% #
  gather('param','value',2:4) %>%
  mutate(param=factor(param,labels=c('Perc (%)','Sugar (mg)','Vol (uL)'))) 

with(filter(mohrData,param=='Vol (uL)'),mean(value[2:7]-value[1:6])/2) #Mean nectar production - 0.32 uL/hr
with(filter(mohrData,param=='Perc (%)'),mean(value[2:7]+value[1:6])/2) #Mean % sugar - 48.9 %
brix2mgul(48.9)*0.32 #Mean sugar production - 0.23 mg/hr

#Nectar production using "periodic removal"
mohrDataPer=data.frame(Time=as.POSIXlt(paste(seq(8,20,2),rep(':00',7),sep=''),format='%H:%M'),Vol=c(0.94,0.86,0.74,0.53,0.57,0.32,0.24),
                       Perc=c(36.7,48.3,56.5,63.1,66.5,66.3,64.5)) %>%
  mutate(Sugar=Vol*brix2mgul(Perc)) %>%
  gather('param','value',2:4) %>%
  mutate(param=factor(param,labels=c('Perc (%)','Sugar (mg)','Vol (uL)'))) 

ggplot(mohrData,aes(Time,value))+geom_point()+geom_line()+facet_wrap(~param,scales='free_y',ncol=1)+
  labs(y=NULL,x='Time of Day',title='Mohr and Jay (1990), Continous Removal')+
  theme(title=element_text(size=15))

ggplot(mohrDataPer,aes(Time,value))+geom_point()+geom_line()+facet_wrap(~param,scales='free_y',ncol=1)+
  labs(y=NULL,x='Time of Day',title='Mohr and Jay (1990), Periodic Removal')+
  theme(title=element_text(size=15))

bind_rows(mutate(mohrData,Removal='Continuous'),mutate(mohrDataPer,Removal='Periodic')) %>%
  ggplot(aes(Time,value,col=Removal))+geom_point()+geom_line()+facet_wrap(~param,scales='free_y',ncol=1)+
  labs(y=NULL,x='Time of Day')+scale_colour_manual(values=c('red','blue'))+
  theme(title=element_text(size=15))

#Data from Mesquida and Renard 1988, using Oleifera variety

#Mean sugar production
(1/brix2mgul(48.05))*0.29 #1.44ul/mg sugar*0.29mg sugar produced over 24 hrs = 0.42uL
#15 Apr vs 20 Apr production
(1/brix2mgul(c(41.08,48.05)))*c(.37,.3) #0.65 vs 0.43uL produced over 24 hrs

#Szabo 1985 - 0.402 ul/24 hr period (0.46 for Regent variety)

# Inferring nectar - Possingham's eqns -------------------------

flowersAll %>%
  filter(Irrigated!='Irrigated') %>%
  select(Year,Area,Honeybee,Time,Vol,FlDens) %>%
  na.omit %>%
  ggplot(aes((Honeybee/(FlDens*Time/60)),Vol))+geom_point()+
#  geom_smooth(method='gam',formula=y~s(x))+
  facet_grid(Year~Area)+
  scale_colour_manual(values=c('blue','red'))+
  labs(x='Honeybee visits/Flower*hr',y='Nectar volume')

#Eqn 6 from Possingham (equilibrium model)
#Mean Volume = Vol_max/(DepletionRate*Vol_max+1)

temp <- filter(flowersAll) %>%
  mutate(VisRate=Total/(FlDens*Time/60)) %>% #Visits per flower per hour (all visitors)
  mutate(BeeVisRate=Honeybee/(FlDens*Time/60)) %>% #Honeybees only
  unite(group,Year,Area) %>% mutate(plot=paste(group,Field,Distance,sep='_')) %>% 
  dplyr::select(Vol,VisRate,BeeVisRate,group,plot) %>% 
  na.omit

#All visitors model
nectarMod <- nlsList(Vol~(Vol_max/((VisRate/Lambda)*Vol_max+1))|group, #Separate models for each year/location
                  data=temp,start=list(Lambda=0.3,Vol_max=1))
nectarMod_all <- nls(Vol~(Vol_max/((VisRate/Lambda)*Vol_max+1)), #Using all visitors
               data=temp,start=list(Lambda=0.3,Vol_max=1))
summary(nectarMod)
summary(nectarMod_all) #Vol_max = 0.82, se:0.022, lambda =0.14, se:0.025 

var(predict(nectarMod_all))/var(temp$Vol) #Poor R^2 values: only 0.064

pred <- predict(nectarMod,newdata=data.frame(VisRate=seq(0,1,0.05)))
pred <- data.frame(VisRate=seq(0,1,0.05),
                group=names(pred),
                predVol=unname(pred)) %>%
  arrange(group,VisRate) %>%
  separate(group,c('Year','Area'),sep='_') %>%
  unique()

annotations <- data.frame(coef(nectarMod)) %>%
  round(.,2) %>%
  mutate(group=rownames(.)) %>%
  separate(group,c('Year','Area'),sep='_') 

p1 <- separate(temp,group,c('Year','Area'),sep='_') %>%
  ggplot(aes(VisRate,Vol))+geom_point(size=1)+
  facet_grid(Year~Area)+
  geom_line(data=pred,aes(x=VisRate,y=predVol),col='red')+
  xlim(0,1.1)+ylim(0,3)+
  labs(y='Nectar (uL)',x='Total Visits/Flower*Hr')+
  geom_label(data=annotations,aes(label=paste('Lambda =',Lambda,'; Vol_max =',Vol_max)),x=0.7,y=2,size=3)

#Bee-only model
nectarMod=nlsList(Vol~(Vol_max/((BeeVisRate/Lambda)*Vol_max+1))|group,
                  data=temp,start=list(Lambda=0.3,Vol_max=1))

nectarMod_all=nls(Vol~(Vol_max/((BeeVisRate/Lambda)*Vol_max+1)),
               data=temp,start=list(Lambda=0.3,Vol_max=1))

summary(nectarMod)
summary(nectarMod_all)

pred=predict(nectarMod,newdata=data.frame(BeeVisRate=seq(0,1,0.05)))
pred=data.frame(BeeVisRate=seq(0,1,0.05),
                   group=names(pred),
                   predVol=unname(pred)) %>%
  arrange(group,BeeVisRate) %>%
  separate(group,c('Year','Area'),sep='_') %>%
  unique()

annotations=data.frame(coef(nectarMod)) %>%
  round(.,2) %>%
  mutate(group=rownames(.)) %>%
  separate(group,c('Year','Area'),sep='_') 

p2=separate(temp,group,c('Year','Area'),sep='_') %>%
  ggplot(aes(BeeVisRate,Vol))+geom_point(size=1)+
  facet_grid(Year~Area)+
  geom_line(data=pred,aes(x=BeeVisRate,y=predVol),col='red')+
  xlim(0,1.1)+ylim(0,3)+
  labs(y='Nectar (uL)',x='Bee Visits/Flower*Hr')+
  geom_label(data=annotations,aes(label=paste('Lambda =',Lambda,'; Vol_max =',Vol_max)),x=0.7,y=2,size=3)


grid::grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size="last"))

#Aggregate models
nectarModAll_bees=nls(Vol~(Vol_max/((BeeVisRate/Lambda)*Vol_max+1)),
                 data=temp,start=list(Lambda=0.3,Vol_max=1))
summary(nectarModAll_bees) #Perhaps it would be better to use this, given large amount of uncertainty

nectarModAll=nls(Vol~(Vol_max/((VisRate/Lambda)*Vol_max+1)),
                 data=temp,start=list(Lambda=0.3,Vol_max=1))
summary(nectarModAll) 


pred=data.frame(VisRate=seq(0,1,0.01),
                predVol=predict(nectarModAll,newdata=data.frame(VisRate=seq(0,1,0.01))))

annotations=round(coef(nectarModAll),3)

p1=temp %>% 
  ggplot(aes(VisRate,Vol))+geom_point(size=1)+
  geom_line(data=pred,aes(x=VisRate,y=predVol),col='red',size=1)+
  xlim(0,1.1)+ylim(0,3)+
  labs(y='Nectar (uL)',x='Total visits per Flower per Hour')

ggsave(paste(folder,'nectVolVisAll.png',sep='\\'),p1,width=8,height=6)

#Alternate form using the mean of each plot - roughly the same answer as before
temp2 <- group_by(temp,plot) %>% 
  summarize(Vol=mean(Vol),VisRate=first(VisRate),BeeVisRate=first(BeeVisRate)) 

nectarVis_bees <- nls(Vol~(Vol_max/((BeeVisRate/Lambda)*Vol_max+1)),
                      data=temp2,start=list(Lambda=0.3,Vol_max=1))
nectarVis_all <- nls(Vol~(Vol_max/((VisRate/Lambda)*Vol_max+1)),
                     data=temp2,start=list(Lambda=0.3,Vol_max=1))

summary(nectarVis_all) #Lambda=0.14, Vol_max=0.82
summary(nectarVis_bees) #Lambda=0.09, Vol_max=0.77


#Alternate form using probability from Eqn 4

#At t = Inf, p(x) = D_lambda(1 - x/x_m)^(D_lambda*x_m-1)

#Function that returns log likelihood of pars (nectar production, max volume), given Eqn 4 
necProb <- function(pars,x,visrate,LL=F){
  necProd <- pars[1]
  x_m <- pars[2]
  D_lambda <- visrate/necProd
  result <- D_lambda*(1 - (x/x_m))^((D_lambda*x_m)-1)
  result <- ifelse(is.nan(result)|result<1e-10,1e-10,result)
  if(LL) return(-sum(log(result))) else return(result)
}

necProb(c(0.15,100),temp$Vol,temp$VisRate,T)

#Comes up with a lower nectar nectar production (0.04) and higher max volume (1.9), but this is sensitive to how volumes that are greater than max nectar are penalized (currently -4.6 on log scale). I'll likely go with the NLS solution from above.
optim(par=c(0.15,1),fn=necProb,
         x=temp$Vol,visrate=temp$VisRate,LL=T,
         method = "L-BFGS-B",lower=c(0,0),upper=c(10,10))

optim(par=c(0.15,1),fn=necProb, #Similar for bees only: production = 0.05, max nect = 1.21
      x=temp$Vol,visrate=temp$BeeVisRate,LL=T,
      method = "L-BFGS-B",lower=c(0,0),upper=c(10,10))

optim(par=c(0.15,1),fn=necProb, #Similar results, nectar production(0.04) and max volume (3.87)
      x=temp2$Vol,visrate=temp2$VisRate,LL=T,
      method = "L-BFGS-B",lower=c(0,0),upper=c(10,10),hessian=T)

optim(par=c(0.15,1),fn=necProb, #Not properly converging
      x=temp2$Vol,visrate=temp2$BeeVisRate,LL=T,
      method = "L-BFGS-B",lower=c(1e-10,1e-10),upper=c(10,10),hessian=T)






# Pollen deposition plots --------------------------------------------------

temp=filter(flowersAll,!is.na(Pollen),Distance!=1)

pr=list(R=list(V=1,nu=0.02),
      G=list(G1=list(V = diag(1), nu = 1,alpha.mu=c(0),alpha.V=diag(1)*2000)))
pol_b=MCMCglmm(Pollen~Distance*factor(Year),random=~Field,data=temp,family='poisson',nitt=45000,burnin=15000,thin=30,prior=pr) 

predB=with(temp,data.frame(round(predict(pol_b,interval='confidence'),2),Distance,Year)) %>%
  distinct()

p=ggplot(predB,aes(Distance,fit))+geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.3)+
  geom_line()+
  geom_point(data=temp,aes(x=Distance,y=Pollen),position=position_jitter(width=2),size=1,alpha=0.5)+
  facet_wrap(~Year,ncol=1)+
  ylim(0,2000)+labs(y='Pollen grains/stigma',x='Distance from honeybees(m)')
ggsave(paste(folder,'pollenDist.png',sep='\\'),p,width=8,height=6)

# Visitation and Nectar Plots ---------------------------------------------------
ylab=expression(paste('Nectar volume (',mu,'L)',sep=''))

ggplot(flowersAll,aes(60*Honeybee/Time,Sugar,colour=as.factor(NumHives)))+geom_point()+
  geom_smooth(method='lm',se=F)+labs(x='Honeybee Visits/hr',y=ylab,col='# of\nHives')
  #Visitation vs Nectar volume

ggplot(flowersAll,aes(60*Fly/Time,Sugar))+geom_point()+geom_smooth(method='lm',se=F)+labs(x='Fly Visits/hr',y=ylab) #Fly Visitation vs Nectar volume 

ggplot(flowersAll,aes(60*Hoverfly/Time,Sugar))+geom_point()+geom_smooth(method='lm',se=F)+labs(x='Hoverfly Visits/hr',y=ylab) #Fly Visitati #Hoverfly visitation vs Nectar volume 

summary(lme(log(Vol)~Honeybee+offset(Time),random=~1|Field,data=flowers,na.action=na.omit))

# Visitation and Nectar Models ---------------------------------------------
temp=filter(flowers,!is.na(Total),!is.na(Vol))

# pr_var=diag(15)*1e+6
# pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
# pr=list(B=list(mu=matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),15),V=pr_var),
pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = 1, nu = 1,alpha.mu=0,alpha.V=2000),
               G2=list(V = 1, nu = 1,alpha.mu=0,alpha.V=2000)))

visnec_1=MCMCglmm(log(Vol)~Total,
                     random=~Field+Variety,
                  data=temp,prior=pr,nitt=50000,thin=30)

visnec_bees=MCMCglmm(log(Vol)~Honeybee,
                     random=~Field+Variety,
                     data=temp,prior=pr,nitt=50000,thin=30)

visnec_flies=MCMCglmm(log(Vol)~Fly,
                      random=~Field+Variety,
                      data=temp,prior=pr,nitt=50000,thin=30)

pred=mutate(data.frame(predict(visnec_bees,interval='confidence'),cbind(temp[,c('Time','Honeybee')])),Rate=Honeybee*(60/Time))%>%
  select(-Time,-Honeybee)%>%
  mutate(fit=exp(fit),lwr=exp(lwr),upr=exp(upr))%>% #Back-transform from log space
  unique()
#write.csv(pred,"C:\\Users\\Samuel\\Desktop\\vol_visits.csv") #Summary for Shelley

p=ggplot(pred,aes(x=Rate,y=fit,ymax=upr,ymin=lwr))+ #Nectar volume with bee visits
  geom_pointrange(position=position_dodge(width=5),size=1)+
  geom_smooth(se=F,size=1)+
  labs(x='Honeybee Visits/hr',y=ylab)
ggsave(paste(folder,'nectarVolBees.png',sep='\\'),p,width=8,height=6)

# Visitation and Pollen Plots ---------------------------------------------------

ylab='Pollen Grains/Stigma'

ggplot(flowers,aes(Total*60/Time,Pollen))+
  geom_point(alpha=(1/5))+
  geom_smooth(method='lm',se=F)+
  labs(x='Total Visits/hr',y=ylab)+theme(legend.position='none') #Total Visitation vs Pollen Deposition
  #geom_smooth(method='nls',formula='y~(a^2)+(2*a*b*x)+(b^2)*x^2',start=list(a=1,b=1),se=F,size=1,colour='black') #Total Visitation vs Pollen Deposition, with a square root transformed regression

ggplot(flowers,aes(Honeybee*60/Time,Pollen))+geom_point(alpha=(1/5))+
  labs(x='Honeybee Visits/hr',y=ylab)+
  geom_smooth(method='lm',se=F) #Honeybee Visitation vs Pollen Deposition 

ggplot(flowers,aes(Fly*60/Time,Pollen))+geom_point(alpha=(1/5))+
  labs(x='Fly Visits/hr',y=ylab)+theme(legend.position='none')+
  geom_smooth(method='lm',se=F,colour='black',size=1.2) #Fly Visitation vs Pollen Deposition 

ggplot(flowers,aes(HoverflyRate,Pollen,colour=Field))+geom_point(alpha=(1/5))+
  labs(x='Hoverfly Visits/hr',y=ylab)+theme(legend.position='none')+
  geom_smooth(method='lm',se=F,colour='black',size=1.2) #Hoverfly Visitation vs Pollen Deposition

# Visitation and Pollen Models ---------------------------------------------

temp=filter(flowers,!is.na(Pollen),!is.na(Distance))
  
summary(lme(log(Pollen+1)~offset(log(Time))+NumHives,random=~1|Field,data=temp)) #No real relationship b/w visitation and pollen load

pr_var=diag(4)*1e+6
pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
pr_mu=matrix(rep(0,4))
pr_mu[2]=1
pr=list(B=list(mu=pr_mu,V=pr_var),
        R=list(V=1,nu=0.02),
        G=list(G1=list(V = 1, nu = 1))) #,alpha.mu=0,alpha.V=2000

polvis_all=MCMCglmm(Pollen~Total,random=~Field,family='poisson',data=temp)
polvis_bee=MCMCglmm(Pollen~Honeybee,random=~Field,family='poisson',data=temp)

# Seed count models and plots (pod-level) ----------------------
#NO DISTANCE, VISITATION, OR NUMHIVES EFFECT. PLANT MASS APPEARS TO BE IMPORTANT, BUT THERE ARE A FEW LARGE PLANT OUTLIERS THAT MAY DRIVE THE PATTERN

temp=select(seeds,Field,Plot,Plant,VegMass,SeedMass,PodCount,PodMass,Distance,BeeYard,NumHives,Honeybee,Flies,Total)%>%
  mutate(PlotField=paste(Field,Plot,sep='.'),PlantField=paste(Field,Plot,Plant,sep='.')) %>%
  mutate(centDist=Distance-125.2) %>% #Centers Distance variable
  na.omit()

ggplot(temp,aes(Total,VegMass))+geom_jitter()+
  geom_smooth(method='glm',family='poisson',se=F)

pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1, alpha.mu=0,alpha.V=2000),
               G2=list(V = diag(1), nu = 1,alpha.mu=0,alpha.V=2000),
               G3=list(V = diag(1), nu = 1,alpha.mu=0,alpha.V=2000))) #Priors

seedNum2=MCMCglmm(PodCount~Distance*NumHives,random=~Field+PlotField+PlantField,family='poisson',
                  nitt=30000,data=temp,prior=pr)

seedNum3=MCMCglmm(PodCount~VegMass,random=~Field+PlotField+PlantField,family='poisson',
                  nitt=30000,data=temp,prior=pr) #Vegetative mass model
seedNum3_reduced=MCMCglmm(PodCount~VegMass,random=~Field+PlotField+PlantField,family='poisson',
                  nitt=30000,data=filter(temp,VegMass<20),prior=pr) #Vegetative mass model with very big plants removed

pred=mutate(cbind(temp[,c('Distance','NumHives')],predict(seedNum2,interval='confidence')))%>%
  unique() 
p=ggplot(pred,aes(x=Distance,y=fit,ymax=upr,ymin=lwr,colour=as.factor(NumHives),group=NumHives))+ 
  geom_pointrange(position=position_dodge(width=5),size=1)+
  geom_smooth(method='lm',se=F,size=1,linetype='dashed')+
  ylab('Seeds per Pod')+labs(colour='Stocking Rate')+  
  scale_colour_manual(values=c('green','orange','red'))
ggsave(paste(folder,'seedsperpodDistance.png',sep='\\'),p,width=8,height=6)

# Seed weight models and plots (pod-level)--------------------------------------------
#NO DISTANCE, VISITATION, OR NUMHIVES EFFECT. PLANT MASS APPEARS TO BE IMPORTANT, BUT THERE ARE A FEW LARGE PLANT OUTLIERS THAT MAY DRIVE THE PATTERN
ggplot(temp,aes(Distance,PodMass,colour=BeeYard))+geom_point()+
  geom_smooth(method='lm',se=F)

podMass=MCMCglmm(PodMass~Distance*NumHives,random=~Field+PlotField+PlantField,
                  data=temp,prior=pr,nitt=30000,thin=20)

pred=mutate(cbind(temp[,c('Distance','NumHives')],predict(podMass,interval='confidence')))%>%
  unique() 
p=ggplot(pred,aes(x=Distance,y=fit*1000,ymax=upr*1000,ymin=lwr*1000,colour=as.factor(NumHives),group=NumHives))+ 
  geom_pointrange(position=position_dodge(width=5),size=1)+
  geom_smooth(method='lm',se=F,size=1,linetype='dashed')+
  ylab('Pod Mass (mg)')+labs(colour='Stocking Rate')+  
  scale_colour_manual(values=c('green','orange','red'))
ggsave(paste(folder,'podmassDistance.png',sep='\\'),p,width=8,height=6)

# Seed count:weight models and plots (pod-level) --------------------------
temp=filter(seeds,!is.na(PodMass),!is.na(PodCount))

ggplot(temp,aes(PodMass,PodCount,colour=Honeybee))+geom_point()+xlim(0,0.2)+ylim(0,40)

seedCountWeight1=MCMCglmm(PodMass~PodCount*NumHives,data=temp)

# Seed weight models and plots (plant-level)--------------------------------------------
#STOCKING RATE HAS SOME KIND OF INFLUENCE ON SEED MASS, BUT NO TIE-IN WITH DISTANCE OR VISITATION RATE.
# ggplot(plants,aes(Distance,SeedMass,colour=as.factor(NumHives),group=NumHives))+ 
#   geom_point()+
#   geom_smooth(method='lm',se=F,size=1,linetype='dashed')+
#   ylab('Seed weight per plant')+labs(colour='Stocking Rate')+  
#   scale_colour_manual(values=c('green','orange','red'))

temp=filter(plants,!is.na(VegMass),!is.na(SeedMass),!is.na(Variety)) %>%
  mutate(PlotField=paste(Field,Plot,sep='.')) %>%
  mutate(centDist=Distance-125.2) #Centers Distance variable

pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1, alpha.mu=0,alpha.V=2000),
               G2=list(V = diag(1), nu = 1,alpha.mu=0,alpha.V=2000))) #Priors

seedMass1=MCMCglmm(sqrt(SeedMass)~VegMass*NumHives,random=~Field+PlotField,prior=pr,data=temp)
seedMass2=MCMCglmm(sqrt(SeedMass)~VegMass+NumHives,random=~Field+PlotField,prior=pr,data=temp)
seedMass3=MCMCglmm(sqrt(SeedMass)~VegMass,random=~Field+PlotField,prior=pr,data=temp)
DIC(seedMass1,seedMass2,seedMass3)

pred=mutate(cbind(temp[,c('VegMass','NumHives')],predict(seedMass1,interval='confidence'))) %>%
 unique()

ggplot(temp,aes(Distance,SeedMass))+geom_point()

ggplot(temp,aes(VegMass,SeedMass,colour=as.factor(NumHives)))+geom_point()+
  scale_colour_manual(values=c('green','orange','red'))+xlim(0,50)+ylim(0,25)

# Pod success models and plots (plant-level)--------------------------------------------
# NO EFFECT OF DISTANCE, NUMHIVES, IRRIGATION, OR VEGMASS. SEEDMASS APPEARS TO BE SOMEWHAT IMPORTANT, BUT THIS IS PROBABLY THE WRONG DIRECTION OF CORRELATION (I.E. SEEDMASS~PODSUCCESS)
temp=select(plantsAll,Field,Distance,Plant,Area,Year,Variety,Distance,PlDens,Irrigated,
            BeeYard,NumHives,VegMass,SeedMass,Pods,Missing) %>%
  mutate(Pods=round(Pods),TotalFls=round(Pods+Missing),FieldPlot=paste(Field,Distance,sep='.'),
         Year=factor(Year,levels=c('2014','2015'))) %>%
  filter(!is.na(Pods),!is.na(Missing),Missing>=0,!is.na(Distance),!is.na(VegMass),Distance!=1)

pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1))) #Priors
flWins=MCMCglmm(cbind(Pods,Missing)~VegMass*BeeYard*Distance*Irrigated-(VegMass:BeeYard:Distance:Irrigated),
                random=~Field,family='multinomial2',nitt=25000,burnin=5000,thin=20,
                data=temp,prior=pr,verbose=F) #Original: VegMass*BeeYard*Distance
summary(flWins) #No effect of Area (p=0.20), Year (0.97), or Area:Year(0.40)
varComp(flWins) #33% variance at field level,15% at Plot, 51% at residual level  

#Marginal plots
#Copies of original model
lowVeg=flWins
massLow=5 #At 5g (0.1 percentile of plant weight)
lowVeg$X[,'VegMass']=massLow
lowVeg$X[,'VegMass:Distance']=massLow*temp$Distance
lowVeg$X[,'VegMass:BeeYardUnstocked']=massLow*as.numeric(temp$BeeYard)
lowVeg$X[,'VegMass:BeeYardUnstocked']=massLow*as.numeric(temp$BeeYard)
lowVeg$X[,'VegMass:BeeYardUnstocked:Distance']=massLow*as.numeric(temp$BeeYard)*temp$Distance
medVeg=flWins
massMed=quantile(temp$VegMass,0.5)
medVeg$X[,'VegMass']=massMed
medVeg$X[,'VegMass:Distance']=massMed*temp$Distance
medVeg$X[,'VegMass:BeeYardUnstocked']=massMed*as.numeric(temp$BeeYard)
medVeg$X[,'VegMass:BeeYardUnstocked:Distance']=massLow*as.numeric(temp$BeeYard)*temp$Distance
highVeg=flWins
massHigh=quantile(temp$VegMass,0.9)
highVeg$X[,'VegMass']=massHigh
highVeg$X[,'VegMass:Distance']=massHigh*temp$Distance
highVeg$X[,'VegMass:BeeYardUnstocked']=massHigh*as.numeric(temp$BeeYard)
highVeg$X[,'VegMass:BeeYardUnstocked:Distance']=massLow*as.numeric(temp$BeeYard)*temp$Distance

#Ran this for an hour, and it still hasn't finished...
predlow=data.frame(Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='Small (5g)',predict(lowVeg,interval='confidence'))
predmed=data.frame(Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='Medium (14g)',predict(medVeg,interval='confidence'))
predhigh=data.frame(Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='High (35g)',predict(highVeg,interval='confidence')) 

predlow$total=predmed$total=predhigh$total=c(temp$Pods+temp$Missing) #Total number of pods

predlow=mutate(predlow,total=c(temp$Pods+temp$Missing),fitprob=fit/total,lwrprob=lwr/total,uprprob=upr/total)
predmed=mutate(predmed,total=c(temp$Pods+temp$Missing),fitprob=fit/total,lwrprob=lwr/total,uprprob=upr/total)
predhigh=mutate(predhigh,total=c(temp$Pods+temp$Missing),fitprob=fit/total,lwrprob=lwr/total,uprprob=upr/total)

pred=rbind(predlow,predmed,predhigh) %>%
  select(-fit:-total) %>%
  mutate_each(funs(round(.,2)),fitprob:uprprob) %>%
  distinct() 

ggplot(pred,aes(Distance,fitprob,linetype=BeeYard))+geom_line(size=1)+
  facet_wrap(~VegMass,ncol=3)+
  geom_ribbon(aes(ymax=uprprob,ymin=lwrprob),alpha=0.3)+
  labs(y='Proportion pod success',linetype='Bee yard')




# Plant-level to pod-level comparison -------------------------------------
# Do we really need pod-level measurements, or are plant-level ones OK?

temp=filter(plantsAll,!is.na(SeedCount),!is.na(Pods),!is.na(AvPodCount),!is.na(AvPodMass))

#Seeds/pod - pod level is consistently higher, therefore oversampling of large pods is likely
ggplot(temp,aes(SeedCount/Pods,AvPodCount))+
  geom_pointrange(aes(ymax=AvPodCount+SEPodCount,ymin=AvPodCount-SEPodCount))+
  geom_smooth(method='rlm')+
  geom_abline(intercept=0,slope=1,linetype='dashed')+
  labs(title='Seeds/Pod',x='Plant Level (#)',y='Pod Level (#)')+ #Hofer 3.1.4 is an outlier
  lims(x=c(5,35),y=c(5,35))

#Weight/pod - ditto
ggplot(temp,aes(SeedMass/Pods,AvPodMass))+
  geom_pointrange(aes(ymax=AvPodMass+SEPodMass,ymin=AvPodMass-SEPodMass))+
  geom_smooth(method='lm')+
  geom_abline(intercept=0,slope=1,linetype='dashed')+
  labs(title='Pod weight',x='Plant Level (g)',y='Pod Level (g)')+ #White 6.4.5 is an outlier
  lims(x=c(0.02,0.1),y=c(0.02,0.1))

#Single seed weight/pod - looks fairly good
ggplot(temp,aes(SeedMass/SeedCount,AvPodMass/AvPodCount))+geom_point()+
  geom_abline(intercept=0,slope=1,linetype='dashed')+
  geom_smooth(method='lm')+
  #geom_pointrange(aes(ymax=AvPodMass/AvPodCount+SEPodMass/SEPodCount,ymin=AvPodMass/AvPodCount-SEPodMass/SEPodCount))+
  labs(title='Per-seed weight',x='Plant Level (g)',y='Pod Level (g)')+
  lims(x=c(0.001,0.005),y=c(0.001,0.005))

#Seed weight/seed number
temp$AvSeedWeight=with(temp,AvPodMass/AvPodCount) #Pod-level
temp$SeedWeight=with(temp,SeedMass/SeedCount) #Plant-level
a=lm(AvSeedWeight~SeedWeight,data=temp,subset=-c(70,27))

# Plant yield with distance -----------------------------------------------

#Plant-level yield
temp=filter(plantsAll,!is.na(Yield),Distance!=1) %>%
  filter(!is.na(SeedMass),!is.na(VegMass),Distance!=1) %>%
  mutate(Stocking=cut(NumHives,breaks=c(-1,20,31,41),labels=c('Zero','Low','High'))) %>%
  unite(FieldPlot,Field,Distance,sep='_',remove=F)
pr=list(R=list(V=1,nu=0.02),#Priors
          G=list(G1=list(V = diag(1), nu = .1),G2=list(V = diag(1), nu = .1))) 
yieldDist=MCMCglmm(Yield~VegMass*Distance*BeeYard+Area,random=~Field+FieldPlot,family='gaussian', #Plant-level yield
                data=temp,verbose=F,nitt=20000+100000,burnin=20000,thin=100,prior=pr) 
summary(yieldDist) #Plots actually yield slightly more with distance
plot(yieldDist) 
varComp(yieldDist) #Variance component: 15% field, 5% plot, 79% residual

#Copies of original model
lowVeg=yieldDist
massLow=quantile(temp$VegMass,0.1)
lowVeg$X[,'VegMass']=massLow
lowVeg$X[,'VegMass:Distance']=massLow*temp$Distance
lowVeg$X[,'VegMass:BeeYardUnstocked']=massLow*as.numeric(temp$BeeYard)
lowVeg$X[,'VegMass:Distance:BeeYardUnstocked']=massLow*as.numeric(temp$BeeYard)*temp$Distance
medVeg=yieldDist
massMed=quantile(temp$VegMass,0.5)
medVeg$X[,'VegMass']=massMed
medVeg$X[,'VegMass:Distance']=massMed*temp$Distance
medVeg$X[,'VegMass:BeeYardUnstocked']=massMed*as.numeric(temp$BeeYard)
medVeg$X[,'VegMass:Distance:BeeYardUnstocked']=massMed*as.numeric(temp$BeeYard)*temp$Distance
highVeg=yieldDist
massHigh=quantile(temp$VegMass,0.9)
highVeg$X[,'VegMass']=massHigh
highVeg$X[,'VegMass:Distance']=massHigh*temp$Distance
highVeg$X[,'VegMass:BeeYardUnstocked']=massHigh*as.numeric(temp$BeeYard)
highVeg$X[,'VegMass:Distance:BeeYardUnstocked']=massHigh*as.numeric(temp$BeeYard)*temp$Distance

pred=rbind(data.frame(Area=temp$Area,Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='Small (5g)',predict(lowVeg,interval='confidence')),
           data.frame(Area=temp$Area,Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='Medium (14g)',predict(medVeg,interval='confidence')),
           data.frame(Area=temp$Area,Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='High (35g)',predict(highVeg,interval='confidence'))) %>%
  mutate_each(funs(round(.,2)),-Area,-Distance,-BeeYard,-VegMass) %>%
  distinct() 

ggplot(pred,aes(Distance,fit,col=VegMass,linetype=BeeYard))+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=VegMass,col=NULL),alpha=0.3)+
  geom_line(size=1)+
  facet_wrap(~Area,ncol=2)+
  labs(y='Yield (bu/ac)',fill='Mass of Plant',col='Mass of Plant',linetype='Bee yard')+
  scale_colour_manual(values=c('blue','red','darkorange'))+
  scale_fill_manual(values=c('blue','red','darkorange'))

#Veg mass dependent on Plant Density
temp=filter(plantsAll,!is.na(SeedMass),!is.na(VegMass),Distance!=1) %>%
  mutate(Stocking=cut(NumHives,breaks=c(-1,20,31,41),labels=c('Zero','Low','High')),Year=factor(Year,levels=c('2014','2015'))) %>%
  filter(Distance!=325) %>% filter(!is.na(PlDens)) %>%
  unite(FieldPlot,Field,Distance,sep='_',remove=F)
  
vegMassDens=MCMCglmm(VegMass~PlDens+Distance,random=~Field,family='gaussian',data=temp,verbose=F)
summary(vegMassDens)

#Seed mass per plant
temp=filter(plantsAll,!is.na(SeedMass),!is.na(VegMass),Distance!=1) %>%
  mutate(Stocking=cut(NumHives,breaks=c(-1,20,31,41),labels=c('Zero','Low','High')),Year=factor(Year,levels=c('2014','2015'))) %>%
  filter(Distance!=325) %>% filter(!is.na(PlDens)) %>%
  unite(FieldPlot,Field,Distance,sep='_',remove=F)
pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(1), nu = .1),G2=list(V = diag(1), nu = .1))) #Priors
seedMassDist=MCMCglmm(SeedMass~VegMass*Distance*BeeYard+Year+PlDens+VegMass:PlDens,random=~Field+FieldPlot,family='gaussian', #Plant-level yield
                        data=temp,verbose=F,nitt=20000+100000,burnin=20000,thin=100,prior=pr) 
summary(seedMassDist) #No effect of Area (p=0.45)
plot(seedMassDist)
varComp(seedMassDist) #Field is 14%, Plot is 2%, residual is 82%

#Copies of original model
lowVeg=seedMassDist
massLow=quantile(temp$VegMass,0.1)
lowVeg$X[,'VegMass']=massLow
lowVeg$X[,'VegMass:Distance']=massLow*temp$Distance
lowVeg$X[,'VegMass:BeeYardUnstocked']=massLow*as.numeric(temp$BeeYard)
lowVeg$X[,'VegMass:Distance:BeeYardUnstocked']=massLow*as.numeric(temp$BeeYard)*temp$Distance
medVeg=seedMassDist
massMed=quantile(temp$VegMass,0.5)
medVeg$X[,'VegMass']=massMed
medVeg$X[,'VegMass:Distance']=massMed*temp$Distance
medVeg$X[,'VegMass:BeeYardUnstocked']=massMed*as.numeric(temp$BeeYard)
medVeg$X[,'VegMass:Distance:BeeYardUnstocked']=massMed*as.numeric(temp$BeeYard)*temp$Distance
highVeg=seedMassDist
massHigh=quantile(temp$VegMass,0.9)
highVeg$X[,'VegMass']=massHigh
highVeg$X[,'VegMass:Distance']=massHigh*temp$Distance
highVeg$X[,'VegMass:BeeYardUnstocked']=massHigh*as.numeric(temp$BeeYard)
highVeg$X[,'VegMass:Distance:BeeYardUnstocked']=massHigh*as.numeric(temp$BeeYard)*temp$Distance

pred=rbind(data.frame(Year=temp$Year,Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='Small (5g)',predict(lowVeg,interval='confidence')),
      data.frame(Year=temp$Year,Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='Medium (14g)',predict(medVeg,interval='confidence')),
      data.frame(Year=temp$Year,Distance=temp$Distance,BeeYard=temp$BeeYard,VegMass='High (35g)',predict(highVeg,interval='confidence'))) %>%
  mutate_each(funs(round(.,2)),-Year,-Distance,-BeeYard,-VegMass) %>%
  distinct() %>%
  mutate(Year=factor(Year,levels=c('2014','2015')))

ggplot(pred,aes(Distance,fit,col=VegMass,linetype=BeeYard))+geom_line(size=1)+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=VegMass,col=NULL),alpha=0.3)+
  facet_wrap(~Year,ncol=2)+
  labs(y='Seed Mass/Plant (g)',fill='Mass of Plant',col='Mass of Plant',linetype='Bee yard')+
  scale_colour_manual(values=c('blue','red','darkorange'))+
  scale_fill_manual(values=c('blue','red','darkorange'))

#Yield per plot
temp=filter(plantsAll,!is.na(SeedMass),!is.na(VegMass),Distance!=1) %>%
  mutate(Stocking=cut(NumHives,breaks=c(-1,20,31,41),labels=c('Zero','Low','High')),Year=factor(Year,levels=c('2014','2015'))) %>%
  unite(FieldPlot,Field,Distance,sep='_',remove=F) %>%
  group_by(FieldPlot) %>%
  summarize(Field=first(Field),Distance=first(Distance),VegMass=mean(VegMass,na.rm=T),
            Area=first(Area),BeeYard=first(BeeYard),Year=first(Year),Yield=mean(Yield)) %>%
  ungroup() %>%
  arrange(Year,Area,Field,Distance) %>%
  select(-FieldPlot)
pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(1), nu = .1))) #Priors
yieldDistPlot=MCMCglmm(Yield~VegMass+Distance*BeeYard-VegMass:Distance:BeeYard,random=~Field,family='gaussian', #Plot-level yield
                      data=temp,verbose=F,nitt=10000+50000,burnin=10000,thin=50,prior=pr) 
summary(yieldDistPlot)
plot(yieldDistPlot)
varComp(yieldDistPlot) #Field accounts for ~55% of variance

#Copies of original model
lowVeg=yieldDistPlot
massLow=5
lowVeg$X[,'VegMass']=massLow
medVeg=yieldDistPlot
massMed=14
medVeg$X[,'VegMass']=massMed
highVeg=yieldDistPlot
massHigh=35
highVeg$X[,'VegMass']=massHigh
pred=rbind(data.frame(Distance=temp$Distance,BeeYard=temp$BeeYard,
                      VegMass='Small (5g)',predict(lowVeg,interval='confidence')),
           data.frame(Distance=temp$Distance,BeeYard=temp$BeeYard,
                      VegMass='Medium (14g)',predict(medVeg,interval='confidence')),
           data.frame(Distance=temp$Distance,BeeYard=temp$BeeYard,
                      VegMass='High (35g)',predict(highVeg,interval='confidence'))) %>%
  mutate_each(funs(round(.,2)),fit:upr) %>%
  distinct()

ggplot(pred,aes(Distance,fit,col=VegMass,linetype=BeeYard))+geom_line(size=1)+
  geom_ribbon(aes(ymax=upr,ymin=lwr,fill=VegMass,col=NULL),alpha=0.3)+
  labs(y='Yield (bu/ac)',fill='Mean Plant Size',col='Mass of Plant',linetype='Bee yard')+
  scale_colour_manual(values=c('blue','red','darkorange'))+
  scale_fill_manual(values=c('blue','red','darkorange'))
  
