#CODE USED IN ANALYSIS OF INSECT VISITATION, NECTAR, AND POLLEN DEPOSITION IN COMMODITY CANOLA FIELDS (2014)

# Load everything in the console ---------------------------------------------------------
library(ggplot2)
library(MCMCglmm)
library(nlme)
library(mgcv)
library(MASS)
library(dplyr)
library(tidyr)
library(reshape2)

prestheme=theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15), 
                axis.title=element_text(size=20),
                title=element_text(size=20),
                panel.grid.major=element_line(size=0.5,colour='black',linetype='dotted'),
                panel.border=element_rect(size=1,colour='black'))
theme_set(theme_bw()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

load("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Commodity field analysis\\commodityfieldData2014.RData")

# fields=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\field info 2014.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# 
# survey=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\survey data 2014.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# 
# nectar=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\nectar calibration.csv",stringsAsFactors=F,strip.white=T,na.strings = c("",'NA'))
# 
# plants=read.csv("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\csv files\\plant info 2014.csv",stringsAsFactors=T,strip.white=T,na.strings = c("",'NA'))
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
#   ggplot(dat,aes(x=params,y=post.mean,ymax=upr,ymin=lwr))+geom_pointrange()+labs(x='Fixed Effects',y='Posterior Mean')+geom_hline(yintercept=0,colour='red')
# }
# 
# se=function(x) sd(x)/sqrt(length(x)) #Convenience SE function
# 
# # Organize data -----------------------------------------------------------
# 
# #General field data
# fields=transmute(fields,Field=paste(Grower,X.), #Concatenates grower name and field number
#        Area=as.factor(Area), #Converts Area to factor
#        Lat=Deg.N+Min.N/60,Lon=-(Deg.W+Min.W/60), #Lat/Lon
#        Variety=as.factor(Variety), #Converts Area to factor
#        Irrigated=factor(Irrigated,labels=c('Unirrigated','Irrigated')),
#        Irrigated=factor(Irrigated,levels=c('Irrigated','Unirrigated')),
#        #Converts Irrigated from logical to factor
#        BeeYard=factor(Number.of.Hives>0,labels=c('Unstocked','Stocked')),
#        BeeYard=factor(BeeYard,levels=c('Stocked','Unstocked')),
#        NumHives=Number.of.Hives, #Converts BeeYard from logical to factor
#        Surveyed=as.POSIXct(Surveyed,format='%b %d, %Y'),
#        StartTime=as.POSIXct(paste(Surveyed,Start.Time),format='%b %d, %Y %H:%M'),
#        EndTime=as.POSIXct(paste(Surveyed,End.Time),format='%b %d, %Y %H:%M'),
#        Temp=Temperature,WindSp=Wind.Speed.m.s.,RH=Humidity,Cloud=Cloud.Cover) %>% #Weather
#   filter(!is.na(Surveyed)) #Strips unsurveyed fields
# 
# matches=match(survey$Field,fields$Field)
# 
# survey=transmute(survey,Field,Distance,Lat=DegN+MinN/60,Lon=-(DegW+MinW/60),
#                  Irrigated=fields$Irrigated[matches], #Matches Irrigation/Stocking from fields to survey
#                  BeeYard=fields$BeeYard[matches],
#                  NumHives=fields$NumHives[matches],
#                  Area=fields$Area[matches],
#                  Variety=fields$Variety[matches], #Matches crop varieties
#                  Temp=fields$Temp[matches], #Air Temperature
#                  WindSp=fields$WindSp[matches], #WindSp Speed
#                  RH=fields$RH[matches], #Relative Humidity
#                  FlDens=Fls.50cm.2*4,  #Flower Density, converting fls/50cm2 to fls/m2
#                  Date=fields$Surveyed[matches],
#                  StartTime=as.POSIXct(paste(Date,StartTime),format='%Y-%m-%d %H:%M'),
#                  EndTime=as.POSIXct(paste(Date,EndTime),format='%Y-%m-%d %H:%M'),
#                  TotalTime=as.numeric(EndTime-StartTime,units='mins'), #Total time spent on survey in minutes
#                  Honeybee,Leafcutterbee,Bumblebee,Otherbee,Hoverfly,Butterfly,
#                  Flies=LargerFly+SmallerFly+CalliphoridFly, #Combines all "fly" categories (-hoverfly))
#                  Total=Honeybee+Leafcutterbee+Bumblebee+Otherbee+Hoverfly+Butterfly+Flies, #Total visits category
#                  Len1,Len2,Len3,Len4,Len5,Perc1,Perc2,Perc3,Perc4,Perc5,Microcap.volume,
#                  Pollen1,Pollen2,Pollen3,Pollen4,Pollen5,MountingOK,PlDens=Plants.m.2)
# #Insect data
# visitors=rowwise(survey) %>%
#   mutate(AvgNec=mean(c(Len1,Len2,Len3,Len4,Len5),na.rm=T)/Microcap.volume,
#          AvgPerc=mean(c(Perc1,Perc2,Perc3,Perc4,Perc5),na.rm=T),
#          AvgPol=mean(c(Pollen1,Pollen2,Pollen3,Pollen4,Pollen5),na.rm=T))%>%
#   select(Field,Area,Distance,Irrigated,BeeYard,NumHives,Temp,WindSp,RH, #Selects survey data
#          TotalTime,Honeybee:Butterfly,Flies,Total, #Flower visitors (minus "other")
#          FlDens,AvgNec,AvgPerc,AvgPol) %>%
#   filter(!is.na(TotalTime)) %>% #Removes NA rows (from extra nectar measurements)
#   gather('Clade','Count',Honeybee:Total) %>% #Reshapes dataframe into long form, with pollinator class as a factor.
#   mutate(Rate=Count*60/as.numeric(TotalTime), #Visits/hr
#        Interaction=Rate/FlDens, #Visits/Flower*hr
#        Field=as.factor(Field),Clade=as.factor(Clade)) #Field and Clade as factors
# 
# #Reshapes nectar length, %, and stigma pollen count
# flowers=select(survey,Field,Area,Variety,Distance,Irrigated,BeeYard,NumHives,StartTime,Len1:MountingOK,FlDens) %>%
#   unite(Samp1,Len1,Perc1,Pollen1) %>% unite(Samp2,Len2,Perc2,Pollen2) %>% unite(Samp3,Len3,Perc3,Pollen3) %>%
#   unite(Samp4,Len4,Perc4,Pollen4) %>% unite(Samp5,Len5,Perc5,Pollen5) %>% #Pairs observations
#   gather(sample,values,Samp1:Samp5) %>% #Gathers samples into 1 column
#   separate(values,c('Len','Perc','Pollen'),sep='_',convert=T) %>% #Separates out observations
#   separate(sample,c('sample','Flower'),sep=4,convert=T) %>% #Separates observation number
#   mutate(Vol=Len*Microcap.volume/32) %>%
#   select(-sample,-Len) #Cleanup
# 
# # Microcap calibration
# nectar=gather(nectar,'Microcap','Len',2:3,convert=T) %>% #Gathers values
#   mutate(Microcap=factor(Microcap,labels=c('2','5')),Microcap=as.numeric(as.character(Microcap)),
#          Measured=(Len/32)*Microcap,Microcap=as.factor(Microcap)) %>% #Converts to numeric, calculates measured volume, then converts back to factor
#   filter(!is.na(Len)) #Strips NAs
# 
# #Conversion formula for %brix (g/100g solution) to mg/uL sugar
# brix2mgul=function(B){ #B is corrected %brix read from refractometer
#   1.177*((0.9988*B) +(0.389*(B^2)/100) +(0.113*(B^3)/10000) +(0.091*(B^4)/1000000) -(0.039*(B^5)/100000000))/100 #Pyke's formula (from Ralph)
# }
# 
# a=lm(Vol~Measured+Measured:Microcap,nectar) #Model actual volume based on measured volume
# flowers=mutate(flowers,Vol=predict(a,data.frame(Measured=Vol,Microcap=as.factor(Microcap.volume))),
#                #Correct the nectar volume for field measurements ("flowers" dataframe)
#                Field=as.factor(flowers$Field), #Convert field to factor
#                Sugar=Vol*brix2mgul(Perc), #Calculates mg sugar per uL 
#                FullPollen=Pollen>159) %>% #Is flower fully pollinated? (Mequida & Renard 1984 estimate)
#   filter(!(is.na(Pollen)&is.na(Perc)&is.na(Sugar))) #Removes extra nectar measurements
# 
# matches=match(paste(flowers$Field,flowers$Distance),paste(visitors$Field, visitors$Distance)) #Matches visitation rates to distance
# flowers=mutate(flowers,Time=as.numeric(visitors$TotalTime)[matches],
#                Honeybee=filter(visitors,Clade=='Honeybee')$Count[matches],
#                Fly=filter(visitors,Clade=='Flies')$Count[matches],
#                Hoverfly=filter(visitors,Clade=='Hoverfly')$Count[matches],
#                Total=filter(visitors,Clade=='Total')$Count[matches]) %>%
#   select(-c(Microcap.volume)) #Cleanup field measurements
# 
# survey=filter(survey,!is.na(StartTime)) #Removes extra plots
# rm(a,matches,brix2mgul,nectar)
# 
# #Plant data
# plants=arrange(plants,FieldName,Plot,Plant) #Re-orders plant measurements
# seeds=transmute(plants,Field=FieldName,Plot,Plant,
#                 VegMass=TotalMass-SeedMass,SeedMass,Branch,
#                 Pods=Pods+Bag+Bag_Tscar,Missing=Missing-Bag,
#                 Pod1,Pod2,Pod3,Pod4,Pod5,
#                 Weigh1,Weigh2,Weigh3,Weigh4,Weigh5,SeedCount) %>% #Selects correct columns
#   gather('Pod','Value',9:18) %>% #Melts pod count and pod weight
#   separate(Pod,into=c('Param','PodNum'),sep=-2) %>%
#   mutate(Param=as.factor(Param),PodNum=as.numeric(PodNum)) %>%
#   spread(Param,Value) %>% #Casts pod count and pod weight back
#   rename(PodCount=Pod,PodMass=Weigh,Pod=PodNum) %>%
#   group_by(Field,Plot,Plant,Pod) %>%
#   mutate(FieldPlot=paste(Field,Plot,sep='_'))%>%
#   filter(paste(Field,Plot,Plant,sep='.')!='McKee 5.1.5') #Removes McKee5.1.5 (DISPUTE IN LABELLING)
# 
# plants=summarise_each(group_by(seeds,Field,Plot,Plant),funs(mean,se),vars=-FieldPlot)%>% #Mean and SE of metrics for each plant (all mean values except for pod measurements will equal plant-level metrics)
#   select(Field,Plot,Plant,VegMass=VegMass_mean,SeedMass=SeedMass_mean,
#          Branch=Branch_mean,Pods=Pods_mean,Missing=Missing_mean,SeedCount=SeedCount_mean,
#          AvPodCount=PodCount_mean,SEPodCount=PodCount_se,AvPodMass=PodMass_mean,SEPodMass=PodMass_se) %>%
#   mutate(FieldPlot=as.factor(paste(Field,Plot,sep='_'))) #Creates FieldPlot factor
# 
# #Merge change Plot to Distance, merging from some other dataframe. Merge other info (variety, stocking, etc.) as well
# dist=select(survey,Field,Distance,Honeybee,Flies,Total,PlDens,Irrigated,BeeYard,NumHives,Area,Variety) %>%
#   distinct(.) %>%
#   group_by(Field) %>%
#   mutate(Plot=rank(Distance)) %>%
#   ungroup(.) %>%
#   mutate(FieldPlot=paste(Field,Plot,sep='_')) %>%
#   select(-Field,-Plot)
# 
# #Merges Distances and other data into seed and plant dataframes
# seeds=left_join(seeds,dist,by='FieldPlot') %>%
#   select(-Plot,-FieldPlot) %>%
#   arrange(Field,Distance,Plant)
# plants=left_join(plants,dist,by='FieldPlot') %>%
#   select(-Plot,-FieldPlot) %>%
#   mutate(PropMissing=Missing/(Pods+Missing)) %>%
#   arrange(Field,Distance,Plant)
# rm(dist) #Remove distance
# 
# # Save workspace
# save.image("C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Seed field analysis\\commodityfieldData2014.RData")

folder="C:\\Users\\Samuel\\Documents\\Projects\\Presentations\\ESA Meeting 2015"
# Other folders
# "C:\\Users\\Samuel\\Documents\\Projects\\UofC\\Canola survey paper\\R figures"

#Data frame descriptions:
#FIELDS: general field information (not really used)
#SURVEY: plot-level info taken from data sheets (not really used)
#FLOWERS: flower-level nectar and pollen counts
#VISITORS: plot-level visitation, with mean nectar, %sugar and pollen counts taken from flowers. Broken down by Clade.
#PLANTS: plant-level yield metrics, matched with Honeybee, Fly, and Total visit
#SEEDS: flower-level seed mass & seed count, matched with field- and plot-level measurements

# Visitation plots --------------------------------------------------------

#OVERALL VISITORS

# Average Visitation Rate
# p=ggplot(summarize(group_by(visitors,Area,Clade),VisitRate=mean(Rate),se=sd(Rate)/sqrt(length(Rate))),
#          aes(Clade,VisitRate,colour=Area))+
#   geom_pointrange(aes(ymin=VisitRate-se,ymax=VisitRate+se),size=1.1,position=position_jitter(width=0.1))+
#   theme(axis.text.x = element_text(size=15,angle = 90,hjust = 1,vjust=0.5))+labs(y='Visits/hr')+ 
#   scale_colour_manual(values=c('darkcyan','darkorange'))
# ggsave(paste(folder,'visitsClade.png',sep='\\'),width=9,height=6)

# #Total Visitation ~ WindSp ... negative relationship
ggplot(filter(visitors,Clade=='Total'),aes(WindSp*3.6,Count))+
  geom_jitter()+labs(x='WindSp Speed (km/hr)')+
  geom_smooth(method='glm.nb')
summary(glm.nb(Count~WindSp,data=filter(visitors,Clade=='Total')))

# #Total Visitation ~ Temperature ... no relationship
# ggplot(filter(visitors,Clade=='Total'),aes(Temp,Count,colour=Area))+geom_jitter()+
#   labs(x='Air Temperature (C)')
# summary(glm.nb(Count~Temp,data=filter(visitors,Clade=='Total')))

# #Total Visitation ~ RH ... no relationship
# ggplot(filter(visitors,Clade=='Total'),aes(RH,Count,colour=Area))+geom_jitter()+
#   labs(x='% RH')
# summary(glm.nb(Count~RH,data=filter(visitors,Clade=='Total')))


# Average Visitation Rate broken down by stocking - apparently fly visitation rates are higher when field is stocked (??)
ggplot(summarize(group_by(visitors,Area,Clade,BeeYard),VisitRate=mean(Rate),se=sd(Rate)/sqrt(length(Rate))),
       aes(Clade,VisitRate,colour=Area))+
  geom_pointrange(aes(ymin=VisitRate-se,ymax=VisitRate+se),size=1.1,position=position_jitter(width=0.1))+
  theme(axis.text.x = element_text(size=15,angle = 90,hjust = 1,vjust=0.5))+labs(y='Visits/hr')+ 
  scale_colour_manual(values=c('darkcyan','darkorange'))+facet_wrap(~BeeYard,ncol=1)

# Average Interaction Rate
# ggplot(visitors,aes(Clade,Rate))+ 
#   geom_jitter()+facet_wrap(~Area,ncol=1)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+labs(y='Visits/(flower*hr)') 
 
#Jittered version of bar plot
# ggplot(visitors,aes(Clade,Interaction))+geom_jitter()+
#   facet_wrap(~Area,ncol=1)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+labs(y='Visits/(flower*hr)') 

#VISITATION VS DISTANCE

#Overall visitation rate with distance
ggplot(filter(visitors,Clade=='Total'&!is.na(Distance)&!is.na(Rate)),
       aes(Distance,round(Rate),colour=Area))+
  geom_smooth(method='glm.nb',se=T,size=1.5)+
  labs(y='Visits/hr')+
  theme(legend.position='none')+ 
  scale_colour_manual(values=c('darkcyan','darkorange'))+facet_wrap(~BeeYard,ncol=1)

#Honeybee visitation rate with distance
ggplot(filter(visitors,Clade=='Honeybee'&!is.na(Distance)&!is.na(Rate)),
       aes(Distance,round(Rate),colour=Area))+
  geom_smooth(method='glm.nb',se=T,size=1.5)+
  labs(y='Visits/hr')+
  #theme(legend.position='none')+ 
  scale_colour_manual(values=c('darkcyan','darkorange'))+facet_wrap(~BeeYard,ncol=1)


#Interaction rate with distance
# ggplot(filter(visitors,Clade=='Total'&!is.na(Distance)&!is.na(Rate)),
#        aes(Distance,Rate/FlDens,group=Field))+
#   geom_point()+
#   geom_smooth(method='lm',se=F,size=.5,linetype='dashed',colour='black')+
#   geom_smooth(method='lm',se=F,size=1.5,colour='black',aes(group=NULL))+
#   #geom_smooth(method='nls',formula='y~a*x^b',start=list(a=1,b=-0.5),se=F,size=1,colour='black')+  
#   labs(y='Visits/flower*hr')+ylim(0,0.5)+scale_x_log10()
# 
# #Visitation rate with distance, split by stocking, colours designating fields
# ggplot(subset(visitors,Clade=='Total'),aes(Distance,Rate,group=Field))+geom_point()+
#   labs(y='Visits/hr')+facet_wrap(~BeeYard,ncol=2)+scale_x_log10()+
#   geom_smooth(method='nls',formula='y~a*x+b',start=list(a=1,b=-0.5),se=F,size=.5,linetype='dashed',colour='black')+
#   geom_smooth(method='nls',formula='y~a*x+b',start=list(a=1,b=-0.5),se=F,size=1.5,colour='black',aes(group=NULL))+
#   ylim(0,100)


#Visitation rate at each distance, split by stocking
# ggplot(droplevels(subset(visitors,Clade=='Total'&Distance!=325)),aes(as.factor(Distance),Rate))+
#   geom_boxplot()+
#   labs(y='Visits/hr')+facet_wrap(~BeeYard,ncol=2)+
#   theme(legend.position='none')
# 
# ggplot(droplevels(subset(visitors,(Clade=='Total'|Clade=='HoneyBee')&Distance!=325)),aes(as.factor(Distance),Rate))+
#   geom_boxplot()+
#   labs(y='Visits/hr')+facet_wrap(~BeeYard,ncol=2)+
#   theme(legend.position='none')
# 
# 
# #Visitation rate with distance, colours designating stocking
# p=ggplot(subset(visitors,Clade=='Total'),aes(Distance,Rate,colour=BeeYard))+geom_point()+
#   geom_smooth(method='lm',se=F)+
#   labs(y='Visits/hr')+
#   ylim(0,100)+scale_colour_manual(values=c('orange','darkgreen'))+
#   geom_smooth(method='nls',formula='y~a*x+b',start=list(a=1,b=-0.5),se=F,size=1)
# ggsave(paste(folder,'visitsDistStocking.png',sep='\\'),width=9,height=6)
# 
# 
# #Visitation rate with distance,  by stocking, colours designating fields
# ggplot(subset(visitors,Clade!='Total'),aes(Distance,Rate,colour=Clade))+
#   geom_point()+geom_smooth(method='lm',se=F)+labs(y='Visits/hr')+
#   facet_wrap(~Area,ncol=2)+scale_x_log10()+ylim(0,100)+
#   theme(strip.text.x = element_text(size = 15))+geom_smooth(method='lm',se=F,size=1,colour='black') #Interaction rate with distance, split by area, and colours designating different types of visitors

# ggplot(subset(visitors,Clade=='Honeybee'),aes(Distance,Interaction))+geom_point()+geom_smooth(method='lm',se=F)+facet_wrap(~Area+BeeYard)+scale_x_log10()+labs(y='Visits/(flower*hr)')+theme(strip.text.x = element_text(size = 15)) #Honey Bees
# 
# ggplot(subset(visitors,Clade!='Honeybee'& Clade!='Total'),aes(Distance,Interaction))+geom_point()+geom_smooth(method='lm',se=F)+facet_wrap(~Area,ncol=1)+scale_x_log10()+labs(y='Visits/(flower*hr)')+theme(strip.text.x = element_text(size = 15)) #Everything-Honey Bees
# 
# ggplot(subset(visitors,Clade=='Flies'),aes(Distance,Interaction))+geom_point()+geom_smooth(method='lm',se=F)+facet_wrap(~Area)+scale_x_log10()+labs(y='Visits/(flower*hr)')+theme(strip.text.x = element_text(size = 15)) #Flies 

# Visitation Models -------------------------------------------------------
# temp=select(visitors,Field,Distance,Clade,Rate) %>%
#   filter(Clade!='Total') %>%
#   dcast(Field+Distance~Clade,value.var='Rate') #Might there be issues with multicollinearity?

filter(visitors,Clade!='Total') %>%
  group_by(Area,Clade) %>%
  summarize(Count=sum(Count),Time=sum(as.numeric(Total.Time))) #Counts for various categories

#OVERALL VISITORS - greater number of flies in Lethbridge, greater number of other wild pollinators in GP. DIC supports Clade*Area model.
temp=filter(visitors,Clade!='Total')

pr_var=diag(15)*1e+6
pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
pr=list(B=list(mu=matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),15),V=pr_var),
        R=list(V=1,nu=0.02),
        G=list(G1=list(V = 1, nu = 1,alpha.mu=0,alpha.V=2000)))

vis_a=MCMCglmm(Count~as.numeric(Total.Time,units='hours')+Clade*Area,
               random=~Field,
               data=temp,family='poisson',prior=pr)

pred=mutate(cbind(temp[,c('Clade','Area','Total.Time')],
                  predict(vis_a,marginal=~Field,interval='confidence'))
            ,Total.Time=as.numeric(Total.Time))%>%
  filter(Total.Time==10) %>%
  unique() %>%
  select(-Total.Time) %>%
  mutate(fit=fit*6,lwr=lwr*6,upr=upr*6) #Visits per hour

p=ggplot(pred,aes(Clade,fit,colour=Area))+ #Way more flies in Lethbridge
  geom_pointrange(aes(ymin=lwr,ymax=upr),size=1.1,position=position_dodge(width=0.3))+
  theme(axis.text.x = element_text(size=20,angle = 90,hjust = 1,vjust=0.5))+
  labs(y='Visits/hr')+ 
  scale_colour_manual(values=c('darkcyan','darkorange'))
ggsave(paste(folder,'visitsClade.png',sep='\\'),p,width=9,height=6)

#DOMINANT VISITORS - only Hoverflies, Flies, & Honeybees
temp=filter(visitors,Clade=='Hoverfly'|Clade=='Flies'|Clade=='Honeybee')

pr_var=diag(13)*1e+6
pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
pr=list(B=list(mu=matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0),13),V=pr_var),
        R=list(V=1,nu=0.02),
        G=list(G1=list(V = 1, nu = 1,alpha.mu=0,alpha.V=2000)))

vis_b=MCMCglmm(Count~as.numeric(Total.Time,units='hours')+Clade*Area*BeeYard,
               random=~Field,data=temp,family='poisson',prior=pr,nitt=30000,thin=30)

pred=mutate(cbind(temp[,c('Clade','Area','BeeYard','Total.Time')],
                  predict(vis_b,marginal=~Field,interval='confidence'))
            ,Total.Time=as.numeric(Total.Time))%>%
  filter(Total.Time==10) %>%
  unique() %>%
  select(-Total.Time) %>%
  mutate(fit=fit*6,lwr=lwr*6,upr=upr*6) #Visits per hour

p=ggplot(pred,aes(Clade,fit,colour=Area))+ 
  geom_pointrange(aes(ymin=lwr,ymax=upr),size=1.1,position=position_dodge(width=0.3))+
  theme(axis.text.x = element_text(size=20,angle = 90,hjust = 1,vjust=0.5),strip.text=element_text(size=15))+
  labs(y='Visits/hr')+facet_wrap(~BeeYard,ncol=2)+
  scale_colour_manual(values=c('darkcyan','darkorange'))
ggsave(paste(folder,'visitsCladeStocking.png',sep='\\'),p,width=9,height=6)

#VISITATION VS DISTANCE - DIC supports Distance*Area*Clade model, also supports idh(1+Distance):Field term

temp=droplevels(filter(visitors,Clade!='Total',Clade!='Leafcutterbee',Clade!='Bumblebee',Clade!='Otherbee')) %>%
  select(Area,Field,Distance,Clade,Count,Total.Time) %>%
  mutate(Time=as.numeric(Total.Time))
  
pr_var=diag(17)*1e+6
pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
pr_mu=matrix(rep(0,17))
pr_mu[2]=1
pr=list(B=list(mu=pr_mu,V=pr_var),
        R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(2), nu = 1,alpha.mu=c(0,0),alpha.V=diag(2)*2000)))
vis_c=MCMCglmm(Count~as.numeric(Total.Time,units='hours')+Distance*Area*Clade,
               random=~idh(1+Distance):Field,
               data=temp,family='poisson',prior=pr,nitt=50000,thin=100)

pred=mutate(cbind(temp[,c('Clade','Distance','Area','Total.Time')],
                  predict(vis_c,interval='confidence')) #Marginalizing across random effects
            ,Total.Time=as.numeric(Total.Time))%>%
  filter(Total.Time==10,Distance!=325) %>%
  unique() %>%
  select(-Total.Time) %>%
  mutate(fit=fit*6,lwr=lwr*6,upr=upr*6) #Per hour estimates
#write.csv(pred,"C:\\Users\\Samuel\\Desktop\\visits_dist_clade.csv") #Summary for Shelley

p=ggplot(pred,aes(Distance,y=fit,ymax=upr,ymin=lwr,colour=Area))+ #
  geom_pointrange(position=position_dodge(width=10),size=1)+
  geom_smooth(method='glm',formula=round(y)~x,family='poisson',se=F,size=1.1)+
  facet_wrap(~Clade)+
  theme(strip.text=element_text(size=15))+
  scale_colour_manual(values=c('darkcyan','darkorange'))+
  ylim(0,45)+ylab('Visits/hr') 
ggsave(paste(folder,'visitsDistClade.png',sep='\\'),p,width=9,height=6)

#VISITATION VS DISTANCE (HONEYBEES ONLY) - No difference b/w Areas, no distance:numhive interaction, weak WindSp:numhive interaction
temp=droplevels(filter(visitors,Clade=='Honeybee')) %>%
  select(Area,Field,Distance,Clade,Count,BeeYard,NumHives,TotalTime,WindSp) %>%
  mutate(Time=as.numeric(TotalTime)) %>%
  mutate(NumHives=ifelse(NumHives>19 & NumHives!=40,20,NumHives)) %>%
  filter(Distance!=325)

ggplot(temp,aes(Distance,Count,colour=as.factor(NumHives)))+
  geom_smooth(method='glm',method.args=list(family='poisson'),se=F)

nfix=3
nrand=2
pr_var=diag(nfix)*1e+6
pr_var[2,2]=1e-06 #Sets Total.Time mean and variance priors to 1 and ~0 ("offset")
pr_mu=matrix(rep(0,nfix))
pr_mu[2]=1
# pr=list(B=list(mu=pr_mu,V=pr_var),
pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = 1,nu =0.002,alpha.mu=0,alpha.V=2000)))
vis_d=MCMCglmm(Count~Distance+NumHives, 
               random=~Field,
               data=temp,family='poisson',nitt=50000,thin=30,prior=pr)

pred=mutate(cbind(temp[,c('Distance','NumHives')],
                  predict(vis_d,interval='confidence')))%>% #Marginalizing across random effects
  unique() %>%
  mutate(fit=fit*6,lwr=lwr*6,upr=upr*6,NumHives=factor(NumHives))  #Per hour estimates
#write.csv(pred,"C:\\Users\\Samuel\\Desktop\\visits_dist_stocking.csv") #Summary for Shelley

p=ggplot(pred,aes(x=Distance,y=fit,fill=NumHives))+ #Honeybee visits with distance and stocking
  geom_line(aes(colour=NumHives),size=1.5)+
  geom_ribbon(aes(ymax=upr,ymin=lwr),alpha=0.4,size=1)+
  geom_point(data=temp,aes(x=Distance,y=Count*6,fill=NULL,colour=as.factor(NumHives)))+
  labs(y='Visits/hr',x='Distance (m)',fill='Hive \nNumber',colour='Hive \nNumber')+  
  scale_colour_manual(values=c('blue','purple','red'))+
  scale_fill_manual(values=c('blue','purple','red'))

folder="C:\\Users\\Samuel\\Documents\\Projects\\Presentations\\Ent Soc Alberta 2016"
ggsave(paste(folder,'beeDistStocking.png',sep='\\'),p,width=8,height=6)


# Nectar volume plots -----------------------------------------------------
# ylab=expression(paste('Nectar volume (',mu,'L)',sep=''))
# 
# temp=unique(data.frame(Area=flowers$Area,Field=flowers$Field)) #Unique field/area
# flowers$Field=factor(flowers$Field,levels=temp[order(temp$Area,temp$Field),]$Field) #Sorts
# 
# ggplot(flowers,aes(Field,Vol,fill=Area))+
#   geom_boxplot()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
#   #facet_wrap(~Area,ncol=1)+
#   labs(y=ylab) #Some large field-to-field differences
# 
# ggplot(flowers,aes(Variety,Vol))+
#   geom_boxplot()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
#   labs(y=ylab) #Cultivar differences

# Nectar volume models ----------------------------------------------------

temp=filter(flowers,!is.na(Vol)) %>%
  mutate(cDist=scale(Distance,T,F))
ylab=expression(paste('Nectar volume (',mu,'L)',sep=''))
# test=group_by(temp,Field,Distance) %>%
#   summarize(Vol=mean(Vol))
# ggplot(test,aes(Distance,Vol))+geom_point()+geom_smooth(method='lm',formula='y~poly(x,2)',se=F)+facet_wrap(~Field)

#NECTAR VOLUME VS DISTANCE - once random effects are properly specified, irrigation and stocking effects are weak, stocking effect strong
pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(2),nu =0.002,alpha.mu=c(0,0),alpha.V=diag(2)*2000),
               G2=list(V = diag(1),nu =0.002,alpha.mu=c(0),alpha.V=diag(1)*2000)))

vol_a=MCMCglmm(sqrt(Vol)~Distance+NumHives,random=~idh(1+Distance):Field+Variety,data=temp,prior=pr,nitt=50000,thin=30) 

pred=mutate(cbind(temp[,c('Distance','NumHives')],
                  predict(vol_a,interval='confidence')))%>% #Marginalizing across random effects
  filter(Distance!=325) %>%
  unique() 
p=ggplot(pred,aes(x=Distance,y=fit^2,ymax=upr^2,ymin=lwr^2,colour=as.factor(NumHives),group=NumHives))+ #Nectar volumes with distance and stocking
#  geom_point(data=temp,aes(Distance,Vol,ymax=NULL,ymin=NULL,colour=NULL))+
  geom_pointrange(position=position_dodge(width=10),size=1)+
  geom_smooth(method='lm',se=F,size=1)+
  ylab(ylab)+labs(colour='Stocking Rate')+ylim(0.25,1.5)+  
  scale_colour_manual(values=c('green','orange','orange','orange','red'))
ggsave(paste(folder,'nectarVolStocking.png',sep='\\'),p,width=9,height=6)

# Nectar concentration plots ----------------------------------------------
ylab='Nectar concentration (%)'

ggplot(flowers,aes(Field,Perc))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+labs(y=ylab) #Nectar concentration split by field. Doesn't appear to be any large differences

ggplot(flowers,aes(Distance,Perc,colour=Irrigated))+geom_point()+labs(y=ylab)+geom_smooth(method='lm',se=F)+scale_colour_manual(values=c('blue','red'))+ylim(c(50,70)) #Nectar concentration vs Distance. 

ggplot(flowers,aes(Vol,Perc))+geom_point()+labs(y=ylab,x='Nectar Volume(uL)')+geom_smooth(method='loess')+facet_wrap(~Irrigated) #Nectar concentration vs Volume split by Irrigation. Positive relationship? Possibly some kind of piecewise model? Weird...

# Nectar concentration models ---------------------------------------------
library(Kendall)

#Perc vs Distance
Perc_a=lm(Perc~log(Distance),data=flowers) 
qqnorm(resid(Perc_a))#Residuals non-normal

#Mixed effects
Perc_b=lme(Perc~Distance,random=~Distance|Field,data=flowers,na.action=na.omit) #SHOULD CONTAIN SOME KIND OF TEST OF AIR TEMP/HEAT INDEX
qqnorm(resid(Perc_b))
qqline(resid(Perc_b)) #Not perfect residuals, but better than lm()

plot(lm(Perc~Distance+Vol+Irrigated,data=flowers,na.action=na.omit),which=2) 

#Look into this more... appears to be a significant positive relationship, but residuals are non-normal. Apply centering, play around with mixed models, GAM, etc.

Perc_c=MCMCglmm(Perc~Irrigated,random=~Field,data=flowers)


# Nectar production -------------------------------------------------------
ylab=expression(paste('Nectar(',mu,'L)',sep='')) #Nectar volume
ylabhr=expression(paste('Nectar(',mu,'L)/hr',sep='')) #Nectar volume/hr
temp= flowers %>%
  filter(Distance==500) %>%
  select(Field,Area,Irrigated,Variety,Distance,StartTime,Honeybee,Total,TotalTime=Time,Flower,Vol) %>%
  left_join(select(fields,Field,Temp:Cloud),by='Field') %>%
  mutate(Time=format(StartTime,format='%H:%M'),
         Time=as.POSIXct(Time,format='%H:%M')-as.POSIXct("8:00",format='%H:%M'),
         Time=as.numeric(Time,units='hours')) %>% #Time in hours after 7AM (sunrise ~5:45AM during mid-July)
  mutate(Cloud=ifelse(grepl('overcast',Cloud),'8/8',Cloud),
         Cloud=as.numeric(substr(Cloud,0,1))/8) 

ggplot(filter(temp,Field!='Hofer 3',Total<2,Irrigated=='Unirrigated'),aes(Time,Vol))+geom_point()+facet_wrap(~Area)+
  labs(y=ylab,x='Hours since 7:00AM')+geom_smooth(method='lm',se=F) #Vol vs Time for low-visited plots
ggplot(filter(temp,Total<2,Field!='Hofer 3'),aes(Time,Vol,colour=Irrigated))+geom_point()+facet_wrap(~Area)+
  scale_colour_manual(values=c('blue','red'))+
  labs(y=ylab,x='Hours since 7:00AM')+geom_smooth(method='lm',se=F) #Vol vs Time with Irrigation

ggplot(filter(temp,Field!='Hofer 3',Irrigated=='Unirrigated'),aes(Time,Vol,colour=factor(Honeybee)))+geom_point()+facet_wrap(~Area)+
  scale_colour_manual(values=c('black','darkgrey'))+
  labs(y=ylab,x='Hours since 7:00AM',colour='Honeybees \n /10mins')+geom_smooth(method='lm',se=F) #Vol vs Time showing hbee visits
ggplot(temp,aes(Temp,Vol,Irrigated=='Unirrigated'))+geom_point()+facet_wrap(~Area) #Vol vs Temperature
ggplot(temp,aes(RH,Vol,Irrigated=='Unirrigated'))+geom_point()+facet_wrap(~Area) #Vol vs RH
ggplot(temp,aes(Cloud,Vol,Irrigated=='Unirrigated'))+geom_point()+facet_wrap(~Area) #Vol vs Cloud cover


ggplot(filter(temp,Field!='Hofer 3',Irrigated=='Unirrigated'),aes(Total,Vol))+geom_point()+facet_wrap(~Area)+
  labs(y=ylab,x='Total visits/10mins')+geom_smooth(method='lm',se=F) #Vol vs Visits


pr=list(R=list(V=1,nu=0.2),
        G=list(G1=list(V = diag(1), nu = 1,alpha.mu=c(0),alpha.V=diag(1)*2000)))

a=MCMCglmm(Vol~Time,random=~Field,data=filter(temp,Area=='Lethbridge',Field!='Hofer 3',Irrigated=='Unirrigated'),
            prior=pr,thin=100,burnin=5000,nitt=5000+100000,verbose=F)
#Time is significant (p=0.02, p=0.10 if Hofer3 included), but Total Visits (p=0.21) doesn't matter (no honeybees were present) 

b=MCMCglmm(Vol~Total,random=~Field,data=filter(temp,Area=='Grand Prairie',Irrigated=='Unirrigated'),
            prior=pr,thin=100,burnin=5000,nitt=5000+100000,verbose=F) 
#Neither Time (p=0.45), Honeybee visits (p=0.15), nor Total visits (p=0.84) are important for Vol in the  G.P.500m plots 

#Data from Mohr and Jay(1990), using Regent variety
mohrData=data.frame(Time=as.POSIXlt(paste(seq(8,20,2),rep(':00',7),sep=''),format='%H:%M'),Vol=c(0.69,1.32,2.01,2.73,3.40,4.01,4.58))
mean(mohrData$Vol[2:7]-mohrData$Vol[1:6])/2 #Mean nectar production estimate from the "continuous removal" method

ggplot(mohrData,aes(Time,Vol))+geom_point()+geom_line()+
  labs(y=ylab,x='Time of Day',title='Mohr and Jay (1990), Continous Removal')+
  theme(title=element_text(size=15))

# Pollen deposition plots --------------------------------------------------
# ylab='Pollen Grains/Stigma'
# flowers$DistFac=as.factor(flowers$Distance)
# 
# #Pollen deposition vs Distance
# ggplot(flowers,aes(Distance,Pollen))+geom_point(position=position_jitter(w=0.05))+scale_x_log10()+
#   geom_smooth(method='lm',se=F,size=1.5,colour='black')+labs(y=ylab)
# 
# ggplot(flowers,aes(Distance,Pollen,colour=Field))+geom_point()+geom_smooth(method='lm',se=F)+geom_smooth(method='lm',se=F,size=1.5,colour='black')+theme(legend.position='none')+labs(y=ylab) #Pollen deposition vs Distance for each field
# 
# #Pollen deposition vs Distance for each field, with colour for stocking
# ggplot(flowers,aes(Distance,Pollen,colour=BeeYard,size=))+geom_point()+geom_smooth(method='lm',se=F,size=1.5)+
#   labs(y=ylab)+scale_colour_manual(values=c('orange','darkgreen')) 
# ggsave(paste(folder,'pollenDistStocking.png',sep='\\'),width=9,height=6)
# 
# ggplot(flowers,aes(Pollen))+geom_histogram(aes(y=..count..),binwidth=50)+facet_grid(DistFac~.)+xlim(0,4000) #Histogram of pollen deposition split by field distance. Possibly a pattern??
# 
# ggplot(flowers,aes(DistFac,Pollen))+geom_violin() #Same as above, but with violin plots.

# Pollen deposition models ------------------------------------------------
temp=filter(flowers,!is.na(Pollen))
pol_a=glm.nb(Pollen~Distance,data=flowers) #NO RELATIONSHIP WITH DISTANCE OR STOCKING

pr=list(R=list(V=1,nu=0.02),
      G=list(G1=list(V = diag(2), nu = 1,alpha.mu=c(0,0),alpha.V=diag(2)*2000)))
pol_b=MCMCglmm(Pollen~Distance*NumHives,random=~idh(1+Distance):Field,data=temp,family='poisson',nitt=30000,thin=30,prior=pr) #Slope, BeeYard and Slope*BeeYard don't appear important

pred=mutate(cbind(temp[,c('Distance','NumHives')],
                  predict(pol_b,interval='confidence')))%>% #Marginalizing across random effects
  filter(Distance!=325) %>%
  unique() 
p=ggplot(pred,aes(x=Distance,y=fit,ymax=upr,ymin=lwr,colour=as.factor(NumHives),group=NumHives))+ #Pollen with distance and stocking
  geom_pointrange(position=position_dodge(width=5),size=1)+
  geom_smooth(method='glm',formula=round(y)~x,family='poisson',se=F,size=1,linetype='dashed')+
  ylab('Pollen Grains per Stigma')+labs(colour='Stocking Rate')+  
  scale_colour_manual(values=c('green','orange','orange','orange','red'))
ggsave(paste(folder,'pollenStocking.png',sep='\\'),p,width=9,height=6)

# Visitation and Nectar Plots ---------------------------------------------------
ylab=expression(paste('Nectar volume (',mu,'L)',sep=''))

ggplot(flowers,aes(60*Honeybee/Time,Sugar,colour=as.factor(NumHives)))+geom_point()+
  geom_smooth(method='lm',se=F)+labs(x='Honeybee Visits/hr',y=ylab)
  #Visitation vs Nectar volume, with a square root transformed linear model 

ggplot(flowers,aes(FlyRate,Sugar,colour=Field))+geom_point()+geom_smooth(method='lm',se=F)+labs(x='Fly Visits/hr',y=ylab) #Fly Visitation vs Nectar volume 

ggplot(flowers,aes(HoverflyRate,Sugar,colour=Field))+geom_point()+geom_smooth(method='lm',se=F)+labs(x='Hoverfly Visits/hr',y=ylab) #Hoverfly visitation vs Nectar volume 

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
ggsave(paste(folder,'nectarVolBees.png',sep='\\'),p,width=9,height=6)

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
ggsave(paste(folder,'seedsperpodDistance.png',sep='\\'),p,width=9,height=6)

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
ggsave(paste(folder,'podmassDistance.png',sep='\\'),p,width=9,height=6)

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

# Flower success models and plots (plant-level)--------------------------------------------
# NO EFFECT OF DISTANCE, NUMHIVES, IRRIGATION, OR VEGMASS. SEEDMASS APPEARS TO BE SOMEWHAT IMPORTANT, BUT THIS IS PROBABLY THE WRONG DIRECTION OF CORRELATION (I.E. SEEDMASS~PODsUCCESS)
temp=select(plants,Field,Plot,Plant,Distance,PlDens,Irrigated,
            BeeYard,NumHives,VegMass,SeedMass,Pods,Missing) %>%
  mutate(Pods=round(Pods),TotalFls=round(Pods+Missing),FieldPlot=paste(Field,Plot,sep='.')) %>%
  na.omit()

ggplot(temp,aes(Distance,Pods/TotalFls,group=Field,colour=BeeYard))+geom_point()+geom_smooth(method='glm',family='binomial',se=F) #Possible decrease?
pr=list(R=list(V=1,nu=0.02),
        G=list(G1=list(V = diag(1), nu = 1),
               G2=list(V = diag(1), nu = 1))) #Priors
flWins=MCMCglmm(cbind(Pods,Missing)~Irrigated,random=~Field+FieldPlot,family='`2',
                data=temp,prior=pr) 
summary(flWins) #No distance effect

predict(flWins)


HPDinterval(flWins$VCV) #Variance components: Field>Plot>units
HPDinterval(flWins$VCV[, 1]/rowSums(flWins$VCV)) #21-78% of var at Field level
HPDinterval(flWins$VCV[, 2]/rowSums(flWins$VCV)) #8-46% of var at Plot level
HPDinterval(flWins$VCV[, 3]/rowSums(flWins$VCV)) #8-23% of var at unit level


#Flower number ~ Distance ?
ggplot(temp,aes(Distance,TotalFls,group=Field,colour=BeeYard))+geom_point()+geom_smooth(method='glm',family='poisson',se=F)

#Pods ~ Distance ?
ggplot(temp,aes(Distance,Pods,group=Field,colour=BeeYard))+geom_point()+geom_smooth(method='glm',family='poisson',se=F)

# Plant-level to pod-level comparison -------------------------------------
# Do we really need pod-level measurements, or are plant-level ones OK?

temp=filter(plants,!is.na(SeedCount),!is.na(Pods),!is.na(AvPodCount),!is.na(AvPodMass))

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
