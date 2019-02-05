library(ggplot2)
library(tidyr)
library(dplyr)
theme_set(theme_classic())

#Data from Clarke 1979
clarke1 <- data.frame(seeds=c(24,28,24,23,18),weight=c(3.3,3.45,3.26,3.23,3.14),
                      type=c('Main raceme (top)','Main raceme (bottom)','Branch 1','Branch 2','Branch 3'),
                      year='1976')
clarke2 <- data.frame(seeds=c(28,24,23,22),weight=c(4.27,3.81,3.59,3.39),
                      type=c('Main raceme (all)','Branch 1','Branch 2','Branch 3'),
                      year='1977')
clarke <- rbind(clarke1,clarke2)
rm(clarke1,clarke2)
ggplot(clarke,aes(seeds,weight,col=year))+
  geom_smooth(method='lm',se=F,size=.2,alpha=0.3)+
  geom_point()+
  geom_text(aes(x=seeds,y=weight,label=type),col='black',size=3,nudge_y=0.02,hjust='inward')+
  labs(x='Seeds per pod',y='Weight per seed (mg)',col='Year',title='Data from Clarke 1979, Table 1')+
  # xlim(17.5,30.5)+
  scale_colour_manual(values=c('red','blue'))

clarke2 <- data.frame(PlantPart=c('Main raceme','Branch 1','Branch 2','Branch 3','Main raceme','Branch 1','Branch 2','Branch 3'),
                      Meas=c(rep('Seeds',4),rep('AbSeeds',4)),
                      Aug3=c(29.5,28.5,26.8,30,1.2,1.3,3.3,1.0),Aug09=c(26.8,24.2,24.5,24.7,4.5,5.3,5.7,6.3),
                      Aug16=c(29,27,25.3,23.8,2.8,3.8,5.8,7),Aug23=c(28.5,25.5,24,27.2,3.3,5.2,7,3.8),
                      Aug31=c(28.8,27.3,25.5,21.8,4.7,4.7,5.5,10),Sep16=c(28,23.8,23,21.5,4,7.7,7.7,10.3))
clarke2 %>% gather('day','count',Aug3:Sep16) %>% spread(Meas,count) %>% 
  gather('Type','Count',AbSeeds:Seeds) %>% 
  mutate(day=as.Date(day,format='%b%d')) %>% 
  mutate(Type=factor(Type,levels=c('Seeds','AbSeeds'),labels=c('Seeds per pod','Aborted seeds per pod'))) %>% 
  ggplot(aes(day,Count,col=PlantPart))+geom_line()+geom_point()+
  facet_wrap(~Type,ncol=1,scales='free_y')+
  scale_colour_manual(values=c('red','blue','purple','darkorange'))+
  labs(x='Day of year',y='Count',col='Plant\nPart',title='Data from Clarke 1979, Table 2')
  




#Data from Angadi et al 2003
angadi <- data.frame(seeds=c(22.3,23.0,23.2,24.1,25.4,24.8,23.9,22.9),weight=c(2.7,2.95,3.06,2.89,3.08,3.07,3.13,3.05),year='1999')
angadi2 <- data.frame(seeds=c(22.4,23.7,24.7,25.8,21.8,16.0,25.1),weight=c(3.37,2.88,3.03,3.06,2.91,3.58,3.41),year='2000')
angadi3 <- data.frame(seeds=c(23.1,20.1,23.6,23.5,20.1,19.8,21.3,19.8),weight=c(3.04,2.77,2.75,3,2.71,2.96,3.08,2.87),year='2001')
angadi <- rbind(angadi,angadi2,angadi3)
rm(angadi,angadi2,angadi3)

ggplot(angadi,aes(seeds,weight,col=year))+geom_point()




