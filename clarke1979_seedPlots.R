library(ggplot2)
theme_set(theme_classic())

#Data from Clarke 1979
clarke1 <- data.frame(seeds=c(24,28,24,23,18),weight=c(3.3,3.45,3.26,3.23,3.14),year='1976')
clarke2 <- data.frame(seeds=c(28,24,23,22),weight=c(4.27,3.81,3.59,3.39),year='1977')
clarke <- rbind(clarke1,clarke2)
rm(clarke1,clarke2)
ggplot(clarke,aes(seeds,weight,col=year))+geom_point()+
  labs(x='Seeds per pod',y='Weight per seed (mg)')+
  geom_smooth(method='lm',se=F)

#Data from Angadi et al 2003
angadi <- data.frame(seeds=c(22.3,23.0,23.2,24.1,25.4,24.8,23.9,22.9),weight=c(2.7,2.95,3.06,2.89,3.08,3.07,3.13,3.05),year='1999')
angadi2 <- data.frame(seeds=c(22.4,23.7,24.7,25.8,21.8,16.0,25.1),weight=c(3.37,2.88,3.03,3.06,2.91,3.58,3.41),year='2000')
angadi3 <- data.frame(seeds=c(23.1,20.1,23.6,23.5,20.1,19.8,21.3,19.8),weight=c(3.04,2.77,2.75,3,2.71,2.96,3.08,2.87),year='2001')
angadi <- rbind(angadi,angadi2,angadi3)
rm(angadi,angadi2,angadi3)

ggplot(angadi,aes(seeds,weight,col=year))+geom_point()




