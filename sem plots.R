# MAKES SEM PLOTS FOR COMMODITY AND SEED MODELS

# Load everything ---------------------------------------------------------

library(tidyverse)
library(beepr)
library(xtable)
library(ggdag)
library(ggpubr)

#Big-text theme, no grid lines (used for Bayer 2016 presentation)
prestheme <- theme(legend.position='right',
                legend.text=element_text(size=15),
                axis.text=element_text(size=15),
                axis.title=element_text(size=20),
                title=element_text(size=20),
                panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                panel.border=element_rect(size=1,colour='black'),
                strip.text=element_text(size=15))
theme_set(theme_bw()+prestheme) #Sets graph theme to B/Ws + prestheme
rm(prestheme)

# setwd('~/Projects/UofC/canola_yield_project') #Multivac path
setwd('~/Documents/canola_yield_project') #Galpern machine path

source('helperFunctions.R')

load('./Models/datalist_commodity.Rdata') #Commodity data
commData <- datalist
load('./Models/datalist_seed.Rdata') #Seed data
seedData <- datalist; rm(datalist)

#Model summaries
load('./Models/modSummaries_commodity.Rdata') 
load('./Models/modSummaries_seed.Rdata') 

# Commodity SEM plot ----------------------------------------------------------

nodeCoords <- data.frame(name=c('numHives','hbeeDist','hbeeVis','pollen',
                                'plSize','plDens','flDens',
                                'flwCount','flwSurv','seedCount','seedWeight'),
                         labs=c('Number\nof Hives','Field Edge\nDistance','Honey bee\nVisits','Pollen\nCount',
                                'Plant\nSize','Plant\nDensity','Flower\nDensity',
                                'Flowers\nper Plant','Pod Set\n(%)','Seeds\nper Pod','Seed\nSize'),
                         x=c(0,1,0.5,0.5,
                             2.5,0,1,
                             3.5,4,2.5,3.5),
                         y=c(4,4,3,2,
                             0,1,1,
                             0,1,3,2))


#Specify model
commDAG <- dagify(plDens ~ hbeeDist,
                  plSize ~ plDens + hbeeDist,
                  flDens ~ plSize + hbeeDist + plDens,
                  hbeeVis ~ hbeeDist + numHives + flDens,
                  pollen ~ hbeeVis + hbeeDist,
                  flwCount ~ plSize + flwSurv,
                  flwSurv ~ hbeeVis + plSize + pollen,
                  seedCount ~ hbeeVis + pollen + plSize + flwSurv + flwCount,
                  seedWeight ~ hbeeVis + pollen + plSize + seedCount + plDens + hbeeDist + flwSurv + flwCount + flDens,
                  coords= list(x = setNames(nodeCoords$x,nodeCoords$name),
                               y = setNames(nodeCoords$y,nodeCoords$name)),
                  labels=setNames(nodeCoords$labs,nodeCoords$name)
)

#Get path coefficients (Z-scores) and match to edges
pathCoefs <- lapply(modSummaries_commodity,function(x) x$summary) %>%  #Get path coefficients
  bind_rows(.id='to') %>% 
  transmute(to,name=param,Z,pval) %>% 
  filter(to!='yield',!grepl('(int|sigma|lambda|Phi|phi|rho)',name)) %>% 
  mutate(to=case_when(
    to=='flDens' & grepl('PlDens$',name) ~ 'plDens',
    to=='flDens' & grepl('PlSize$',name) ~ 'plSize',
    TRUE ~ gsub('avg','seed',to)
  )) %>% 
  mutate(name=gsub('slope','',name)) %>% 
  filter(mapply(grepl,capFirst(to),name)) %>% 
  mutate(name=mapply(gsub,capFirst(to),'',name)) %>% 
  mutate(name=capFirst(name,TRUE)) 

commDAG <- commDAG %>% #Create tidy dagitty set
  tidy_dagitty() 

# ggplot(commDAG,aes(x = x, y = y, xend = xend, yend = yend))+
#   geom_dag_edges() +
#   geom_dag_text(col='black') +
#   theme_dag_blank()

commDAG$data <- commDAG$data %>% #Match coefs to dagitty set
  rename(xstart=x,ystart=y) %>% 
  full_join(x=pathCoefs,y=.,by=c('to','name')) %>% 
  mutate(isNeg=Z<0,isSig=pval<0.05) %>% 
  mutate(L=sqrt(abs(Z)),1) %>% 
  # mutate(edgeLab=ifelse(isSig,as.character(sign(Z)*round(L,1)),'')) %>%
  mutate(edgeLab=ifelse(isSig,as.character(formatC(Z,format='f',digits=1)),'')) %>%
  mutate(C=ifelse(isNeg,'red','black')) %>% 
  mutate(A=ifelse(isSig,1,0.1))

(commSEM <- ggplot(commDAG$data,aes(x = xstart, y = ystart, xend = xend, yend = yend))+
  # geom_dag_point(size=20,colour='black',shape='square')+
  # geom_dag_text(aes(label=label),col='white') +
  annotate('text',x=0.5,y=4.5+0.1,label='Plot Level',size=5)+
  annotate('rect',xmin=-0.5,ymin=0.5,xmax=1.5,ymax=4.5,fill=NA,col='black',
           linetype='dashed',linejoin='round',size=1)+
  
  annotate('text',x=3.5,y=3.5+0.1,label='Plant Level',size=5)+
  annotate('rect',xmin=2,ymin=-0.5,xmax=4.5,ymax=3.5,fill=NA,col='black',
           linetype='dashed',linejoin='round',size=1)+
  annotate('label',x=3.25,y=4.25,label='Commodity fields',size=10)+
  geom_dag_edges(aes(edge_width=L,edge_colour=C,edge_alpha=A),
                 arrow_directed=arrow(angle=20,type='open')) +
  geom_dag_edges(aes(label=edgeLab),label_pos=0.45,edge_alpha=0,label_size=6,fontface='bold',label_colour ='black') +
  geom_dag_edges(aes(label=edgeLab),label_pos=0.45,edge_alpha=0,label_size=6,fontface='plain',label_colour='white') +
  geom_dag_label_repel(aes(label=label),col='black',force=0) +
  theme_dag_blank())

# Seed SEM plot -----------------------------------------------------------

nodeCoords <- data.frame(name=c('hbeeDist','lbeeDist','hbeeVis','lbeeVis','pollen',
                                'plSize','plDens','flDens','cent',
                                'flwCount','flwSurv','seedCount','seedWeight'),
                         labs=c('Field Edge\nDistance','Shelter\nDistance',
                                'Honey bee\nVisits','Leafcutter\nVisits','Pollen\nCount',
                                'Plant\nSize','Plant\nDensity','Flower\nDensity','Bay Centre',
                                'Flowers\nper Plant','Pod Set\n(%)','Seeds\nper Pod','Seed\nSize'),
                         x=c(0,1,0,1,0.5,2.5,0,1,0.5,3.5,4,2.5,3.5),
                         y=c(4,4,3,3,2,0,1,1,0.5,0,1,3,2))


#Specify model
seedDAG <- dagify(plDens ~ hbeeDist,
                  plSize ~ plDens + hbeeDist,
                  flDens ~ hbeeDist,
                  hbeeVis ~ flDens + hbeeDist + lbeeDist + cent,
                  lbeeVis ~ flDens + hbeeDist + lbeeDist + cent,
                  pollen ~ hbeeVis  + lbeeVis + cent + hbeeDist + lbeeDist + flDens,
                  flwCount ~ plSize + cent + flwSurv,
                  flwSurv ~ pollen + plSize + cent + hbeeDist + lbeeDist + flDens,
                  seedCount ~ pollen + plSize + cent + hbeeDist + flDens + flwSurv + flwCount,
                  seedWeight ~ pollen + seedCount + plSize + plDens + lbeeDist,
                  coords= list(x = setNames(nodeCoords$x,nodeCoords$name),
                               y = setNames(nodeCoords$y,nodeCoords$name)),
                  labels=setNames(nodeCoords$labs,nodeCoords$name)
)

pathCoefs <- lapply(modSummaries_seed,function(x) x$summary) %>% #Get path coefficients
  bind_rows(.id='to') %>% 
  transmute(to,name=param,Z,pval) %>% 
  filter(to!='yield',!grepl('(int|sigma|lambda|Phi|phi|rho|nu|theta)',name)) %>% 
  mutate(to=case_when(
    to=='flDens' & grepl('PlDens$',name) ~ 'plDens',
    to=='flDens' & grepl('PlSize$',name) ~ 'plSize',
    TRUE ~ gsub('avg','seed',to)
  )) %>% 
  mutate(name=gsub('slope','',name)) %>% 
  filter(mapply(grepl,capFirst(to),name)) %>% 
  mutate(name=mapply(gsub,capFirst(to),'',name)) %>% 
  mutate(name=capFirst(name,TRUE)) 

tidySeedDAG <- seedDAG %>% #Create tidy dagitty set
  tidy_dagitty() 

tidySeedDAG$data <- tidySeedDAG$data %>% #Match coefs to dagitty set
  rename(xstart=x,ystart=y) %>% 
  full_join(x=pathCoefs,y=.,by=c('to','name')) %>% 
  # filter(!is.na(to)) %>% 
  mutate(isNeg=Z<0,isSig=pval<0.05) %>% 
  mutate(L=sqrt(abs(Z)),1) %>% 
  # mutate(edgeLab=ifelse(isSig,as.character(sign(Z)*round(L,1)),'')) %>%
  mutate(edgeLab=ifelse(isSig,as.character(formatC(Z,format='f',digits=1)),'')) %>%
  mutate(C=ifelse(isNeg,'red','black')) %>% 
  mutate(A=ifelse(isSig,1,0.1))

(seedSEM <- ggplot(tidySeedDAG$data,aes(x = xstart, y = ystart, xend = xend, yend = yend))+
  annotate('text',x=0.5,y=4.5+0.1,label='Plot Level',size=5)+
  annotate('rect',xmin=-0.5,ymin=0,xmax=1.5,ymax=4.5,fill=NA,col='black',
           linetype='dashed',linejoin='round',size=1)+
  annotate('text',x=3.5,y=3.5+0.1,label='Plant Level',size=5)+
  annotate('rect',xmin=2,ymin=-0.5,xmax=4.5,ymax=3.5,fill=NA,col='black',
           linetype='dashed',linejoin='round',size=1)+
  annotate('label',x=3.25,y=4.25,label='Seed fields',size=10)+  
  geom_dag_edges(aes(edge_width=L,edge_colour=C,edge_alpha=A),
                 arrow_directed=arrow(angle=20,type='open')) +
  geom_dag_edges(aes(label=edgeLab),label_pos=0.45,edge_alpha=0,label_size=6,fontface='bold',label_colour ='black') +
  geom_dag_edges(aes(label=edgeLab),label_pos=0.45,edge_alpha=0,label_size=6,fontface='plain',label_colour='white') +
  geom_dag_label_repel(aes(label=label),col='black',force=0) +
  theme_dag_blank()
)

# Combine plots -----------------------------------------------------------

(p <- ggarrange(commSEM,seedSEM,ncol=1,nrow=2))
ggsave('./Figures/allSEM.png',p,height = 14,width=12,bg='white')


