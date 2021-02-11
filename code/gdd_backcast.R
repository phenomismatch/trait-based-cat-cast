###Macrosystems: Trait-based caterpillar cast
###Estimating caterpillar phenology from adults
#February 2021, Elise Larsen
#Code for GDD model of back-casting caterpillar phenology PDF from adult phenology PDF
#In pheno.x, doyl and doyp rescale the adult pdf along the doy axis using GDD

##libraries
library(tidyverse) # tidy
library(lubridate) # for dates
library(ggplot) # for visualizations
#############################################################################
## Scope & Data Imports

#set GDD model quantities
gddP<-154
gddL<-330

##Set parameters
my.hexes<-c(595, 676)
my.years<-c(2017)

##Import Formatted GDD DATA
#Import (these day only include DOY 60-250. Elise has asked Naresh for a version with the whole year.)
load("data/gdd_data2012-2017.RData")

## Import adult Lep phenology pdf data
phenoraw<-read_csv("data/simpleton_pheno_pdfs.csv")
#Filter to RL data in target hexes & years
phenoDF<-filter(phenoraw, HEXcell%in%my.hexes, year%in%my.years, code=="RL")
names(phenoDF)<-c("DOY","pdf","year","hex","code")

##MERGE TABLES
pheno.x<-merge(phenoDF, gddDF,by=intersect(names(phenoDF),names(gddDF)))
pheno.x<-pheno.x %>%
  mutate(gddp=accumGDD-gddP, doyp=0, gddl=accumGDD-gddL, doyl=0)

for(i in 1:nrow(pheno.x)) {
  temp<-filter(gddDF, year==pheno.x$year[i], hex==pheno.x$hex[i])
  #Because we don't have DOY0-59 GDD and we have negative gddp and gddl values:
  day1<-min(pheno.x$DOY[pheno.x$year==pheno.x$year[i] & pheno.x$hex==pheno.x$hex[i] & pheno.x$gddl>0])
  if(pheno.x$gddl[i]<1) { pheno.x$doyl[i]<- floor(pheno.x$DOY[i]/day1*60)
  } else {   pheno.x$doyl[i]<-temp$DOY[which(temp$accumGDD>pheno.x$gddl[i])[1]] }
  if(pheno.x$gddp[i]<1) { pheno.x$doyp[i]<- floor(pheno.x$DOY[i]/day1*60)
  } else { pheno.x$doyp[i]<-temp$DOY[which(temp$accumGDD>pheno.x$gddp[i])[1]] }
}

summary(pheno.x)
####In pheno.x, doyl and doyp rescale the adult pdf along the doy axis using GDD

######################################################
#Visualizations

#plot results
figPhenoX<-ggplot(data=pheno.x, aes(x=doyl, y=pdf, group=interaction(year, hex), color=as.factor(hex))) + 
  geom_line() + xlim(c(0,330)) + theme_classic() + theme(legend.position='none',axis.text.y = element_blank()) + 
  scale_color_manual(values=fig4acolors) + labs(x="DOY", y="") + facet_wrap(~year) + ggtitle("Caterpillar biomass estimates")
figPhenoX


#Plot gdd data
figcol1<-c("purple","blue","darkgreen")
figGddHex<-ggplot(data=filter(gddDF, hex %in% my.hexes), aes(x=DOY,y=accumGDD, group=interaction(as.factor(year),hex), color=hex)) + 
  geom_line() + theme_classic() + theme(legend.position='none') + 
  scale_color_manual(values=figcol1) + 
  labs(x="Ordinal day of year (DOY)", y="Accum. growing degree days") + 
  annotate(geom="text", x=50, y=c(1600,1200), label=c("Chicago", "Boston"), color=figcol1[c(2,1)], hjust=0, size=5) + 
  ggtitle("GDD accumulation 2012-2017")
figGddHex

figGddYears<-ggplot(data=filter(gddDF, hex==618, year%in%c(2012,2017)), aes(x=DOY,y=accumGDD, group=interaction(as.factor(year),hex), color=as.factor(year))) + 
  geom_line(size=2) + theme_classic() + theme(legend.position='none') + 
  scale_color_manual(values=c("seagreen4","darkslategray")) + xlim(100,250) +
  labs(x="DOY", y="Accum. growing degree days") + 
  annotate(geom="text", x=100, y=c(1600,1400), label=c("Warm Year","Cold Year"), color=c("seagreen4","darkslategray"), hjust=0, size=5)
figGddYears

figAdultPhenoHexes<-ggplot(data=filter(phenoDF, HEXcell%in%my.hexes, year%in%c(2015:2017), code=="RL"), aes(x=x, y=y, group=interaction(year, HEXcell), color=as.factor(HEXcell))) + 
  geom_line() + theme_classic() + theme(legend.position='none',axis.text.y = element_blank()) + 
  scale_color_manual(values=fig4acolors) + labs(x="DOY", y="") + facet_wrap(~year) + ggtitle("Adult phenology estimates")
figAdultPhenoHexes

figAdultPhenoYears<-ggplot(data=filter(phenoDF, HEXcell==595, year%in%c(2012,2017), code=="RL"), aes(x=x, y=y, group=interaction(year, HEXcell), color=as.factor(year))) + 
  geom_line(size=2) + theme_classic() + theme(legend.position='none',axis.text.y = element_blank()) + coord_flip()   + 
  scale_color_manual(values=c("darkslategray","seagreen4")) + labs(x="DOY", y="")
figAdultPhenoYears


##Michael's output:
phenoraw<-read_csv("adult_data/simpleton_pheno_pdfs.csv")
#add index to filter data

plotPhen1<- function(
  data) { #data tibble
  title<-data$HEXcell[1]
  t1<-ggplot(data=data, aes(x=x, y=y,color=code)) + 
    geom_line() + xlim(60,330) + 
    labs(x="DOY",y="Density",title=paste(title)) + theme_classic()
  return(t1)
}

for(i in my.hexes)  {
  print(plotPhen1(filter(phenoraw, HEXcell==i, year %in% my.years)) + 
          facet_wrap(~year))
}
