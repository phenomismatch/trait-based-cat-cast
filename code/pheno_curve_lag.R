###Macrosystems: Trait-based caterpillar cast
###Identifying species-specific lags between curves

##libraries
library(tidyverse) #tidy
library(readxl) #for importing excel data
library(rbms) # for gam models
library(lubridate) # for dates
library(MuMIn) # for r-squared estimation

#import other code
source("code/gam_fx.R")

#############################################################################
## Scope & Data Imports
#Splist is currently in species codes but will probably generally be scientific names
splist<-c("AMCO","AMLA","BLSW","CAWH","CORI","CWNY","EACO","ETBL","GSFR","MOCL","MONA","PECR","PISW","READ","SPAZ","SSSK","ZASK")
regions<-read_excel("data/regions.xlsx")
year.scope<-c(2016:2019)
doy.scope<-c(1:366)


##Adult data import
adult.data<-read.csv("data/ma-bfly.csv")
adult.data$Date<-as.Date(adult.data$Date, format="%Y-%m-%d")
adult.data<-adult.data %>%
  select(SITE_ID, DATE = Date, Year, sp = Code, English.Name, count = abund) %>%
  arrange(SITE_ID, DATE)
adult.format<-"abundance"


## Larval data import & formatting
larvae<-read.csv("data/caterpillars-of-eastern-north-america.csv")
larvae$observed_on<-as.Date(larvae$observed_on)
larval.data<-larvae %>%
  filter(between(latitude, regions[3,2],regions[3,3]) & between(longitude, regions[3,4],regions[3,5])) %>%
  select(scientific_name, common_name, observed_on, latitude, longitude, quality_grade, num_identification_agreements, num_identification_disagreements) %>%
  mutate(year=year(observed_on), doy=yday(observed_on)) 
larval.format<-"occurrence"


################################################################
##  ESTIMATE LAGS BETWEEN LARVAL AND ADULT PHENOLOGY
#Create adult phenology curve
phen.lag<-list()
phen.lag.sp<-NULL

## Loop through species list 
for(spi in splist) {

  ## Get ADULT PHENOLOGY
  #get adult data for species, including 0s for days not observed
  if(adult.format=="abundance") { 
    adult.spi<-adult.data %>%
      mutate(CT = ifelse(sp==spi,count,0)) %>%
      group_by(SITE_ID, DATE) %>%
      arrange(DATE) %>%
      summarize(SPECIES=spi,COUNT=sum(CT))
  
    #for MA data this gets the species common name
    spname<-adult.data$English.Name[adult.data$sp==spi][1]
  
  ################################################################
    #GAM estimate of species flight curves for years selected in scope
    adult.spi.phen<-list()
    for(i in 1:length(year.scope) ){
      adult.spi.phen[[i]]<-flightcurve(filter(adult.spi, year(DATE)==year.scope[i]))
    }
  } #end adult GAM
  
  #######################################################################
  for(i in 1:length(year.scope) ){
    
    # larval phenology
    larval.i<-filter(larval.data,common_name==spname, year==year.scope[i])
    ## Larval data is currently occurrence data
    if(larval.format=="occurrence") {
      #verify at least 3 larval records, then: calculate smoothed data density 
      if(nrow(larval.i)>2) {
        temp.larval<-density(larval.i$doy,n = length(adult.spi.phen[[i]]), from=1, to=length(adult.spi.phen[[i]]))      
        temp.larval<-data.frame(doy=doy.scope[1:length(adult.spi.phen[[i]])], larv=round(temp.larval$y[1:length(adult.spi.phen[[i]])]*100,3))
      }
    }
    
    ## Estimate Lag
        temp.ad<-data.frame(doy=doy.scope[1:length(adult.spi.phen[[i]])], adult=adult.spi.phen[[i]])
        lagresult<-tibble(region=character(0),species=character(0),year=numeric(0),lag=numeric(0),r2=numeric(0))
        lagr2<-NULL
        #loop lags
        lag.range<-c(-60:60)
        for(lagi in lag.range) {
          #larvae before adults
          temp.lalag<-temp.larval %>% add_row(doy=(1-lagi),larv=0) %>% mutate(doy=doy+lagi)
          temp.datalag<-merge(temp.ad,temp.lalag, by=intersect(names(temp.ad), names(temp.lalag)))
          lag.lm<-lm(larv~adult, data=temp.datalag)
          lagr2<-c(lagr2,round(summary(lag.lm)[[9]],2))
          } #close lag loop
        bestlags<-list(lag.range[which(lagr2==max(lagr2))])[[1]]
        lagresult[c(1:length(bestlags))+nrow(lagresult),]<-tibble(region="MA",species=spname,year=year.scope[i],lag=bestlags,r2=max(lagr2))
      } #close years
  } #close species
    
    ### formatted to here
    ##rsq_compar[nrow(rsq_compar)+1,]<-c(spnames[spi],years[i],lagi,round(summary(temp.lm)[[9]],2))
    ##lag_result<-lag result %>% arrange(-r.squared)  

  
  
  ############################################################

 ################   END 11-2-2020 ############################

