#NABA data for phenometrics & species composition

library(tidyverse)
library(RODBC)
library(sp)
library(sf)
library(maptools)

#Hex cells
hex_sf <- read_sf("data/maps/hex_grid_crop.shp")
hex_psf<-st_transform(hex_sf, 3857) #3857 is projected pseudo-mercator, 4326 is lat/long

latlongproj<-("+proj=longlat +datum=WGS84 +no_defs")
myproj<-("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

hex1<-st_as_sf(hex_sf)
#hex_sp<-st_read("data/maps/hex_grid_crop.shp")
hex.proj<-spTransform(hex1, proj4string(spatial.proj))







#NABA data - raw data inputs
con1 <- odbcConnect("NABA")
sqlTables(con1)

sightings2010.13<-sqlFetch(con1, "NABA_OBS_FlatFile") %>%
  filter(Year<2014, between(Latitude,20,60), between(Longitude,-100,-60)) %>%
  select(Year, Month, Day, SurveyDate, Latitude, Longitude, UMD_Species_Code, Total_Count)

con2 <- odbcConnect("2019BFLY")
sqlTables(con2)

data4mike<-sqlFetch(con2, "Simpleton NABA raw") %>%
  filter(ObsYear>2013)

#Data spatial join to hex cells & formatting
spat.naba<-SpatialPointsDataFrame(cbind(sightings2010.13$Longitude, sightings2010.13$Latitude), sightings2010.13, coords.nrs = numeric(0), 
                                  proj4string = CRS(as.character(latlongproj)),  bbox = NULL)
spat.naba<-st_as_sf(spat.naba, coords = c("Longitude", "Latitude"), crs = latlongproj, agr = "constant")
hex.data<-st_read('C:/Users/eal109/Documents/Bird_Phenology/Data/hex_grid_crop/hex_grid_crop_data.shp')
spat.proj<-st_transform(spat.naba,crs=myproj)
spat.hex<-st_transform(hex.data,crs=myproj)

data.hex<-st_within(spat.proj, spat.hex)
spat.hex_cells <- as.numeric(as.character(data.hex))
hex1<-spat.hex$cell[spat.hex_cells]
sightings<-bind_cols(sightings2010.13, hex1) 

points(x=sightings$Longitude, y=sightings$Latitude, type="p", col="red")





t3<-sqlFetch(con2, "Simpleton species codes") %>%
  dplyr::select(`UMD-CODE`,Family,Genus, Species, Code)

try1<-merge(sightings, t3, by.x="UMD_Species_Code", by.y="UMD-CODE", all.x=T)%>%
  filter(!is.na(Total_Count & Total_Count>0)
  
try1a<-try1  %>%
  filter(Code %in% c("RE","RL","RP","Oth")) %>%
  group_by(hex,Year,Month,Day,Code) %>%
  summarize(ntot=sum(Total_Count, na.rm=T))





try2<-merge(data4mike, t3, by="UMD-CODE", all.x=T) %>% 
  filter(!is.na(SumOfNumSeen1 & SumOfNumSeen1>0)
try2a<-try2  %>%
  filter(Code %in% c("RE","RL","RP","Oth"), TempEventType%in%c("Opportunistic","Trip"), !is.na(HEXcell))  %>%
  group_by(HEXcell, ObsYear,ObsMonth,ObsDay,Code) %>%
  summarize(ntot=sum(SumOfNumSeen1))
names(try2a)<-names(try1a)

naba_opp_data<-bind_rows(try1a, try2a) %>% filter(Year>1999)

table(naba_opp_data$Year, naba_opp_data$Code)

write.csv(naba_opp_data, file="data/derived_data/naba.data.pheno.csv")



### Family Comp
try1b<-try1  %>%
  select(Family, hex, Year, ntot=Total_Count) 
try2b<-try2 %>% select(Family, hex=HEXcell, Year=ObsYear, ntot=SumOfNumSeen1) %>% mutate(ntot=ifelse(ntot==0,1,ntot))
naba_comp_data<-bind_rows(try1b, try2b) %>% filter(Year>1999) %>%
  filter(!is.na(hex) & !is.na(Family) & !is.na(hex) & !is.na(ntot)) %>%
  group_by(Year, hex, Family) %>%
  summarize(tot_count=sum(ntot, na.rm=T),n.recs=n())
  

summary(naba_comp_data)

write.csv(naba_comp_data, file="data/derived_data/naba.taxonomy.csv")




dup<-try1 %>%  #this works
  group_by(HEXcell, ObsYear,ObsMonth,ObsDay,Code) %>%
  summarize(ntot=sum(SumOfNumSeen1), bph=sum(SumOfNumSeen1/SumOfPartyHours))



justopp<-try1 %>%  #this works
  filter(TempEventType=="Opportunistic") %>%
  group_by(HEXcell, ObsYear,ObsMonth,ObsDay,Code) %>%
  summarize(ntot=sum(SumOfNumSeen1))

justft<-try1 %>%  #this works
  filter(TempEventType=="Trip") %>%
  group_by(HEXcell, ObsYear,ObsMonth,ObsDay,Code) %>%
  summarize(ntot=sum(SumOfNumSeen1))

write.csv(justopp, file="NABA.opp.data.csv")
write.csv(justft, file="NABA.trip.data.csv")

spcomp<-try1 %>%  
  group_by(HEXcell, ObsYear,Code, Family) %>%
  summarize(nobs=n(),ntot=sum(SumOfNumSeen1))

spcomp.opp<-try1 %>%  
  filter(TempEventType=="Opportunistic") %>%
  group_by(HEXcell, ObsYear,Code, Family) %>%
  summarize(nobs=n(),ntot=sum(SumOfNumSeen1))

spcomp.trip<-try1 %>%  
  filter(TempEventType=="Trip") %>%
  group_by(HEXcell, ObsYear,Code, Family) %>%
  summarize(nobs=n(),ntot=sum(SumOfNumSeen1))

write.csv(spcomp.opp, file="NABA.composition.oppdata.csv")
write.csv(spcomp.trip, file="NABA.composition.tripdata.csv")

write.csv(spcomp, file="NABA.species.composition.csv")

