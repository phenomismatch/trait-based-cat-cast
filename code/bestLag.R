# Function for calculating best fit lag between adult and larval phenology

# Need to allow for either abundance (e.g. MA Butterfly Survey) or occurrence (e.g. iNaturalist) data for adults,
# which will determine whether a GAM or kernel density estimator is fit to the data, respectively.

# For caterpillars, also need to allow for abundance (e.g. Hubbard Brook, Coweeta, Caterpillars Count!) or occurrence
# (e.g. iNaturalist) data, but in either case, kernel density estimator will be used.



# Temporary code for creating dataframes for testing:
# libraries
library(dplyr)
library(lubridate)
library(readxl)
library(rbms) # for gam models
library(MuMIn) # for r-squared estimation

#import other code
source("code/gam_fx.R")

# Until scientific names are merged in, plan on using Monarch, Danaus plexippus for testing
adultDF = read.csv('data/ma-bfly.csv') %>%
  mutate(DATE = as.Date(Date, format = "%m/%d/%Y"),
         scientific_name = English.Name) %>%
  select(SITE_ID, DATE, Year, sp = Code, scientific_name, count = abund) %>%
  arrange(SITE_ID, DATE)

adultDF$scientific_name = gsub("Monarch", "Danaus plexippus", adultDF$scientific_name)

mothDF = read.csv('data/discoverlife_3sp_ma.csv', header = T) %>%
  mutate(scientific_name = gsub("\\.", " ", species),
         DATE = as.Date(paste(year, month, day, sep = "-"), format = "%Y-%m-%d")) %>%
  select(scientific_name, year, DATE, count) %>%
  filter(count >= 0) %>% # filters out records with count = -1 where no sampling was conducted
  mutate(scientific_name = gsub("Halysidota harrisii--tessellaris", "Halysidota harrisii", scientific_name)) # fixing test species

regions<-read_excel("data/regions.xlsx")

larvalDF = read.csv("data/caterpillars-of-eastern-north-america.csv") %>%
  mutate(DATE = as.Date(observed_on),
         year = year(DATE), 
         doy = yday(DATE)) %>% 
  filter(between(latitude, regions[3,2],regions[3,3]) & between(longitude, regions[3,4],regions[3,5])) %>%
  select(scientific_name, common_name, DATE, doy) %>%
  na.omit() # remove records with NA's in date or name
  


# 1) Get adultDF and larvalDF in proper format with scientific_name, year, DATE, count, and (optionally) SITE_ID
#    --this may involve merging in scientific name based on species code or English name

# 2) Restrict geographically as needed ('regions' has lat-long bounding boxes for each region; could instead use hex cells)

# 3) Create species-year loops

# 4) Use the lagFits() function available in gam_fx.R

