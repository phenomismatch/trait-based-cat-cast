### Relative and absolute comparisons of lep phenology

library(tidyverse)
library(phenesse)
library(lubridate)
library(sf)
library(tmap)

## ggplot theme

theme_set(theme_classic(base_size = 15))

## Read in CatCount data (weekly phenology at hex cells, sites with at least 6 good weeks)

cc_site_data <- read_csv("data/derived_data/cc_subset_trait_based_pheno.csv")

### Estimate 10% and 50% phenometrics for Caterpillars Count! data

cc_pres_only <- cc_site_data %>%
  filter(totalCount > 0)

cc_pheno <- cc_pres_only %>%
  group_by(Name, cell, Latitude, Longitude, Year) %>%
  nest() %>%
  mutate(n_dates = map_dbl(data, ~n_distinct(.$julianweek))) %>%
  filter(n_dates >= 4) %>%
  mutate(pres_dates = map(data, ~{
    df <- .
    
    rep(df$julianweek, df$totalCount)
  })) %>%
  mutate(cc_10 = map_dbl(pres_dates, ~weib_percentile(., percentile = 0.1, iterations = 500)),
         cc_50 = map_dbl(pres_dates, ~weib_percentile(., percentile = 0.5, iterations = 500)))

# cc_pheno_write <- cc_pheno %>%
#   select(-data, -pres_dates)
# write.csv(cc_pheno_write, "data/derived_data/cc_weibull_10_50.csv", row.names = F)

cc_pheno_unnest <- cc_pheno %>%
  select(-pres_dates) %>%
  unnest(cols = c("data"))

### Deviations in 10% and 50% phenometrics

## Use only sites within hex with >= 2 years, meet min year threshold