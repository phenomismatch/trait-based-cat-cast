### Comparing adult lep phenocurves to Caterpillars Count! phenology

#### Libraries ####

library(tidyverse)
library(lubridate)
library(purrr)
library(sf)
library(tmap)

#### Functions ####

# fn: lag lep data by x weeks, find R2 between y weekly caterpillar phenometric and 3 adult lep trait groups
lag_r2 <- function(weeks, cat_data, lep_data) {
  l <- weeks*7
  
  lep_data$lag_x <- lep_data$x - l
  
  lep_codes <- unique(lep_data$code)
  
  cc_lep_pheno <- cat_data %>%
    left_join(lep_data, by = c("julianweek" = "lag_x")) %>%
    filter(!is.na(code)) %>%
    pivot_wider(names_from = "code", values_from = "y")
  
  # Some hex cell/years missing some adult lep groups, only calculate R2 if those groups are present
  
  if("RE" %in% lep_codes) {
    re_r2 <- summary(lm(val ~ RE, data = cc_lep_pheno))$r.squared
  } else {re_r2 <- NA}
  
  if("RL" %in% lep_codes) {
    rl_r2 <- summary(lm(val ~ RL, data = cc_lep_pheno))$r.squared
  } else {rl_r2 <- NA}
  
  if("RP" %in% lep_codes) {
    rp_r2 <- summary(lm(val ~ RP, data = cc_lep_pheno))$r.squared
  } else {rp_r2 <- NA}
  
  return(data.frame(code = c("RE", "RL", "RP"), r2 = c(re_r2, rl_r2, rp_r2)))
}

# fn to compute lags for any hex, year, cat_data combination
output_lags <- function(cell, Year, data, ...) {
  cats <- data
  
  leps <- lep_pheno %>%
    filter(year == Year, HEXcell == cell)
  
  res <- data.frame(lag = c(), code = c(), r2 = c())
  for(w in -8:8){
    
    r2_df <- lag_r2(w, cats, leps)
    
    r2_df$lag <- w
    
    res <- rbind(res, r2_df)
    
  }
  
  res
  
}

## Backcasting fn - from code/gdd_backcast.R

backcast_cats <- function(my.hexes, my.years) {
  #set GDD model quantities
  gddP<-154
  gddL<-330
  
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
  
  return(pheno.x)
}

# fn to compute lags with backcasted lep data for any hex, year, cat_data combination
output_lags_backcast <- function(cell, Year, data, ...) {
  cats <- data
  
  leps <- backcast_pheno %>%
    filter(year == Year, hex == cell) %>%
    rename(x = "doyl", y = "pdf") # match with lep_pheno colnames
  
  res <- data.frame(lag = c(), code = c(), r2 = c())
  for(w in -8:8){
    
    r2_df <- lag_r2(w, cats, leps)
    
    r2_df$lag <- w
    
    res <- rbind(res, r2_df)
    
  }
  
  res
  
}


#### Set up ####

## ggplot theme

theme_set(theme_classic(base_size = 15))

## Read in adult lep curves

lep_pheno <- read_csv("data/simpleton_pheno_pdfs.csv")

## Read in CatCount data (weekly phenology at hex cells, sites with at least 6 good weeks)

cc_site_data <- read_csv("data/cc_subset_trait_based_pheno.csv")

## CatCount data availability

cc_sample_size <- cc_site_data %>%
  group_by(cell, Year) %>%
  summarize(n_sites = n_distinct(Name))

ggplot(cc_sample_size, aes(x = Year, y = as.factor(cell), size = n_sites)) + 
  geom_point() +
  scale_x_continuous(breaks = c(seq(2010, 2020, by = 2))) +
  labs(x = "Year", y = "Hex cell", size = "CC sites")
# ggsave("figures/cc_hex_data_avail.pdf")

# Map of hex cells represented

hex_sf <- read_sf("data/maps/hex_grid_crop.shp") %>%
  mutate(cell_num = as.numeric(cell)) %>%
  st_transform("+proj=ortho +lon_0=-75 +lat_0=40")

na_map <- read_sf("data/maps/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 %in% c("USA", "CAN")) %>%
  st_crop(c(xmin = -100, ymin = 20, xmax = -59, ymax = 51)) %>%
  st_transform("+proj=ortho +lon_0=-75 +lat_0=40")

sites_sf <- hex_sf %>%
  right_join(cc_sample_size, by = c("cell_num" = "cell")) %>%
  group_by(cell) %>%
  summarize(n_sites = max(n_sites),
            first_year = min(Year))

cc_site_map <- tm_shape(na_map) + tm_polygons() +
  tm_shape(sites_sf) + tm_polygons(col = "n_sites", palette = "YlGnBu", title = "CC! sites", alpha = 0.65) +
  tm_shape(sites_sf) + tm_text(text = "first_year") +
  tm_layout(legend.text.size = 1, legend.title.size = 1.5)
# tmap_save(cc_site_map, "figures/cc_site_map.pdf")

#### Correlations of lep adult curves and CC pheno curves ####

## CC hex phenology

cc_pheno <- cc_site_data %>%
  group_by(cell, Year, julianweek) %>%
  summarize(mean_dens= mean(meanDensity),
            mean_fracsurv = mean(fracSurveys),
            mean_biomass = mean(meanBiomass))

## Indiv phenology for Coweeta, Prairie Ridge, Bot Garden sites

cc_site_pheno <- cc_site_data %>%
  filter(Name %in% c("Coweeta - BB", "Coweeta - BS", "NC Botanical Garden", "Prairie Ridge Ecostation"))

## Check lags from -8 to 8 weeks, calculate R2 by hex-year

## Lags between adult lep pheno curves and caterpillar biomass, density, and fraction of surveys

cc_lep_lags <- cc_pheno %>%
  pivot_longer(names_to = "cat_metric", values_to = "val", mean_dens:mean_biomass) %>%
  group_by(cat_metric, cell, Year) %>%
  nest() %>%
  mutate(r2 = pmap(list(cell, Year, data), output_lags)) %>%
  select(-data) %>%
  unnest(cols = c("r2"))

## Plot R2 results

pdf(paste0(getwd(), "/figures/cc_adult_correlation.pdf"), height = 8, width = 10)

for(c in unique(cc_lep_lags$cell)) {
  
  plot_df <- cc_lep_lags %>%
    filter(cell == c)
  
  plot <- ggplot(plot_df, aes(x = lag, y = r2, col = as.factor(Year), group = Year)) + geom_line() + facet_grid(code ~ cat_metric)+
    labs(x = "Lag (weeks)", y = expression(R^2), title = paste0("Hex cell ", c), col = "Year") + theme_bw(base_size = 15)
  
  print(plot)
}

dev.off()

## Single site caterpillar to adult lep hex comparisons

cc_site_lags <- cc_site_pheno %>%
  pivot_longer(names_to = "cat_metric", values_to = "val", meanDensity:meanBiomass) %>%
  select(Name, cell, Year, julianweek, cat_metric, val) %>%
  group_by(cat_metric, Name, cell, Year) %>%
  nest() %>%
  mutate(r2 = pmap(list(cell, Year, data), output_lags)) %>%
  select(-data) %>%
  unnest(cols = c("r2"))

## Plot R2 results

pdf(paste0(getwd(), "/figures/cc_adult_correlation_NCsites.pdf"), height = 8, width = 10)

for(n in unique(cc_site_lags$Name)) {
  
  plot_df <- cc_site_lags %>%
    filter(Name == n)
  
  plot <- ggplot(plot_df, aes(x = lag, y = r2, col = as.factor(Year), group = Year)) + geom_line() + facet_grid(code ~ cat_metric)+
    labs(x = "Lag (weeks)", y = expression(R^2), title = paste0("Site: ", n), col = "Year") + theme_bw(base_size = 15)
  
  print(plot)
}

dev.off()

#### Correlations of backcasted cat curves and CC pheno curves ####

## Backcast pheno PDFs for hex-years with CC data

cc_pheno_1217 <- cc_pheno %>%
  filter(Year >= 2012 & Year <= 2017)

hexes <- unique(cc_pheno_1217$cell)
years <- unique(cc_pheno_1217$Year)

backcast_pheno <- backcast_cats(hexes, years)

## Test lags with CC hexes

cc_cat_lags <- cc_pheno_1217 %>%
  pivot_longer(names_to = "cat_metric", values_to = "val", mean_dens:mean_biomass) %>%
  group_by(cat_metric, cell, Year) %>%
  nest() %>%
  mutate(r2 = pmap(list(cell, Year, data), output_lags_backcast)) %>%
  select(-data) %>%
  unnest(cols = c("r2")) %>%
  filter(code == "RL")

## Plot R2 results

pdf(paste0(getwd(), "/figures/cc_catcast_correlation.pdf"), height = 4, width = 10)

for(c in unique(cc_cat_lags$cell)) {
  
  plot_df <- cc_cat_lags %>%
    filter(cell == c)
  
  plot <- ggplot(plot_df, aes(x = lag, y = r2, col = as.factor(Year), group = Year)) + geom_line() + facet_wrap(~cat_metric) +
    labs(x = "Lag (weeks)", y = expression(R^2), title = paste0("Hex cell ", c), col = "Year") + theme_bw(base_size = 15)
  
  print(plot)
}

dev.off()

## Test lags with CC NC sites

cc_site_cat_lags <- cc_site_pheno %>%
  pivot_longer(names_to = "cat_metric", values_to = "val", meanDensity:meanBiomass) %>%
  select(Name, cell, Year, julianweek, cat_metric, val) %>%
  filter(Year >= 2012 & Year <= 2017) %>%
  group_by(cat_metric, Name, cell, Year) %>%
  nest() %>%
  mutate(r2 = pmap(list(cell, Year, data), output_lags_backcast)) %>%
  select(-data) %>%
  unnest(cols = c("r2")) %>%
  filter(code == "RL")

## Plot R2 results

pdf(paste0(getwd(), "/figures/cc_catcast_correlation_NCsites.pdf"), height = 4, width = 10)

for(n in unique(cc_site_cat_lags$Name)) {
  
  plot_df <- cc_site_cat_lags %>%
    filter(Name == n)
  
  plot <- ggplot(plot_df, aes(x = lag, y = r2, col = as.factor(Year), group = Year)) + geom_line() + facet_wrap(~cat_metric)+
    labs(x = "Lag (weeks)", y = expression(R^2), title = paste0("Site: ", n), col = "Year") + theme_bw(base_size = 15)
  
  print(plot)
}

dev.off()

##### Phenometrics ####

## Compare early/late years in catcast and CC sites using deviations from mean 
## Calculate curve centroids for cat backcast and CC sites

lep_join <- lep_pheno %>%
  filter(code == "RL", year <= 2017 & year >= 2012, HEXcell %in% c(702, 703))

pheno_centroids <- backcast_pheno %>%
  left_join(cc_pheno_1217, by = c("year" = "Year", "hex" = "cell", "doyl" = "julianweek")) %>%
  group_by(hex, year) %>%
  mutate(catcount_centr = sum(doyl*mean_fracsurv, na.rm = T)/sum(mean_fracsurv, na.rm = T),
         catcast_centr = sum(doyl*pdf, na.rm = T)/sum(pdf, na.rm = T)) %>%
  right_join(lep_join, by = c("year", "code", "doyl" = "x", "hex" = "HEXcell")) %>%
  mutate(adult_peak = doyl[y == max(y, na.rm = T)])

## Plot phenocurves for catcast

pdf(paste0(getwd(), "/figures/catcast_phenocurves.pdf"), height = 8, width = 12)

for(c in unique(pheno_centroids$hex)) {
  
  plot_df <- pheno_centroids %>%
    filter(hex == c)
  
  plot <- ggplot(plot_df, aes(x = doyl)) + 
    geom_smooth(aes(y = pdf, col = "Larvae"), se = F) +
    geom_line(aes(y = y, col = "Adults")) + 
    geom_vline(aes(xintercept = catcast_centr, col = "Larvae"), lty = 2) +
    facet_wrap(~year) +
    labs(x = "Day of year", y = "PDF", col = "Life stage",
         title = paste0("Hex cell ", c)) + theme_bw(base_size = 15)
  
  print(plot)
}

dev.off()

## Correlation between early/late years 

pheno_dev <- pheno_centroids %>%
  distinct(year, hex, catcount_centr, catcast_centr) %>%
  group_by(hex) %>%
  mutate(mean_catcast_centr = mean(catcast_centr, na.rm = T),
         mean_catcount_centr = mean(catcount_centr, na.rm = T),
         catcast_dev = catcast_centr - mean_catcast_centr,
         catcount_dev = catcount_centr - mean_catcount_centr)

r_dev <- cor(pheno_dev$catcast_dev, pheno_dev$catcount_dev, use = "pairwise.complete.obs")

ggplot(pheno_dev, aes(x = catcount_dev, y = catcast_dev)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  annotate(geom = "text", x = 5, y = 30, label = paste0("r  = ", round(r_dev, 2)), size = 6) +
  labs(x = "Caterpillars Count! centroid deviation", y = "CatCast centroid deviation")
ggsave("figures/catcast_catcount_deviation_1to1.pdf")

## Early/late years in catcast and CC sites with temperature 

hex_temps <- read.csv("data/hex_mean_temps.csv", stringsAsFactors = F)

pheno_temps <- pheno_centroids %>%
  ungroup() %>%
  distinct(year, hex, code, catcount_centr, catcast_centr, adult_peak) %>%
  filter(!is.na(catcount_centr)) %>%
  left_join(hex_temps, by = c("hex" = "cell", "year")) %>%
  pivot_longer(names_to = "data", values_to = "pheno", catcount_centr:adult_peak) %>%
  mutate(plot_labels = case_when(data == "catcast_centr" ~ "CatCast",
                                 data == "catcount_centr" ~ "CatCount",
                                 data == "adult_peak" ~ "Adult curve"))

catcount_mod <- summary(lm(pheno ~ mean_temp, data = filter(pheno_temps, data == "catcount_centr")))

catcast_mod <- summary(lm(pheno ~ mean_temp, data = filter(pheno_temps, data == "catcast_centr")))

adult_mod <- summary(lm(pheno ~ mean_temp, data = filter(pheno_temps, data == "adult_peak")))

# ggplot cols
cols <- scales::hue_pal()(3)

ggplot(pheno_temps, aes(x = mean_temp, y = pheno, col = plot_labels)) + 
  geom_point(size = 2, alpha = 0.5) + 
  geom_smooth(method = "lm", se = F) +
  annotate(geom = "text", x = 15.5, y = 125, 
           label = paste0("p = ", round(catcount_mod$coefficients[2,4], 2), "; R2 = ", round(catcount_mod$r.squared, 2)), 
           col = cols[3]) +
  annotate(geom = "text", x = 15.5, y = 130, 
           label = paste0("p = ", round(catcast_mod$coefficients[2,4], 2), "; R2 = ", round(catcast_mod$r.squared, 2)), 
           col = cols[2]) +
  annotate(geom = "text", x = 15.5, y = 135, 
           label = paste0("Slope p = ", round(adult_mod$coefficients[2,4], 2), "; R2 = ", round(adult_mod$r.squared, 2)), 
           col = cols[1]) +
  labs(x = "Avg spring temperature (March-June)", y = "Peak/centroid date", col = "")
ggsave("figures/catcast_catcount_temp.pdf", units = "in", height = 6, width = 8)
