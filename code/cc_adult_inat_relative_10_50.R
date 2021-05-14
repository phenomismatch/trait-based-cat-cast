### Relative and absolute comparisons of lep phenology

library(tidyverse)
library(phenesse)
library(lubridate)
library(sf)
library(tmap)
library(cowplot)

## ggplot theme

theme_set(theme_classic(base_size = 15))

## Hex grid and maps

hex_sf <- read_sf("data/maps/hex_grid_crop.shp")

hex_sf <- hex_sf %>%
  mutate(centroid = st_centroid(hex_sf)$geometry,
         latitude = map_dbl(centroid, ~{
           c <- .
           c[[2]]
         }))

nam_sf <- read_sf("data/maps/ne_50m_admin_1_states_provinces_lakes.shp") %>%
  filter(sr_adm0_a3 %in% c("CAN", "USA"), iso_3166_2 != "US-AK", iso_3166_2 != "US-HI") %>%
  st_crop(xmin = -100, ymin = 25, xmax = -50, ymax = 55)

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

cc_pheno_unnest <- read_csv("data/derived_data/cc_weibull_10_50.csv")

## Visualize CC Weibull fits

# phenesse::create_predict_df
create_predict_df <- function(observations){
  # previous create_cdf_ends()
  weib <- fitdistrplus::fitdist(observations, distr = "weibull",
                                method = "mle")
  cdf0 <- as.numeric(weib$estimate['scale']*
                       (-log(1-0.01))^(1/weib$estimate['shape']))
  cdf100 <- as.numeric(weib$estimate['scale']*
                         (-log(1-0.99))^(1/weib$estimate['shape']))
  added_vec <- sort(append(observations, values = c(cdf0, cdf100)),
                    decreasing = FALSE)
  
  new_vec <- seq(from = min(added_vec), to = max(added_vec), by = 0.5)
  
  cdfadded <- 1 - exp(-(new_vec/weib$estimate['scale'])^weib$estimate['shape'])
  
  cdf_df <- data.frame(x = new_vec, y = cdfadded)
  ends <- data.frame(x = c(min(added_vec - 1), max(added_vec + 1)),
                     y = c(-0.001,1.001))
  cdf_df <- rbind(cdf_df, ends)
  cdf_df <- cdf_df[order(cdf_df$x, decreasing = FALSE),]
  
  return(cdf_df)
}

cc_weibull <- cc_pheno %>%
  mutate(cdf = map(pres_dates, ~create_predict_df(.)),
         pdf = map(cdf, ~{
           df <- .
           
           n <- nrow(df)
           
           res <- data.frame(x = df$x[1:(n-1)], y = diff(df$y)/diff(df$x))
           
           res[2:(n-2), ]
         }))

pdf("figures/cc_weibull_cdf.pdf", height = 8, width = 10)
for(n in unique(cc_weibull$Name)) {
  df <- cc_weibull %>%
    filter(Name == n) %>%
    select(-data)
  
  cdfs <- df %>%
    select(-pres_dates) %>%
    unnest(cols = c("cdf"))
  
  pts <- df %>%
    select(-cdf) %>%
    unnest(cols = c("pres_dates"))
  
  p <- ggplot(cdfs, aes(x = x, y = y)) + geom_line() + facet_wrap(~Year) +
    geom_vline(aes(xintercept = cc_10, col = "10%"), lty = 2) +
    geom_vline(aes(xintercept = cc_50, col = "50%"), lty = 2) +    
    geom_point(data = pts, aes(x = pres_dates, y = 0), alpha = 0.1) +
    labs(title = n, x = "Day of year", y = "CDF", col = "Percentile")

  print(p)
  
}
dev.off()

pdf("figures/cc_weibull_pdfs.pdf", height = 8, width = 10)
for(n in unique(cc_weibull$Name)) {
  df <- cc_weibull %>%
    filter(Name == n) %>%
    select(-cdf)
  
  pdfs <- df %>%
    select(-pres_dates) %>%
    unnest(cols = c("pdf"))
  
  pts <- df %>%
    select(-pdf) %>%
    unnest(cols = c("data"))
  
  p <- ggplot(pdfs, aes(x = x, y = y)) + geom_line() + facet_wrap(~Year) +
    geom_vline(aes(xintercept = cc_10, col = "10%"), lty = 2) +
    geom_vline(aes(xintercept = cc_50, col = "50%"), lty = 2) +    
    geom_line(data = pts, aes(x = julianweek, y = fracSurveys/1000), col = "darkgray") +
    labs(title = n, x = "Day of year", y = "PDF", col = "Percentile")
  
  print(p)
  
}
dev.off()


### Deviations in 10% and 50% phenometrics
## Use only sites within hex with >= 2 years, make sure in cells w/ multiple years, not flickering betw sites

good_sites <- cc_pheno_unnest %>%
  group_by(cell, Name) %>%
  summarize(n_year = n_distinct(Year)) %>%
  filter(n_year >= 2) %>%
  filter(Name != "UNC Chapel Hill Campus")

## Compare deviations in cc 10/50 with inat cats and adult bfly

inat_cats <- read_csv("data/derived_data/allCaterpillars_phenometrics.csv")

inat_cats_dev <- inat_cats %>%
  group_by(HEXcell) %>%
  mutate(mean10 = mean(w10, na.rm = T),
         mean50 = mean(w50, na.rm = T),
         dev10 = w10 - mean10,
         dev50 = w50 - mean50)

adult_bfly <- read_csv("data/derived_data/adult_bfly_phenometrics_phenesse.csv")

adult_bfly_dev <- adult_bfly %>%
  group_by(HEXcell, code) %>%
  mutate(mean10 = mean(w10, na.rm = T),
         mean50 = mean(w50, na.rm = T),
         dev10 = w10 - mean10,
         dev50 = w50 - mean50)

#### Figure 1 #####
# Data availability hex map for three datasets

cc_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(good_sites)

cc_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(cc_sf) + 
  tm_polygons(col = "n_year", alpha = 0.5, title = "Years", palette = "YlGnBu", breaks = c(1, 5, 10, 15, 20)) +
  tm_layout(main.title = "A) Caterpillars Count!", scale = 1.5)

inat_yrs <- inat_cats %>%
  group_by(HEXcell) %>%
  summarize(n_year = n_distinct(year))

inat_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(inat_yrs, by = c("cell" = "HEXcell"))

inat_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(inat_sf) + 
  tm_polygons(col = "n_year", alpha = 0.5, title = "Years", palette = "YlGnBu", breaks = c(1, 5, 10, 15, 20), legend.show = F) +
  tm_layout(main.title = "B) iNaturalist Caterpillars", scale = 1.5)

bfly_yrs <- adult_bfly %>%
  group_by(HEXcell) %>%
  summarize(n_year = n_distinct(year))

bfly_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(bfly_yrs, by = c("cell" = "HEXcell"))

bfly_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(bfly_sf) + 
  tm_polygons(col = "n_year", alpha = 0.5, title = "Years", palette = "YlGnBu", breaks = c(1, 5, 10, 15, 20), legend.show = F) +
  tm_layout(main.title = "C) Adult Butterflies", scale = 1.5)

data_map <- tmap_arrange(cc_map, inat_map, bfly_map, ncol = 2)
tmap_save(data_map, "figures/fig1_data_availability.pdf", units = "in", height = 8, width = 10)

#### Figure 2 ####
# 10% date maps, 2018

cc_2018 <- cc_pheno_unnest %>%
  filter(Name %in% good_sites$Name) %>%
  group_by(cell, Year) %>%
  summarize(mean10 = mean(cc_10, na.rm = T),
            mean50 = mean(cc_50, na.rm = T)) %>%
  filter(Year == 2018)

cc18_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(cc_2018)

cc18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(cc18_sf) + 
  tm_polygons(col = "mean10", alpha = 0.5, title = "10% date (DOY)", palette = "YlGn", breaks = seq(110, 230, by = 20)) +
  tm_layout(main.title = "A) Caterpillars Count!", scale = 1.5)

inat_2018 <- inat_cats %>%
  filter(year == 2018)

inat18_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(inat_2018, by = c("cell" = "HEXcell"))

inat18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(inat18_sf) + 
  tm_polygons(col = "w10", alpha = 0.5, title = "10% date (DOY)", palette = "YlGn", breaks = seq(110, 230, by = 20), legend.show = F) +
  tm_layout(main.title = "B) iNaturalist Caterpillars", scale = 1.5)

bfly_2018 <- adult_bfly %>%
  filter(year == 2018, code == "RL")

bfly18_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(bfly_2018, by = c("cell" = "HEXcell"))

bfly18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(bfly18_sf) + 
  tm_polygons(col = "w10", alpha = 0.5, title = "10% date (DOY)", palette = "YlGn", breaks = seq(110, 230, by = 20), legend.show = F) +
  tm_layout(main.title = "C) Adult Butterflies (overwinter larvae)", scale = 1.5)

data_map <- tmap_arrange(cc18_map, inat18_map, bfly18_map, ncol = 2)
tmap_save(data_map, "figures/fig2_10pct_dates_map.pdf", units = "in", height = 8, width = 10)

# 50% date maps, 2018

cc18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(cc18_sf) + 
  tm_polygons(col = "mean50", alpha = 0.5, title = "50% date (DOY)", palette = "YlGn", breaks = seq(165, 255, by = 15)) +
  tm_layout(main.title = "A) Caterpillars Count!", scale = 1.5)

inat18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(inat18_sf) + 
  tm_polygons(col = "w50", alpha = 0.5, title = "50% date (DOY)", palette = "YlGn", breaks = seq(165, 255, by = 15), legend.show = F) +
  tm_layout(main.title = "B) iNaturalist Caterpillars", scale = 1.5)

bfly18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(bfly18_sf) + 
  tm_polygons(col = "w50", alpha = 0.5, title = "50% date (DOY)", palette = "YlGn", breaks = seq(165, 255, by = 15), legend.show = F) +
  tm_layout(main.title = "C) Adult Butterflies (overwinter larvae)", scale = 1.5)

data_map <- tmap_arrange(cc18_map, inat18_map, bfly18_map, ncol = 2)
tmap_save(data_map, "figures/pct50_dates_map.pdf", units = "in", height = 8, width = 10)

#### Figure 3 #####
## Plot correlations of 10, 50% deviances

cc_dev <- cc_pheno_unnest %>%
  filter(Name %in% good_sites$Name) %>%
  group_by(cell, Name) %>%
  mutate(mean10 = mean(cc_10, na.rm = T),
         mean50 = mean(cc_50, na.rm = T), 
         dev10 = cc_10 - mean10,
         dev50 = cc_50 - mean50) %>%
  group_by(cell, Year) %>%
  summarize(dev10 = mean(dev10),
            dev50 = mean(dev50))

quant_dev <- cc_dev %>%
  left_join(select(adult_bfly_dev, year, HEXcell, code, dev10, dev50), by = c("Year" = "year", "cell" = "HEXcell"), suffix = c("_cc", "_adult")) %>%
  left_join(select(inat_cats_dev, year, HEXcell, dev10, dev50), by = c("Year" = "year", "cell" = "HEXcell"))

inat_bfly_dev <- select(adult_bfly_dev, year, HEXcell, code, dev10, dev50) %>%
  left_join(inat_cats_dev, by = c("year", "HEXcell"), suffix = c("_bfly", "_inat"))

# 10%
inat_cc10 <- ggplot(quant_dev, aes(x = dev10_cc, y = dev10)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 10% Caterpillars Count!", y = "Deviance 10% iNaturalist caterpillars") +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "in"))

inat_adult10 <- ggplot(filter(inat_bfly_dev, !is.na(code)), aes(x = dev10_bfly, y = dev10_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(method = "lm", se = F) +
  xlim(-40, 40) +
  labs(x = "Deviance 10% Adult butterflies", y = "Deviance 10% iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF")) +
  theme(legend.position = "none", plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "in"))

cc_adult10 <- ggplot(filter(quant_dev, !is.na(code)), aes(y = dev10_adult, x = dev10_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 10% Adult butterflies", x = "Deviance 10% Caterpillars Count!", col = "Overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.15, 0.8), plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "in"), legend.background = element_rect(fill = "transparent"))

plot_grid(inat_cc10, inat_adult10, cc_adult10, ncol = 2, labels = c("A", "B", "C"), label_size = 15, hjust = c(-2))
ggsave("figures/fig3_relative_adult_inat_cc_10.pdf", units = "in", height = 8, width = 10)

# 50%
inat_cc50 <- ggplot(quant_dev, aes(x = dev50_cc, y = dev50)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 50% Caterpillars Count!", y = "Deviance 50% iNaturalist caterpillars") +
  theme(plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "in"))

inat_adult50 <- ggplot(filter(inat_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-30, 30) +
  labs(x = "Deviance 50% Adult butterflies", y = "Deviance 50% iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF")) +
  theme(legend.position = "none", plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "in"))

cc_adult50 <- ggplot(filter(quant_dev, !is.na(code)), aes(y = dev50_adult, x = dev50_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 50% Adult butterflies", x = "Deviance 50% Caterpillars Count!", col = "Overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.15, 0.8), plot.margin = unit(c(0.5, 0.1, 0.1, 0.1), "in"), legend.background = element_rect(fill = "transparent"))

plot_grid(inat_cc50, inat_adult50, cc_adult50, ncol = 2, labels = c("A", "B", "C"))
ggsave("figures/relative_adult_inat_cc_50.pdf", units = "in", height = 8, width = 10)

#### Figure 2 sensitivity: different subsets of data ####
## Forest only subset

for_inat_cats <- read_csv("data/derived_data/FORESTONLY_allCaterpillars_phenometrics.csv")

for_inat_cats_dev <- for_inat_cats %>%
  group_by(HEXcell) %>%
  mutate(mean10 = mean(w10, na.rm = T),
         mean50 = mean(w50, na.rm = T),
         dev10 = w10 - mean10,
         dev50 = w50 - mean50)

for_adult_bfly <- read_csv("data/derived_data/ForestOnly_adult_bfly_phenometrics_phenesse.csv")

for_adult_bfly_dev <- for_adult_bfly %>%
  group_by(HEXcell, code) %>%
  mutate(mean10 = mean(w10, na.rm = T),
         mean50 = mean(w50, na.rm = T),
         dev10 = w10 - mean10,
         dev50 = w50 - mean50)

for_quant_dev <- cc_dev %>%
  left_join(select(for_adult_bfly_dev, year, HEXcell, code, dev10, dev50), by = c("Year" = "year", "cell" = "HEXcell"), suffix = c("_cc", "_adult")) %>%
  left_join(select(for_inat_cats_dev, year, HEXcell, dev10, dev50), by = c("Year" = "year", "cell" = "HEXcell"))

for_inat_bfly_dev <- select(for_adult_bfly_dev, year, HEXcell, code, dev10, dev50) %>%
  left_join(for_inat_cats_dev, by = c("year", "HEXcell"), suffix = c("_bfly", "_inat"))

# 10%
for_inat_cc10 <- ggplot(for_quant_dev, aes(x = dev10_cc, y = dev10)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 10% Caterpillars Count!", y = "Deviance 10% iNaturalist caterpillars")

for_inat_adult10 <- ggplot(filter(for_inat_bfly_dev, !is.na(code)), aes(x = dev10_bfly, y = dev10_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-25,50)+
  labs(x = "Deviance 10% Adult butterflies", y = "Deviance 10% iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = "none")

for_cc_adult10 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(y = dev10_adult, x = dev10_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 10% Adult butterflies", x = "Deviance 10% Caterpillars Count!", col = "Overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.85, 0.2))

plot_grid(for_inat_cc10, for_inat_adult10, for_cc_adult10, ncol = 2, labels = c("Forest only"))
ggsave("figures/relative_adult_inat_cc_10_forest.pdf", units = "in", height = 8, width = 10)

# 50%
for_inat_cc50 <- ggplot(for_quant_dev, aes(x = dev50_cc, y = dev50)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 50% Caterpillars Count!", y = "Deviance 50% iNaturalist caterpillars")
 
for_inat_adult50 <- ggplot(filter(for_inat_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-20, 30) +
  labs(x = "Deviance 50% Adult butterflies", y = "Deviance 50% iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = "none")

for_cc_adult50 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(y = dev50_adult, x = dev50_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 50% Adult butterflies", x = "Deviance 50% Caterpillars Count!", col = "Overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.85, 0.15), legend.background = element_rect(fill = "transparent"))

plot_grid(for_inat_cc50, for_inat_adult50, for_cc_adult50, ncol = 2)
ggsave("figures/relative_adult_inat_cc_50_forest.pdf", units = "in", height = 8, width = 10)

## iNat only data for adults

adult_bfly_inat <- read_csv("data/derived_data/iNatONLY_adult_bfly_phenometrics_phenesse.csv")

adult_inat_dev <- adult_bfly_inat %>%
  group_by(HEXcell, code) %>%
  mutate(mean10 = mean(w10, na.rm = T),
         mean50 = mean(w50, na.rm = T),
         dev10 = w10 - mean10,
         dev50 = w50 - mean50)

inat_only_bfly_dev <- select(adult_inat_dev, year, HEXcell, code, dev10, dev50) %>%
  left_join(inat_cats_dev, by = c("year", "HEXcell"), suffix = c("_bfly", "_inat"))

inat_only_10 <- ggplot(filter(inat_only_bfly_dev, !is.na(code)), aes(x = dev10_bfly, y = dev10_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-25, 25) +
  labs(x = "Deviance 10% Adult butterflies", y = "Deviance 10% iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = "none")

inat_only_50 <- ggplot(filter(inat_only_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-20, 20) +
  labs(x = "Deviance 50% Adult butterflies", y = "Deviance 50% iNaturalist caterpillars", col = "Overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.85, 0.15))

plot_grid(inat_only_10, inat_only_50)
ggsave("figures/relative_10_50_inat_only_.pdf", units = "in", height = 5, width = 10)

### Caterpillars Count! Obs excluded from iNat Cats

### Absolute comparisons: lag predicted by GDD and time

# GDD data
hex_gdd <- read_csv("data/derived_data/hex_gdd_2000-2020.csv")

accum_gdd <- hex_gdd %>%
  group_by(hex_cell, year) %>%
  mutate(accumGDD = cumsum(meanGDD)) %>%
  select(hex_cell, doy, year, accumGDD)

# 10 and 50 quantiles

inat_cats
adult_bfly

cc_quant <- cc_pheno_unnest %>%
  filter(Name %in% good_sites$Name) %>%
  group_by(cell, Year) %>%
  summarize(mean10 = round(mean(cc_10, na.rm = T)),
         mean50 = round(mean(cc_50, na.rm = T))) 

# Difference between 10/50 across datasets (units days)
# Difference between 10/50 across datasets (units GDD)

inat_cc_diff <- cc_quant %>%
  left_join(inat_cats, by = c("Year" = "year", "cell" = "HEXcell")) %>%
  mutate_at(c("w10", "w50"), ~round(.)) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "mean10" = "doy")) %>%
  rename("cc_gdd_10" = accumGDD) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "mean50" = "doy")) %>%
  rename("cc_gdd_50" = accumGDD) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "w10" = "doy")) %>%
  rename("inat_gdd_10" = accumGDD) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "w50" = "doy")) %>%
  rename("inat_gdd_50" = accumGDD) %>%
  mutate(diff_10_days = w10 - mean10,
         diff_50_days = w50 - mean50,
         diff_10_gdd = inat_gdd_10 - cc_gdd_10,
         diff_50_gdd = inat_gdd_50 - cc_gdd_50)

inat_bfly_diff <- inat_cats %>%
  left_join(adult_bfly, by = c("year", "HEXcell"), suffix = c("_inat", "_bfly")) %>%
  mutate_at(c("w10_inat", "w50_inat", "w10_bfly", "w50_bfly"), ~round(.)) %>%
  left_join(accum_gdd, by = c("year", "HEXcell" = "hex_cell", "w10_bfly" = "doy")) %>%
  rename("bfly_gdd_10" = accumGDD) %>%
  left_join(accum_gdd, by = c("year", "HEXcell" = "hex_cell", "w50_bfly" = "doy")) %>%
  rename("bfly_gdd_50" = accumGDD) %>%
  left_join(accum_gdd, by = c("year", "HEXcell" = "hex_cell", "w10_inat" = "doy")) %>%
  rename("inat_gdd_10" = accumGDD) %>%
  left_join(accum_gdd, by = c("year", "HEXcell" = "hex_cell", "w50_inat" = "doy")) %>%
  rename("inat_gdd_50" = accumGDD) %>%
  mutate(diff_10_days = w10_inat - w10_bfly,
         diff_50_days = w50_inat - w50_bfly,
         diff_10_gdd = inat_gdd_10 - bfly_gdd_10,
         diff_50_gdd = inat_gdd_50 - bfly_gdd_50)

cc_bfly_diff <- cc_quant %>%
  left_join(adult_bfly, by = c("Year" = "year", "cell" = "HEXcell")) %>%
  mutate_at(c("w10", "w50"), ~round(.)) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "mean10" = "doy")) %>%
  rename("cc_gdd_10" = accumGDD) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "mean50" = "doy")) %>%
  rename("cc_gdd_50" = accumGDD) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "w10" = "doy")) %>%
  rename("bfly_gdd_10" = accumGDD) %>%
  left_join(accum_gdd, by = c("Year" = "year", "cell" = "hex_cell", "w50" = "doy")) %>%
  rename("bfly_gdd_50" = accumGDD) %>%
  mutate(diff_10_days = mean10 - w10,
         diff_50_days = mean50 - w50,
         diff_10_gdd =  cc_gdd_10 - bfly_gdd_10,
         diff_50_gdd =  cc_gdd_50 - bfly_gdd_50)

# 10% lags, days and gdd for each dataset

inat_cc_plot <- ggplot(inat_cc_diff, aes(x = diff_10_days, y = diff_10_gdd)) +
  geom_point() + labs(x = "Lag iNat - CC! (days)", y = "Lag iNat - CC! GDD")

inat_adult_plot <- ggplot(inat_bfly_diff, aes(x = diff_10_days, y = diff_10_gdd, col = code)) +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  geom_point() + labs(x = "Lag iNat - Bfly (days)", y = "Lag iNat - Bfly GDD", col = "Overwinter")

adult_cc_plot <- ggplot(cc_bfly_diff, aes(x = diff_10_days, y = diff_10_gdd, col = code)) +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  geom_point() + labs(x = "Lag CC! - Bfly (days)", y = "Lag CC! - Bfly GDD", col = "Overwinter")

plot_grid(inat_cc_plot, inat_adult_plot, adult_cc_plot, nrow = 2)
ggsave("figures/lag10_days_gdd.pdf", units = "in", height = 8, width = 10)

# 50% lags
inat_cc_plot <- ggplot(inat_cc_diff, aes(x = diff_50_days, y = diff_50_gdd)) +
  geom_point() + labs(x = "Lag iNat - CC! (days)", y = "Lag iNat - CC! GDD")

inat_adult_plot <- ggplot(inat_bfly_diff, aes(x = diff_50_days, y = diff_50_gdd, col = code)) +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  geom_point() + labs(x = "Lag iNat - Bfly (days)", y = "Lag iNat - Bfly GDD", col = "Overwinter")

adult_cc_plot <- ggplot(cc_bfly_diff, aes(x = diff_50_days, y = diff_50_gdd, col = code)) +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  geom_point() + labs(x = "Lag CC! - Bfly (days)", y = "Lag CC! - Bfly GDD", col = "Overwinter")

plot_grid(inat_cc_plot, inat_adult_plot, adult_cc_plot, nrow = 2)
ggsave("figures/lag50_days_gdd.pdf", units = "in", height = 8, width = 10)

#### Figure 4 #####
## Density plots for lags in days - 10%

inat_cc_plot <- ggplot(inat_cc_diff, aes(x = diff_10_days)) +
  geom_density(fill = "gray") + labs(x = "10% Lag iNat - CC! (days)",  y = "Density")

inat_adult_plot <- ggplot(inat_bfly_diff, aes(x = diff_10_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "10% Lag iNat - Bfly (days)", fill = "Overwinter", y = "") +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae"))

adult_cc_plot <- ggplot(cc_bfly_diff, aes(x = diff_10_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "10% Lag CC! - Bfly (days)", y = "Density") +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae"))

legend <- get_legend(inat_adult_plot)

plot_grid(inat_cc_plot, inat_adult_plot + theme(legend.position = "none"), adult_cc_plot + theme(legend.position = "none"), legend,
          nrow = 2, labels = c("A", "B", "C"), label_size = 15)
ggsave("figures/fig4_lag10_days_density.pdf", units = "in", height = 8, width = 10)

## Density plots for lags in days - 50%
inat_cc_plot <- ggplot(inat_cc_diff, aes(x = diff_50_days)) +
  geom_density(fill = "gray") + labs(x = "50% Lag iNat - CC! (days)", y = "Density")

inat_adult_plot <- ggplot(inat_bfly_diff, aes(x = diff_50_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "50% Lag iNat - Bfly (days)", y = "", fill = "Overwinter") +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae"))

adult_cc_plot <- ggplot(cc_bfly_diff, aes(x = diff_50_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "50% Lag CC! - Bfly (days)", y = "Density") +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae"))

legend <- get_legend(inat_adult_plot)

plot_grid(inat_cc_plot, inat_adult_plot + theme(legend.position = "none"), adult_cc_plot + theme(legend.position = "none"), legend,
          nrow = 2, labels = c("A", "B", "C"), label_size = 15)
ggsave("figures/lag50_days_density.pdf", units = "in", height = 8, width = 10)

# Variation (standard devation) across years w/in cell
# 4 maps, 1 per var/unit combo, shaded hexes by SD, marginal histograms

### Model variation in lags

## Map: recent year, lag values

inat_cc_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(inat_cc_diff)

inat_bfly_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(inat_bfly_diff, by = c("cell" = "HEXcell"))

cc_bfly_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(cc_bfly_diff)

# 2018, 10%
inat_cc_map <- tm_shape(nam_sf) + tm_polygons() +
  tm_shape(inat_cc_sf %>%
             filter(Year == 2018)) + tm_polygons(col = "diff_10_days", alpha = 0.5, palette = "YlGnBu", title = "10% iNat - CC!")

inat_bfly_map <-  tm_shape(nam_sf) + tm_polygons() +
  tm_shape(inat_bfly_sf %>%
             filter(year == 2018, code == "RL")) + tm_polygons(col = "diff_10_days", alpha = 0.5, palette = "YlGnBu", title = "10% iNat - Bfly")

cc_bfly_map <- tm_shape(nam_sf) + tm_polygons() +
  tm_shape(cc_bfly_sf %>%
             filter(Year == 2018, code == "RL")) + tm_polygons(col = "diff_10_days", alpha = 0.5, palette = "YlGnBu", title = "10% CC! - Bfly")

panel_10 <- tmap_arrange(inat_cc_map, inat_bfly_map, cc_bfly_map, nrow = 2)
tmap_save(panel_10, "figures/lag10_2018_map.pdf", units = "in", height = 8, width = 10)

# 2018, 50%
inat_cc_map <- tm_shape(nam_sf) + tm_polygons() +
  tm_shape(inat_cc_sf %>%
             filter(Year == 2018)) + tm_polygons(col = "diff_50_days", alpha = 0.5, palette = "YlGnBu", title = "50% iNat - CC!")

inat_bfly_map <-  tm_shape(nam_sf) + tm_polygons() +
  tm_shape(inat_bfly_sf %>%
             filter(year == 2018, code == "RL")) + tm_polygons(col = "diff_50_days", alpha = 0.5, palette = "YlGnBu", title = "50% iNat - Bfly")

cc_bfly_map <- tm_shape(nam_sf) + tm_polygons() +
  tm_shape(cc_bfly_sf %>%
             filter(Year == 2018, code == "RL")) + tm_polygons(col = "diff_50_days", alpha = 0.5, palette = "YlGnBu", title = "50% CC! - Bfly")

panel_50 <- tmap_arrange(inat_cc_map, inat_bfly_map, cc_bfly_map, nrow = 2)
tmap_save(panel_50, "figures/lag50_2018_map.pdf", units = "in", height = 8, width = 10)

#### Figure 5 ####
## Model: lag ~ latitude + temp + temp:lat + (overwintering?)

hex_temps <- read_csv("data/derived_data/hex_mean_temps.csv")

inat_cc_mod <- inat_cc_sf %>%
  left_join(hex_temps, by = c("cell", "Year" = "year")) %>%
  mutate(dataset = "iNat - CC!") %>%
  group_by(dataset) %>%
  nest()

inat_bfly_mod <- inat_bfly_sf %>%
  left_join(hex_temps, by = c("cell", "year")) %>%
  mutate(dataset = "iNat - Bfly") %>%
  group_by(dataset) %>%
  nest()

cc_bfly_mod <- cc_bfly_sf %>%
  left_join(hex_temps, by = c("cell", "Year" = "year")) %>%
  mutate(dataset = "CC! - Bfly") %>%
  group_by(dataset) %>%
  nest()

mod_all <- inat_cc_mod %>%
  rbind(inat_bfly_mod, cc_bfly_mod) %>%
  mutate(mod10 = map(data, ~lm(diff_10_days ~ latitude + mean_temp + latitude:mean_temp, data = .)),
         mod50 = map(data, ~lm(diff_50_days ~ latitude + mean_temp + latitude:mean_temp, data = .)),
         tidy10 = map(mod10, ~broom::tidy(.)),
         tidy50 = map(mod50, ~broom::tidy(.)))

mod_ests_10 <- mod_all %>%
  select(dataset, tidy10) %>%
  unnest(cols = c("tidy10")) %>%
  filter(term != "(Intercept)")

ggplot(mod_ests_10, aes(x = term, y = estimate, col = dataset)) + 
  geom_point(cex = 2, position = position_dodge(width = 0.3)) + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), 
                                                          width = 0.1, cex = 1, position = position_dodge(width = 0.3)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = "Estimate", col = "Lag 10% date") +
  scale_x_discrete(labels = c("mean_temp" = "Spring temp", "latitude:mean_temp" = "Latitude:Spring temp", "latitude" = "Latitude")) +
  scale_color_manual(values = c("#3E4A89FF", "#1F9E89FF", "#6DCD59FF")) +
  theme(legend.position = c(0.2, 0.2)) +
  coord_flip()
ggsave("figures/fig5_lag10_mod_ests.pdf")

mod_ests_50 <- mod_all %>%
  select(dataset, tidy50) %>%
  unnest(cols = c("tidy50")) %>%
  filter(term != "(Intercept)")

ggplot(mod_ests_50, aes(x = term, y = estimate, col = dataset)) + 
  geom_point(cex = 2, position = position_dodge(width = 0.3)) + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), 
                width = 0.1, cex = 1, position = position_dodge(width = 0.3)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = "Estimate", col = "Lag 50% date") +
  scale_x_discrete(labels = c("mean_temp" = "Spring temp", "latitude:mean_temp" = "Latitude:Spring temp", "latitude" = "Latitude")) +
  scale_color_manual(values = c("#3E4A89FF", "#1F9E89FF", "#6DCD59FF")) +
  theme(legend.position = c(0.2, 0.2)) +
  coord_flip()
ggsave("figures/lag50_mod_ests.pdf")

sjPlot::plot_model(mod_all$mod10[[2]], type = "int")
sjPlot::plot_model(mod_all$mod50[[2]], type = "int")
