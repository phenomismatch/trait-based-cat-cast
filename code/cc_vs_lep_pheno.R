### Comparing adult lep phenocurves to Caterpillars Count! phenology

library(tidyverse)
library(purrr)
library(sf)
library(tmap)

## Read in adult lep curves

lep_pheno <- read_csv("data/simpleton_pheno_pdfs.csv")

## Read in CatCount data (weekly phenology at hex cells, sites with at least 6 good weeks)

cc_site_data <- read_csv("data/cc_subset_trait_based_pheno.csv")

## CatCount data availability

cc_sample_size <- cc_site_data %>%
  group_by(cell, Year) %>%
  summarize(n_sites = n_distinct(Name))

theme_set(theme_classic(base_size = 15))
ggplot(cc_sample_size, aes(x = Year, y = as.factor(cell), size = n_sites)) + 
  geom_point() +
  scale_x_continuous(breaks = c(seq(2010, 2020, by = 2))) +
  labs(x = "Year", y = "Hex cell", size = "CC sites")
ggsave("figures/cc_hex_data_avail.pdf")

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
tmap_save(cc_site_map, "figures/cc_site_map.pdf")

## CC hex phenology

cc_pheno <- cc_site_data %>%
  group_by(cell, Year, julianweek) %>%
  summarize(mean_dens= mean(meanDensity),
            mean_fracsurv = mean(fracSurveys),
            mean_biomass = mean(meanBiomass))

## Correlation with 3 diff wintering stage lep phenocurves, no lag

cc_lep_pheno <- cc_pheno %>%
  left_join(lep_pheno, by = c("julianweek" = "x", "Year" = "year", "cell" = "HEXcell")) %>%
  filter(!is.na(code)) %>%
  pivot_wider(names_from = "code", values_from = "y")

cor(select(ungroup(cc_lep_pheno), 
           mean_dens, mean_fracsurv, mean_biomass, RE, RL, RP),
    use = "pairwise.complete.obs")

## Check lags from -8 to 8 weeks

# fn: lag lep data by x weeks, find R2 between y caterpillar metric and 3 adult lep trait groups
lag_r2 <- function(weeks, cat_data) {
  l <- weeks*7
  
  df <- lep_pheno
  df$lag_x <- df$x - l
  
  cc_lep_pheno <- cat_data %>%
    left_join(df, by = c("julianweek" = "lag_x", "Year" = "year", "cell" = "HEXcell")) %>%
    filter(!is.na(code)) %>%
    pivot_wider(names_from = "code", values_from = "y")
  
  re_r2 <- summary(lm(val ~ RE, data = cc_lep_pheno))$r.squared
  rl_r2 <- summary(lm(val ~ RL, data = cc_lep_pheno))$r.squared
  rp_r2 <- summary(lm(val ~ RP, data = cc_lep_pheno))$r.squared
  
  return(data.frame(code = c("RE", "RL", "RP"), r2 = c(re_r2, rl_r2, rp_r2)))
}

## Lags between adult lep pheno curves and caterpillar biomass, density, and fraction of surveys

cc_lep_lags <- cc_pheno %>%
  pivot_longer(names_to = "cat_metric", values_to = "val", mean_dens:mean_biomass) %>%
  group_by(cat_metric) %>%
  nest() %>%
  mutate(r2 = map(data, ~{
    cats <- .
    
    res <- data.frame(lag = c(), code = c(), r2 = c())
    for(w in -8:8){
      
      r2_df <- lag_r2(w, cats)
      
      r2_df$lag <- w
      
      res <- rbind(res, r2_df)
      
    }
    
    res

  })) %>%
  select(-data) %>%
  unnest(cols = c("r2"))

## Plot R2 results

ggplot(cc_lep_lags, aes(x = lag, y = r2)) + geom_line() + facet_grid(code ~ cat_metric)+
  labs(x = "Lag (weeks)", y = expression(R^2))
ggsave("figures/cc_adult_correlation.pdf")
