### Relative and absolute comparisons of lep phenology

library(tidyverse)
library(phenesse)
library(lubridate)
library(sf)
library(tmap)
library(cowplot)
library(grid)

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
# 
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

adult_bfly <- read_csv("data/derived_data/adult_bfly_phenometrics_iNatOnlyData.csv")

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
  tm_polygons(col = "n_year", alpha = 0.5, title = "Years", palette = "YlGnBu", breaks = c(2, 5, 10, 15, 20, 25)) +
  tm_shape(cc_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "A) Caterpillars Count!", scale = 1.5)

inat_yrs <- inat_cats %>%
  group_by(HEXcell) %>%
  summarize(n_year = n_distinct(year))

inat_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(inat_yrs, by = c("cell" = "HEXcell"))

inat_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(inat_sf) + 
  tm_polygons(col = "n_year", alpha = 0.5, title = "Years", palette = "YlGnBu", breaks = c(2, 5, 10, 15, 20, 25), legend.show = F) +
  tm_shape(cc_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "B) iNaturalist Caterpillars", scale = 1.5)

bfly_yrs <- adult_bfly %>%
  group_by(HEXcell) %>%
  summarize(n_year = n_distinct(year))

bfly_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(bfly_yrs, by = c("cell" = "HEXcell"))

bfly_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(bfly_sf) + 
  tm_polygons(col = "n_year", alpha = 0.5, title = "Years", palette = "YlGnBu", breaks = c(2, 5, 10, 15, 20, 25), legend.show = F) +
  tm_shape(cc_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "C) Adult Butterflies", scale = 1.5)

data_map <- tmap_arrange(cc_map, inat_map, bfly_map, ncol = 2)
tmap_save(data_map, "figures/fig1_data_availability.pdf", units = "in", height = 8, width = 10)

#### Figure 2 ####
# 10% date maps, 2018

cols <- RColorBrewer::brewer.pal(9, "YlGn")
fig2_palette <- cols[c(1,3,5,7,9)]

cc_2018 <- cc_pheno_unnest %>%
  filter(Name %in% good_sites$Name) %>%
  group_by(cell, Year) %>%
  summarize(mean10 = mean(cc_10, na.rm = T),
            mean50 = mean(cc_50, na.rm = T)) %>%
  filter(Year == 2018)

cc18_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(cc_2018)

cc18_map <- tm_shape(nam_sf) + tm_polygons(col = "gray95") + tm_shape(cc18_sf) + 
  tm_polygons(col = "mean10", alpha = 0.75, title = "10% date (DOY)", palette = "YlGn", breaks = seq(50, 250, by = 50)) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "A) Caterpillars Count!", scale = 1.5)

inat_2018 <- inat_cats %>%
  filter(year == 2018)

inat18_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(inat_2018, by = c("cell" = "HEXcell"))

inat18_map <- tm_shape(nam_sf) + tm_polygons(col = "gray95") + tm_shape(inat18_sf) + 
  tm_polygons(col = "w10", alpha = 0.75, title = "10% date (DOY)", palette = "YlGn", breaks = seq(50, 250, by = 50), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "B) iNaturalist Caterpillars", scale = 1.5)

bfly_2018 <- adult_bfly %>%
  filter(year == 2018, code == "RL")

bfly18_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(bfly_2018, by = c("cell" = "HEXcell"))

bfly18_map <- tm_shape(nam_sf) + tm_polygons(col = "gray95") + tm_shape(bfly18_sf) + 
  tm_polygons(col = "w10", alpha = 0.75, title = "10% date (DOY)", palette = "YlGn", breaks = seq(50, 250, by = 50), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "C) Adult Butterflies (overwinter as larvae)", scale = 1.5)

egg_2018 <- adult_bfly %>%
  filter(year == 2018, code == "RE")

egg_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(egg_2018, by = c("cell" = "HEXcell"))

egg_map <- tm_shape(nam_sf) + tm_polygons(col = "gray95") + tm_shape(egg_sf) + 
  tm_polygons(col = "w10", alpha = 0.75, title = "10% date (DOY)", palette = "YlGn", breaks = seq(50, 250, by = 50), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "D) Adult Butterflies (overwinter as eggs)", scale = 1.5)

pup_2018 <- adult_bfly %>%
  filter(year == 2018, code == "RP")

pup_sf <- hex_sf %>%
  mutate_at(c("cell"), ~as.numeric(.)) %>%
  right_join(pup_2018, by = c("cell" = "HEXcell"))

pup_map <- tm_shape(nam_sf) + tm_polygons(col = "gray95") + tm_shape(pup_sf) + 
  tm_polygons(col = "w10", alpha = 0.75, title = "10% date (DOY)", palette = "YlGn", breaks = seq(50, 250, by = 50), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "E) Adult Butterflies (overwinter as pupae)", scale = 1.5,
            inner.margins = c(0.12, 0.02, 0.02, 0.02), 
            outer.margins = c(0.01,0.01,0.01,0.01))

inat18_df <- inat18_sf %>%
  st_set_geometry(NULL)
cc18_df <- cc18_sf %>%
  st_set_geometry(NULL) 
bfly18_df <- bfly18_sf %>%
  st_set_geometry(NULL) 
egg_df <- egg_sf %>%
  st_set_geometry(NULL)
pup_df <- pup_sf %>%
  st_set_geometry(NULL)

lat_regression <- ggplot() + 
  geom_point(data = inat18_df, aes(x = latitude, y = w10, col = "iNaturalist caterpillars")) + 
  geom_smooth(data = inat18_df,aes(x = latitude, y = w10, col = "iNaturalist caterpillars", fill = "iNaturalist caterpillars"), method = "lm", show.legend = F) +
  geom_point(data = cc18_df, aes(x = latitude, y = mean10, col = 'Caterpillars Count!')) +
  geom_smooth(data = cc18_df, aes(x = latitude, y = mean10, col = 'Caterpillars Count!', fill = "Caterpillars Count!"), method = "lm", show.legend = F) +
  geom_point(data = bfly18_df, aes(x = latitude, y = w10, col = 'Adults overwinter as larvae')) +
  geom_smooth(data = bfly18_df, aes(x = latitude, y = w10, col = 'Adults overwinter as larvae', fill = 'Adults overwinter as larvae'), method = "lm", show.legend = F) +
  geom_point(data = egg_df, aes(x = latitude, y = w10, col = 'Adults overwinter as eggs')) +
  geom_smooth(data = egg_df, aes(x = latitude, y = w10, col = 'Adults overwinter as eggs', fill = 'Adults overwinter as eggs'), method = "lm", show.legend = F) +
  geom_point(data = pup_df, aes(x = latitude, y = w10, col = 'Adults overwinter as pupae')) +
  geom_smooth(data = pup_df, aes(x = latitude, y = w10, col = 'Adults overwinter as pupae', fill = 'Adults overwinter as pupae'), method = "lm", show.legend = F)+
  labs(x = "Latitude (centerpoint of hex cell)", y = "10% date (DOY)", col = "", title = "F)") +
  ylim(30, 220) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.7, 0.25), legend.background = element_rect(fill = "transparent"),
        title = element_text(size = 16)) +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF", "#D7BDE2","dodgerblue4")) +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF", "#D7BDE2","dodgerblue4"))

pdf(paste0(getwd(), "/figures/fig2_10pct_dates_map.pdf"), height = 15, width = 12)
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
print(cc18_map, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(inat18_map, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(bfly18_map, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(egg_map, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(pup_map, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(lat_regression, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()

# 50% date maps, 2018

cc18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(cc18_sf) + 
  tm_polygons(col = "mean50", alpha = 0.75, title = "50% date (DOY)", palette = "YlGn", breaks = seq(165, 255, by = 15)) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "A) Caterpillars Count!", scale = 1.5)

inat18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(inat18_sf) + 
  tm_polygons(col = "w50", alpha = 0.75, title = "50% date (DOY)", palette = "YlGn", breaks = seq(165, 255, by = 15), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "B) iNaturalist Caterpillars", scale = 1.5)

bfly18_map <- tm_shape(nam_sf) + tm_polygons() + tm_shape(bfly18_sf) + 
  tm_polygons(col = "w50", alpha = 0.75, title = "50% date (DOY)", palette = "YlGn", breaks = seq(165, 255, by = 15), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "C) Adult Butterflies (overwinter larvae)", scale = 1.5)

egg_map <- tm_shape(nam_sf) + tm_polygons(col = "gray95") + tm_shape(egg_sf) + 
  tm_polygons(col = "w50", alpha = 0.75, title = "10% date (DOY)", palette = "YlGn", breaks = seq(50, 250, by = 50), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "D) Adult Butterflies (overwinter as eggs)", scale = 1.5)

pup_map <- tm_shape(nam_sf) + tm_polygons(col = "gray95") + tm_shape(pup_sf) + 
  tm_polygons(col = "w50", alpha = 0.75, title = "10% date (DOY)", palette = "YlGn", breaks = seq(50, 250, by = 50), legend.show = F) +
  tm_shape(cc18_sf) + tm_borders(lwd = 2) +
  tm_layout(main.title = "E) Adult Butterflies (overwinter as pupae)", scale = 1.5)


inat18_df <- inat18_sf %>%
  st_set_geometry(NULL)
cc18_df <- cc18_sf %>%
  st_set_geometry(NULL) 
bfly18_df <- bfly18_sf %>%
  st_set_geometry(NULL) 
egg_df <- egg_sf %>%
  st_set_geometry(NULL)
pup_df <- pup_sf %>%
  st_set_geometry(NULL)

lat_regression <- ggplot() + 
  geom_point(data = inat18_df, aes(x = latitude, y = w50, col = "iNaturalist caterpillars")) + 
  geom_smooth(data = inat18_df,aes(x = latitude, y = w50, col = "iNaturalist caterpillars", fill = "iNaturalist caterpillars"), method = "lm", show.legend = F) +
  geom_point(data = cc18_df, aes(x = latitude, y = mean50, col = 'Caterpillars Count!')) +
  geom_smooth(data = cc18_df, aes(x = latitude, y = mean50, col = 'Caterpillars Count!', fill = "Caterpillars Count!"), method = "lm", show.legend = F) +
  geom_point(data = bfly18_df, aes(x = latitude, y = w50, col = 'Adults overwinter as larvae')) +
  geom_smooth(data = bfly18_df, aes(x = latitude, y = w50, col = 'Adults overwinter as larvae', fill = 'Adults overwinter as larvae'), method = "lm", show.legend = F) +
  geom_point(data = egg_df, aes(x = latitude, y = w50, col = 'Adults overwinter as eggs')) +
  geom_smooth(data = egg_df, aes(x = latitude, y = w50, col = 'Adults overwinter as eggs', fill = 'Adults overwinter as eggs'), method = "lm", show.legend = F) +
  geom_point(data = pup_df, aes(x = latitude, y = w50, col = 'Adults overwinter as pupae')) +
  geom_smooth(data = pup_df, aes(x = latitude, y = w50, col = 'Adults overwinter as pupae', fill = 'Adults overwinter as pupae'), method = "lm", show.legend = F)+
  labs(x = "Latitude (centerpoint of hex cell)", y = "50% date (DOY)", col = "", title = "F)") +
  ylim(130, 210) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.7, 0.25), legend.background = element_rect(fill = "transparent")) +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF", "#D7BDE2","dodgerblue4")) +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF", "#D7BDE2","dodgerblue4"))

pdf(paste0(getwd(), "/figures/pct50_dates_map.pdf"), height = 15, width = 12)
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 2)))
print(cc18_map, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(inat18_map, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(bfly18_map, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(egg_map, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(pup_map, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(lat_regression, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()

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
  annotate(geom = "text", x = -3, y = 50, label = "Deviance in 10% date", size = 7) +
  annotate(geom = "text", x = 6, y = -55, label = "r = 0.56***", size = 5) +
  labs(x = expression(paste(Delta, " Caterpillars Count!")),
       y = expression(paste(Delta, " iNaturalist caterpillars"))) +
  theme_set(theme_classic(base_size  = 16))

inat_adult10 <- ggplot(filter(inat_bfly_dev, !is.na(code)), aes(x = dev10_bfly, y = dev10_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_set(theme_classic(base_size  = 16)) +
  geom_smooth(method = "lm", se = F) +
  xlim(-40, 40) +
  ylim(-80, 100) +
  labs(x = expression(paste(Delta, " Adult butterflies")), y = expression(paste(Delta, " iNaturalist caterpillars")), col = "Adults overwinter as") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), 
                     labels = c("RE" = "eggs, r = 0.26**", "RL" = "larvae, r = 0.35***", "RP" = "pupae, r = 0.39***")) +
  theme(legend.position = c(0.25, 0.8), legend.background = element_rect(fill = "transparent"))

cc_adult10 <- ggplot(filter(quant_dev, !is.na(code)), aes(y = dev10_adult, x = dev10_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_set(theme_classic(base_size = 16)) +
  geom_smooth(method = "lm", se = F) +
  ylim(-40, 40) +
  labs(y = expression(paste(Delta, " Adult butterflies")), x = expression(paste(Delta, " Caterpillars Count!")), col = "Adults overwinter as") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), 
                     labels = c("RE" = "eggs, r = 0.70", "RL" = "larvae, r = 0.21", "RP" = "pupae, r = 61**")) +
  theme(legend.position = c(0.25, 0.8), legend.background = element_rect(fill = "transparent")) +
  annotate(geom= "text", x = -8, y =-38, label = "early", size = 5) +
  annotate(geom= "text", x = 8, y = -38, label = "late", size = 5)

plot_grid(inat_cc10, inat_adult10, cc_adult10, ncol = 2, labels = c("A", "B", "C"), label_size = 15)
ggsave("figures/fig3_relative_adult_inat_cc_10.pdf", units = "in", height = 8, width = 10)

# 50%
inat_cc50 <- ggplot(quant_dev, aes(x = dev10_cc, y = dev50)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_smooth(method = "lm", se = F) +
  annotate(geom = "text", x = -3, y = 50, label = "Deviance in 50% date", size = 7) +
  annotate(geom = "text", x = 6, y = -55, label = "r = 0.32*", size = 5) +
  labs(x = expression(paste(Delta, " Caterpillars Count!")),
       y = expression(paste(Delta, " iNaturalist caterpillars"))) +
  theme_set(theme_classic(base_size  = 16))

inat_adult50 <- ggplot(filter(inat_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_set(theme_classic(base_size  = 16)) +
  geom_smooth(method = "lm", se = F) +
  labs(x = expression(paste(Delta, " Adult butterflies")), y = expression(paste(Delta, " iNaturalist caterpillars")), col = "Adults overwinter as") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), 
                     labels = c("RE" = "eggs, r = 0.33***", "RL" = "larvae, r = 0.41***", "RP" = "pupae, r = 0.45***")) +
  theme(legend.position = c(0.25, 0.8), legend.background = element_rect(fill = "transparent"))

cc_adult50 <- ggplot(filter(quant_dev, !is.na(code)), aes(y = dev50_adult, x = dev50_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  theme_set(theme_classic(base_size = 16)) +
  geom_smooth(method = "lm", se = F) +
  labs(y = expression(paste(Delta, " Adult butterflies")), x = expression(paste(Delta, " Caterpillars Count!")), col = "Adults overwinter as") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), 
                     labels = c("RE" = "eggs, r = 0.06", "RL" = "larvae, r = 0.06", "RP" = "pupae, r = 0.31")) +
  theme(legend.position = c(0.25, 0.8), legend.background = element_rect(fill = "transparent")) +
  annotate(geom= "text", x = -8, y =-28, label = "early", size = 5) +
  annotate(geom= "text", x = 4, y = -28, label = "late", size = 5)

plot_grid(inat_cc50, inat_adult50, cc_adult50, ncol = 2, labels = c("A", "B", "C"))
ggsave("figures/relative_adult_inat_cc_50.pdf", units = "in", height = 8, width = 10)

## Correlation table

cor_res <- vector(mode = "list", length = 14)

cor_res[[1]] <- cor.test(quant_dev$dev10_cc, quant_dev$dev10)
cor_res[[2]] <- cor.test(quant_dev$dev50_cc, quant_dev$dev50)

cor_res[[3]] <- cor.test(quant_dev$dev10_cc[quant_dev$code == "RE"], quant_dev$dev10_adult[quant_dev$code == "RE"])
cor_res[[4]] <- cor.test(quant_dev$dev50_cc[quant_dev$code == "RE"], quant_dev$dev50_adult[quant_dev$code == "RE"])

cor_res[[5]] <- cor.test(quant_dev$dev10_cc[quant_dev$code == "RL"], quant_dev$dev10_adult[quant_dev$code == "RL"])
cor_res[[6]] <- cor.test(quant_dev$dev50_cc[quant_dev$code == "RL"], quant_dev$dev50_adult[quant_dev$code == "RL"])

cor_res[[7]] <- cor.test(quant_dev$dev10_cc[quant_dev$code == "RP"], quant_dev$dev10_adult[quant_dev$code == "RP"])
cor_res[[8]] <- cor.test(quant_dev$dev50_cc[quant_dev$code == "RP"], quant_dev$dev50_adult[quant_dev$code == "RP"])

cor_res[[9]] <- cor.test(inat_bfly_dev$dev10_inat[inat_bfly_dev$code == "RE"], inat_bfly_dev$dev10_bfly[inat_bfly_dev$code == "RE"])
cor_res[[10]] <- cor.test(inat_bfly_dev$dev50_inat[inat_bfly_dev$code == "RE"], inat_bfly_dev$dev50_bfly[inat_bfly_dev$code == "RE"])

cor_res[[11]] <- cor.test(inat_bfly_dev$dev10_inat[inat_bfly_dev$code == "RL"], inat_bfly_dev$dev10_bfly[inat_bfly_dev$code == "RL"])
cor_res[[12]] <- cor.test(inat_bfly_dev$dev50_inat[inat_bfly_dev$code == "RL"], inat_bfly_dev$dev50_bfly[inat_bfly_dev$code == "RL"])

cor_res[[13]] <- cor.test(inat_bfly_dev$dev10_inat[inat_bfly_dev$code == "RP"], inat_bfly_dev$dev10_bfly[inat_bfly_dev$code == "RP"])
cor_res[[14]] <- cor.test(inat_bfly_dev$dev50_inat[inat_bfly_dev$code == "RP"], inat_bfly_dev$dev50_bfly[inat_bfly_dev$code == "RP"])

cor_df <- map_dfr(cor_res, ~{
  res <- .
  
  data.frame(data = res$data.name,
             r_est = res$estimate,
             conf_lo = res$conf.int[1],
             conf_hi = res$conf.int[2],
             p_val = res$p.value)
})

write.csv(cor_df, "data/derived_data/relative_pheno_corr_table.csv", row.names = F)

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
  annotate(geom = "text", x = -3, y = 50, label = "Deviance in 10% date", size = 7) +
  labs(x = "Caterpillars Count!", y = "iNaturalist caterpillars")

for_inat_adult10 <- ggplot(filter(for_inat_bfly_dev, !is.na(code)), aes(x = dev10_bfly, y = dev10_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-25,50)+
  labs(x = "Adult butterflies", y = "iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = "none")

for_cc_adult10 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(y = dev10_adult, x = dev10_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Adult butterflies", x = "Caterpillars Count!", col = "Adult overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.85, 0.2))

plot_grid(for_inat_cc10, for_inat_adult10, for_cc_adult10, ncol = 2)
ggsave("figures/relative_adult_inat_cc_10_forest.pdf", units = "in", height = 8, width = 10)

# 50%
for_inat_cc50 <- ggplot(for_quant_dev, aes(x = dev50_cc, y = dev50)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  annotate(geom = "text", x = -3, y = 50, label = "Deviance in 50% date", size = 7) +
  labs(x = "Caterpillars Count!", y = "iNaturalist caterpillars")
 
for_inat_adult50 <- ggplot(filter(for_inat_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-20, 30) +
  labs(x = "Adult butterflies", y = "iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = "none")

for_cc_adult50 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(y = dev50_adult, x = dev50_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Adult butterflies", x = "Caterpillars Count!", col = "Adult overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.85, 0.15), legend.background = element_rect(fill = "transparent"))

plot_grid(for_inat_cc50, for_inat_adult50, for_cc_adult50, ncol = 2)
ggsave("figures/relative_adult_inat_cc_50_forest.pdf", units = "in", height = 8, width = 10)

## iNat only data for adults

adult_bfly_inat <- read_csv("data/derived_data/adult_bfly_phenometrics_iNatOnlyData.csv")

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
  annotate(geom = "text", x = -5, y = 50, label = "Deviance in 10% date", size = 7) +
  labs(x = "Adult butterflies", y = "iNaturalist caterpillars") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = "none")

inat_only_50 <- ggplot(filter(inat_only_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  annotate(geom = "text", x = -5, y = 50, label = "Deviance in 50% date", size = 7) +
  xlim(-20, 20) +
  labs(x = "Adult butterflies", y = "iNaturalist caterpillars", col = "Overwinter") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.85, 0.15))

plot_grid(inat_only_10, inat_only_50)
ggsave("figures/relative_10_50_inat_only_.pdf", units = "in", height = 5, width = 10)

### Caterpillars Count! Obs excluded from iNat Cats

#### Absolute comparisons: lag predicted by GDD and time ####

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
  annotate(geom = "text", x = -10, y = 0.017, label = "Lag in 10% date", size = 7) +
  geom_density(fill = "gray96") + labs(x = "iNaturalist caterpillars - Caterpillars Count!",  y = "Density")

inat_adult_plot <- ggplot(inat_bfly_diff, aes(x = diff_10_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "iNaturalist caterpillars - Adult butterflies", fill = "Adults overwinter as", y = "") +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "eggs", "RL" = "larvae", "RP" = "pupae"))

adult_cc_plot <- ggplot(cc_bfly_diff, aes(x = diff_10_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "Caterpillars Count! - Adult butterflies", y = "Density") +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae"))

legend <- get_legend(inat_adult_plot)

plot_grid(inat_cc_plot, inat_adult_plot + theme(legend.position = "none"), adult_cc_plot + theme(legend.position = "none"), legend,
          nrow = 2, labels = c("A", "B", "C"), label_size = 15)
ggsave("figures/fig4_lag10_days_density.pdf", units = "in", height = 8, width = 10)

## Density plots for lags in days - 50%
inat_cc_plot <- ggplot(inat_cc_diff, aes(x = diff_50_days)) +
  annotate(geom = "text", x = 40, y = 0.027, label = "Lag in 50% date", size = 7) +
  geom_density(fill = "gray96") + labs(x = "iNaturalist caterpillars - Caterpillars Count!",  y = "Density")

inat_adult_plot <- ggplot(inat_bfly_diff, aes(x = diff_50_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "iNaturalist caterpillars - Adult butterflies", fill = "Adults overwinter as", y = "") +
  scale_fill_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae"))

adult_cc_plot <- ggplot(cc_bfly_diff, aes(x = diff_50_days, fill = code)) +
  geom_density(alpha = 0.5) + labs(x = "Caterpillars Count! - Adult butterflies", y = "Density") +
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

#### Figure 6 ####
## Model: lag ~ latitude + temp + (overwintering?)

load("data/cpc.tempC.2010-2020.RData")

hex_temps <- cpc %>%
  filter(doy >=172 & doy <= 266) %>%
  group_by(hex_cell, year) %>%
  summarize(mean_tmin = mean(t.min, na.rm = T),
            mean_tmax = mean(t.max, na.rm = T),
            mean_temp = mean(mean_tmin, mean_tmax, na.rm = T))

temp_dev <- hex_temps %>%
  group_by(hex_cell) %>%
  mutate(hex_mean = mean(mean_temp, na.rm = T),
         temp_dev = hex_mean - mean_temp)

inat_cc_mod <- inat_cc_sf %>%
  left_join(temp_dev, by = c("cell" = "hex_cell", "Year" = "year")) %>%
  mutate(dataset = "iNat - CC!") %>%
  group_by(dataset) %>%
  nest()

inat_bfly_mod <- inat_bfly_sf %>%
  left_join(temp_dev, by = c("cell" = "hex_cell", "year")) %>%
  mutate(dataset = "iNat - Bfly") %>%
  group_by(dataset, code) %>%
  filter(!is.na(code)) %>%
  nest()

cc_bfly_mod <- cc_bfly_sf %>%
  left_join(temp_dev, by = c("cell" = "hex_cell", "Year" = "year")) %>%
  mutate(dataset = "CC! - Bfly") %>%
  group_by(dataset, code) %>%
  filter(!is.na(code)) %>%
  nest()

mod_all <- inat_cc_mod %>%
  rbind(inat_bfly_mod, cc_bfly_mod) %>%
  mutate(mod10 = map(data, ~lm(abs(diff_10_days) ~ latitude + temp_dev, data = .)),
         mod50 = map(data, ~lm(abs(diff_50_days) ~ latitude + temp_dev, data = .)),
         tidy10 = map(mod10, ~broom::tidy(.)),
         tidy50 = map(mod50, ~broom::tidy(.)))

mod_ests_10 <- mod_all %>%
  select(dataset, code, tidy10) %>%
  unnest(cols = c("tidy10")) %>%
  mutate_at(c("code"), ~ifelse(is.na(.), "None", .)) %>%
  filter(!(dataset == "CC! - Bfly" & code == "RE")) %>%
  filter(term != "(Intercept)")

temp10 <- ggplot(filter(mod_ests_10, term == "temp_dev"), aes(x = term, y = estimate, col = dataset, shape = code)) + 
  geom_point(cex = 4, position = position_dodge(width = 0.6)) + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), 
                                                          width = 0.3, cex = 1, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = "Estimate", col = "Lag in 10% date", shape = "Adults overwinter as") +
  scale_x_discrete(labels = c("temp_dev" = "Deviance in\nspring temp", "latitude" = "Latitude")) +
  scale_color_manual(values = c("#3E4A89FF", 
                                "#1F9E89FF",
                                "#6DCD59FF"),
                     labels = c("CC! - Bfly" = "Caterpillars Count! - Adult butterflies",
                                "iNat - Bfly" = "iNaturalist caterpillars - adult butterflies",
                                "iNat - CC!" = "iNaturalist caterpillars - Caterpillars Count!")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  coord_flip()

lat10 <- ggplot(filter(mod_ests_10, term == "latitude"), aes(x = term, y = estimate, col = dataset, shape = code)) + 
  geom_point(cex = 4, position = position_dodge(width = 0.6)) + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), 
                width = 0.3, cex = 1, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = "Estimate", col = "Lag 10% date model", shape = "Adults overwinter as") +
  scale_x_discrete(labels = c("mean_temp" = "Spring temp", "latitude" = "      Latitude")) +
  scale_color_manual(values = c("#3E4A89FF", 
                                "#1F9E89FF",
                                "#6DCD59FF"),
                     labels = c("CC! - Bfly" = "Caterpillars Count! - Adult butterflies",
                                "iNat - Bfly" = "iNaturalist caterpillars - adult butterflies",
                                "iNat - CC!" = "iNaturalist caterpillars - Caterpillars Count!")) +
  scale_shape_manual(values = c(15, 16, 17, 18), labels = c("None" = "none", "RL" = "larvae",
                                                            "RE" = "eggs", "RP" = "pupae")) +
  coord_flip()

ests_legend <- get_legend(lat10)

ests_multi <- plot_grid(temp10 + theme(legend.position = "none"), 
                        lat10 + theme(legend.position = "none"), nrow = 2)
plot_grid(ests_multi, ests_legend, ncol = 2, rel_widths = c(0.5, 0.5))

ggsave("figures/fig5_lag10_mod_ests.pdf", units = "in", height = 5, width = 10)

mod_ests_50 <- mod_all %>%
  select(dataset, code, tidy50) %>%
  unnest(cols = c("tidy50")) %>%
  mutate_at(c("code"), ~ifelse(is.na(.), "None", .)) %>%
  filter(term != "(Intercept)")

temp50 <- ggplot(filter(mod_ests_50, term == "mean_temp"), aes(x = term, y = estimate, col = dataset, shape = code)) + 
  geom_point(cex = 4, position = position_dodge(width = 0.6)) + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), 
                width = 0.3, cex = 1, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = "Estimate", col = "Lag in 50% date", shape = "Adults overwinter as") +
  scale_x_discrete(labels = c("mean_temp" = "Spring temp", "latitude" = "Latitude")) +
  scale_color_manual(values = c("#3E4A89FF", 
                                "#1F9E89FF",
                                "#6DCD59FF"),
                     labels = c("CC! - Bfly" = "Caterpillars Count! - Adult butterflies",
                                "iNat - Bfly" = "iNaturalist caterpillars - adult butterflies",
                                "iNat - CC!" = "iNaturalist caterpillars - Caterpillars Count!")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  coord_flip()

lat50 <- ggplot(filter(mod_ests_50, term == "latitude"), aes(x = term, y = estimate, col = dataset, shape = code)) + 
  geom_point(cex = 4, position = position_dodge(width = 0.6)) + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error), 
                width = 0.3, cex = 1, position = position_dodge(width = 0.6)) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "", y = "Estimate", col = "Lag 50% date model", shape = "Adults overwinter as") +
  scale_x_discrete(labels = c("mean_temp" = "Spring temp", "latitude" = "      Latitude")) +
  scale_color_manual(values = c("#3E4A89FF", 
                                "#1F9E89FF",
                                "#6DCD59FF"),
                     labels = c("CC! - Bfly" = "Caterpillars Count! - Adult butterflies",
                                "iNat - Bfly" = "iNaturalist caterpillars - adult butterflies",
                                "iNat - CC!" = "iNaturalist caterpillars - Caterpillars Count!")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  coord_flip()

ests_legend <- get_legend(lat50)

ests_multi <- plot_grid(temp50 + theme(legend.position = "none"), 
                        lat50 + theme(legend.position = "none"), nrow = 2)
plot_grid(ests_multi, ests_legend, ncol = 2, rel_widths = c(0.5, 0.5))

ggsave("figures/lag50_mod_ests.pdf", units = "in", height = 5, width = 10)

## Summary model table

mod_summary <- mod_all %>%
  mutate(r2_10 = map_dbl(mod10, ~summary(.)$r.squared),
         r2_50 = map_dbl(mod50, ~summary(.)$r.squared)) %>%
  dplyr::select(-data, -mod10, -mod50) %>%
  pivot_longer(tidy10:tidy50, names_to = "pct", values_to = "data") %>%
  pivot_longer(r2_10:r2_50, names_to = "pct_r2", values_to = "r2") %>%
  filter(pct == "tidy10" & pct_r2 == "r2_10" | pct == "tidy50" & pct_r2 == "r2_50") %>%
  unnest(cols = c("data")) %>%
  select(-pct_r2) %>%
  distinct()
write.csv(mod_summary, "data/derived_data/lag_mods.csv", row.names = F)

#### Figure 5: Temperature deviance relationships ####

temp_dev <- hex_temps %>%
  group_by(hex_cell) %>%
  mutate(hex_mean = mean(mean_temp, na.rm = T),
         temp_dev = hex_mean - mean_temp)

pheno_dev_temp <- quant_dev %>%
  left_join(temp_dev, by = c("cell" = "hex_cell", "Year" = "year"))

inat_dev_temp <- inat_cats_dev %>%
  left_join(temp_dev, by = c("HEXcell" = "hex_cell", "year"))

adult_bfly_dev_temp <- adult_bfly_dev %>%
  left_join(temp_dev, by = c("HEXcell" = "hex_cell", "year"))

cc_temp <- ggplot(pheno_dev_temp, aes(x = temp_dev, y = dev10_cc)) + geom_point() + geom_smooth(method = "lm", se = F) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Spring temperature deviance (°C)", y = "Caterpillars Count!", title = "Deviance in 10% date")

inat_temp <- ggplot(inat_dev_temp, aes(x = temp_dev, y = dev10)) + geom_point() + geom_smooth(method = "lm", se = F) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Spring temperature deviance (°C)", y = "iNaturalist caterpillars")

adult_temp <- ggplot(adult_bfly_dev_temp, aes(x = temp_dev, y = dev10, col = code)) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() + geom_smooth(method = "lm", se = F) + 
  labs(x = "Spring temperature deviance (°C)", y = "Adult butterflies", col = "Adults overwinter as") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "eggs", "RL" = "larvae", "RP" = "pupae")) +
  theme(legend.position = c(0.25, 0.85), legend.background = element_rect(fill = "transparent"))

plot_grid(cc_temp, inat_temp, adult_temp, nrow = 2, labels = c("A", "B", "C"))
ggsave("figures/fig6_temp_deviance.pdf", height = 8, width = 10)

## Summary of models

summary(lm(dev10_cc ~ temp_dev, pheno_dev_temp))

summary(lm(dev10 ~ temp_dev, inat_dev_temp))

summary(lm(dev10 ~ temp_dev, data = filter(adult_bfly_dev_temp, code == "RP")))
summary(lm(dev10 ~ temp_dev, data = filter(adult_bfly_dev_temp, code == "RL")))

## Temp deviance relationships, 50% date
cc_temp <- ggplot(pheno_dev_temp, aes(x = temp_dev, y = dev50_cc)) + geom_point() + geom_smooth(method = "lm", se = F) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Spring temperature deviance (°C)", y = "Caterpillars Count!", title = "Deviance in 50% date")

inat_temp <- ggplot(inat_dev_temp, aes(x = temp_dev, y = dev50)) + geom_point() + geom_smooth(method = "lm", se = F) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Spring temperature deviance (°C)", y = "iNaturalist caterpillars")

adult_temp <- ggplot(adult_bfly_dev_temp, aes(x = temp_dev, y = dev50, col = code)) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() + geom_smooth(method = "lm", se = F) + 
  labs(x = "Spring temperature deviance (°C)", y = "Adult butterflies", col = "Adults overwinter as") +
  scale_color_manual(values = c("#482878FF", "#26828EFF", "#B4DE2CFF"), labels = c("RE" = "Eggs", "RL" = "Larvae", "RP" = "Pupae")) +
  theme(legend.position = c(0.25, 0.25), legend.background = element_rect(fill = "transparent"))

plot_grid(cc_temp, inat_temp, adult_temp, nrow = 2, labels = c("A", "B", "C"))
ggsave("figures/temp_deviance_50.pdf", height = 8, width = 10)

