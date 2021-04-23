### Relative and absolute comparisons of lep phenology

library(tidyverse)
library(phenesse)
library(lubridate)
library(sf)
library(tmap)
library(cowplot)

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

cc_pheno_unnest <- read_csv("data/derived_data/cc_weibull_10_50.csv")

### Deviations in 10% and 50% phenometrics
## Use only sites within hex with >= 2 years, make sure in cells w/ multiple years, not flickering betw sites

good_sites <- cc_pheno_unnest %>%
  group_by(cell, Name) %>%
  summarize(n_year = n_distinct(Year)) %>%
  filter(n_year >= 2) %>%
  filter(Name != "UNC Chapel Hill Campus")

## Compare deviations in cc 10/50 with inat cats and adult bfly

inat_cats <- read_csv("data/derived_data/caterpillar_inat_phenometrics_phenesse.csv")

inat_cats_dev <- inat_cats %>%
  group_by(HEXcell, code) %>%
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

## Plot correlations of 10, 50% deviances

quant_dev <- cc_dev %>%
  left_join(select(adult_bfly_dev, year, HEXcell, code, dev10, dev50), by = c("Year" = "year", "cell" = "HEXcell"), suffix = c("_cc", "_adult")) %>%
  left_join(select(inat_cats_dev, year, HEXcell, code, dev10, dev50), by = c("Year" = "year", "cell" = "HEXcell", "code"))

inat_bfly_dev <- select(adult_bfly_dev, year, HEXcell, code, dev10, dev50) %>%
  left_join(inat_cats_dev, by = c("year", "HEXcell", "code"), suffix = c("_bfly", "_inat"))

# 10%
inat_cc10 <- ggplot(filter(quant_dev, !is.na(code)), aes(x = dev10_cc, y = dev10, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 10% Caterpillars Count!", y = "Deviance 10% iNaturalist caterpillars", col = "Overwinter") +
  theme(legend.position = c(0.8, 0.2))

inat_adult10 <- ggplot(filter(inat_bfly_dev, !is.na(code)), aes(x = dev10_bfly, y = dev10_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-40, 40) +
  labs(x = "Deviance 10% Adult butterflies", y = "Deviance 10% iNaturalist caterpillars") +
  theme(legend.position = "none")

cc_adult10 <- ggplot(filter(quant_dev, !is.na(code)), aes(y = dev10_adult, x = dev10_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 10% Adult butterflies", x = "Deviance 10% Caterpillars Count!") +
  theme(legend.position = "none")

plot_grid(inat_cc10, inat_adult10, cc_adult10, ncol = 2)
ggsave("figures/relative_adult_inat_cc_10.pdf", units = "in", height = 8, width = 10)

# 50%
inat_cc50 <- ggplot(filter(quant_dev, !is.na(code)), aes(x = dev50_cc, y = dev50, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 50% Caterpillars Count!", y = "Deviance 50% iNaturalist caterpillars", col = "Overwinter") +
  theme(legend.position = c(0.85, 0.15), legend.background = element_rect(fill = "transparent"))

inat_adult50 <- ggplot(filter(inat_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-30, 30) +
  labs(x = "Deviance 50% Adult butterflies", y = "Deviance 50% iNaturalist caterpillars") +
  theme(legend.position = "none")

cc_adult50 <- ggplot(filter(quant_dev, !is.na(code)), aes(y = dev50_adult, x = dev50_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 50% Adult butterflies", x = "Deviance 50% Caterpillars Count!") +
  theme(legend.position = "none")

plot_grid(inat_cc50, inat_adult50, cc_adult50, ncol = 2)
ggsave("figures/relative_adult_inat_cc_50.pdf", units = "in", height = 8, width = 10)

## Forest only subset

for_inat_cats <- read_csv("data/derived_data/ForestOnly_caterpillar_inat_phenometrics_phenesse.csv")

for_inat_cats_dev <- for_inat_cats %>%
  group_by(HEXcell, code) %>%
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
  left_join(select(for_inat_cats_dev, year, HEXcell, code, dev10, dev50), by = c("Year" = "year", "cell" = "HEXcell", "code"))

for_inat_bfly_dev <- select(for_adult_bfly_dev, year, HEXcell, code, dev10, dev50) %>%
  left_join(for_inat_cats_dev, by = c("year", "HEXcell", "code"), suffix = c("_bfly", "_inat"))

# 10%
for_inat_cc10 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(x = dev10_cc, y = dev10, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 10% Caterpillars Count!", y = "Deviance 10% iNaturalist caterpillars", col = "Overwinter") +
  theme(legend.position = c(0.85, 0.2))

for_inat_adult10 <- ggplot(filter(for_inat_bfly_dev, !is.na(code)), aes(x = dev10_bfly, y = dev10_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-25,50)+
  labs(x = "Deviance 10% Adult butterflies", y = "Deviance 10% iNaturalist caterpillars") +
  theme(legend.position = "none")

for_cc_adult10 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(y = dev10_adult, x = dev10_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 10% Adult butterflies", x = "Deviance 10% Caterpillars Count!") +
  theme(legend.position = "none")

plot_grid(for_inat_cc10, for_inat_adult10, for_cc_adult10, ncol = 2, labels = c("Forest only"))
ggsave("figures/relative_adult_inat_cc_10_forest.pdf", units = "in", height = 8, width = 10)

# 50%
for_inat_cc50 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(x = dev50_cc, y = dev50, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Deviance 50% Caterpillars Count!", y = "Deviance 50% iNaturalist caterpillars", col = "Overwinter") +
  theme(legend.position = c(0.85, 0.15), legend.background = element_rect(fill = "transparent"))

for_inat_adult50 <- ggplot(filter(for_inat_bfly_dev, !is.na(code)), aes(x = dev50_bfly, y = dev50_inat, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  xlim(-20, 30) +
  labs(x = "Deviance 50% Adult butterflies", y = "Deviance 50% iNaturalist caterpillars") +
  theme(legend.position = "none")

for_cc_adult50 <- ggplot(filter(for_quant_dev, !is.na(code)), aes(y = dev50_adult, x = dev50_cc, col = code)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "lm", se = F) +
  labs(y = "Deviance 50% Adult butterflies", x = "Deviance 50% Caterpillars Count!") +
  theme(legend.position = "none")

plot_grid(for_inat_cc50, for_inat_adult50, for_cc_adult50, ncol = 2)
ggsave("figures/relative_adult_inat_cc_50_forest.pdf", units = "in", height = 8, width = 10)

### Absolute comparisons: lag predicted by GDD and time
