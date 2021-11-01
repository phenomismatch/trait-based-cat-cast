## Jaccard similarity of caterpillar and adult lep datasets

library(tidyverse)

## iNat caterpillar species
inat_spp <- read_csv("data/derived_data/total_caterpillar_obs_byHEX.csv") %>%
  filter(year == 2018) %>%
  group_by(taxon_family_name) %>%
  summarize(nObs = n()) %>%
  filter(!(is.na(taxon_family_name)))

## Adult lep species
naba_opp <- read_csv("data/derived_data/inat_adultFamilies.csv") %>%
  filter(year == 2018) %>%
  group_by(family) %>%
  summarize(nObs = sum(observations)) %>%
  filter(!is.na(family))

## Caterpillars Count! species
cc_spp <- read_csv("data/derived_data/total_caterpillar_obs_byHEX.csv") %>%
  filter(user_login == "caterpillarscount") %>%
  filter(cell %in% c(593, 595, 702, 703)) %>%
  filter(year == 2018) %>%
  group_by(taxon_family_name) %>%
  summarize(nObs = n()) %>%
  filter(!(is.na(taxon_family_name)))

## Jaccard indices of species lists

## iNat & adult leps
sum(inat_spp$taxon_family_name %in% naba_opp$family)/length(unique(c(naba_opp$family, inat_spp$taxon_family_name)))
#0.1
# 5 overlap

## iNat & cc
sum(inat_spp$taxon_family_name %in% cc_spp$taxon_family_name)/length(unique(c(cc_spp$taxon_family_name, inat_spp$taxon_family_name)))
# 0.204
# 10 overlap

## CC & adult lep
sum(cc_spp$taxon_family_name %in% naba_opp$family)/length(unique(c(naba_opp$family, cc_spp$taxon_family_name)))
# 0.067
# 1 overlap - Papilionidae

# Total families
length(unique(c(naba_opp$family, cc_spp$taxon_family_name, inat_spp$taxon_family_name)))
# 50

## Summary table of families x dataset
fam_table <- data.frame(family = unique(c(naba_opp$family, cc_spp$taxon_family_name, inat_spp$taxon_family_name))) %>%
  mutate(cc = ifelse(family %in% cc_spp$taxon_family_name, "X", " "),
         inat_cat = ifelse(family %in% inat_spp$taxon_family_name, "X", " "),
         inat_adult = ifelse(family %in% naba_opp$family, "X", " "))
