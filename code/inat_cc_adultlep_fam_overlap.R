## Jaccard similarity of caterpillar and adult lep datasets

library(tidyverse)

## iNat caterpillar species
inat_spp <- read_csv("data/derived_data/total_caterpillar_obs_byHEX.csv") %>%
  filter(year == 2018) %>%
  group_by(taxon_family_name) %>%
  summarize(nObs = n()) %>%
  filter(!(is.na(taxon_family_name)))

## Adult lep species
naba_opp <- read_csv("data/derived_data/NABA.taxonomy.csv") %>%
  filter(ObsYear == 2018) %>%
  group_by(Family) %>%
  summarize(nObs = sum(nobs)) %>%
  filter(!is.na(Family))

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
sum(inat_spp$taxon_family_name %in% naba_opp$Family)/length(unique(c(naba_opp$Family, inat_spp$taxon_family_name)))
#0.1
# 5 overlap

## iNat & cc
sum(inat_spp$taxon_family_name %in% cc_spp$taxon_family_name)/length(unique(c(cc_spp$taxon_family_name, inat_spp$taxon_family_name)))
# 0.204
# 10 overlap

## CC & adult lep
sum(cc_spp$taxon_family_name %in% naba_opp$Family)/length(unique(c(naba_opp$Family, cc_spp$taxon_family_name)))
# 0.067
# 1 overlap - Papilionidae

# Total families
length(unique(c(naba_opp$Family, cc_spp$taxon_family_name, inat_spp$taxon_family_name)))
# 50

