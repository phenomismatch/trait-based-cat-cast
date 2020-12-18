# trait-based-cat-cast project
This repo is part of the Macrosystems Phenological Mismatch project. Here we are using species-specific data and traits to identify appropriate proxies for forest caterpillar abundance across space and time, to examine relationships with bird demography.

Repository structure:

data/ - Datasets relevant for project

CatCountSpeciesList.csv - species names in Caterpillars Count! data
DiscLifeTraits.csv - Species names in discover life moth data and associated traits
caterpillars-of-eastern-north-america.csv - opportunistic occurrence data file for iNaturalist caterpillars of ENA projecy
discoverlife_3sp_ma.csv - semi-standard survey occurrence data for 3 moth species in Massachusetts, discoverf life adult data
lep.names.query.csv - Standardized (Genus_Species) names and variations (altname) seen in datasets  with boolean ID for whether data are available (with.data)
lep.trait.csv - File using standarized species names, for taxonomic heirarchy and traits
ma-bfly.csv - Field Trip occurrence data from the Massachusetts butterfly club
regions.xlsx - information on preliminary regions being used to examine datasets

code/  - code (currently all R) 
bestLag.R - not needed
overwintering_lists.R - managineg trait data for overwinter stage
pheno_curve_lag.R - R code aligning phenology curves from GAMs or smoothed density surves and identifying the best lags between adults and larvae
gam_fx.R - File for R functions for GAM phenology estimation and adult-larval Lag comparison
