# trait-based-cat-cast project
This repo is part of the Macrosystems Phenological Mismatch project. Here we are using species-specific data and traits to identify appropriate proxies for forest caterpillar abundance across space and time, to examine relationships with bird demography.

Repository structure:
<br>
data/ - Datasets relevant for project
<ul>
  <li>CatCountSpeciesList.csv - species names in Caterpillars Count! data</li>
  <li>DiscLifeTraits.csv - Species names in discover life moth data and associated traits</li>
  <li>caterpillars-of-eastern-north-america.csv - opportunistic occurrence data file for iNaturalist caterpillars of ENA projecy</li>
  <li>discoverlife_3sp_ma.csv - semi-standard survey occurrence data for 3 moth species in Massachusetts, discover life adult data</li>
  <li>lep.names.query.csv - Standardized (Genus_Species) names and variations (altname) seen in datasets  with boolean ID for whether data are available (with.data)</li>
  <li>lep.trait.csv - File using standarized species names, for taxonomic heirarchy and traits</li>
  <li>ma-bfly.csv - Field Trip occurrence data from the Massachusetts butterfly club</li>
  <li>regions.xlsx - information on preliminary regions being used to examine datasets</li>
</ul>

code/  - code (currently all R) 
<ul>
  <li>bestLag.R - not needed</li>
  <li>overwintering_lists.R - managineg trait data for overwinter stage</li>
  <li>pheno_curve_lag.R - R code aligning phenology curves from GAMs or smoothed density surves and identifying the best lags between adults and larvae</li>
  <li>gam_fx.R - File for R functions for GAM phenology estimation and adult-larval Lag comparison</li>
</ul>
