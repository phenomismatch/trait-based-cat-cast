# Function for calculating best fit lag between adult and larval phenology

# Need to allow for either abundance (e.g. MA Butterfly Survey) or occurrence (e.g. iNaturalist) data for adults,
# which will determine whether a GAM or kernel density estimator is fit to the data, respectively.

# For caterpillars, also need to allow for abundance (e.g. Hubbard Brook, Coweeta, Caterpillars Count!) or occurrence
# (e.g. iNaturalist) data, but in either case, kernel density estimator will be used.

bestLag = function(sci_name,              # scientific name
                   year,                  # year of comparison
                   doy.scope = c(1, 366), # day of year scope for estimating phenology
                   larvaDF,               # dataframe of larval occurrences with columns: (code, common_name, or scientific_name), year, doy, count
                   adultDF,               # dataframe of adult occurrences with columns: (code, common_name, or scientific_name), year, doy, count,
                                          # and optionally a SITE_ID field (as in MA Butterfly data)
                   minNumLarvalRecs = 3,  # minimum number of larval records
                   lagRange = -60:60,     # vector of lags to evaluate
                   multipleSites = FALSE  # are count data from multiple sites over which phenology is to be estimated? 
                                          # If TRUE, then expects SITE_ID field in adultDF which is used in hierarchical model
                   ) {
  
  if (!multipleSites) {
    adultDF$SITE_ID = 1
  }
  
  if (!'SITE_ID' %in% names(adultDF)) {
    warning("The dataframe 'adultDF' must have a SITE_ID field if multipleSites = TRUE.")
    return(NULL)
  }
  
  # Calculate adult counts by day
  #     Need to use specific columns in ALL CAPS for the flight_curve function
  
  adult.spi <- adultDF %>%
    filter(year == year) %>%
    mutate(CT = ifelse(scientific_name == sci_name, count, 0)) %>%
    group_by(SITE_ID, doy) %>%
    arrange(doy) %>%
    summarize(SPECIES = sci_name, COUNT = sum(CT))
  
  adult.spi.phen <- flight_curve(adult.spi)
  
  # Calculate larval phenology
  larval.i <- filter(larval.data, scientific_name == sci_name, year == year)
  #verify at least 3 larval records, then: calculate smoothed data density 
  if(nrow(larval.i) >= minNumLarvalRecs) {
    temp.larval<-density(larval.i$doy, n = length(adult.spi.phen), from = 1, to = length(adult.spi.phen))      
    temp.larval<-data.frame(doy = doy.scope[1:length(adult.spi.phen)], 
                            larv = round(temp.larval$y[1:length(adult.spi.phen)]*100,3))
  } else {
    warning("Not enough larval records.")
    return(NULL)
  }
  
}