# Function for calculating best fit lag between adult and larval phenology

bestLag = function(speciesID,             # species identifier
                   IDtype = 'code',       # possible values are 'code' (4-letter butterfly code), 'common_name', or 'scientific_name'
                   year,                  # year of comparison
                   doy.scope = c(1, 366), # day of year scope for estimating phenology
                   larvaDF,               # dataframe of larval occurrences with columns: (code, common_name, or scientific_name), year, doy, count
                   adultDF,               # dataframe of adult occurrences with columns: (code, common_name, or scientific_name), year, doy, count
                   minNumLarvalRecs = 3   # minimum number of larval records
                   ) {
  
  
  
  # Calculate adult counts by day
  adult.spi <- adultDF %>%
    filter(year == year) %>%
    mutate(CT = ifelse(get(IDtype) == speciesID, count, 0)) %>%
    group_by(SITE_ID, doy) %>%
    arrange(doy) %>%
    summarize(Species = spi, Count = sum(CT))
  
  adult.spi.phen <- flightcurve(adult.spi)
  
  # Calculate larval phenology
  larval.i <- filter(larval.data, get(IDtype) == speciesID, year == year)
  ## Larval data is currently occurrence data
  if(larval.format=="occurrence") {
    #verify at least 3 larval records, then: calculate smoothed data density 
    if(nrow(larval.i) >= minNumLarvalRecs) {
      temp.larval<-density(larval.i$doy,n = length(adult.spi.phen), from=1, to=length(adult.spi.phen))      
      temp.larval<-data.frame(doy = doy.scope[1:length(adult.spi.phen)], 
                              larv = round(temp.larval$y[1:length(adult.spi.phen)]*100,3))
    }
  }
  
}