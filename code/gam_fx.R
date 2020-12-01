#phenology & lag functions for survey and incidental data

#if(!requireNamespace("devtools")) install.packages("devtools")
#devtools::install_github("RetoSchmucki/rbms")
require(rbms)

#butterfly_day function estimates annual abundance at specific site(s) 
#butterfly_day function input is site_year_sp_count from GAM flight curve phenology
butterfly_day <- function(site_year_sp_count, WeekCount = TRUE){
  if (WeekCount == TRUE){
    b_day <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0 & WEEK_DAY == 4, FITTED, by = .(SITE_ID, M_YEAR, WEEK)][,sum(FITTED), by = .(SITE_ID, M_YEAR)]
  } else {
    b_day <- site_year_sp_count$impute_count[COMPLT_SEASON == 1 & M_SEASON != 0, FITTED, by = .(SITE_ID, M_YEAR, WEEK)][,sum(FITTED), by = .(SITE_ID, M_YEAR)]
  }
  data.table::setnames(b_day,"V1","BUTTERFLY_DAY")
  return(b_day)
}


#phenoFlight function estimates flight curve phenology & returns specific phenometrics
#Phenometrics currently hardcoded but should be functionalized
#Assumes abundance data from repeated surveys to specific site(s)
#Assumes phenology is consistent across sites, abundance varies
fcPhenometrics<-function(sites, visit, cnt, indices) {
  #format visit and count data
  sites<-sites[indices]
  visit<-visit[visit$DATE%in%sites,]
  cnt<-cnt[cnt$DATE%in%sites,] #cnt<-count1
  if(nrow(cnt)>5) {
  ts_date <- rbms::ts_dwmy_table(InitYear = as.numeric(format(visit$DATE[1],"%Y")), LastYear = as.numeric(format(visit$DATE[1],"%Y")), WeekDay1 = 'monday')

  ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 2, EndMonth = 10, StartDay = 1, EndDay = NULL,
                               CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 7)
  ts_season_visit <- rbms::ts_monit_site(visit, ts_season)

  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, cnt, sp = cnt$SPECIES[1])

  ## compute the flight curve, using the
  ## regionalGAM method
  ##=========================================
  dataset_y <- ts_season_count[, .(SPECIES, SITE_ID, DATE, WEEK, WEEK_DAY, DAY_SINCE, M_YEAR, M_SEASON, COUNT, ANCHOR)]
  dataset_y[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1]
  try(
  system.time(
    #ts_flight_curve <- fit_gam(dataset_y, NbrSample = 200, GamFamily = 'nb', MaxTrial = 3)

    ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 200, MinVisit = 4, MinOccur = 4, MinNbrSite = 1,
                                  MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)
  )   )

  #check if successful
  if(is.na(max(ts_flight_curve$pheno[, NM]))){ }else {

   test1<-ts_flight_curve$pheno %>%
      group_by(M_YEAR) %>%
      mutate(bd=cumsum(NM))

    ##onset
    #test1[test1$bd>0,c(3,4,15,16,19)][1,]
    #10% curve
    p10<-test1[test1$bd>0.09999,]$trimDAYNO[1]
    #50% curve
    p50<-test1[test1$bd>0.49999,]$trimDAYNO[1]
    return(c(p10,p50))
  } #end if regional GAM was successful
  } else {
    p10<-0
    p50<-0
    return(c(p10,p50))
  }
} # end if data density sufficient


#fcPhenoboot function also estimates flight curve phenometrics from abundance data
#assumes repeated surveys at specific site(s)
#allows bootstrapping
fcPhenoboot<-function(index, visit, cnt, indices) {
  #format visit and count data
  ind<-index[indices]
  visit<-visit[visit$DATE%in%ind,]
  cnt<-cnt[cnt$DATE%in%ind,] #cnt<-count1
  if(nrow(cnt)>5) {
    ts_date <- ts_dwmy_table(InitYear = as.numeric(format(visit$DATE[1],"%Y")), LastYear = as.numeric(format(visit$DATE[1],"%Y")), WeekDay1 = 'monday')

    ts_season <- ts_monit_season(ts_date, StartMonth = 2, EndMonth = 10, StartDay = 1, EndDay = NULL,
                                 CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 7)
    ts_season_visit <- ts_monit_site(visit, ts_season)

    ts_season_count <- ts_monit_count_site(ts_season_visit, cnt, sp = cnt$SPECIES[1])

    ## compute the flight curve, using the
    ## regionalGAM method
    ##=========================================
    dataset_y <- ts_season_count[, .(SPECIES, SITE_ID, DATE, WEEK, WEEK_DAY, DAY_SINCE, M_YEAR, M_SEASON, COUNT, ANCHOR)]
    dataset_y[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1]
    try(
      system.time(
        ts_flight_curve <- flight_curve(ts_season_count, NbrSample = 200, MinVisit = 4, MinOccur = 3, MinNbrSite = 1,
                                        MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)
      )
    )

    #check if successful
    if(is.na(max(ts_flight_curve$pheno[, NM]))){ }else {

      test1<-ts_flight_curve$pheno %>%
        group_by(M_YEAR) %>%
        mutate(bd=cumsum(NM))

      ##onset
      #test1[test1$bd>0,c(3,4,15,16,19)][1,]
      #10% curve
      p10<-test1[test1$bd>0.09999,]$trimDAYNO[1]
      #50% curve
      p50<-test1[test1$bd>0.49999,]$trimDAYNO[1]
      return(c(p10,p50))
    } #end if regional GAM was successful
  } else {
    p10<-0
    p50<-0
    return(c(p10,p50))
  }
} # end if data density sufficient


#flightcurve funtion estimates phenology/abundance flight curve from abundance data
#flight curve assumes repeated visit to specific site(s)
#flightcurve uses a GAM model to smooth abundances over time
#phenology is assumed consistent across sites, abundance varies across sites
#flight curve returns relative abundance over time 
flightcurve<-function(mydata) {
  visit<-mydata %>%
    group_by(SITE_ID, DATE) %>%
    select(SITE_ID, DATE)
  count<-mydata %>%
    filter(COUNT>0)
  ts_date <- rbms::ts_dwmy_table(InitYear = as.numeric(format(visit$DATE[1],"%Y")), LastYear = as.numeric(format(visit$DATE[1],"%Y")), WeekDay1 = 'monday')

  ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 2, EndMonth = 10, StartDay = 1, EndDay = NULL,
                                     CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 14)
  ts_season_visit <- rbms::ts_monit_site(visit, ts_season)

  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, count, sp = count$SPECIES[1])

  ## compute the flight curve, using the
  ## regionalGAM method
  ##=========================================
  dataset_y <- ts_season_count[, .(SPECIES, SITE_ID, DATE, WEEK, WEEK_DAY, DAY_SINCE, M_YEAR, M_SEASON, COUNT, ANCHOR)]
  dataset_y[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1]

  system.time(
    ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 200, MinVisit = 2, MinOccur = 1, MinNbrSite = 1,
                                          MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)
  )

return(ts_flight_curve$pheno$NM)
}

#flightcurveBoot allows the flight curve funtion to be bootstrapped
flightcurveBoot<-function(mydata, indices) {
  mydata<-mydata[indices,]
  visit<-mydata %>%
    group_by(SITE_ID, DATE) %>%
    select(SITE_ID, DATE)
  count<-mydata %>%
    filter(COUNT>0)
  ts_date <- rbms::ts_dwmy_table(InitYear = as.numeric(format(visit$DATE[1],"%Y")), LastYear = as.numeric(format(visit$DATE[1],"%Y")), WeekDay1 = 'monday')

  ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 2, EndMonth = 10, StartDay = 1, EndDay = NULL,
                                     CompltSeason = TRUE, Anchor = TRUE, AnchorLength = 7, AnchorLag = 14)
  ts_season_visit <- rbms::ts_monit_site(visit, ts_season)

  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, count, sp = count$SPECIES[1])

  ## compute the flight curve, using the
  ## regionalGAM method
  ##=========================================
  dataset_y <- ts_season_count[, .(SPECIES, SITE_ID, DATE, WEEK, WEEK_DAY, DAY_SINCE, M_YEAR, M_SEASON, COUNT, ANCHOR)]
  dataset_y[, trimDAYNO := DAY_SINCE - min(DAY_SINCE) + 1]

  system.time(
    ts_flight_curve <- rbms::flight_curve(ts_season_count, NbrSample = 200, MinVisit = 2, MinOccur = 1, MinNbrSite = 1,
                                          MaxTrial = 3, FcMethod = 'regionalGAM', GamFamily = 'nb', SpeedGam = FALSE, CompltSeason = TRUE)
  )
  return(ts_flight_curve$pheno$NM)
}


#lagFits function calculates the adult-larval lag for a species/group by:
#1 - estimates adult phenology using a GAM 
#2 - estimates larval phenology using smoothed data density
#3 - shifts larval phenology by # days in lagRange
#4 - estimates how well phenology curves match @ each lag using a linear model
lagFits = function(sci_name,              # scientific name
                   Year,                  # year of comparison
                   doy.scope = 1:366,     # day of year scope for estimating phenology
                   larvaDF,               # dataframe of larval occurrences with columns: scientific_name, year, DATE, count
                   adultDF,               # dataframe of adult occurrences with columns: scientific_name, year, DATE, count,
                   # and optionally a SITE_ID field (as in MA Butterfly data)
                   minNumRecs = 5,        # minimum number of larval or adult records
                   lagRange = -60:60,     # vector of lags to evaluate
                   multipleSites = FALSE, # are adult count data from multiple sites over which phenology is to be estimated? 
                   # If TRUE, then expects SITE_ID field in adultDF which is used in hierarchical model
                   bestFitOnly = FALSE    # If TRUE, only return lag with the highest R2, otherwise return R2's of all lags
) {
  
  # If the data do not span multiple sites, create a constant column, SITE_ID = 1
  if (!multipleSites) {
    adultDF$SITE_ID = 1
  }
  
  if (!'SITE_ID' %in% names(adultDF)) {
    warning("The dataframe 'adultDF' must have a SITE_ID field if multipleSites = TRUE.")
    return(NULL)
  }
  
  # Calculate adult counts by day
  #   (Need to use specific columns in ALL CAPS for the flight_curve function)
  
  adult.spi <- adultDF %>%
    filter(year == Year) %>%
    mutate(CT = ifelse(scientific_name == sci_name, count, 0)) %>%
    group_by(SITE_ID, DATE) %>%
    arrange(DATE) %>%
    summarize(SPECIES = sci_name, COUNT = sum(CT))
  
  # Check that there are enough adult records (i.e. dates with count > 0)
  if (sum(adult.spi$COUNT > 0) >= minNumRecs) {
    
    adult.spi.phen <- flightcurve(adult.spi)
    temp.ad <- data.frame(doy = doy.scope[1:length(adult.spi.phen)], adult = adult.spi.phen)
    
    # Calculate larval phenology
    larval.i <- filter(larvalDF, scientific_name == sci_name, year == year)
    #verify at least 3 larval records, then: calculate smoothed data density 
    if(nrow(larval.i) >= minNumRecs) {
      
      temp.larval <- density(larval.i$doy, n = length(adult.spi.phen), from = 1, to = length(adult.spi.phen))      
      temp.larval <- data.frame(doy = doy.scope[1:length(adult.spi.phen)], 
                                larva = round(temp.larval$y[1:length(adult.spi.phen)]*100,3))
    } else {
      warning("Not enough larval records for this species in this year.")
      return(NULL)
    }
    
    ## Estimate Lag
    
    lagresult = tibble(species = sci_name, year = year, lag = lagRange, r2 = NA)
    
    #loop lags
    for(lagi in lagRange) {
      #larvae before adults
      temp.lalag <- temp.larval %>% 
        add_row(doy = (1 - lagi), larva = 0) %>% 
        mutate(doy = doy + lagi)
      
      temp.datalag <- merge(temp.ad, temp.lalag, by = intersect(names(temp.ad), names(temp.lalag)))
      
      lag.lm <- lm(larva ~ adult, data = temp.datalag)
      
      lagresult$r2[lagresult$lag == lagi] = summary(lag.lm)$r.squared
      
    } #close lag loop
    
    # Only return best fit lag if bestFitOnly = TRUE
    if (bestFitOnly) {
      lagresult = lagresult[lagresult$r2 == max(lagresult$r2), ]
    }
    
    return(lagresult)
    
  } else {
    warning("Not enough adult records for this species in this year.")
  }
  
}
