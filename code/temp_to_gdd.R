## Process daily temperature data and calculate GDD at 50 km and hex scales
## Input climate data: CPC 50 km min and max temperature data
## Link to CPC climate downloads: https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html#plot
## Grace Di Cecco

library(raster)
library(ncdf4)
library(tidyverse)
library(purrr)
library(dggridR)
library(lubridate)

# Directory with CPC climate raw files
cpc_dir <- "data/cpc/"

cpc_files <- list.files(cpc_dir)

# Place to save output GDD files
output_dir <- "data/derived_data/gdd/"

# Function to calculate growing degree days given tmin and tmax
source("code/gdd_calc.R")

# Quiet degreedays function - no warning messages
quiet_degreedays <- quietly(degreedays)

# Function: read CPC nc file and format into raster brick
# Inputs: directory with CPC files, file name
# Outputs: raster brick with daily temperature data

cpc_raster <- function(cpc_dir, file) {
  
  year <- word(file, 2, sep = fixed('.'))
  var <- word(file, 1, sep = fixed('.'))
  
  # open temp NC file
  nc_data <- nc_open(paste0(cpc_dir, file))
  
  # get lon, lat, time dimensions
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat")
  t <- ncvar_get(nc_data, "time") # hours since 1900-01-01 00:00:00
  
  # get temp array
  temp <- ncvar_get(nc_data, var)
  
  # get missing value
  fillvalue <- ncatt_get(nc_data, var, "missing_value")
  
  nc_close(nc_data)
  
  # Replace missing values with NAs
  temp[temp == fillvalue$value] <- NA
  
  # Raster of temp data
  r <- brick(temp, 
             ymn=min(lon), ymx=max(lon), xmn=min(lat), xmx=max(lat), 
             crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  # Transpose to correct dimensions
  r_transp <- t(r)
  
  # Rotate to convert 0-360 longitude to -180-180
  r_rot <- rotate(r_transp)
  
  # Add day names to layers
  names(r_rot) <- t
  
  return(r_rot)
}

# Function: calculate GDD from CPC daily tmin/tmax data
# Input: CPC tmin and tmax raster bricks
# Other inputs: spatial window for GDD calculations
# Output: data frame containing daily gdd for the year

cpc_gdd <- function(tmin_raster, tmax_raster,
                    lat_min = 15, lat_max = 60, 
                    lon_min = -100, lon_max = -40) {
  
  # Convert tmin brick to df
  tmin_df <- rasterToPoints(tmin_raster) %>%
    as.data.frame() %>%
    pivot_longer(names_to = "day", values_to = "tmin", 3:(nlayers(tmin_raster) + 2)) %>%
    mutate(day = substr(day, 2, length(day)))
  
  # Convert tmax brick to df
  tmax_df <- rasterToPoints(tmax_raster) %>%
    as.data.frame() %>%
    pivot_longer(names_to = "day", values_to = "tmax", 3:(nlayers(tmin_raster) + 2)) %>%
    mutate(day = substr(day, 2, length(day)))
  
  # Set origin of CPC time
  origin <- as.Date("1900-01-01")
  
  # Join tmax and tmin by lat/lon/day of year, calculate GDD for each cell/day
  # Filter to just eastern North America lat-lon (GDD calc step is rate limiting)
  # Convert hours since 1900-01-01 to day of year
  temp_df <- tmin_df %>%
    left_join(tmax_df, by = c("x", "y", "day")) %>%
    rename(lon = "x", lat = "y") %>%
    filter(lat > lat_min & lat < lat_max, lon > lon_min & lon < lon_max) %>%
    mutate(gdd = map2_dbl(tmin, tmax, ~quiet_degreedays(.x, .y)$result),
           day.num = as.numeric(day)) %>%
    mutate(date = origin + day.num/24, 
           doy = yday(date),
           year = year(date))
  
  return(temp_df)  
} 

## Process GDD from CPC data for study years

years <- c(2000:2020)

Sys.time()
for(y in years) {
  
  # file names of tmin and tmax data
  f_tmin <- cpc_files[grepl(y, cpc_files) & grepl("tmin", cpc_files)]
  f_tmax <- cpc_files[grepl(y, cpc_files) & grepl("tmax", cpc_files)]
  
  # create tmin and tmax raster bricks
  tmin_raster <- cpc_raster(cpc_dir, file = f_tmin)
  tmax_raster <- cpc_raster(cpc_dir, file = f_tmax)
  
  # process tmin and tmax daily temps into GDD
  gdd_df <- cpc_gdd(tmin_raster, tmax_raster)
  
  # write GDD data file to output directory
  write.csv(gdd_df, paste0(output_dir, "cpc_gdd_", y, ".csv"), row.names = F)
  
  print(paste(Sys.time(), "Completed year", y))
  
}

## Aggregate GDD at hex scale

# make hex grid
hex <- dggridR::dgconstruct(res = 6)

# read in yearly GDD files
gdd_files <- list.files(output_dir)[grepl("cpc_gdd_", list.files(output_dir))]

gdd_allyears <- tibble(file = gdd_files) %>%
  mutate(data = purrr::map(file, ~read_csv(paste0(output_dir, .))))

# add hex ID, average GDD by hex cell
hex_gdd <- gdd_allyears %>%
  mutate(hex_df = purrr::map(data, ~{
    df <- .
           df$hex_cell <- dgGEO_to_SEQNUM(hex, df$lon, df$lat)$seqnum
           df
         }),
         hex_gdd = purrr::map(hex_df, ~{
           df <- .
           df %>%
             select(-day) %>%
             group_by(hex_cell, day.num, date, doy, year) %>%
             summarize(meanGDD = mean(gdd, na.rm = T))
         }))

hex_unnest <- hex_gdd %>%
  select(-data, -hex_df) %>%
  unnest(cols = c("hex_gdd"))
write.csv(hex_unnest, "data/derived_data/hex_gdd_2000-2020.csv", row.names = F)
