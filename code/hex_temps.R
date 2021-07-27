### Get Daymet temperature means for hex cells

### Libraries

library(sf)
library(rgdal)
library(rgeos)
library(daymetr)
library(tidyverse)
library(raster)
library(ncdf4)
library(lubridate)
library(units)
library(purrr)

### Hex cells

hex <- st_read("data/maps/hex_grid_crop.shp", stringsAsFactors= F) %>%
  mutate(cell.num = as.numeric(cell)) %>%
  dplyr::select(-cell) %>%
  rename(cell = "cell.num")

## 2012-2019
## Mar-Jun average

# Function: get average temperature for one hex cell polygon

daymetMean <- function(id) {
  hex <- filter(hex_transf, cell == id)
  
  daymet_crop <- raster::crop(daymet_mean, hex)
  daymet_mask <- raster::mask(daymet_crop, hex)
  
  local_extract <- raster::extract(daymet_mask, hex, fun = mean, na.rm = T, df = T)
  
  return(mean(local_extract$layer))
}

# Function errors if some routes are outside the extent of DAYMET raster - use possibly function to ignore those routes

possibly_daymetMean <- possibly(daymetMean, otherwise = NA)

# Calculate mean spring temp at hex cells
# Make sure you have folders "daymet" and "daymet_out" in directory where you will download Daymet files

wd <- "/Users/gracedicecco/Desktop/"

daymet_out_wd <- "/Users/gracedicecco/Desktop/daymet_out"

years <- c(2000:2020)

for(y in years) {
  download_daymet_ncss(location = c(60, -100, 20, -52),
                       start = y,
                       end = y,
                       param = c("tmin", "tmax"), 
                       frequency = "daily",
                       path = paste0(wd, "daymet/"))
  
  # Read in data
  files <- list.files(paste0(wd, "daymet/"))
  
  for(f in files) {
    daymet_nc <- nc_open(paste0(wd, "daymet/", f))
    daymet_raster <- brick(paste0(wd, "daymet/", f))
    crs(daymet_raster) <- "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0"
    
    daymet_spring <- daymet_raster[[3:6]]
    
    daymet_mean <- mean(daymet_spring, na.rm = T)
    
    crs_daymet <- crs(daymet_spring)
    
    hex_transf <- st_transform(hex_subset, "+proj=lcc +datum=WGS84 +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +units=km +lat_1=25 +lat_2=60 +ellps=WGS84 +towgs84=0,0,0")
    
    hexclim <- data.frame(cell = hex_subset$cell) %>%
      mutate(mean_temp = purrr::map_dbl(cell, possibly_daymetMean))
    
    write.csv(hexclim, paste0(daymet_out_wd, "/", f, ".csv"), row.names = F)
    print(f)
    nc_close(daymet_nc)
  }
  
  print(y)
  sapply(paste0(wd, "daymet/", files), unlink)
  
}

## Combine files to one dataset

hexDAYMET <- data.frame(cell = c(), year = c(), mean_tmax = c(), mean_tmin = c())

for(y in years) {
  files <- list.files(daymet_out_wd)
  files_y <- files[grepl(y, files)]
  
  tmax <- read.csv(paste0(daymet_out_wd, "/", files_y[grepl("tmax", files_y)])) %>%
    rename(tmax = "mean_temp")
  tmin <- read.csv(paste0(daymet_out_wd, "/", files_y[grepl("tmin", files_y)])) %>%
    rename(tmin = "mean_temp")
  
  tmp <- tmax %>%
    left_join(tmin) %>%
    mutate(year = y)
  print(nrow(tmp))
  
  hexDAYMET <- rbind(hexDAYMET, tmp)
}

hexDAYMET$mean_temp <- rowMeans(hexDAYMET[, 2:3]) ## Update with spring equinox to summer solstice dates

write.csv(hexDAYMET, "data/derived_data/hex_mean_temps.csv", row.names = F)
