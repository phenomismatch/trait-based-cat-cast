library(tidyverse)
library(splitstackshape)

# read adult bfly pdfs
pdfs <- read.csv("data/derived_data/simpleton_pheno_pdfs-OutlierDetection_CIs.csv")

head(pdfs)

# split into list by year, overwintering code, and hexcell
pdf_list <- split(pdfs, 
                  f = list(pdfs$year, 
                           pdfs$code,
                           pdfs$HEXcell), 
                  drop = TRUE)

# write function to calculate 10% & 50% phenometrics
pheno_fun <- function(x){
  
  pdf_df <- pdf_list[[x]] %>% 
    mutate(
      cum_prob = cumsum(meanPDF),
      cum_perc = cum_prob / max(cum_prob),
      tenth = DOY[which.max(cum_perc >= 0.1)],
      fiftieth = DOY[which.max(cum_perc >= 0.5)],
      
      cum_prob_lowCI = cumsum(PDF2.5),
      cum_perc_lowCI = cum_prob_lowCI / max(cum_prob_lowCI),
      tenth_highCI = DOY[which.max(cum_perc_lowCI >= 0.1)],
      fiftieth_lowCI = DOY[which.max(cum_perc_lowCI >= 0.5)],
      
      cum_prob_highCI = cumsum(PDF97.5),
      cum_perc_hichCI = cum_prob_highCI / max(cum_prob_highCI),
      tenth_lowCI = DOY[which.max(cum_perc_hichCI >= 0.1)],
      fiftieth_highCI = DOY[which.max(cum_perc_hichCI >= 0.5)],
      
      fiftieth_CI_days = fiftieth_highCI - fiftieth_lowCI,
      tenth_CI_days = tenth_highCI - tenth_lowCI
    ) %>% 
    dplyr::select(year, HEXcell, code, obsDays, abundance,
                  tenth, tenth_lowCI, tenth_highCI, 
                  fiftieth, fiftieth_lowCI, fiftieth_highCI,
                  fiftieth_CI_days, tenth_CI_days) 
  
  pdf_df2 <- pdf_df[1,]
  
  return(pdf_df2)
}

# calculate phenometrics for each grouping
pheno_list <- lapply(X = 1:length(pdf_list), FUN = pheno_fun)
# combine outputs into df
pheno_df_output <- do.call(bind_rows, pheno_list)

write.csv(x = pheno_df_output, "outputs/adult_bfly_phenometrics.csv", row.names = F)


# data exploration, examining the effect of observation days and total number of 
# observations on the number of days the CI spans

#ggplot(pheno_df_output) +
#  geom_point(mapping = aes(x = abundance, y = abs(tenth_CI_days), color = code)) +
#  geom_smooth(mapping = aes(x = abundance, y = abs(tenth_CI_days), color = code)) +
#  geom_hline(yintercept = 0) +
#  scale_x_log10() +
#  theme_bw()
#
#ggplot(pheno_df_output) +
#  geom_point(mapping = aes(x = obsDays, y = abs(tenth_CI_days), color = code)) +
#  geom_smooth(mapping = aes(x = obsDays, y = abs(tenth_CI_days), color = code)) +
#  geom_hline(yintercept = 0) +
#  scale_x_log10() +
#  theme_bw()

# read in data
bfly_pheno <- pheno_df_output %>% 
  mutate(cell_code = paste(code, HEXcell, sep = "_")) # add hexcell for joining

#filter to 2016 - 2020
pheno_2016 <- filter(bfly_pheno, year >= 2016)

consq_yrs2 <- pheno_2016 %>% 
  group_by(HEXcell, code) %>% 
  summarise(years = n(), 
            mean10 = mean(tenth),
            mean50 = mean(fiftieth))

tot_2016.2020 <- filter(consq_yrs2, years == 5) %>% 
  mutate(cell_code = paste(code, HEXcell, sep = "_")) # add hexcell for joining/antijoining

head(tot_2016.2020)

#### Filter to only keep cells with at least 5 years
bfly_pheno_prun <- bfly_pheno %>% 
  filter(cell_code %in% tot_2016.2020$cell_code)

bfly_pheno_mean <- left_join(bfly_pheno_prun, tot_2016.2020)

# calculate deviations from mean
bfly_pheno_mean <- bfly_pheno_mean %>% 
  mutate(dev_10 = tenth - mean10,
         dev_50 = fiftieth - mean50)

write.csv(bfly_pheno_mean,
          "outputs/adult_bfly_phenometrics_deviations.csv", row.names = F)

## do some real quick visualizations
#ggplot(bfly_pheno_mean) + 
#  geom_point(mapping = aes(x = year, y = dev_10, color = code)) +
#  geom_smooth(mapping = aes(x = year,y = dev_10, color = code), method = "lm") + 
#  facet_wrap(~HEXcell)
#
#ggplot(bfly_pheno_mean) + 
#  geom_point(mapping = aes(x = year, y = dev_50, color = code)) +
#  geom_smooth(mapping = aes(x = year,y = dev_50, color = code), method = "lm") + 
#  facet_wrap(~HEXcell)
#

# There are some major outliers. These should probably be removed before analysis.
# i.e., remove any years that have have a deviation >60 (?) days 