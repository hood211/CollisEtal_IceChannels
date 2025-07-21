# calculating AFDM
# Oct 2022, JMH
# Oct 2022, LMC

# libraries
library(tidyverse)

# data ----
## 2016 ----
# bm16 doesn't contain N-fixation AFDM-- we're bringing that data in further down in the script b/c St. Kates didn't give us raw AFDM values-- they gave us "Total Sample AFDM (g)" which represents the AFDM in each chamber during the N-fix measurements
# units here are g AFDM/m2
bm16 <- read.csv("00_BestDataFromLyndsie/IceChan_AFDM_2016_updated.csv", row.names = 1) %>% 
                select(-c(DryWt_g, afdm_g_m2)) %>% # since we're re-calculating these values in this script, I wanted to remove these (LMC)
                mutate(
                  DryWt_g = PanFilterSample_Dried_g - PanWt_g - FilterWt_g,
                  gDM_m2 = DryWt_g * (TotalVol_mL / FilteredVol_mL) / (TileLength_m * TileLength_m * num_tiles),
                  gAFDM = (PanFilterSample_Dried_g - PanFilterSample_Ashed_g),
                  gAFDM_m2 = (PanFilterSample_Dried_g - PanFilterSample_Ashed_g) * (TotalVol_mL/FilteredVol_mL) / (TileLength_m*TileLength_m * num_tiles),
                  perOM = gAFDM_m2 / gDM_m2 * 100) %>%
                filter(!grepl("Biomass", sample_event))

# one has neg gDM_m2 
plot(bm16$gAFDM ~ bm16$DryWt_g)
hist(bm16[bm16$perOM >0,]$perOM, breaks = 20)

# Some AFDM samples seem to have 'weird' %OM values (i.e., >100% or <0%)-- to correct for this we're going to identify samples w/ %OM in the >90th percentile and <10th percentile and re-estimate the AFDM value using the 50th percentile %OM value
summary(bm16[bm16$perOM >0,]$perOM)
per90.16 <- quantile(bm16[bm16$perOM > 0,]$perOM, prob = 0.9, na.rm = T)
per10.16 <- quantile(bm16[bm16$perOM > 0,]$perOM, prob = 0.1, na.rm = T)
per50.16 <- quantile(bm16[bm16$perOM > 0,]$perOM, prob = 0.5, na.rm = T)

ggplot(bm16 %>% 
         filter(gDM_m2 >0) %>% 
         mutate(OverUnder = ifelse(perOM < per10.16 | perOM > per90.16, "OverUnder","Not" )), aes(y = gAFDM_m2, x = gDM_m2, color = OverUnder)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

bm16_2 <- bm16 %>% 
          # if per OM is < or > 10th or 90th percentile than impute using median per OM
          mutate(gAFDM_m2_jh = ifelse(perOM < per10.16 | perOM > per90.16, gDM_m2 * (per50.16 / 100), gAFDM_m2)) %>% 
          select(-gAFDM_m2) %>% 
          # overwrite with imputed
          rename(gAFDM_m2 = gAFDM_m2_jh) %>% 
          select(SampleDate:channel, gAFDM_m2) %>%
          mutate(sample_event2 = ifelse(grepl("Metabolism", sample_event), "Met",
                                         ifelse(grepl("N Uptake", sample_event), "Nup",
                                            ifelse(grepl("P Uptake", sample_event), "Pup", "BLAH"))),
                       ExpDay2 = ifelse(ExpDay %in% c(54,55), "D1",
                                    ifelse(ExpDay %in% c( 60,61), "D2", "BLAH")),
                 gAFDM_m2 = replace(gAFDM_m2, which(gAFDM_m2 < 0), NA)) %>% 
          pivot_wider(id_cols = c(Year, ExpDay2, channel), names_from = sample_event2, values_from = gAFDM_m2) %>% 
          rowwise() %>% 
          mutate(Met = ifelse(is.na(Met), mean(c(Met, Nup, Pup), na.rm = T), Met),
                         Nup = ifelse(is.na(Nup), mean(c(Met, Nup, Pup), na.rm = T), Nup),
                         Pup = ifelse(is.na(Pup), mean(c(Met, Nup, Pup), na.rm = T), Pup))

# join N-fix AFDM with rest of the AFDM data
# LMC: We have 2015 raw data to confirm that "Total_AFDM_sample" was calculated as: Dry Sample - Ash Sample * (Total Vol Filtered / AFDM Vol Filtered)-- so "Total_AFDM_sample" represents the AFDM in each chamber during the N-fix measurements.
bm16nfix <- read.csv("00_BestDataFromLyndsie/2015_2017 Channel N-fixation Data_JMH.csv") %>% 
              mutate(gAFDM_m2 = Total_AFDM_sample / (0.025 * 0.025 * 4)) # Dividing by tile length and width (0.025m and number of tiles (4) gets us the areal rate

# missing C17 AFDM and n fixation data (row missing-- will become an "NA" when we join with metab/uptake AFDM data)
bm16_3 <- bm16_2 %>% 
          left_join(bm16nfix %>% 
                      select(Channel, Smp_event, gAFDM_m2), by = c("channel" = "Channel", "ExpDay2" = "Smp_event"))

names(bm16_3)[7] <- "Nfix"
names(bm16_3)[4:7] <- paste0(names(bm16_3)[4:7],"_gAFDM_m2")


write.csv(bm16_3, "00_BestDataFromLyndsie/2016ChannelAFDM_cleaned.csv")


################################################################3
# 2017
# LMC Oct 2022
# bm17 doesn't contain N-fixation AFDM-- we're bringing that data in further down in the script b/c St. Kates didn't give us raw AFDM values-- they gave us "Total Sample AFDM (g)" which represents the AFDM in each chamber during the N-fix measurements
# units here are g AFDM/m2

# also need N fixation for this
bm17 <- read.csv("00_BestDataFromLyndsie/IceChan_AFDM_2017_updated.csv", row.names = 1) %>% 
  mutate(
    DryWt_g = SmpFilter_Dry_wt_g - FilterWt_g,
    gDM_m2 = DryWt_g * (TotalVol_mL / FilteredVol_mL) / (TileLength_m * TileLength_m * num_tiles),
    gAFDM = SmpFilter_Dry_wt_g - SmpFilter_Ash_wt_g,
    AFDM_g_m2JMH = (SmpFilter_Dry_wt_g - SmpFilter_Ash_wt_g) * (TotalVol_mL / FilteredVol_mL) / (TileLength_m * TileLength_m * num_tiles), # just checking that the afdm_g_m2 calculations are correct
    perOM = AFDM_g_m2JMH / gDM_m2 * 100) %>% 
  filter(!grepl("BMS", sample_event)) # remove all AFDM time series data from the experiment

# one negative dry weight
plot(bm17$gAFDM ~ bm17$DryWt_g)
hist(bm17[bm16$perOM >0,]$perOM, breaks = 20)

# confirming the afdm caculations are correct-- they are
ggplot(bm17, aes(y = AFDM_g_m2JMH, x = afdm_g_m2)) +
  geom_point() +
  geom_abline()

# Some AFDM samples seem to have 'weird' %OM values (i.e., >100% or <0%)-- to correct for this we're going to identify samples w/ %OM in the >90th percentile and <10th percentile and re-estimate the AFDM value using the 50th percentile %OM value
summary(bm17[bm17$perOM > 0,]$perOM)
per90.17 <- quantile(bm17[bm17$perOM > 0,]$perOM, prob = 0.9, na.rm = T)
per10.17 <- quantile(bm17[bm17$perOM > 0,]$perOM, prob = 0.1, na.rm = T)
per50.17 <- quantile(bm17[bm17$perOM > 0,]$perOM, prob = 0.5, na.rm = T)

bm17_2 <- bm17 %>%
  # if per OM is < or > 10th or 90th percentile than impute using median per OM
  mutate(gAFDM_m2_jh = ifelse(perOM < per10.17 | perOM > per90.17, gDM_m2 * (per50.17 / 100), afdm_g_m2),
         gAFDM_m2_lc = ifelse(perOM < 0, afdm_g_m2, gAFDM_m2_jh)) # LMC- I'm adding this b/c C19 on Met D2 had a neg %OM, but it seems to be due to a weird/incorrect filter weight, and not because of the dry - ash calculation. To avoid losing this data point I'm allowing the dry - ash calculation to be the value we use for this sample's AFDM.
  
# confirming the afdm changes are correct-- yupp it only changed the one point I wanted changed
ggplot(bm17_2, aes(y = gAFDM_m2_jh, x = gAFDM_m2_lc)) +
  geom_point() +
  geom_abline()

bm17_3 <- bm17_2 %>% 
  # overwrite with imputed
  rename(gAFDM_m2 = gAFDM_m2_lc) %>% 
  select(SampleDate:channel, gAFDM_m2) %>% 
  mutate(sample_event2 = ifelse(grepl("Metab", sample_event), "Met",
                                ifelse(grepl("N Uptake", sample_event), "Nup",
                                       ifelse(grepl("P Uptake", sample_event), "Pup", "BLAH"))),
         ExpDay2 = ifelse(ExpDay %in% c(48, 49), "D1",
                          ifelse(ExpDay %in% c(59, 60), "D2", "BLAH"))) %>% 
  pivot_wider(id_cols = c(Year, ExpDay2, channel), names_from = sample_event2, values_from = gAFDM_m2) %>% 
  rowwise() %>% 
  mutate(Met = ifelse(is.na(Met), mean(c(Met, Nup, Pup), na.rm = T), Met),
         Nup = ifelse(is.na(Nup), mean(c(Met, Nup, Pup), na.rm = T), Nup),
         Pup = ifelse(is.na(Pup), mean(c(Met, Nup, Pup), na.rm = T), Pup))

names(bm17_3)[4:6] <- paste0(names(bm17_3)[4:6],"_gAFDM_m2")

# N-fix AFDM data
# LMC: SOMETHING IS WRONG WITH THESE DATA-- values are off by a factor of ~100 compared to 2015 and 2016
# LMC: We have 2015 raw data to confirm that "Total_AFDM_Sample" was calculated as: Dry Sample - Ash Sample * (Total Vol Filtered / AFDM Vol Filtered)-- so "Total_AFDM_sample" is supposed to represent the AFDM in each chamber during the N-fix measurements.
# However, it's clear something is wrong with the 2017 data-- since we don't have raw values we can't figure out what's wrong
bm17nfix <- read.csv("00_BestDataFromLyndsie/2017_Channel_Nfix_AFDM_data_LMC.csv") %>% 
  mutate(gAFDM_m2 = Total_AFDM_sample_g / (0.025 * 0.025 * 4)) # Dividing by tile length and width (0.025m and number of tiles (4) gets us the areal rate

# CAN'T JOIN N-FIX SINCE IT'S NOT CORRECT
# bm17_4 <- bm17_3 %>% 
#   left_join(bm17nfix %>% 
#               select(Channel, Smp_event, gAFDM_m2), by = c("channel" = "Channel", "ExpDay2" = "Smp_event"))

# names(bm17_4)[7] <- "Nfix"
# names(bm17_4)[4:7] <- paste0(names(bm17_4)[4:7],"_gAFDM_m2")


write.csv(bm17_3, "00_BestDataFromLyndsie/2017ChannelAFDM_cleaned.csv")
