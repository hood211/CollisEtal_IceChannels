# 2015 AFDM data
# JMH Oct 2022
# LMC OCt 2022


# library
library(tidyverse)

#data ---
BM15 <- read.csv("00_BestDataFromLyndsie/Data_forR/2015ChannelAFDM_rawJMH.csv") %>% 
    mutate(PanPlusdryFilterWeightG_2 = as.numeric(ifelse(PanPlusdryFilterWeightG == "PAN NOT ASHED", AshedPanWeight_g + DrySamplePlusFilterWeight_g, PanPlusdryFilterWeightG)),
           DMg_jh =   DrySamplePlusFilterWeight_g - AFDMFilterWeight_g,
           AFDMg_jh = PanPlusdryFilterWeightG_2 - PanAshedFilterWeightG,
           perOM = AFDMg_jh / DMg_jh * 100) 

  
  plot(BM15$AFDMg_jh ~BM15$DMg_jh)
  hist(BM15[BM15$perOM >0 & BM15$perOM <100,]$perOM, breaks = 30) # removing some extreme outliers for visualization
  # if per OM exceeds the inner quartile use DM and median to calculate
  # removing the unreasonable values: < 0 and > 100
  summary(BM15[BM15$perOM >0 & BM15$perOM <100,]$perOM)
  perOM15_10 <- quantile(BM15[BM15$perOM >0 & BM15$perOM <100,]$perOM, prob = 0.1, na.rm = T) # just checking how different it is if you use 10th and 90th percentile
  perOM15_90 <- quantile(BM15[BM15$perOM >0 & BM15$perOM <100,]$perOM, prob = 0.9, na.rm = T)
  perOM_25 <- quantile(BM15[BM15$perOM >0 & BM15$perOM <100,]$perOM, prob = 0.25, na.rm = T)
  perOM_75 <- quantile(BM15[BM15$perOM >0 & BM15$perOM <100,]$perOM, prob = 0.75, na.rm = T)
  perOM_50 <- quantile(BM15[BM15$perOM >0 & BM15$perOM <100,]$perOM, prob = 0.5, na.rm = T)
  
  # On D1, Nfix was measured over 3 days, 
  # Based on ChannelN2FixationDAta2015Jim.xls it seems that the AFDM data for 7/16/2015 is consistent 
  # with what Delore thought was the best Nfix measurements, even though dates do not make sense
  SampleDatesMonographs <- c("7/14/2015", "7/15/2015", "7/16/2015",  "7/27/2015", "7/28/2015",  "7/29/2015", "7/30/2015")
  SampleDatesMonographsD1 <- c("7/14/2015", "7/15/2015", "7/16/2015")
  SampleDatesMonographsD2 <- c("7/27/2015", "7/28/2015",  "7/29/2015", "7/30/2015")
  
  
BM15_2 <- BM15 %>% 
            # exclude unneeded dates
            filter(SampleDate %in% SampleDatesMonographs) %>% 
            # if per OM exceeds the inner quartile use DM and median to calculate - see above.
            mutate(AFDMg_jh2 = ifelse(perOM < perOM_25 | perOM > perOM_75, DMg_jh * (perOM_50 / 100), AFDMg_jh),
                   SubSmpPer = AFDMSampleVolMl / RecordedTotalSampleVolMl,
                   AFDMg_jh2_tot = AFDMg_jh2 / SubSmpPer,
                   AFDMgm2 = AFDMg_jh2_tot / (TileCM2 / 10000)) %>% 
            select(SampleDate, PairedAnalysis, channel, AFDMgm2) %>% 
            mutate(PairedAnalysis = ifelse(PairedAnalysis == "N FIX'N", "Nfix", PairedAnalysis),
                    MetDate = case_when(SampleDate %in% SampleDatesMonographsD1 ~ "D1",
                                       SampleDate %in% SampleDatesMonographsD2 ~ "D2"))
  
BM15s <- BM15_2 %>% 
          pivot_wider(id_cols = c(MetDate, channel), names_from = PairedAnalysis, values_from = AFDMgm2) %>% 
          # deal with NA's
          rowwise() %>% 
          mutate(AFDMgm2_mean = mean(c(Met, Nup, Pup, Nfix), na.rm = T)) %>% 
          # if missing data fill in with average
          # if zero fill with average - issue with C24, day 2
          mutate(Met = ifelse(is.na(Met) | Met == 0, AFDMgm2_mean, Met),
                 Nup = ifelse(is.na(Nup) | Met == 0, AFDMgm2_mean, Nup),
                 Pup = ifelse(is.na(Pup) | Met == 0, AFDMgm2_mean, Pup),
                 Nfix = ifelse(is.na(Nfix) | Met == 0, AFDMgm2_mean, Nfix)) %>% 
          select(-AFDMgm2_mean)
  
names(BM15s)[3:6] <- paste0(names(BM15s)[3:6],"_gAFDM_m2")
  
write.csv(BM15s, "00_BestDataFromLyndsie/2015ChannelAFDM_cleaned.csv")