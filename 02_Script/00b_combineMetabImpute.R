# Combine dataframes of imputed values from all experiment years
# Create a "master dataframe" of all values from the experiment
# LMC Nov 2022

# LOAD LIBRARIES
library(tidyverse)
library(lubridate)
library(ggthemes)

## THEME
theme_set(theme_few(base_size = 16, base_family = "Helvetica"))
windowsFonts(Helvetica=windowsFont("TT Helvetica")) 

# read in data
# removing the stoich data from each df 
impute15 <- read.csv("01_Data/00a_IceChan_2015_Metab_imputed.csv") %>% 
  select(-c(P_uM.stoich:Pred_P_dem))
impute16 <- read.csv("01_Data/00a_IceChan_2016_Metab_imputed.csv") %>% 
  select(-Year.y) %>% 
  rename(Year = Year.x) %>% 
  select(-c(P_uM.stoich:Pred_P_dem))
impute17 <- read.csv("01_Data/00a_IceChan_2017_Metab_imputed.csv") %>% 
  select(-Year.y) %>% 
  rename(Year = Year.x) %>% 
  select(-c(P_uM.stoich:Pred_P_dem))

allyrs <- bind_rows(impute15, impute16, impute17)

allyrs2 <- allyrs %>% 
  select(-c(NEP_uM_C_m2h, R_uM_C_m2h, GPP_uM_C_m2h, GPP_uM_CperL)) %>% # removing the un-imputed metab variables, to replace with the imputed values
  rename(NEP_uM_C_m2h = NEP_impute_uMC_m2h, # rename the imputed metab columns
         R_uM_C_m2h = ER_impute_uMC_m2h,
         GPP_uM_C_m2h = GPP_impute_uMC_m2h,
         lightL_Met = lightL) %>%
  mutate(BMS_N2fixers_g_m2 = (RelBV_allN2fixers / 100) * Met_gAFDM_m2,
         BMS_nonN2fixers_g_m2 = (RelBV_allnonN2fixers / 100) * Met_gAFDM_m2) %>% 
  select(c(MetDate, Year, sample_event, ID2:channel, chan, Block, repF:NPratio, N_cat, P_cat, StartDate, lightL_Met, NEP_uM_C_m2h, NEP_is_imputed, R_uM_C_m2h, ER_is_imputed, 
           GPP_uM_C_m2h, GPP_is_imputed, Met_gAFDM_m2, UpDate:PUp_uM_P_m2_hr, Nup_gAFDM_m2, Pup_gAFDM_m2, NfixDate:Nfix_uM_N_m2h, Nfix_gAFDM_m2, Comments, 
           RelBV_allN2fixers:RelBV_otheralgalgroups, BMS_N2fixers_g_m2, BMS_nonN2fixers_g_m2, propNfixer, logitNfixer)) 

write.csv(allyrs2, "01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv")
