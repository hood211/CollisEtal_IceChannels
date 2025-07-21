# Data prep for SEM analysis
# LMC Nov 2022

## LOAD LIBRARIES
library(tidyverse)
library(ggthemes)
library(gridExtra)
library(GGally)

## SET THEME ##################################################
theme_set(theme_few(base_size = 16, base_family = "Helvetica"))
windowsFonts(Helvetica=windowsFont("TT Helvetica"))

## LOAD DATA
dataSEM <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv") %>% 
  select(-X)

### Log-transform data and rename variables
## Subsetting 2017 data b/c that's the only year used for the SEM
## Not using relative abundances because they are perfectly neg cor- using BMS of Nfixers and nonNfixers instead (also makes more sense linking to GPP)
dataSEM1 = dataSEM %>% 
  filter(Year == "2017") %>%
  mutate(ln_GPP_uMCm2h = log(GPP_uM_C_m2h), 
         ln_NEP_uMCm2h = log(NEP_uM_C_m2h), 
         ln_ER_uMCm2h = log(R_uM_C_m2h), 
         ln_Nfix_uMNm2h = log(Nfix_uM_N_m2h),
         ln_BMS_N2fixers_g_m2 = log(BMS_N2fixers_g_m2 + 0.01791971), # contains zeros, so adding minimum value for 2017 data
         ln_BMS_nonN2fixers_g_m2 = log(BMS_nonN2fixers_g_m2 + 0.1136800), # contains zeros, so adding minimum value for 2017 data
         ln_Nup_uMNm2h = log(NUp_uM_N_m2_hr),
         ln_Pup_uMNm2h = log(PUp_uM_P_m2_hr),
         ln_AFDM_gm2 = log(Met_gAFDM_m2),
         ln_temp = log(MeanPre2wksTemp),
         Ncat = as.factor(N_uM),
         Ratiocat = case_when(NPratio <= 0.93 ~ "low",
                              NPratio > 0.93 ~ "high"),
         Ratiocat = as.factor(Ratiocat)) %>% 
  select(c(Year, ln_temp, N_uM, P_uM, Ncat, NPratio, Ratiocat, ln_GPP_uMCm2h, ln_NEP_uMCm2h, ln_ER_uMCm2h, ln_AFDM_gm2, ln_Nup_uMNm2h, ln_Pup_uMNm2h, ln_Nfix_uMNm2h, 
           ln_BMS_N2fixers_g_m2, ln_BMS_nonN2fixers_g_m2)) %>% 
  rename(temp = ln_temp,
         N = N_uM,
         P = P_uM,
         GPP = ln_GPP_uMCm2h,
         NEP = ln_NEP_uMCm2h,
         ER = ln_ER_uMCm2h,
         AFDM = ln_AFDM_gm2,
         Nup = ln_Nup_uMNm2h,
         Pup = ln_Pup_uMNm2h,
         Nfix = ln_Nfix_uMNm2h,
         Nfixers = ln_BMS_N2fixers_g_m2,
         nonNfixers = ln_BMS_nonN2fixers_g_m2)

# Subset data: variables to potentially use in SEM
dat2017 <- dataSEM1 %>% 
  mutate(Ncat = droplevels(Ncat)) %>%
  select(temp, N, P, NPratio, Ncat, Ratiocat, GPP, NEP, ER, Nup, Nfix, Nfixers, nonNfixers, AFDM) 

# Check relationships between variables
dataSEM1 %>% 
  select(temp, N, P, NPratio, GPP, NEP, ER, Nup, Nfix, Nfixers, nonNfixers, AFDM) %>% 
  ggpairs(upper = list(continuous = wrap("cor", size = 2.75, alignPercent = 1)),
          lower = list(continuous = wrap("points", size = 0.75, alpha = 0.3))) +
    theme_few(base_size = 8) 

write.csv(dat2017, "01_Data/06_IceChan_2017_dataforSEM.csv") 