# Impute 2016 Metabolism values
# LMC Oct 2022

# LOAD LIBRARIES
library(tidyverse)
library(lubridate)
library(ggthemes)

## THEME
theme_set(theme_few(base_size = 16, base_family = "Helvetica"))
windowsFonts(Helvetica=windowsFont("TT Helvetica")) 

## Read in data
chan <- read.csv("01_Data/Chan_data_all_091520.csv") %>% 
  select(-X)
AFDM17 <- read.csv("01_Data/2017ChannelAFDM_cleaned.csv") %>% 
  select(-X) %>% 
  rename(sample_event = ExpDay2,
         chan = channel)
str(chan)

chan1 <- chan %>% 
  mutate(StartDate = case_when(Year == "2015" ~ "06/05/2015",
                               Year == "2016" ~ "06/03/2016",
                               Year == "2017" ~ "05/24/2017"),  # Specify experiment start date for each year
         StartDate = as.POSIXct(StartDate, format = "%m/%d/%Y"), # Convert to POSIXct
         sample_event = ifelse(grepl("Metab1", ID2), "D1",
                                       ifelse(grepl("Metab2", ID2), "D2", "BLAH")),
         N_cat = case_when(N_uM == 0.11 ~ "AmbN", # creating categorical variables for the N treatments
                           N_uM == 1.9 | N_uM == 3.68 ~ "MedN",
                           N_uM == 7.25 | N_uM == 10.82 | N_uM == 14.4 ~ "HighN",
                           TRUE ~ "NA"),
         P_cat = case_when(P_uM == 0.36 ~ "AmbP", # creating categorial variables for the P treatments
                           P_uM == 1.17 | P_uM == 1.97 ~ "MedP",
                           P_uM == 3.59 | P_uM == 3.94 | P_uM == 5.2 | P_uM == 6.81 ~ "HighP",
                           TRUE ~ "NA"),
         chan = gsub("C", "", as.character(channel)),
         chan = as.numeric(chan),
         Block = case_when(chan <= 15 ~ 1, # categorizing the channel blocks, since in 2017 rep 1 was in chan 1-15 and rep 2 was in chan 16-30
                           chan > 15 ~ 2),
         propNfixer = RelBV_allN2fixers / 100,
         logitNfixer = log((propNfixer + 0.001) / (1 - propNfixer + 0.001))) %>% # using for rel. biovolume models-- adding minimum proportion N-fixers
  select(-c(AFDM_mg, afdm_g_m2, AFDM_Nfixn_g, NupAFDM_mg, PupAFDM_mg, Nfix_uM_N_mgAFDM)) 
# Note for the "select()": AFDM_mg was the AFDM in the CHAMBER during Metabolism measurements, not the areal AFDM. Removing so we don't accidentally use that as the areal AFDM value. Doing the same for NupAFDM_mg and PupAFDM_mg, since it's also a chamber measurement
# Also removing all areal AFDM values, since we corrected AFDM values in "00_NewAFDM_2015.R" and "00_NewAFDM_2016_2017.R" and are now using those values

# Join corrected AFDM data 
chan17 <- chan1 %>% 
  filter(Year == 2017) %>% 
  left_join(AFDM17, by = c("sample_event", "chan"))

#### IMPUTE MISSING NEP VALUES ####
# one missing NEP value-- needs imputation
## Imputing by fitting a linear regression to the Areal NEP vs. Areal ADFM and using that relationship to impute the missing NEP values
# not removing any points to fit the regression

# Plot NEP vs. AFDM -- unlogged
ggplot(chan17, aes(y = (NEP_uM_C_m2h), x = (Met_gAFDM_m2), color = as.factor(N_uM), size = MeanPre2wksTemp)) + 
  geom_point() +
  ylab(expression(paste("NEP ("~mu~"M C ",m^-2," ", h^-1, ")"))) +
  labs(size = "Temp", color = "Nitrogen") 

# Plot NEP vs. AFDM -- log-transformed
ggplot(chan17, aes(y = log(NEP_uM_C_m2h), x = log(Met_gAFDM_m2), size = MeanPre2wksTemp)) + 
  geom_point()

# Run linear model: log(NEP) ~ log(AFDM)
NEP2017_lm <- lm(log(NEP_uM_C_m2h) ~ log(Met_gAFDM_m2), data = chan17); summary(NEP2017_lm) # all points inicluded
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)        6.65819    0.21745  30.620  < 2e-16 ***
#   log(Met_gAFDM_m2)  0.79107    0.08418   9.397 3.49e-13 ***

# Impute missing NEP value, using regression coefficients from NEP2017_lm
chan17_impute <- chan17 %>% 
  mutate(NEP_impute_uMC_m2h = ifelse(is.na(NEP_uM_C_m2h), exp((0.79107 * log(Met_gAFDM_m2)) + 6.65819), NEP_uM_C_m2h),
         NEP_is_imputed = ifelse(is.na(NEP_uM_C_m2h), TRUE, FALSE)) 

# Check imputation-- looks good
ggplot(chan17_impute, aes(y = log(NEP_impute_uMC_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), shape = as.factor(NEP_is_imputed))) +
  geom_point(size = 3)

##### IMPUTE ER VALUES #######
# 1 missing ER measurement
# keeping all points for the linear regression

# Plot ER vs. AFDM
ggplot(chan17, aes(y = log(R_uM_C_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), size = MeanPre2wksTemp), shape = as.factor(is_imputed)) + 
  geom_point()

# Fit linear model, including all points
ER2017_lm <- lm(log(R_uM_C_m2h) ~ log(Met_gAFDM_m2), data = chan17); summary(ER2017_lm)

# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         5.6585     0.2549  22.198  < 2e-16 ***
# log(Met_gAFDM_m2)   0.9008     0.1002   8.993 1.59e-12 ***

# Impute missing ER value, using regression coefficients from ER2017_lm
chan17_impute2 <- chan17_impute %>% 
  mutate(ER_impute_uMC_m2h = ifelse(is.na(R_uM_C_m2h), exp((0.9008 * log(Met_gAFDM_m2)) + 5.6585), R_uM_C_m2h),
         ER_is_imputed = ifelse(is.na(R_uM_C_m2h), TRUE, FALSE)) 

# Check imputation-- looks good
ggplot(chan17_impute2, aes(y = log(ER_impute_uMC_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), shape = as.factor(ER_is_imputed))) +
  geom_point(size = 3)

##### GPP #######
# Now calculate GPP using imputed NEP and ER
chan17_impute3 <- chan17_impute2 %>% 
  mutate(GPP_impute_uMC_m2h = NEP_impute_uMC_m2h + ER_impute_uMC_m2h,
         GPP_is_imputed = case_when(NEP_is_imputed == "TRUE" | ER_is_imputed == "TRUE" ~ "TRUE",
                                    TRUE ~ "FALSE"))

# Check GPP imputation-- looks good
ggplot(chan17_impute3, aes(y = log(GPP_impute_uMC_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), shape = as.factor(GPP_is_imputed))) +
  geom_point(size = 3)

ggplot(chan17_impute3, aes(y = log(GPP_impute_uMC_m2h), x = MeanPre2wksTemp, color = as.factor(N_uM), shape = as.factor(GPP_is_imputed))) +
  geom_point(size = 3) 

# make dataframe for export
write.csv(chan17_impute3, "IcelandChannels_MS/Lyndsie_Rcode_files/output/00a_IceChan_2017_Metab_imputed.csv")
