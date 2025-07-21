# Impute 2015 Metabolism values
# LMC Oct 2022

# LOAD LIBRARIES
library(tidyverse)
library(lubridate)
library(ggthemes)

## SET THEME
theme_set(theme_few(base_size = 16, base_family = "Helvetica"))
windowsFonts(Helvetica=windowsFont("TT Helvetica")) 

## LOAD DATA
chan <- read.csv("01_Data/Chan_data_all_091520.csv") %>% 
  select(-X)
AFDM15 <- read.csv("01_Data/2015ChannelAFDM_cleaned.csv") %>% 
  select(-X) %>% 
  rename(sample_event = MetDate,
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
# Note for the "select()": AFDM_mg was the AFDM in the CHAMBER during Metabolism measurements, not the areal AFDM. 
# Removing so we don't accidentally use that as the areal AFDM value. Doing the same for NupAFDM_mg and PupAFDM_mg, since it's also a chamber measurement
# Also removing all areal AFDM values, since we corrected AFDM values in "00_NewAFDM_2015.R" and "00_NewAFDM_2016_2017.R" and are now using those values
# (new datasheets are YEARChannelAFDM_cleaned.csv) in the 00_BestDatafromLyndsie folder on OneDrive

#### IMPUTE MISSING NEP, ER AND GPP VALUES ####
## Looking at Areal Metabolism (GPP, NEP, ER) vs. Areal ADFM (afdm_g_m2) and use that relationship to impute the missing metabolism values

# Join corrected AFDM data 
chan15 <- chan1 %>% 
  filter(Year == 2015) %>% 
  left_join(AFDM15, by = c("sample_event", "chan"))

# Plot NEP vs. ER for each year
ggplot(chan1, aes(y = (NEP_uM_C_m2h), x = (R_uM_C_m2h), color = as.factor(Year))) + 
  geom_point(size = 2) +
  scale_color_brewer(palette = "Oranges") +
  geom_smooth(method = lm) +
  geom_abline(intercept = 0, slope = 2.73) +
  facet_grid(Year~.)

# Plot NEP vs. AFDM (unlogged)
ggplot(chan15, aes(y = (NEP_uM_C_m2h), x = (Met_gAFDM_m2), color = as.factor(N_uM), size = MeanPre2wksTemp)) + 
  geom_point() +
  scale_color_brewer(palette = "Spectral") +
  ylab(expression(paste("NEP ("~mu~"M C ",m^-2," ", h^-1, ")"))) +
  labs(size = "Temp", color = "Nitrogen") 

# Plot NEP vs. AFDM (logged)
ggplot(chan15, aes(y = log(NEP_uM_C_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), size = MeanPre2wksTemp)) + 
  geom_point() +
  scale_color_brewer(palette = "Spectral")

### Imputing missing NEP values via regression imputation: ####
#   Including all points to fit the regression
#   We fit a linear model between areal NEP and AFDM (log-transformed)
#   Predicting areal NEP for missing values based off meausred AFDM for the fitted model

# Run linear model: log(NEP) ~ log(AFDM)
NEP2015_lm1 <- lm(log(NEP_uM_C_m2h) ~ log(Met_gAFDM_m2), data = chan15); summary(NEP2015_lm1)

# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)         7.7097     0.7157  10.772 9.63e-15 ***
#   log(Met_gAFDM_m2)   0.4689     0.2177   2.154    0.036 * 
  
# Impute missing NEP values, using regression coefficients from NEP2015_lm models
chan15_impute <- chan15 %>% 
  mutate(NEPimpute_lm1 = ifelse(is.na(NEP_uM_C_m2h), exp((0.4689 * log(Met_gAFDM_m2)) + 7.7097), NEP_uM_C_m2h), # coefficients from lm
         NEP_is_imputed = ifelse(is.na(NEP_uM_C_m2h), TRUE, FALSE)) 

##### IMPUTE ER VALUES #######
# Only missing 1 ER measurement

# Plot ER vs. AFDM
ggplot(chan15, aes(y = log(R_uM_C_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), size = MeanPre2wksTemp), shape = as.factor(is_imputed)) + 
  geom_point() +
  scale_color_brewer(palette = "Spectral")

# See what lm looks like with all points included
ggplot(chan15, aes(y = log(R_uM_C_m2h), x = log(Met_gAFDM_m2))) + 
  geom_point(size = 2) +
  scale_color_brewer(palette = "Spectral") +
  geom_smooth(method = "lm", se = FALSE)

# Fit linear model, including all points
ER2015_lm <- lm(log(R_uM_C_m2h) ~ log(Met_gAFDM_m2), data = chan15); summary(ER2015_lm)

# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)        6.61203    0.30209  21.888  < 2e-16 ***
#   log(Met_gAFDM_m2)  0.61273    0.09472   6.469 2.43e-08 ***

# Impute missing ER value, using regression coefficients from ER2015_lm
chan15_impute2 <- chan15_impute %>% 
  mutate(ER_impute_uMC_m2h = ifelse(is.na(R_uM_C_m2h), exp((0.61273 * log(Met_gAFDM_m2)) + 6.61203), R_uM_C_m2h),
         ER_is_imputed = ifelse(is.na(R_uM_C_m2h), TRUE, FALSE)) 

# Check imputation-- looks good
ggplot(chan15_impute2, aes(y = log(ER_impute_uMC_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), shape = as.factor(ER_is_imputed))) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Spectral")

# Now calculate GPP using imputed NEP and ER
chan15_impute3 <- chan15_impute2 %>% 
  # select(-c(NEPimpute_lm2, NEPimpute_lm3)) %>% # remove NEP imputions not used
  rename(NEP_impute_uMC_m2h = NEPimpute_lm1) %>% # rename imputed NEP column
  mutate(GPP_impute_uMC_m2h = NEP_impute_uMC_m2h + ER_impute_uMC_m2h,
         GPP_is_imputed = case_when(NEP_is_imputed == "TRUE" | ER_is_imputed == "TRUE" ~ "TRUE",
                                    TRUE ~ "FALSE"))

# Check GPP imputation-- looks good
ggplot(chan15_impute3, aes(y = log(GPP_impute_uMC_m2h), x = log(Met_gAFDM_m2), color = as.factor(N_uM), shape = as.factor(GPP_is_imputed))) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Spectral")

ggplot(chan15_impute3, aes(y = log(GPP_impute_uMC_m2h), x = MeanPre2wksTemp, color = as.factor(N_uM), shape = as.factor(GPP_is_imputed))) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Spectral")

# make dataframe for export
write.csv(chan15_impute3, "01_Data/00a_IceChan_2015_Metab_imputed.csv")