# Impute 2016 Metabolism values
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
AFDM16 <- read.csv("01_Data/2016ChannelAFDM_cleaned.csv") %>% 
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
# Note for the "select()": AFDM_mg was the AFDM in the CHAMBER during Metabolism measurements, not the areal AFDM. 
# Removing so we don't accidentally use that as the areal AFDM value. Doing the same for NupAFDM_mg and PupAFDM_mg, since it's also a chamber measurement
# Also removing all areal AFDM values, since we corrected AFDM values in "00_NewAFDM_2015.R" and "00_NewAFDM_2016_2017.R" and are now using those values

#### IMPUTE MISSING NEP, ER AND GPP VALUES ####
## Looking at Areal Metabolism (GPP, NEP, ER) vs. Areal ADFM (afdm_g_m2) and use that relationship to impute the missing metabolism values

# Join corrected AFDM data 
chan16 <- chan1 %>% 
  filter(Year == 2016) %>% 
  left_join(AFDM16, by = c("sample_event", "chan"))

# Plot NEP vs. AFDM -- unlogged
ggplot(chan16, aes(y = (NEP_uM_C_m2h), x = (Met_gAFDM_m2), color = as.factor(N_uM), size = MeanPre2wksTemp)) + 
  geom_point() +
  scale_color_brewer(palette = "Spectral") +
  ylab(expression(paste("NEP ("~mu~"M C ",m^-2," ", h^-1, ")"))) +
  labs(size = "Temp", color = "Nitrogen") 

# Plot NEP vs.AFDM -- log-transformed
ggplot(chan16, aes(y = log(NEP_uM_C_m2h), x = log(Met_gAFDM_m2), size = MeanPre2wksTemp)) + 
  geom_point()

# Check what NEP-AFDM relationship looks like if we use linear vs. polynomial regression
# change the "geom_smooth" function to see each version
# We're going to use a second-order polynomial, based on the visual fit
chan16 %>% 
  ggplot(aes(y = log(NEP_uM_C_m2h), x = log(Met_gAFDM_m2))) + 
    geom_point() +
    geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE) # polynomial
    # geom_smooth(method = lm, formula = y ~ x, se = FALSE) # linear


### Imputing missing NEP values via regression imputation:
#   Including all points to fit the regression
#   Then fitting a second-order polynomial model between areal NEP and AFDM
#   Predicting areal NEP for missing values based off meausred AFDM for the fitted model

NEP2016_lm <- lm(log(NEP_uM_C_m2h) ~ log(Met_gAFDM_m2) + I(log(Met_gAFDM_m2)^2), data = chan16); summary(NEP2016_lm)

# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)             4.08943    0.36784  11.117 3.08e-15 ***
#   log(Met_gAFDM_m2)       3.21775    0.34619   9.295 1.47e-12 ***
#   I(log(Met_gAFDM_m2)^2) -0.46638    0.07671  -6.080 1.54e-07 ***

# Impute NEP values  
chan16_impute <- chan16 %>% 
  mutate(NEP_impute_uMC_m2h = ifelse(is.na(NEP_uM_C_m2h), exp((-0.46638 * log(Met_gAFDM_m2)^2) + (3.21775 * log(Met_gAFDM_m2)) + 4.08943), NEP_uM_C_m2h),
         NEP_is_imputed = ifelse(is.na(NEP_uM_C_m2h), TRUE, FALSE)) 

# Check imputation
# Imputed values look good
ggplot(chan16_impute, aes(y = log(NEP_impute_uMC_m2h), x = log(Met_gAFDM_m2), color = as.factor(P_uM), shape = as.factor(NEP_is_imputed))) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Spectral")

ggplot(chan16_impute, aes(y = log(NEP_impute_uMC_m2h), x = MeanPre2wksTemp, color = as.factor(P_uM), shape = as.factor(NEP_is_imputed))) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Spectral")

##### ER VALUES #######
# No missing ER values for 2016
# Just creating an "impute" column to keep dataframes consistent across years
chan16_impute2 <- chan16_impute %>% 
  mutate(ER_impute_uMC_m2h = R_uM_C_m2h, # all values are the same since nothing needs to be imputed
         ER_is_imputed = ifelse(is.na(R_uM_C_m2h), TRUE, FALSE))

##### GPP #######
# Now calculate GPP using imputed NEP (no ER values needed to be imputed)
chan16_impute3 <- chan16_impute2 %>% 
  mutate(GPP_impute_uMC_m2h = NEP_impute_uMC_m2h + ER_impute_uMC_m2h,
         GPP_is_imputed = case_when(NEP_is_imputed == "TRUE" | ER_is_imputed == "TRUE" ~ "TRUE",
                                    TRUE ~ "FALSE"))

# Check GPP imputation-- looks good
ggplot(chan16_impute3, aes(y = log(GPP_impute_uMC_m2h), x = log(Met_gAFDM_m2), color = as.factor(P_uM), shape = as.factor(GPP_is_imputed))) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Spectral")

ggplot(chan16_impute3, aes(y = log(GPP_impute_uMC_m2h), x = MeanPre2wksTemp, color = as.factor(P_uM), shape = as.factor(GPP_is_imputed))) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Spectral")

# make dataframe for export
write.csv(chan16_impute3, "IcelandChannels_MS/Lyndsie_Rcode_files/output/00a_IceChan_2016_Metab_imputed.csv")
