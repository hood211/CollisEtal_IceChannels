# Areal responses for 2017
# JMH, Dec 2022/Jan 2023, Jun 24


# Libraries ----
# general
library(tidyverse)
# gam models
library(mgcv)
library(mgcViz)
library(MuMIn)


# Data ----
chan <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv", row.names = 1) %>% 
  mutate(across(c(MetDate, UpDate, NfixDate), as.POSIXct, format = "%m/%d/%Y")) %>% 
  mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d")) %>% 
  mutate(NEP_uM_C_gAFDM_h = NEP_uM_C_m2h/Met_gAFDM_m2,
         R_uM_C_gAFDM_h = R_uM_C_m2h/Met_gAFDM_m2,
         GPP_uM_C_gAFDM_h = GPP_uM_C_m2h/Met_gAFDM_m2,
         NUp_uM_N_gAFDM_h = (NUp_uM_N_m2_hr)/Nup_gAFDM_m2,
         # N fix AFDM not usable, use average of met and Nup
         Nfix_uM_N_gAFDM_h = Nfix_uM_N_m2h/((Met_gAFDM_m2+Nup_gAFDM_m2)/2),
         TotNAssim_uM_N_gAFDM_h = NUp_uM_N_gAFDM_h + Nfix_uM_N_gAFDM_h) %>% 
  # C12 & C25 on this date have anomalously low AFDM/m2 compared to GPP/ER rates and AFDM/m2 for N uptake
  # Models are VERY SENSITIVE to this
  mutate(GPP_uM_C_gAFDM_h = ifelse(MetDate == "2017-07-11" & channel %in% c("C12", "C25"), 
                                   GPP_uM_C_m2h/Nup_gAFDM_m2,
                                   GPP_uM_C_m2h/Met_gAFDM_m2),
         R_uM_C_gAFDM_h = ifelse(MetDate == "2017-07-11" & channel %in% c("C12", "C25"), 
                                 R_uM_C_m2h/Nup_gAFDM_m2,
                                 R_uM_C_m2h/Met_gAFDM_m2),
         NEP_uM_C_gAFDM_h = ifelse(MetDate == "2017-07-11" & channel %in% c("C12", "C25"), 
                                   NEP_uM_C_m2h/Nup_gAFDM_m2,
                                   NEP_uM_C_m2h/Met_gAFDM_m2))


chanF <- chan %>% 
  filter(Year == "2017") %>% 
  select(channel, MetDate, MeanPre2wksTemp, 
         N_uM, P_uM, NPratio,
         NEP_uM_C_gAFDM_h, R_uM_C_gAFDM_h, GPP_uM_C_gAFDM_h, 
         NUp_uM_N_gAFDM_h, Nfix_uM_N_gAFDM_h, TotNAssim_uM_N_gAFDM_h) %>% 
  # combining lower NP treatments into low
  mutate(NPratio = ifelse(NPratio %in% c("0.31", "0.93"), "Low", "High")) %>% 
  mutate(across(c(NPratio, N_uM, P_uM), factor)) %>% 
  mutate(N_uM = fct_relevel(N_uM, c("0.11", "3.68"))) %>% 
  mutate(P_uM = fct_relevel(P_uM, c("0.36", "3.94"))) %>% 
  mutate(NPratio = fct_relevel(NPratio, c("Low",  "High"))) %>% 
  pivot_longer(cols = NEP_uM_C_gAFDM_h:TotNAssim_uM_N_gAFDM_h, names_to = "Response", values_to = "Values") %>% 
  mutate(MetDate2 = as.factor(MetDate),
         MetDate2 = fct_reorder(MetDate2, MetDate),
         N_uM_date = as.factor(paste0(N_uM, "_", MetDate2)),
         P_uM_date = as.factor(paste0(P_uM,  "_",MetDate2)),
         NPratio_date = as.factor(paste0(NPratio,  "_",MetDate2)))




fun_ms_gams <- function(Response_i, ChanDF, ModFamily){
  # TESTING
  # Response_i <- "R_uM_C_gAFDM_h"
  # ChanDF <- chan15
  # ModFamily <- ERModFam
  # NOTE THIS WORKS (AND IS NEEDED) IN FUNCTION BUT NOT IN ISOLATION
  Response_i <- substitute(Response_i)
  ChanDF_i <- ChanDF %>% 
    filter(Response == Response_i)
  
  # All models have random intercept for date
  # splines for temp for each nut by date combo
  # N_uM_date combined the N concentration and date together into a factor with four levels
  g1_TxNxdate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = N_uM_date) + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  g2_TxPxdate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = P_uM_date) + s(MetDate2, bs = "re"), 
                         family = ModFamily,
                         data = ChanDF_i, 
                         method = "ML")
  
  g3_TxNPxdate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = NPratio_date) + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  # spline for temp for each nut
  g4_TxN <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = N_uM) + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  g5_TxP <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = P_uM) + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  g6_TxNP <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = NPratio) + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  # spline for temp by date + Nxdate as cat variable
  g7_T_date_pN_date <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = MetDate2) + N_uM_date + s(MetDate2, bs = "re"), 
                      family = ModFamily,
                      data = ChanDF_i, 
                      method = "ML")
  
  g8_T_date_pP_date <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = MetDate2) + P_uM_date + s(MetDate2, bs = "re"), 
                      family = ModFamily,
                      data = ChanDF_i, 
                      method = "ML")
  
  g9_T_date_pNP_date <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = MetDate2) + NPratio_date + s(MetDate2, bs = "re"), 
                       family = ModFamily,
                       data = ChanDF_i, 
                       method = "ML")
  
  
  # spline for temp + Nxdate as cat variable
  g10_TpN_date <- mgcv::gam(Values ~ s(MeanPre2wksTemp) + N_uM_date + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  g11_TpP_date <- mgcv::gam(Values ~ s(MeanPre2wksTemp) + P_uM_date + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  g12_TpNP_date <- mgcv::gam(Values ~ s(MeanPre2wksTemp) + NPratio_date + s(MetDate2, bs = "re"), 
                            family = ModFamily,
                            data = ChanDF_i, 
                            method = "ML")
  
  # spline for temp by date + N as cat variable
    g13_T_date_pN <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = MetDate2) + N_uM + s(MetDate2, bs = "re"), 
                       family = ModFamily,
                       data = ChanDF_i, 
                       method = "ML")
  
  g14_T_date_pP <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = MetDate2) + P_uM + s(MetDate2, bs = "re"), 
                       family = ModFamily,
                       data = ChanDF_i, 
                       method = "ML")
  
  g15_T_date_pNP <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = MetDate2) + NPratio + s(MetDate2, bs = "re"), 
                        family = ModFamily,
                        data = ChanDF_i, 
                        method = "ML")  
  
  # spline for temp + N as cat variable
  g16_TpN <- mgcv::gam(Values ~ s(MeanPre2wksTemp) + N_uM + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  g17_TpP <- mgcv::gam(Values ~ s(MeanPre2wksTemp) + P_uM + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  g18_TpNP <- mgcv::gam(Values ~ s(MeanPre2wksTemp) + NPratio + s(MetDate2, bs = "re"), 
                            family = ModFamily,
                            data = ChanDF_i, 
                            method = "ML")
  
  # spline for temp x date NO nut
  g19_T_date <- mgcv::gam(Values ~ s(MeanPre2wksTemp, by = MetDate2) + s(MetDate2, bs = "re"), 
                        family = ModFamily,
                        data = ChanDF_i, 
                        method = "ML")
  
  # spline for temp  NO nut
  g20_T <- mgcv::gam(Values ~ s(MeanPre2wksTemp) + s(MetDate2, bs = "re"), 
                          family = ModFamily,
                          data = ChanDF_i, 
                          method = "ML")
  
  # spline for nutrients NO temp
  g21_N <- mgcv::gam(Values ~ N_uM + s(MetDate2, bs = "re"), 
                     family = ModFamily,
                     data = ChanDF_i, 
                     method = "ML")
  
  
  g22_P <- mgcv::gam(Values ~ P_uM + s(MetDate2, bs = "re"), 
                     family = ModFamily,
                     data = ChanDF_i, 
                     method = "ML")
  
  
  g23_NP <- mgcv::gam(Values ~ NPratio + s(MetDate2, bs = "re"), 
                      family = ModFamily,
                      data = ChanDF_i, 
                      method = "ML")
  
  # only random effects
  g24_IntOnly <- mgcv::gam(Values ~ 1 + s(MetDate2, bs = "re"), 
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  
  
  
  # Create model selection table, rank based on AICc
  ModelSectionTable = as.data.frame(MuMIn::model.sel(g1_TxNxdate, g2_TxPxdate, g3_TxNPxdate, 
                                                     g4_TxN, g5_TxP, g6_TxNP,
                                                     g7_T_date_pN_date, g8_T_date_pP_date, g9_T_date_pNP_date,
                                                     g10_TpN_date, g11_TpP_date, g12_TpNP_date,
                                                     g13_T_date_pN, g14_T_date_pP, g15_T_date_pNP,
                                                     g16_TpN, g17_TpP, g18_TpNP,
                                                     g19_T_date,
                                                     g20_T, g21_N, g22_P, g23_NP, g24_IntOnly,
                                                     rank = AICc,
                                                     extra = c(R2 = function(x) summary(x)$r.sq)))
  
  ModelSectionTable
}

# Models ----


# Issues associated with two very high values for GPP and ER
# chanF[chanF$Response == "GPP_uM_C_gAFDM_h" & chanF$Values > 3000,]
# 1 C12     2017-07-11 00:00:00
# 2 C25     2017-07-11 00:00:00

# These two have very high areal GPP/ER rates, but moderate met AFDM/m2. They have the two highest AFDM/m2 for N uptake on this date
# I dealt with this above by using the N uptake AFDM/m2 for to estiamte ms rates for these two channels on this date.

## MS Respiration ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "R_uM_C_gAFDM_h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
ERModFam <- tw(link = "log") # gaussian(link = "log") 

# model selection
ER_ModSel <- fun_ms_gams(R_uM_C_gAFDM_h, chanF, ERModFam); ER_ModSel

# "Best model
chanF_ER <- chanF %>% 
  filter(Response == "R_uM_C_gAFDM_h")

ER_MostLikely <- gam(
  Values ~ s(MeanPre2wksTemp, by = P_uM) + s(MetDate2, bs = "re"),
                    family = ERModFam,
                    data = chanF_ER, 
                    method = "REML")

summary(ER_MostLikely)
check(getViz(ER_MostLikely))  
print(plot(getViz(ER_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(ER_MostLikely, "03_Model_RDS/MassSpec_2017_ER_mostlikely.rds")

## GPP ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "GPP_uM_C_gAFDM_h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
GPPModFam <- tw(link = "log") 
  # Tweedie(p = 1.55, link = power(0.1))
  #Gamma(link = "log") ## gaussian(link = "log") 

# model selectiona
# Again! Best model has P!
GPP_ModSel <- fun_ms_gams(GPP_uM_C_gAFDM_h, chanF, GPPModFam); GPP_ModSel

# "Best model
chanF_GPP <- chanF %>% 
  filter(Response == "GPP_uM_C_gAFDM_h")

GPP_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + P_uM_date + s(MetDate2, bs = "re"), 
                     family = GPPModFam,
                     data = chanF_GPP, 
                     method = "REML")

summary(GPP_MostLikely)
check(getViz(GPP_MostLikely))  
# adding P lowers GPP on one date, but doesn't do anything on the other
print(plot(getViz(GPP_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(GPP_MostLikely, "03_Model_RDS/MassSpec_2017_GPP_mostlikely.rds")

## NEP ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "NEP_uM_C_gAFDM_h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NEPModFam <- tw(link = "log") 

# model selection
NEP_ModSel <- fun_ms_gams(NEP_uM_C_gAFDM_h, chanF, NEPModFam); NEP_ModSel

# "Best model
chanF_NEP <- chanF %>% 
  filter(Response == "NEP_uM_C_gAFDM_h")

# P again!
NEP_MostLikely <-  gam(Values ~ s(MeanPre2wksTemp) + P_uM_date + s(MetDate2, bs = "re"),
                      family = NEPModFam,
                      data = chanF_NEP, 
                      method = "REML")

summary(NEP_MostLikely)
check(getViz(NEP_MostLikely))  
# most of the differences here are due to date
# there's a P effect on the first date, but not second.
print(plot(getViz(NEP_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(NEP_MostLikely, "03_Model_RDS/MassSpec_2017_NEP_mostlikely.rds")


## Nup ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "NUp_uM_N_gAFDM_h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NupModFam <- tw(link = "log") 

# model selection
Nup_ModSel <- fun_ms_gams(NUp_uM_N_gAFDM_h, chanF, NupModFam); Nup_ModSel

# "Best model
chanF_Nup <- chanF %>% 
  filter(Response == "NUp_uM_N_gAFDM_h")

# N
Nup_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM + s(MetDate2, bs = "re"),
                      family = NupModFam,
                      data = chanF_Nup, 
                      method = "REML")

summary(Nup_MostLikely)
check(getViz(Nup_MostLikely))  
# N up decreases with temp!
# less N uptake at high N:P ratios, but 'threshold' is in wrong place
print(plot(getViz(Nup_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(Nup_MostLikely, "03_Model_RDS/MassSpec_2017_Nup_mostlikely.rds")


## Nfix ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "Nfix_uM_N_gAFDM_h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NfixModFam <- tw(link = "log") 

# model selection
Nfix_ModSel <- fun_ms_gams(Nfix_uM_N_gAFDM_h, chanF, NfixModFam); Nfix_ModSel

# "Best model
chanF_Nfix <- chanF %>% 
  filter(Response == "Nfix_uM_N_gAFDM_h")

# N only here.
Nfix_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM_date + s(MetDate2, bs = "re"),
                      family = NfixModFam,
                      data = chanF_Nfix, 
                      method = "REML")

summary(Nfix_MostLikely)
check(getViz(Nfix_MostLikely))  
# complicated relationship with temp. 
# fixation is lower at high N on both dates. P or N:P doesn't matter
print(plot(getViz(Nfix_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(Nfix_MostLikely, "03_Model_RDS/MassSpec_2017_Nfix_mostlikely.rds")


## N assimilation ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "TotNAssim_uM_N_gAFDM_h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NassModFam <- tw(link = "log") 

# model selection
Nass_ModSel <- fun_ms_gams(TotNAssim_uM_N_gAFDM_h, chanF, NassModFam); Nass_ModSel

# "Best model
chanF_Nass <- chanF %>% 
  filter(Response == "TotNAssim_uM_N_gAFDM_h")

# Best model has N:P
Nass_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM + s(MetDate2, bs = "re"),
                       family = NassModFam,
                       data = chanF_Nass, 
                       method = "REML")

summary(Nass_MostLikely)
check(getViz(Nass_MostLikely))  
# Negative relationship between N assim and Temp - driven by N uptake?
# N assim is highest at lowest N:P. Threshold is in wrong place to support N fixation hypothesis
print(plot(getViz(Nass_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(Nass_MostLikely, "03_Model_RDS/MassSpec_2017_Nass_mostlikely.rds")


# save image ----
# save.image("02b_Script_SavedImages/08_2017_GAMs_MassSpecific.Rdata")
# load("02b_Script_SavedImages/08_2017_GAMs_MassSpecific.Rdata")
