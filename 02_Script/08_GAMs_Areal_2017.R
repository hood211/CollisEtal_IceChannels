# Areal responses for 2017
# JMH, Dec 2022, Jul 2024

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
  mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d"),
         TotAssim_uM_N_m2_hr = NUp_uM_N_m2_hr + Nfix_uM_N_m2h)

chanF <- chan %>% 
  filter(Year == "2017") %>% 
  select(channel, MetDate, MeanPre2wksTemp, 
         N_uM, P_uM, NPratio,
         NEP_uM_C_m2h, R_uM_C_m2h, GPP_uM_C_m2h, 
         Met_gAFDM_m2,
         NUp_uM_N_m2_hr, Nfix_uM_N_m2h, TotAssim_uM_N_m2_hr,
         propNfixer) %>% 
  # combining lower NP treatments into low
  mutate(NPratio = ifelse(NPratio %in% c("0.31", "0.93"), "Low", "High")) %>% 
  mutate(across(c(NPratio, N_uM, P_uM), factor)) %>% 
  mutate(N_uM = fct_relevel(N_uM, c("0.11", "3.68"))) %>% 
  mutate(P_uM = fct_relevel(P_uM, c("0.36", "3.94"))) %>% 
  mutate(NPratio = fct_relevel(NPratio, c("Low",  "High"))) %>% 
  pivot_longer(cols = NEP_uM_C_m2h:propNfixer, names_to = "Response", values_to = "Values") %>% 
  mutate(MetDate2 = as.factor(MetDate),
         MetDate2 = fct_reorder(MetDate2, MetDate),
         N_uM_date = as.factor(paste0(N_uM, "_", MetDate2)),
         P_uM_date = as.factor(paste0(P_uM,  "_",MetDate2)),
         NPratio_date = as.factor(paste0(NPratio,  "_",MetDate2)))

chanAFDM <- chan %>% 
  filter(Year == "2017") %>% 
  select(channel, MetDate, MeanPre2wksTemp, 
         N_uM, P_uM, NPratio,
         Met_gAFDM_m2, Nup_gAFDM_m2, Pup_gAFDM_m2) %>% 
  # combining lower NP treatments into low
  mutate(NPratio = ifelse(NPratio %in% c("0.31", "0.93"), "Low", "High")) %>% 
  mutate(across(c(NPratio, N_uM, P_uM), factor)) %>% 
  mutate(N_uM = fct_relevel(N_uM, c("0.11", "3.68"))) %>% 
  mutate(P_uM = fct_relevel(P_uM, c("0.36", "3.94"))) %>% 
  mutate(NPratio = fct_relevel(NPratio, c("Low",  "High"))) %>% 
  pivot_longer(cols = c(Met_gAFDM_m2:Pup_gAFDM_m2), names_to = "AFDM_smp", values_to = "Values") %>% 
  mutate(AFDM_smp = as.factor(AFDM_smp),
         Response = "gAFDM_m2",
         MetDate2 = as.factor(MetDate),
         MetDate2 = fct_reorder(MetDate2, MetDate),
         N_uM_date = as.factor(paste0(N_uM, "_", MetDate2)),
         P_uM_date = as.factor(paste0(P_uM,  "_",MetDate2)),
         NPratio_date = as.factor(paste0(NPratio,  "_",MetDate2))) %>% 
  # Outlier - C12, 2017-07-11, Nup
  filter(Values < 150)


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
## Respiration ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "R_uM_C_m2h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
ERModFam <- tw(link = "log") 

# model selection
ER_ModSel <- fun_ms_gams(R_uM_C_m2h, chanF, ERModFam); ER_ModSel

# "Best model
chanF_ER <- chanF %>% 
  filter(Response == "R_uM_C_m2h")

ER_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM_date + s(MetDate2, bs = "re"), 
                    family = ERModFam,
                    data = chanF_ER, 
                    method = "REML")

summary(ER_MostLikely)
check(getViz(ER_MostLikely))  
print(plot(getViz(ER_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(ER_MostLikely, "03_Model_RDS/Areal_2017_ER_mostlikely.rds")

## GPP ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "GPP_uM_C_m2h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
GPPModFam <- tw(link = "log") 

# model selection
GPP_ModSel <- fun_ms_gams(GPP_uM_C_m2h, chanF, GPPModFam); GPP_ModSel

# "Best model
chanF_GPP <- chanF %>% 
  filter(Response == "GPP_uM_C_m2h")

GPP_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM_date + s(MetDate2, bs = "re"), 
                     family = GPPModFam,
                     data = chanF_GPP, 
                     method = "REML")

summary(GPP_MostLikely)
check(getViz(GPP_MostLikely))  
print(plot(getViz(GPP_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(GPP_MostLikely, "03_Model_RDS/Areal_2017_GPP_mostlikely.rds")

## NEP ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "NEP_uM_C_m2h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NEPModFam <- tw(link = "log") 

# model selection
NEP_ModSel <- fun_ms_gams(NEP_uM_C_m2h, chanF, NEPModFam); NEP_ModSel

# "Best model
chanF_NEP <- chanF %>% 
  filter(Response == "NEP_uM_C_m2h")

# model with N:P is similar, however, the big difference is between the two lowest sub 1 N:P ratios, 
# not consistent with expectations
NEP_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM_date + s(MetDate2, bs = "re"), 
                      family = NEPModFam,
                      data = chanF_NEP, 
                      method = "REML")

summary(NEP_MostLikely)
check(getViz(NEP_MostLikely))  
print(plot(getViz(NEP_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(NEP_MostLikely, "03_Model_RDS/Areal_2017_NEP_mostlikely.rds")

## Nup ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "NUp_uM_N_m2_hr"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NupModFam <- tw(link = "log") 

# model selection
Nup_ModSel <- fun_ms_gams(NUp_uM_N_m2_hr, chanF, NupModFam); Nup_ModSel

# "Best model
chanF_Nup <- chanF %>% 
  filter(Response == "NUp_uM_N_m2_hr")

# N:P model is most likely, but % dev explained with N is 8% more.
# Dif in AIC due to edf of MetDate?
Nup_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM + s(MetDate2, bs = "re"),
                      family = NupModFam,
                      data = chanF_Nup, 
                      method = "REML")

summary(Nup_MostLikely)
check(getViz(Nup_MostLikely))  
print(plot(getViz(Nup_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(Nup_MostLikely, "03_Model_RDS/Areal_2017_Nup_mostlikely.rds")


## Nfix ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "Nfix_uM_N_m2h"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NfixModFam <- tw(link = "log") 

# model selection
Nfix_ModSel <- fun_ms_gams(Nfix_uM_N_m2h, chanF, NfixModFam); Nfix_ModSel

# "Best model
chanF_Nfix <- chanF %>% 
  filter(Response == "Nfix_uM_N_m2h")

# first model with N:P >2
Nfix_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM + s(MetDate2, bs = "re"),
                      family = NfixModFam,
                      data = chanF_Nfix, 
                      method = "REML")

summary(Nfix_MostLikely)
check(getViz(Nfix_MostLikely))  
print(plot(getViz(Nfix_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(Nfix_MostLikely, "03_Model_RDS/Areal_2017_Nfix_mostlikely.rds")


## N assimilation ----
# qq line plot
ggplot(chanF %>% 
         filter(Response == "TotAssim_uM_N_m2_hr"), aes(sample = log(Values))) +
  stat_qq() +
  stat_qq_line()

# model family
NassModFam <- tw(link = "log") 

# model selection
Nass_ModSel <- fun_ms_gams(TotAssim_uM_N_m2_hr, chanF, NassModFam); Nass_ModSel

# "Best model
chanF_Nass <- chanF %>% 
  filter(Response == "TotAssim_uM_N_m2_hr")

# first model with N:P #2, < 2
# again, big differences is between 0.31 and 0.93 indicating that N drives this

Nass_MostLikely <- gam(Values ~ s(MeanPre2wksTemp) + N_uM + s(MetDate2, bs = "re"),
                       family = NassModFam,
                       data = chanF_Nass, 
                       method = "REML")

summary(Nass_MostLikely)
check(getViz(Nass_MostLikely))  
print(plot(getViz(Nass_MostLikely), allTerms = T), pages = 1)

# output model
saveRDS(Nass_MostLikely, "03_Model_RDS/Areal_2017_Nass_mostlikely.rds")


## Biomass ----
  # qq line plot
  ggplot(chanAFDM %>% 
           filter(Response == "gAFDM_m2"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  AFDM_ModFam <- tw(link = "log") # Gamma(link = "log") #
  
  # model selection
  # removal of outlier changes from g7 to g12
  AFDM_ModSel <- fun_ms_gams(gAFDM_m2, chanAFDM, AFDM_ModFam); AFDM_ModSel
  
  # "Best model
  chanAFDM2 <- chanAFDM %>% 
    filter(Response == "gAFDM_m2")
  
  # most likely
  AFDM_MostLikely = gam(Values ~ s(MeanPre2wksTemp) + N_uM + s(MetDate2, bs = "re"),
                        family = AFDM_ModFam,
                        data = chanAFDM2)
  
  summary(AFDM_MostLikely)
  check(getViz(AFDM_MostLikely))  
  print(plot(getViz(AFDM_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(AFDM_MostLikely, "03_Model_RDS/Areal_2017_AFDM_mostlikely.rds")  


# save image ----
# save.image("02b_Script_SavedImages/08_2017_GAMs_Areal.Rdata")
# load("02b_Script_SavedImages/08_2017_GAMs_Areal.Rdata")
