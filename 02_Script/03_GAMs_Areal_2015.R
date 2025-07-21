# Areal responses for 2015
# JMH, Nov 2022

# Libraries ----
# general
library(tidyverse)
# gam models
library(mgcv)
library(mgcViz)
library(MuMIn)
# tables
library(sjPlot)
library(stargazer)

# Data ----
chan <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv", row.names = 1) %>% 
  mutate(across(c(MetDate, UpDate, NfixDate), as.POSIXct, format = "%m/%d/%Y")) %>% 
  mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d"),
         # Adding 2015 minimum to Nuptake to allow modeling
         NUp_uM_N_m2_hr = NUp_uM_N_m2_hr + 202,
         TotAssim_uM_N_m2_hr = NUp_uM_N_m2_hr + Nfix_uM_N_m2h) 

chanF <- chan %>% 
  filter(Year == "2015") %>% 
  select(channel, MetDate, MeanPre2wksTemp, NutTrt = N_uM, 
         NEP_uM_C_m2h, R_uM_C_m2h, GPP_uM_C_m2h, 
         Met_gAFDM_m2,
         NUp_uM_N_m2_hr, Nfix_uM_N_m2h, TotAssim_uM_N_m2_hr,
         propNfixer) %>% 
  pivot_longer(cols = NEP_uM_C_m2h:propNfixer, names_to = "Response", values_to = "Values") %>% 
  mutate(MetDate2 = as.factor(MetDate),
         MetDate2 = fct_reorder(MetDate2, MetDate))

chanAFDM <- chan %>% 
  filter(Year == "2015") %>% 
  select(channel, MetDate, MeanPre2wksTemp, NutTrt = N_uM,
         Met_gAFDM_m2, Nup_gAFDM_m2, Pup_gAFDM_m2) %>% 
  pivot_longer(cols = c(Met_gAFDM_m2:Pup_gAFDM_m2), names_to = "AFDM_smp", values_to = "Values") %>% 
  mutate(AFDM_smp = as.factor(AFDM_smp),
         Response = "gAFDM_m2",
         MetDate2 = as.factor(MetDate),
         MetDate2 = fct_reorder(MetDate2, MetDate))



# Model section ----

## Fun for model selection
fun_ms_gams <- function(Response_i, ChanDF, NutTrt, ModFamily){
  # TESTING
  # Response_i <- "R_uM_C_gAFDM_h"
  # ChanDF <- chan15
  # ModFamily <- ERModFam
  # NOTE THIS WORKS (AND IS NEEDED) IN FUNCTION BUT NOT IN ISOLATION
  Response_i <- substitute(Response_i)
  ChanDF_i <- ChanDF %>% 
    filter(Response == Response_i)
  
  # All models have random intercept for date
  # Temp x N interaction - for each date
  g1_ResSurf_ByDate <- mgcv::gam(Values ~ te(MeanPre2wksTemp, NutTrt, bs = "ts", by = MetDate2) + s(MetDate2, bs = "re"),
                                 family = ModFamily,
                                 data = ChanDF_i, 
                                 method = "ML")
  
  # Temp x N interaction - dates combined
  g2_ResSurf <-  mgcv::gam(Values ~ te(MeanPre2wksTemp, NutTrt, bs = "ts") + s(MetDate2, bs = "re"),
                           family = ModFamily,
                           data = ChanDF_i, 
                           method = "ML")
  
  # Additive temp and N-- both splines separated by date
  g3_Additive_BothByDate <-  mgcv::gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                                         s(NutTrt, bs = "ts", k = 5, by = MetDate2) + 
                                         s(MetDate2, bs = "re"), 
                                       family = ModFamily,
                                       data = ChanDF_i, 
                                       method = "ML")
  
  # Additive temp and N-- temp spline separated by date
  g4_Additive_TbyDate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                                     s(NutTrt, bs = "ts", k = 5) + 
                                     s(MetDate2, bs = "re"), 
                                   family = ModFamily,
                                   data = ChanDF_i, 
                                   method = "ML")
  
  # Additive temp and N-- N spline separated by date
  g5_Additive_NbyDate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                                     s(NutTrt, bs = "ts", k = 5, by = MetDate2) + 
                                     s(MetDate2, bs = "re"), 
                                   family = ModFamily,
                                   data = ChanDF_i, 
                                   method = "ML")
  
  # Additive temp and N-- single spline for T and N
  g6_Additive_NoDate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                                    s(NutTrt, bs = "ts", k = 5) + 
                                    s(MetDate2, bs = "re"), 
                                  family = ModFamily,
                                  data = ChanDF_i, 
                                  method = "ML")
  
  # Temp only - by date
  g7_TempOnly_ByDate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                                    s(MetDate2, bs = "re"), 
                                  family = ModFamily,
                                  data = ChanDF_i, 
                                  method = "ML")
  
  # Temp only - NoDate
  g8_TempOnly_NoDate <- mgcv::gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                                    s(MetDate2, bs = "re"), 
                                  family = ModFamily,
                                  data = ChanDF_i, 
                                  method = "ML")
  
  # N only - by date
  g9_NitrogenOnly_byDate <- mgcv::gam(Values ~ s(NutTrt, bs = "ts", by = MetDate2, k = 5) + 
                                        s(MetDate2, bs = "re"), 
                                      family = ModFamily,
                                      data = ChanDF_i, 
                                      method = "ML")
  
  # N only - No date
  g10_NitrogenOnly_NoDate <- mgcv::gam(Values ~ s(NutTrt, bs = "ts", k = 5) + 
                                         s(MetDate2, bs = "re"), 
                                       family = ModFamily,
                                       data = ChanDF_i, 
                                       method = "ML")
  
  # Intercept only
  g11_InterceptOnly = mgcv::gam(Values ~ 1 + s(MetDate2, bs = "re"), 
                                family = ModFamily,
                                data = ChanDF_i, 
                                method = "ML")
  
  # Create model selection table, rank based on AICc
  ModelSectionTable = as.data.frame(MuMIn::model.sel(g1_ResSurf_ByDate, g2_ResSurf, 
                                                     g3_Additive_BothByDate, g4_Additive_TbyDate, g5_Additive_NbyDate, g6_Additive_NoDate, 
                                                     g7_TempOnly_ByDate, g8_TempOnly_NoDate,
                                                     g9_NitrogenOnly_byDate, g10_NitrogenOnly_NoDate,
                                                     g11_InterceptOnly, 
                                                     rank = AICc,
                                                     extra = c(R2 = function(x) summary(x)$r.sq)))
  
  ModelSectionTable
}


## Respiration ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "R_uM_C_m2h"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  ERModFam <- tw(link = "log")#gaussian(link = "log")#tw(link = "log")
  
  # model selection
  ER_ModSel <- fun_ms_gams(R_uM_C_m2h, chanF, NutTrt, ERModFam); ER_ModSel
  
  # "Best model
  chanF_ER <- chanF %>% 
    filter(Response == "R_uM_C_m2h")
  
  ER_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                        s(NutTrt, bs = "ts", k = 5, by = MetDate2) + 
                        s(MetDate2, bs = "re"),
                      family = ERModFam,
                      data = chanF_ER)
  
  summary(ER_MostLikely)
  check(getViz(ER_MostLikely))  
  print(plot(getViz(ER_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(ER_MostLikely, "03_Model_RDS/Areal_2015_ER_mostlikely.rds")
  
## GPP ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "GPP_uM_C_m2h"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  GPPModFam <- tw(link = "log")#gaussian(link = "log")#tw(link = "log")
  
  # model selection
  GPP_ModSel <- fun_ms_gams(GPP_uM_C_m2h, chanF, NutTrt, GPPModFam); GPP_ModSel
  
  # "Best model
  chanF_GPP <- chanF %>% 
    filter(Response == "GPP_uM_C_m2h")
  
  GPP_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                         s(NutTrt, bs = "ts", k = 5) + 
                         s(MetDate2, bs = "re"),
                      family = GPPModFam,
                      data = chanF_GPP)
  
  summary(GPP_MostLikely)
  check(getViz(GPP_MostLikely))  
  print(plot(getViz(GPP_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(GPP_MostLikely, "03_Model_RDS/Areal_2015_GPP_mostlikely.rds")
  
## NEP ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "NEP_uM_C_m2h"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  NEPModFam <- tw(link = "log")
  
  # model selection
  NEP_ModSel <- fun_ms_gams(NEP_uM_C_m2h, chanF, NutTrt, NEPModFam); NEP_ModSel
  
  # "Best model
  chanF_NEP <- chanF %>% 
    filter(Response == "NEP_uM_C_m2h")
  
  NEP_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                         s(NutTrt, bs = "ts", k = 5) + 
                         s(MetDate2, bs = "re"),
                       family = NEPModFam,
                       data = chanF_NEP)
  
  summary(NEP_MostLikely)
  check(getViz(NEP_MostLikely))  
  print(plot(getViz(NEP_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(NEP_MostLikely, "03_Model_RDS/Areal_2015_NEP_mostlikely.rds")
  
## N uptake ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "NUp_uM_N_m2_hr"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  NupModFam <- tw(link = "log")
  
  # model selection
  Nup_ModSel <- fun_ms_gams(NUp_uM_N_m2_hr, chanF, NutTrt, NupModFam); Nup_ModSel
  
  # "Best model
  chanF_Nup <- chanF %>% 
    filter(Response == "NUp_uM_N_m2_hr")
  
  Nup_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                         s(NutTrt, bs = "ts", k = 5) + 
                         s(MetDate2, bs = "re"),
                       family = NupModFam,
                       data = chanF_Nup)
  
  summary(Nup_MostLikely)
  check(getViz(Nup_MostLikely))  
  print(plot(getViz(Nup_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(Nup_MostLikely, "03_Model_RDS/Areal_2015_Nup_mostlikely.rds")

## N fixation ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "Nfix_uM_N_m2h"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  NfixModFam <- tw(link = "log")
  
  # model selection
  Nfix_ModSel <- fun_ms_gams(Nfix_uM_N_m2h, chanF, NutTrt, NfixModFam); Nfix_ModSel
  
  # "Best model
  chanF_Nfix <- chanF %>% 
    filter(Response == "Nfix_uM_N_m2h")
  
  Nfix_MostLikely = gam(Values ~ te(MeanPre2wksTemp, NutTrt, bs = "ts") + s(MetDate2, bs = "re"),
                       family = NfixModFam,
                       data = chanF_Nfix)
  
  summary(Nfix_MostLikely)
  check(getViz(Nfix_MostLikely))  
  print(plot(getViz(Nfix_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(Nfix_MostLikely, "03_Model_RDS/Areal_2015_Nfix_mostlikely.rds")

## N assimilation ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "TotAssim_uM_N_m2_hr"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  TotAssimModFam <- tw(link = "log")
  
  # model selection
  TotAssim_ModSel <- fun_ms_gams(TotAssim_uM_N_m2_hr, chanF, NutTrt, TotAssimModFam); TotAssim_ModSel
  
  # "Best model
  chanF_TotAssim <- chanF %>% 
    filter(Response == "TotAssim_uM_N_m2_hr")
  
  TotAssim_MostLikely = gam(Values ~ te(MeanPre2wksTemp, NutTrt, bs = "ts") + s(MetDate2, bs = "re"),
                        family = TotAssimModFam,
                        data = chanF_TotAssim)
  
  summary(TotAssim_MostLikely)
  check(getViz(TotAssim_MostLikely))  
  print(plot(getViz(TotAssim_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(TotAssim_MostLikely, "03_Model_RDS/Areal_2015_NTotAssim_mostlikely.rds")
  
  ## Biomass ----
  # qq line plot
  ggplot(chanAFDM %>% 
           filter(Response == "gAFDM_m2"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  AFDM_ModFam <- tw(link = "log")
  
  # model selection
  AFDM_ModSel <- fun_ms_gams(gAFDM_m2, chanAFDM, NutTrt, AFDM_ModFam); AFDM_ModSel
  
  # "Best model
  chanAFDM2 <- chanAFDM %>% 
    filter(Response == "gAFDM_m2")
  
  AFDM_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                          s(NutTrt, bs = "ts", k = 5) + 
                          s(MetDate2, bs = "re"),
                            family = AFDM_ModFam,
                            data = chanAFDM2)
  
  summary(AFDM_MostLikely)
  check(getViz(AFDM_MostLikely))  
  print(plot(getViz(AFDM_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(AFDM_MostLikely, "03_Model_RDS/Areal_2015_AFDM_mostlikely.rds")
  
  
  
# save.images ----
  # save.image("02b_Script_savedImages/03_2015_GAMs_Areal_Rdat")
  # load("02b_Script_savedImages/03_2015_GAMs_Areal_Rdat")
  