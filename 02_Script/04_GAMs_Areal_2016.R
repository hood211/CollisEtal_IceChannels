# Areal responses for 2016
# JMH, Nov 2022
# NOTE: 2016 lacks uptake data
# proportion N fixer is modeled in 05_PropNfix_2015_2016

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
  mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d")) 

chanF <- chan %>% 
  filter(Year == "2016") %>% 
  select(channel, MetDate, MeanPre2wksTemp, NutTrt = P_uM, 
         NEP_uM_C_m2h, R_uM_C_m2h, GPP_uM_C_m2h, 
         Met_gAFDM_m2,
         Nfix_uM_N_m2h) %>% 
  pivot_longer(cols = NEP_uM_C_m2h:Nfix_uM_N_m2h, names_to = "Response", values_to = "Values") %>% 
  mutate(MetDate2 = as.factor(MetDate),
         MetDate2 = fct_reorder(MetDate2, MetDate))

chanAFDM <- chan %>% 
  filter(Year == "2016") %>% 
  select(channel, MetDate, MeanPre2wksTemp, NutTrt = P_uM,
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
  
  # Random Effects (calling intercept) only
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


## ER ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "R_uM_C_m2h"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  ERModFam <- tw(link = "log")
  
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
  saveRDS(ER_MostLikely, "03_Model_RDS/Areal_2016_ER_mostlikely.rds")
  
## GPP ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "GPP_uM_C_m2h"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  GPPModFam <- tw(link = "log")
  
  # model selection
  GPP_ModSel <- fun_ms_gams(GPP_uM_C_m2h, chanF %>% 
                                          # removing big outlier
                                          filter(Values < 75000), NutTrt, GPPModFam); GPP_ModSel
  
  # "Best model
  chanF_GPP <- chanF %>% 
    filter(Response == "GPP_uM_C_m2h")%>% 
    filter(Values < 75000)
  
  GPP_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                         s(NutTrt, bs = "ts", k = 5, by = MetDate2) + 
                         s(MetDate2, bs = "re"),
                      family = GPPModFam,
                      data = chanF_GPP)
  
  summary(GPP_MostLikely)
  # there's one big outlier (C23 [high T med Nut]; day 2), do I get the same model if I remove it?
  # g5 and g6 switch, N by date becomes better
  # changes the nature of the P curve

  check(getViz(GPP_MostLikely))  
  print(plot(getViz(GPP_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(GPP_MostLikely, "03_Model_RDS/Areal_2016_GPP_mostlikely.rds")
  
## NEP ----
  # qq line plot
  ggplot(chanF %>% 
           filter(Response == "NEP_uM_C_m2h"), aes(sample = log(Values))) +
    stat_qq() +
    stat_qq_line()
  
  # model family
  NEPModFam <- tw(link = "log")
  
  # model selection
  NEP_ModSel <- fun_ms_gams(NEP_uM_C_m2h, chanF%>% 
                              # removing big outlier, same as for GPP
                              filter(Values < 75000), NutTrt, NEPModFam); NEP_ModSel
  
  # "Best model
  chanF_NEP <- chanF %>% 
    filter(Response == "NEP_uM_C_m2h")%>% 
    # removing big outlier, same as for GPP
    filter(Values < 75000)
  
  NEP_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                         s(NutTrt, bs = "ts", k = 5, by = MetDate2) + 
                         s(MetDate2, bs = "re"),
                       family = NEPModFam,
                       data = chanF_NEP)
  
  summary(NEP_MostLikely)
  check(getViz(NEP_MostLikely))  
  print(plot(getViz(NEP_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(NEP_MostLikely, "03_Model_RDS/Areal_2016_NEP_mostlikely.rds")
  
## N uptake ----
  # N uptake data lost
  

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
  # additive no date is similar, P spline is not sig
  chanF_Nfix <- chanF %>% 
    filter(Response == "Nfix_uM_N_m2h")
  
  Nfix_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                          s(MetDate2, bs = "re"),
                       family = NfixModFam,
                       data = chanF_Nfix)
  
  summary(Nfix_MostLikely)
  check(getViz(Nfix_MostLikely))  
  print(plot(getViz(Nfix_MostLikely), allTerms = T), pages = 1)
  
  # output model
  saveRDS(Nfix_MostLikely, "03_Model_RDS/Areal_2016_Nfix_mostlikely.rds")

## N assimilation ----
  # N uptake samples lost in freezer malfunction
  
  

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
  
  # most likely
  AFDM_MostLikely = gam(Values ~ te(MeanPre2wksTemp, NutTrt, bs = "ts", by = MetDate2) + s(MetDate2, bs = "re"),
                        family = AFDM_ModFam,
                        data = chanAFDM2,
                        select = T)
  
  # Temp only by date: not nearly as good based on AIC
  # But I feel like the ts may be overfitting
  
  AFDM_MostLikely2 = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2)+ s(MetDate2, bs = "re"),
                        family = AFDM_ModFam,
                        data = chanAFDM2)
  
  summary(AFDM_MostLikely2)
  check(getViz(AFDM_MostLikely2))  
  print(plot(getViz(AFDM_MostLikely), allTerms = T), pages = 1)
  
  # output model
  # using the temp only model here, see justification in text
  saveRDS(AFDM_MostLikely2, "03_Model_RDS/Areal_2016_AFDM_mostlikely.rds")  
  
  
# save.image ----
  # save.image("02b_Script_SavedImages/04_2016_GAMs_Areal_Rdat")
  # load("02b_Script_SavedImages/04_2016_GAMs_Areal_Rdat")
  