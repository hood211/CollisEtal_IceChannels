# Mass-specific responses for 2015
# JMH, Nov 2022

# Libraries ----
library(tidyverse)
library(mgcv)
library(mgcViz)
library(MuMIn)

# Data ----
  chan <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed.csv", row.names = 1) %>% 
            mutate(across(c(MetDate, UpDate, NfixDate), as.POSIXct, format = "%m/%d/%Y")) %>% 
            mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d"),
                   NEP_uM_C_gAFDM_h = NEP_uM_C_m2h/Met_gAFDM_m2,
                   R_uM_C_gAFDM_h = R_uM_C_m2h/Met_gAFDM_m2,
                   GPP_uM_C_gAFDM_h = GPP_uM_C_m2h/Met_gAFDM_m2,
                   # N uptake has negative values so adding min
                   NUp_uM_N_gAFDM_h = (NUp_uM_N_m2_hr + 202)/Nup_gAFDM_m2,
                   Nfix_uM_N_gAFDM_h = Nfix_uM_N_m2h/Nfix_gAFDM_m2,
                   TotNAssim_uM_N_gAFDM_h = NUp_uM_N_gAFDM_h + Nfix_uM_N_gAFDM_h) 
  
  chan15 <- chan %>% 
    filter(Year == "2015") %>% 
    select(channel, MetDate, MeanPre2wksTemp, NutTrt = N_uM, NEP_uM_C_gAFDM_h:GPP_uM_C_gAFDM_h, 
           NUp_uM_N_gAFDM_h, Nfix_uM_N_gAFDM_h, TotNAssim_uM_N_gAFDM_h) %>% 
    pivot_longer(cols = NEP_uM_C_gAFDM_h:TotNAssim_uM_N_gAFDM_h, names_to = "Response", values_to = "Values") %>% 
    mutate(MetDate2 = as.factor(MetDate),
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
  
  
  ## ER ----
    # qq line plot
    ggplot(chan15 %>% 
             filter(Response == "R_uM_C_gAFDM_h"), aes(sample = log(Values))) +
              stat_qq() +
              stat_qq_line()
    
    # model family
    ERModFam <- Gamma(link = "log") 
    
    # model selection
    ER_ModSel <- fun_ms_gams(R_uM_C_gAFDM_h, chan15, NutTrt, ERModFam); ER_ModSel
    
    # "Best model
    chan15_ER <- chan15 %>% 
      filter(Response == "R_uM_C_gAFDM_h")
    
    # this isn't most likely, but it has del AIC < 2 and an R2 0.05 from top model
    ER_MostLikely = gam(Values ~ 
                          s(NutTrt, bs = "ts", k = 5) + 
                          s(MetDate2, bs = "re"),
                        family = ERModFam,
                        data = chan15_ER)
    
    summary(ER_MostLikely)
    check(getViz(ER_MostLikely))  
    print(plot(getViz(ER_MostLikely), allTerms = T), pages = 1)
    
    # output model
    saveRDS(ER_MostLikely, "03_Model_RDS/MS_2015_ER_mostlikely.rds")
    
  ## GPP ----
    # qq line plot
    ggplot(chan15 %>% 
             filter(Response == "GPP_uM_C_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    GPPModFam <- Gamma(link = "log")
    
    # model selection
    GPP_ModSel <- fun_ms_gams(GPP_uM_C_gAFDM_h, chan15, NutTrt, GPPModFam); GPP_ModSel
    # seems like limited support for the T model
    
    # "Best model
    chan15_GPP <- chan15 %>% 
      filter(Response == "GPP_uM_C_gAFDM_h")
    
    GPP_MostLikely = gam(Values ~ s(NutTrt, bs = "ts", k = 5) + 
                           s(MetDate2, bs = "re"),
                        family = GPPModFam,
                        data = chan15_GPP)
    
    summary(GPP_MostLikely)
    check(getViz(GPP_MostLikely))  
    print(plot(getViz(GPP_MostLikely), allTerms = T), pages = 1)
    
    # output model
    saveRDS(GPP_MostLikely, "03_Model_RDS/MS_2015_GPP_mostlikely.rds")
    
## NEP ----
    # qq line plot
    ggplot(chan15 %>% 
             filter(Response == "NEP_uM_C_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    NEPModFam <- Gamma(link = "log")
    
    # model selection
    NEP_ModSel <- fun_ms_gams(NEP_uM_C_gAFDM_h, chan15, NutTrt, NEPModFam); NEP_ModSel
    # seems like limited support for the T model
    
    # "Best model
    chan15_NEP <- chan15 %>% 
      filter(Response == "NEP_uM_C_gAFDM_h")
    
    NEP_MostLikely = gam(Values ~ s(NutTrt, bs = "ts", k = 5) + 
                           s(MetDate2, bs = "re"),
                         family = NEPModFam,
                         data = chan15_NEP)
    
    summary(NEP_MostLikely)
    check(getViz(NEP_MostLikely))  
    print(plot(getViz(NEP_MostLikely), allTerms = T), pages = 1)
    
    # output model
    saveRDS(NEP_MostLikely, "03_Model_RDS/MS_2015_NEP_mostlikely.rds")
    
 ## N uptake ----
    # qq line plot
    ggplot(chan15 %>% 
             filter(Response == "NUp_uM_N_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    NupModFam <- tw(link = "log")
    
    # model selection
    Nup_ModSel <- fun_ms_gams(NUp_uM_N_gAFDM_h, chan15, NutTrt, NupModFam); Nup_ModSel
    # seems like limited support for the T model
    
    # "Best model
    chan15_Nup <- chan15 %>% 
      filter(Response == "NUp_uM_N_gAFDM_h")
    
    # the neg N uptake value is a bit of an outlier, but the results are the same with or without
    Nup_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                           s(NutTrt, bs = "ts", k = 5) + 
                           s(MetDate2, bs = "re"),
                         family = NupModFam,
                         data = chan15_Nup)
    
    summary(Nup_MostLikely)
    check(getViz(Nup_MostLikely))  
    print(plot(getViz(Nup_MostLikely), allTerms = T), pages = 1)
    
    # output model
    saveRDS(Nup_MostLikely, "03_Model_RDS/MS_2015_Nup_mostlikely.rds")
    
## N fixation ----
    # qq line plot
    ggplot(chan15 %>% 
             filter(Response == "Nfix_uM_N_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    NfixModFam <- tw(link = "log")
    
    # model selection
    Nfix_ModSel <- fun_ms_gams(Nfix_uM_N_gAFDM_h, chan15, NutTrt, NfixModFam); Nfix_ModSel
    # seems like limited support for the T model
    
    # "Best model
    chan15_Nfix <- chan15 %>% 
      filter(Response == "Nfix_uM_N_gAFDM_h")
    
    Nfix_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                            s(NutTrt, bs = "ts", k = 5) + 
                            s(MetDate2, bs = "re"),
                         family = NfixModFam,
                         data = chan15_Nfix)
    
    summary(Nfix_MostLikely)
    check(getViz(Nfix_MostLikely))  
    print(plot(getViz(Nfix_MostLikely), allTerms = T), pages = 1)
    
    
    # output model
    saveRDS(Nfix_MostLikely, "03_Model_RDS/MS_2015_Nfix_mostlikely.rds")

 ## Total N assimilation ----
    # qq line plot
    ggplot(chan15 %>% 
             filter(Response == "TotNAssim_uM_N_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    NAssimModFam <- tw(link = "log")
    
    # model selection
    NAssim_ModSel <- fun_ms_gams(TotNAssim_uM_N_gAFDM_h, chan15, NutTrt, NAssimModFam); NAssim_ModSel
    # seems like limited support for the T model
    
    # "Best model
    chan15_NAssim <- chan15 %>% 
      filter(Response == "TotNAssim_uM_N_gAFDM_h")
    
    NAssim_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts", by = MetDate2) + 
                              s(NutTrt, bs = "ts", k = 5) + 
                              s(MetDate2, bs = "re"),
                          family = NAssimModFam,
                          data = chan15_NAssim)
    
    summary(NAssim_MostLikely)
    check(getViz(NAssim_MostLikely))  
    print(plot(getViz(NAssim_MostLikely), allTerms = T), pages = 1)

    # output model
    saveRDS(NAssim_MostLikely, "03_Model_RDS/MS_2015_NAssim_mostlikely.rds")
    
# Save image
    # save.image("02b_Script_SavedImages/O1_2015_GAMS_massSpecific_Rdat")
    # load("02b_Script_SavedImages/O1_2015_GAMS_massSpecific_Rdat")
    