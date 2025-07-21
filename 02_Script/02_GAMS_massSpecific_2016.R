# Mass-specific responses for 2015
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
  chan <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed.csv", row.names = 1) %>% 
            mutate(across(c(MetDate, UpDate, NfixDate), as.POSIXct, format = "%m/%d/%Y")) %>% 
            mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d"),
                   NEP_uM_C_gAFDM_h = NEP_uM_C_m2h/Met_gAFDM_m2,
                   R_uM_C_gAFDM_h = R_uM_C_m2h/Met_gAFDM_m2,
                   GPP_uM_C_gAFDM_h = GPP_uM_C_m2h/Met_gAFDM_m2,
                   # N uptake data was lost
                   Nfix_uM_N_gAFDM_h = Nfix_uM_N_m2h/Nfix_gAFDM_m2) 
  

  chan16 <- chan %>% 
    filter(Year == "2016") %>% 
    select(channel, MetDate, MeanPre2wksTemp, NutTrt = P_uM, NEP_uM_C_gAFDM_h:GPP_uM_C_gAFDM_h, 
           Nfix_uM_N_gAFDM_h) %>% 
    pivot_longer(cols = NEP_uM_C_gAFDM_h:Nfix_uM_N_gAFDM_h, names_to = "Response", values_to = "Values") %>% 
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
  
  
  ## Respiration ----
    # qq line plot
    ggplot(chan16 %>% 
             filter(Response == "R_uM_C_gAFDM_h"), aes(sample = log(Values))) +
              stat_qq() +
              stat_qq_line()
    
    # model family
    ERModFam <- Gamma(link = "log")
    
    # model selection
    ER_ModSel <- fun_ms_gams(R_uM_C_gAFDM_h, chan16, NutTrt, ERModFam); ER_ModSel
    
    # "Best model
    chan16_ER <- chan16 %>% 
      filter(Response == "R_uM_C_gAFDM_h")
    
    ER_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                          s(MetDate2, bs = "re"),
                        family = ERModFam,
                        data = chan16_ER)
    
    summary(ER_MostLikely)
    check(getViz(ER_MostLikely))  
    print(plot(getViz(ER_MostLikely), allTerms = T), pages = 1)
    
    # output model
    saveRDS(ER_MostLikely, "03_Model_RDS/MS_2016_ER_mostlikely.rds")
    
  ## GPP ----
    # qq line plot
    ggplot(chan16 %>% 
             filter(Response == "GPP_uM_C_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    GPPModFam <- tw(link = "log")
    
    # model selection
    # C23 on day 2 is a huge outlier ~15,000 when next highest point is ~2,500
    # removing doesn't change the best model
    GPP_ModSel <- fun_ms_gams(GPP_uM_C_gAFDM_h, chan16 %>% 
                                                  filter(!(channel == "C23" & MetDate2 == "2016-08-02")), 
                                                  NutTrt, GPPModFam); GPP_ModSel
    
    # "Best model
    chan16_GPP <- chan16 %>% 
      filter(Response == "GPP_uM_C_gAFDM_h") %>% 
      filter(!(channel == "C23" & MetDate2 == "2016-08-02"))
    
    GPP_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                           s(NutTrt, bs = "ts", k = 5, by = MetDate2) + 
                           s(MetDate2, bs = "re"),
                        family = GPPModFam,
                        data = chan16_GPP)
    
    summary(GPP_MostLikely)
    check(getViz(GPP_MostLikely))  
    print(plot(getViz(GPP_MostLikely), allTerms = T), pages = 1)
    
    # output model
    saveRDS(GPP_MostLikely, "03_Model_RDS/MS_2016_GPP_mostlikely.rds")
    
## NEP ----
    # qq line plot
    ggplot(chan16 %>% 
             filter(Response == "NEP_uM_C_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    NEPModFam <- Gamma(link = "log")
    
    # model selection
    # C23 on day 2 is a huge outlier ~15,000 when next highest point is ~2,500

    NEP_ModSel <- fun_ms_gams(NEP_uM_C_gAFDM_h, chan16 %>% 
                                filter(!(channel == "C23" & MetDate2 == "2016-08-02")), 
                              NutTrt, NEPModFam); NEP_ModSel
    
    
    # "Best model
    chan16_NEP <- chan16 %>% 
      filter(Response == "NEP_uM_C_gAFDM_h") %>% 
      filter(!(channel == "C23" & MetDate2 == "2016-08-02"))
    
    NEP_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                           s(NutTrt, bs = "ts", k = 5, by = MetDate2) + 
                           s(MetDate2, bs = "re"),
                         family = NEPModFam,
                         data = chan16_NEP)
    
    summary(NEP_MostLikely)
    check(getViz(NEP_MostLikely))  
    print(plot(getViz(NEP_MostLikely), allTerms = T), pages = 1)
    
    # output model
    saveRDS(NEP_MostLikely, "03_Model_RDS/MS_2016_NEP_mostlikely.rds")
    
 ## N uptake ----
    # LOST
    
## N fixation ----
    # qq line plot
    ggplot(chan16 %>% 
             filter(Response == "Nfix_uM_N_gAFDM_h"), aes(sample = log(Values))) +
      stat_qq() +
      stat_qq_line()
    
    # model family
    NfixModFam <- Gamma(link = "log")
    
    # model selection
    Nfix_ModSel <- fun_ms_gams(Nfix_uM_N_gAFDM_h, chan16, NutTrt, NfixModFam); Nfix_ModSel
    # seems like limited support for the T model
    
    # "Best model
    chan16_Nfix <- chan16 %>% 
      filter(Response == "Nfix_uM_N_gAFDM_h")
    
    Nfix_MostLikely = gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                            s(MetDate2, bs = "re"),
                         family = NfixModFam,
                         data = chan16_Nfix,
                         method = "REML")
    
    blah <- mgcv::gam(Values ~ s(MeanPre2wksTemp, bs = "ts") + 
                        s(NutTrt, bs = "ts", k = 5) + 
                        s(MetDate2, bs = "re"), 
                      family = NfixModFam,
                      data = chan16_Nfix,
                      method = "REML")

    model.sel(Nfix_MostLikely, blah)
    
    summary(Nfix_MostLikely)
    check(getViz(Nfix_MostLikely))  
    print(plot(getViz(Nfix_MostLikely), allTerms = T), pages = 1)

    # output model
    saveRDS(Nfix_MostLikely, "03_Model_RDS/MS_2016_Nfix_mostlikely.rds")

 ## N Assimilation ----
    # N uptake was lost
    
# Model selection Table ----
    ModelSelTable2016 <- rbind(ER_ModSel %>% 
                                 mutate(Response = "ER") %>% 
                                 select(Response, R2, df, AICc, delta, weight),
                               GPP_ModSel %>% 
                                 mutate(Response = "GPP") %>% 
                                 select(Response, R2, df, AICc, delta, weight),
                               NEP_ModSel %>% 
                                 mutate(Response = "NEP") %>% 
                                 select(Response, R2, df, AICc, delta, weight),
                               Nup_ModSel %>% 
                                 mutate(Response = "Nup") %>% 
                                 select(Response, R2, df, AICc, delta, weight),
                               Nfix_ModSel %>% 
                                 mutate(Response = "Nfix") %>% 
                                 select(Response, R2, df, AICc, delta, weight),
                               NAssim_ModSel %>% 
                                 mutate(Response = "Nassim") %>% 
                                 select(Response, R2, df, AICc, delta, weight))
    ModelSelTable2016$Model <- row.names(ModelSelTable2016)
    ModelSelTable2016_2 <- ModelSelTable2016 %>% 
                            select(Response, Model, df, AICc, delta, weight, R2)
    

    
    # https://www.r-bloggers.com/2021/11/publication-ready-tables-with-flextable-and-own-theme-in-r/
    #   https://rfortherestofus.com/2019/11/how-to-make-beautiful-tables-in-r/

    # https://dmyee.files.wordpress.com/2016/03/table_workshop.pdf    
    tab_df(ModelSelTable2016_2,
           alternate.rows = T,
           digits = 2,
           title = "Table S?. Model selection table for mass-specific fluxes from experiment 2.")
    
    # https://towardsdatascience.com/create-flawless-tables-from-your-dataframe-ready-for-publication-7e3fe2d63a52
    library(gt)
    # https://gt.rstudio.com
    # constants ----
    n = 0
    c_col = c("#1e3048", "#274060", "#2f5375", "#4073a0", "#5088b9")
    c_col_light_blue = c("#edf2fb", "#e2eafc", "#d7e3fc", "#ccdbfd", "#c1d3fe")
    c_container_width = px(800)
    c_table_width = px(650)
    c_rn = 30
    c_save = TRUE
    c_format = "html"
    
    
        ModelSelTable2016_2_T <- ModelSelTable2016_2 %>% 
                              gt(groupname_col = "Response",
                                 rowname_col = "Model") %>% 
                              cols_label(R2 = expression(paste(R^2))) %>% 
          fmt_number(columns = c("AICc", "delta", "weight", "R2"),
                     decimals = 2)
        ModelSelTable2016_2_T
    
# Save image
    # save.image("02b_Script_SavedImages/O1_2016_GAMS_massSpecific_Rdat")
    load("02b_Script_SavedImages/O1_2016_GAMS_massSpecific_Rdat")
    