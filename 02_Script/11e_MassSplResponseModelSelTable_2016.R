# Model selection tables for mass specific GAMMs from 2016
# JMH Jan 2023

# Libraries ----
library(tidyverse)
library(mgcv)
library(MuMIn)
library(gt)
library(glue)

# Saved image from modeling script ----
load("02b_Script_savedImages/01_2016_GAMS_massSpecific_Rdat")


# Model selection table ----
  ## General stuff & Functions ----
### Model FANCY Names ----
  # these apply to all models
  ModelFancyNames1 <- c("g1_ResSurf_ByDate",
                        "g2_ResSurf",
                        "g3_Additive_BothByDate",
                        "g4_Additive_TbyDate",
                        "g5_Additive_NbyDate",
                        "g6_Additive_NoDate",
                        "g7_TempOnly_ByDate",
                        "g8_TempOnly_NoDate",
                        "g9_NitrogenOnly_byDate",
                        "g10_NitrogenOnly_NoDate",
                        "g11_InterceptOnly")
  
  ModelFancyNames2 <- c("te(Temperature, N, Date)",
                        "te(Temperature, N)",
                         "s(Temperature, Date) + s(N, Date)",
                         "s(Temperature, Date) + s(N)",
                        "s(Temperature) + s(N, Date)",
                        "s(Temperature) + s(N)",
                        "s(Temperature, Date)",
                        "s(Temperature)",
                        "s(N, Date)",
                        "s(N)",
                        "Intercept-only")
  
  ModelFancyNames <- as.data.frame(cbind(ModelFancyNames1, ModelFancyNames2))
  
### Function for generating nice table ----
  ModelSelTableFun <- function(ModelSelTab, ResponseName, ModelFancyNames){
    # for coding
      # ModelSelTab <- ER_ModSel
      # ResponseName <- "ER"
    ModelSelTab$Models <- row.names(ModelSelTab)
    row.names(ModelSelTab) <- seq(1,dim(ModelSelTab)[1])
    ModelSelTab2 <- ModelSelTab %>% 
      mutate(ResponseName = ResponseName) %>% 
      select(ResponseName, Models, df:weight, R2) %>% 
      mutate(R2 = round(R2,2),
             logLik = round(logLik, 1),
             AICc = round(AICc, 1),
             delta = round(delta, 2),
             weight = round(weight, 2)) %>% 
      full_join(ModelFancyNames, by = c("Models" = "ModelFancyNames1")) %>% 
      select(Response = ResponseName, Model = ModelFancyNames2,
             df, logLik, AICc, 
             "Delta AIC" = delta, Weight = weight, R2)
  }
  
  
## Model Selection tables ---
  

  ModelSelTable_ER_2016 <- ModelSelTableFun(ER_ModSel, "ER (µM C/g AFDM/h)", ModelFancyNames)
  ModelSelTable_GPP_2016 <- ModelSelTableFun(GPP_ModSel, "GPP (µM C/g AFDM/h)", ModelFancyNames)
  ModelSelTable_NEP_2016 <- ModelSelTableFun(NEP_ModSel, "NEP (µM C/g AFDM/h)", ModelFancyNames)
  ModelSelTable_Nfix_2016 <- ModelSelTableFun(Nfix_ModSel, "N fixation (µM N/g AFDM/h)", ModelFancyNames)
  # ModelSelTable_Nup_2016 <- ModelSelTableFun(Nup_ModSel, "NH4 uptake (µM N/g AFDM/h)", ModelFancyNames)
  # ModelSelTable_TotAssim_2016 <- ModelSelTableFun(NAssim_ModSel, "Tot. N Assim (µM N/g AFDM/h)", ModelFancyNames)

  # combine
  ModelSelTable_2016 <- rbind(ModelSelTable_ER_2016,
                              ModelSelTable_GPP_2016,
                              ModelSelTable_NEP_2016,
                              ModelSelTable_Nfix_2016) 
  
  
  ModelSelTable_2016_2 <- ModelSelTable_2016 %>% 
    group_by(Response) %>% 
    gt() %>% 
    cols_align(align = "center",
               columns = c("df", "logLik", "AICc", "Delta AIC", "Weight", "R2")) %>% 
    gtsave("04_Tables4MS/ModelSelTable_2016_MassSpecific.docx")
  
  