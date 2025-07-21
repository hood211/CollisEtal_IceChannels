# Model selection tables for areal GAMMs from 2017
# JMH Jan 2023, June 24

# Libraries ----
library(tidyverse)
library(mgcv)
library(MuMIn)
library(gt)
library(glue)

# Saved image from modeling script ----
load("02b_Script_savedImages/08_2017_GAMs_MassSpecific.Rdata")


# Model selection table ----
  ## General stuff & Functions ----
### Model FANCY Names ----
  # these apply to all models
  ModelFancyNames1 <- c( "g1_TxNxdate",
                         "g2_TxPxdate",
                         "g3_TxNPxdate",
                         "g4_TxN",
                         "g5_TxP",
                         "g6_TxNP",
                         "g7_T_date_pN_date",
                         "g8_T_date_pP_date",
                         "g9_T_date_pNP_date",
                         "g10_TpN_date",
                         "g11_TpP_date",
                         "g12_TpNP_date",
                         "g13_T_date_pN",
                         "g14_T_date_pP",
                         "g15_T_date_pNP",
                         "g16_TpN",
                         "g17_TpP",
                         "g18_TpNP",
                         "g19_T_date",
                         "g20_T",
                         "g21_N",
                         "g22_P",
                         "g23_NP",
                         "g24_IntOnly")
  
  # need to do this next...
  
  ModelFancyNames2 <- c("s(Temperature, Nf x Date)",
                        "ts(Temperature, Pf x Date)",
                         "s(Temperature, NPf x Date)",
                        "s(Temperature, Nf)",
                        "s(Temperature, Pf)",
                        "s(Temperature, NPf)",
                        "s(Temperature, Date) + Nf x Date",
                        "s(Temperature, Date) + Pf x Date",
                        "s(Temperature, Date) + NPf x Date",
                        "s(Temperature) + Nf x Date",
                        "s(Temperature) + Pf x Date",
                        "s(Temperature) + NPf x Date",
                        "s(Temperature, Date) + Nf",
                        "s(Temperature, Date) + Pf",
                        "s(Temperature, Date) + NPf",
                        "s(Temperature) + Nf",
                        "s(Temperature) + Pf",
                        "s(Temperature) + NPf",
                        "s(Temperature, Date)",
                        "s(Temperature)",
                        "Nf",
                        "Np",
                        "NP",
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

  ModelSelTable_ER_2017 <- ModelSelTableFun(ER_ModSel, "ER (µM C/g AFDM/h)", ModelFancyNames)
  ModelSelTable_GPP_2017 <- ModelSelTableFun(GPP_ModSel, "GPP (µM C/g AFDM/h)", ModelFancyNames)
  ModelSelTable_NEP_2017 <- ModelSelTableFun(NEP_ModSel, "NEP (µM C/g AFDM/h)", ModelFancyNames)
  ModelSelTable_Nfix_2017 <- ModelSelTableFun(Nfix_ModSel, "N fixation (µM N/g AFDM/h)", ModelFancyNames)
  ModelSelTable_Nup_2017 <- ModelSelTableFun(Nup_ModSel, "NH4 uptake (µM N/g AFDM/h)", ModelFancyNames)
  ModelSelTable_TotAssim_2017 <- ModelSelTableFun(Nass_ModSel, "Tot. N Assim (µM N/g AFDM/h)", ModelFancyNames)

  # combine
  ModelSelTable_2017 <- rbind(ModelSelTable_ER_2017,
                              ModelSelTable_GPP_2017,
                              ModelSelTable_NEP_2017,
                              ModelSelTable_Nfix_2017,
                              ModelSelTable_Nup_2017,
                              ModelSelTable_TotAssim_2017) 
  
  
  ModelSelTable_2017_2 <- ModelSelTable_2017 %>% 
    group_by(Response) %>% 
    gt() %>% 
    cols_align(align = "center",
               columns = c("df", "logLik", "AICc", "Delta AIC", "Weight", "R2")) %>% 
    gtsave("04_Tables4MS/ModelSelTable_2017_MassSpecific.docx")
  
  