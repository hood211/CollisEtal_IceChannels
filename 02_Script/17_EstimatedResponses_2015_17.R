# responses over standarized intervals
# 11 and 21C
# 0 and 3.6 N or P
# day 1 and 2
# JMH Jun 2023
#JMH Aug 2023, June 24


# Libraries ----
library(tidyverse)
library(mgcv)
# remotes::install_github("stefanocoretta/tidygam@v0.1.0")
library(tidygam)
library(gt)
library(scales)
library(webshot2)


# this gives an error, but all works out :)
ResponseTabFun <- function(ModelRDS, Day1, Day2, ModelName){
  # ModelRDS = ER_2015
  # Day1 = "2015-07-15"
  # Day2 = "2015-07-27"
  ModelPredictions <- predict_gam(ModelRDS, values = list(MetDate2 = c(Day1, Day2),
                                                   NutTrt = c(0,3.6),
                                                   MeanPre2wksTemp = c(11,21)),
                                  tran_fun = exp)%>% 
    mutate(Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
    pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
    mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
           Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
           N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
           N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_21_0*100,
           Response = ModelName)
  
}

# Generate response tables ----
## 2015 ----
### 2015 Areal ----
# most likely models
ER_2015 <- readRDS("03_Model_RDS/Areal_2015_ER_mostlikely.rds")
GPP_2015 <- readRDS("03_Model_RDS/Areal_2015_GPP_mostlikely.rds")
NEP_2015 <- readRDS("03_Model_RDS/Areal_2015_NEP_mostlikely.rds")
Nup_2015 <- readRDS("03_Model_RDS/Areal_2015_Nup_mostlikely.rds")
Nfix_2015 <- readRDS("03_Model_RDS/Areal_2015_Nfix_mostlikely.rds")
Nassim_2015 <- readRDS("03_Model_RDS/Areal_2015_NTotAssim_mostlikely.rds")
AFDM_2015 <- readRDS("03_Model_RDS/Areal_2015_AFDM_mostlikely.rds")

# response tables
# generates error, but not a problem
RT_ER_2015_areal <- ResponseTabFun(ER_2015, "2015-07-15", "2015-07-27", "2015_ER_Areal") %>% mutate(Nutrient = "N")
RT_GPP_2015_areal <- ResponseTabFun(GPP_2015, "2015-07-15", "2015-07-27", "2015_GPP_Areal") %>% mutate(Nutrient = "N")
RT_NEP_2015_areal <- ResponseTabFun(NEP_2015, "2015-07-15", "2015-07-27", "2015_NEP_Areal") %>% mutate(Nutrient = "N")
RT_Nup_2015_areal <- ResponseTabFun(Nup_2015, "2015-07-15", "2015-07-27", "2015_Nup_Areal") %>% mutate(Nutrient = "N")
RT_Nfix_2015_areal <- ResponseTabFun(Nfix_2015, "2015-07-15", "2015-07-27", "2015_Nfix_Areal") %>% mutate(Nutrient = "N") 
RT_Nassim_2015_areal <- ResponseTabFun(Nassim_2015, "2015-07-15", "2015-07-27", "2015_Nassim_Areal") %>% mutate(Nutrient = "N")
RT_AFDM_2015_areal <- ResponseTabFun(AFDM_2015, "2015-07-15", "2015-07-27", "2015_AFDM_Areal") %>% mutate(Nutrient = "N")

### 2015 MS ----
# most likely models
ER_2015_MS <- readRDS("03_Model_RDS/MS_2015_ER_mostlikely.rds")
GPP_2015_MS <- readRDS("03_Model_RDS/MS_2015_GPP_mostlikely.rds")
NEP_2015_MS <- readRDS("03_Model_RDS/MS_2015_NEP_mostlikely.rds")
Nup_2015_MS <- readRDS("03_Model_RDS/MS_2015_Nup_mostlikely.rds")
Nfix_2015_MS <- readRDS("03_Model_RDS/MS_2015_Nfix_mostlikely.rds")
Nassim_2015_MS <- readRDS("03_Model_RDS/MS_2015_NAssim_mostlikely.rds")

# response tables
# No temp, can't use function
# RT_ER_2015_MS <- ResponseTabFun(ER_2015_MS, "2015-07-15", "2015-07-27", "2015_ER_MS") %>% mutate(Nutrient = "N")
RT_ER_2015_MS0 <- predict_gam(ER_2015_MS, values = list(MetDate2 = c("2015-07-15", "2015-07-27"),
                                                          NutTrt = c(0,3.6)),
                               tran_fun = exp) 

RT_ER_2015_MS <- rbind(RT_ER_2015_MS0 %>% 
                          mutate(MeanPre2wksTemp = "11") ,
                       RT_ER_2015_MS0 %>% 
                          mutate(MeanPre2wksTemp = "21")) %>% 
  mutate(Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
  pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
  mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
         Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
         N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
         N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
         Response = "2015_ER_MS",
         Nutrient = "N")

# No temp, can't use function
RT_GPP_2015_MS0 <- predict_gam(GPP_2015_MS, values = list(MetDate2 = c("2015-07-15", "2015-07-27"),
                                                      NutTrt = c(0,3.6)),
                               tran_fun = exp) 

RT_GPP_2015_MS <- rbind(RT_GPP_2015_MS0 %>% 
                         mutate(MeanPre2wksTemp = "11") ,
                        RT_GPP_2015_MS0 %>% 
                          mutate(MeanPre2wksTemp = "21")) %>% 
                      mutate(Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
                      pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                      mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                             Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                             N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                             N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                             Response = "2015_GPP_MS",
                             Nutrient = "N")
  
  
# No temp, can't use function
RT_NEP_2015_MS0 <- predict_gam(NEP_2015_MS, values = list(MetDate2 = c("2015-07-15", "2015-07-27"),
                                                          NutTrt = c(0,3.6)),
                               tran_fun = exp) 
  
RT_NEP_2015_MS <- rbind(RT_NEP_2015_MS0 %>% 
                          mutate(MeanPre2wksTemp = "11") ,
                        RT_NEP_2015_MS0 %>% 
                          mutate(MeanPre2wksTemp = "21")) %>% 
                    mutate(Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
                    pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                    mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                           Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                           N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                           N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                           Response = "2015_NEP_MS",
                           Nutrient = "N")
  
  
RT_Nup_2015_MS <- ResponseTabFun(Nup_2015_MS, "2015-07-15", "2015-07-27", "2015_Nup_MS") %>% mutate(Nutrient = "N")
RT_Nfix_2015_MS <- ResponseTabFun(Nfix_2015_MS, "2015-07-15", "2015-07-27", "2015_Nfix_MS") %>% mutate(Nutrient = "N")
RT_Nassim_2015_MS <- ResponseTabFun(Nassim_2015_MS, "2015-07-15", "2015-07-27", "2015_Nassim_MS") %>% mutate(Nutrient = "N")

## 2016 ----
### 2016 Areal ----
# most likely models
ER_2016 <- readRDS("03_Model_RDS/Areal_2016_ER_mostlikely.rds")
GPP_2016 <- readRDS("03_Model_RDS/Areal_2016_GPP_mostlikely.rds")
NEP_2016 <- readRDS("03_Model_RDS/Areal_2016_NEP_mostlikely.rds")
Nfix_2016 <- readRDS("03_Model_RDS/Areal_2016_Nfix_mostlikely.rds")
AFDM_2016 <- readRDS("03_Model_RDS/Areal_2016_AFDM_mostlikely.rds")


# response tables
RT_ER_2016_areal <- ResponseTabFun(ER_2016, "2016-07-27", "2016-08-02", "2016_ER_Areal") %>% mutate(Nutrient = "P")
RT_GPP_2016_areal <- ResponseTabFun(GPP_2016, "2016-07-27", "2016-08-02", "2016_GPP_Areal") %>% mutate(Nutrient = "P")
RT_NEP_2016_areal <- ResponseTabFun(NEP_2016, "2016-07-27", "2016-08-02", "2016_NEP_Areal") %>% mutate(Nutrient = "P")

# No N needs to be done manually
RT_Nfix_2016_areal0 <-  predict_gam(Nfix_2016, values = list(MetDate2 = c("2016-07-27", "2016-08-02"),
                                                            MeanPre2wksTemp = c(11,21)),
                                    tran_fun = exp) 

RT_Nfix_2016_areal <- rbind(RT_Nfix_2016_areal0 %>% 
                          mutate(NutTrt = "0") ,
                        RT_Nfix_2016_areal0 %>% 
                          mutate(NutTrt = "3.6")) %>% 
                      mutate(Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
                      pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                      mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                             Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                             N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                             N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                             Response = "2016_Nfix_Areal",
                             Nutrient = "P")

# need to do manually because model predicts neg at O P and 11
RT_AFDM_2016_areal0 <- predict_gam(AFDM_2016, values = list(MetDate2 = c("2016-07-27", "2016-08-02"),
                                                          NutTrt = c(0,3.6),
                                                          MeanPre2wksTemp = c(11,21)),
                                  tran_fun = exp) 

RT_AFDM_2016_areal <- rbind(RT_AFDM_2016_areal0 %>% 
                              mutate(NutTrt = "0"), 
                            RT_AFDM_2016_areal0 %>% 
                              mutate(NutTrt = "3.6")) %>% 
                      mutate(Values = ifelse(Values <0, 0.1, Values),
                                             Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
                      pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                      mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                             Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                             N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                             N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                             Response = "2016_AFDM_Areal",
                             Nutrient = "P")


### 2016 MS ----
# most likely models
ER_2016_MS <- readRDS("03_Model_RDS/MS_2016_ER_mostlikely.rds")
GPP_2016_MS <- readRDS("03_Model_RDS/MS_2016_GPP_mostlikely.rds")
NEP_2016_MS <- readRDS("03_Model_RDS/MS_2016_NEP_mostlikely.rds")
Nfix_2016_MS <- readRDS("03_Model_RDS/MS_2016_Nfix_mostlikely.rds")

# response tables
# No P, needs to be done manually
RT_ER_2016_MS0 <- predict_gam(ER_2016_MS, values = list(MetDate2 = c("2016-07-27", "2016-08-02"),
                                                      MeanPre2wksTemp = c(11,21)),
                              tran_fun = exp) 

RT_ER_2016_MS <- rbind(RT_ER_2016_MS0 %>% 
                              mutate(NutTrt = "0") ,
                            RT_ER_2016_MS0 %>% 
                              mutate(NutTrt = "3.6")) %>% 
                mutate(Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
                pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                       Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                       N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                       N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                       Response = "2016_ER_MS",
                       Nutrient = "P")
  
RT_GPP_2016_MS <- ResponseTabFun(GPP_2016_MS, "2016-07-27", "2016-08-02", "2016_GPP_MS") %>% mutate(Nutrient = "P")
RT_NEP_2016_MS <- ResponseTabFun(NEP_2016_MS, "2016-07-27", "2016-08-02", "2016_NEP_MS") %>% mutate(Nutrient = "P")

# No P, needs to be done manually
RT_Nfix_2016_MS0 <- predict_gam(Nfix_2016_MS, values = list(MetDate2 = c("2016-07-27", "2016-08-02"),
                                                         MeanPre2wksTemp = c(11,21)),
                                tran_fun = exp) 

RT_Nfix_2016_MS <- rbind(RT_Nfix_2016_MS0 %>% 
                         mutate(NutTrt = "0") ,
                       RT_Nfix_2016_MS0 %>% 
                         mutate(NutTrt = "3.6")) %>% 
                  mutate(Temp_N = paste0("TN_",MeanPre2wksTemp,"_", NutTrt)) %>% 
                  pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                  mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                         Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                         N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                         N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                         Response = "2016_Nfix_MS",
                         Nutrient = "P")
  
  
  
## 2017 ----
### 2017 Areal ----
# most likely models
ER_2017 <- readRDS("03_Model_RDS/Areal_2017_ER_mostlikely.rds")
GPP_2017 <- readRDS("03_Model_RDS/Areal_2017_GPP_mostlikely.rds")
NEP_2017 <- readRDS("03_Model_RDS/Areal_2017_NEP_mostlikely.rds")
Nup_2017 <- readRDS("03_Model_RDS/Areal_2017_Nup_mostlikely.rds")
Nfix_2017 <- readRDS("03_Model_RDS/Areal_2017_Nfix_mostlikely.rds")
Nassim_2017 <- readRDS("03_Model_RDS/Areal_2017_Nass_mostlikely.rds")
AFDM_2017 <- readRDS("03_Model_RDS/Areal_2017_AFDM_mostlikely.rds")

# response tables
RT_AFDM_2017_areal <- predict_gam(AFDM_2017, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                       MeanPre2wksTemp = c(11,21)),
                                tran_fun = exp) %>% 
                    mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                           Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
                    pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                    mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                           Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                           N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                           N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                           Response = "2017_AFDM_Areal",
                           Nutrient = "N")

RT_ER_2017_areal <- predict_gam(ER_2017, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                        MeanPre2wksTemp = c(11,21)),
                                tran_fun = exp) %>% 
                      separate(N_uM_date, into = c("N_uM", "MetDateModel"), sep = "_") %>% 
                      # remove predictions where Day differs between NxDay and re Day
                      filter(MetDateModel == MetDate2) %>% 
                      mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                        Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
                      pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                      mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                             Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                             N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                             N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                             Response = "2017_ER_Areal",
                             Nutrient = "N")
  
RT_GPP_2017_areal <- predict_gam(GPP_2017, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                        MeanPre2wksTemp = c(11,21)),
                                 tran_fun = exp) %>% 
                      separate(N_uM_date, into = c("N_uM", "MetDateModel"), sep = "_") %>% 
                      # remove predictions where Day differs between NxDay and re Day
                      filter(MetDateModel == MetDate2) %>% 
                      mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                             Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
                      pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                      mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                             Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                             N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                             N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                             Response = "2017_GPP_Areal",
                             Nutrient = "N")
  
  
RT_NEP_2017_areal <- predict_gam(NEP_2017, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                         MeanPre2wksTemp = c(11,21)),
                                 tran_fun = exp) %>% 
                      separate(N_uM_date, into = c("N_uM", "MetDateModel"), sep = "_") %>% 
                      # remove predictions where Day differs between NxDay and re Day
                      filter(MetDateModel == MetDate2) %>% 
                      mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                             Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
                      pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                      mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                             Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                             N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                             N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                             Response = "2017_NEP_Areal",
                             Nutrient = "N")

RT_Nup_2017_areal <- predict_gam(Nup_2017, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                         MeanPre2wksTemp = c(11,21)),
                                 tran_fun = exp) %>% 
                      mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                             Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
                      pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                      mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                             Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                             N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                             N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                             Response = "2017_Nup_Areal",
                             Nutrient = "N")
  
RT_Nfix_2017_areal <- predict_gam(Nfix_2017, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                          MeanPre2wksTemp = c(11,21)),
                                  tran_fun = exp) %>% 
                    mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                           Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
                    pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                    mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                           Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                           N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                           N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                           Response = "2017_Nfix_Areal",
                           Nutrient = "N")
  
RT_Nassim_2017_areal <- predict_gam(Nassim_2017, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                             MeanPre2wksTemp = c(11,21)),
                                    tran_fun = exp) %>% 
                  mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                         Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
                  pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
                  mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                         Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                         N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                         N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                         Response = "2017_Nassim_Areal",
                         Nutrient = "N")



### 2017 MS ----
# most likely models
ER_2017_MS <- readRDS("03_Model_RDS/MassSpec_2017_ER_mostlikely.rds")
GPP_2017_MS <- readRDS("03_Model_RDS/MassSpec_2017_GPP_mostlikely.rds")
NEP_2017_MS <- readRDS("03_Model_RDS/MassSpec_2017_NEP_mostlikely.rds")
Nup_2017_MS <- readRDS("03_Model_RDS/MassSpec_2017_Nup_mostlikely.rds") # uses NP, now N
Nfix_2017_MS <- readRDS("03_Model_RDS/MassSpec_2017_Nfix_mostlikely.rds")
Nassim_2017_MS <- readRDS("03_Model_RDS/MassSpec_2017_Nass_mostlikely.rds") # uses NP, now N

# response tables

RT_ER_2017_MS <- predict_gam(ER_2017_MS, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                       MeanPre2wksTemp = c(11,21)),
                             tran_fun = exp) %>% 
          # these were coded as the actual P, not added P
          mutate(P_uM = ifelse(P_uM == 0.36, 0, 3.6),
                 Temp_N = paste0("TN_",MeanPre2wksTemp,"_", P_uM)) %>% 
          pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
          mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                 Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                 N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                 N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                 Response = "2017_ER_MS",
                 Nutrient = "P")

RT_GPP_2017_MS <- predict_gam(GPP_2017_MS, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                          MeanPre2wksTemp = c(11,21)),
                              tran_fun = exp) %>% 
  separate(P_uM_date, into = c("P_uM", "MetDateModel"), sep = "_") %>% 
  # remove predictions where Day differs between NxDay and re Day
  filter(MetDateModel == MetDate2) %>% 
  # these were coded as the actual P, not added P
  mutate(P_uM = ifelse(P_uM == 0.36, 0, 3.6),
         Temp_N = paste0("TN_",MeanPre2wksTemp,"_", P_uM)) %>% 
  pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
  mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
         Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
         N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
         N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
         Response = "2017_GPP_MS",
         Nutrient = "P")


RT_NEP_2017_MS <- predict_gam(NEP_2017_MS, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                            MeanPre2wksTemp = c(11,21)),
                              tran_fun = exp) %>% 
  separate(P_uM_date, into = c("P_uM", "MetDateModel"), sep = "_") %>% 
  # remove predictions where Day differs between NxDay and re Day
  filter(MetDateModel == MetDate2) %>% 
  # these were coded as the actual P, not added P
  mutate(P_uM = ifelse(P_uM == 0.36, 0, 3.6),
         Temp_N = paste0("TN_",MeanPre2wksTemp,"_", P_uM)) %>% 
  pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
  mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
         Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
         N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
         N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
         Response = "2017_NEP_MS",
         Nutrient = "P")


RT_Nfix_2017_MS <- predict_gam(Nfix_2017_MS, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                            MeanPre2wksTemp = c(11,21)),
                               tran_fun = exp) %>% 
        separate(N_uM_date, into = c("N_uM", "MetDateModel"), sep = "_") %>% 
        # remove predictions where Day differs between NxDay and re Day
        filter(MetDateModel == MetDate2) %>% 
        # these were coded as the actual P, not added P
        mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
               Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
        pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
        mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
               Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
               N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
               N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
               Response = "2017_Nfix_MS",
               Nutrient = "N")

RT_Nup_2017_MS <- predict_gam(Nup_2017_MS, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                           MeanPre2wksTemp = c(11,21)),
                               tran_fun = exp) %>% 
          mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
                 Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
          pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
          mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
                 Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
                 N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
                 N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
                 Response = "2017_Nup_MS",
                 Nutrient = "N")


RT_Nassim_2017_MS <- predict_gam(Nassim_2017_MS, values = list(MetDate2 = c("2017-07-11", "2017-07-22"),
                                                         MeanPre2wksTemp = c(11,21)),
                              tran_fun = exp) %>% 
        mutate(N_uM = ifelse(N_uM == 0.11, 0, 3.6),
               Temp_N = paste0("TN_",MeanPre2wksTemp,"_", N_uM)) %>% 
        pivot_wider(id_cols = MetDate2, names_from = Temp_N, values_from = Values) %>% 
        mutate(Temp_LowN_PerChange = (TN_21_0-TN_11_0)/TN_11_0*100,
               Temp_HighN_PerChange = (TN_21_3.6-TN_11_3.6)/TN_11_3.6*100,
               N_Cold_PerChange = (TN_11_3.6- TN_11_0)/TN_11_0*100,
               N_Warm_PerChange = (TN_21_3.6- TN_21_0)/TN_11_0*100,
               Response = "2017_Nassim_MS",
               Nutrient = "N")


# Make tables ----

# get list of response tables
objects(pattern = "RT")

ResponseTable <- rbind(
                        #2015
                        RT_AFDM_2015_areal, RT_ER_2015_areal, RT_GPP_2015_areal, RT_NEP_2015_areal,  RT_Nassim_2015_areal,RT_Nfix_2015_areal,RT_Nup_2015_areal,
                       
                       RT_ER_2015_MS, RT_GPP_2015_MS, RT_NEP_2015_MS, RT_Nassim_2015_MS,RT_Nfix_2015_MS,RT_Nup_2015_MS,
                       
                       # 2016
                       RT_AFDM_2016_areal, RT_ER_2016_areal,RT_GPP_2016_areal,RT_NEP_2016_areal,RT_Nfix_2016_areal,
                       
                       RT_ER_2016_MS, RT_GPP_2016_MS, RT_NEP_2016_MS, RT_Nfix_2016_MS,
                       
                       # 2017
                       RT_AFDM_2017_areal, RT_ER_2017_areal, RT_GPP_2017_areal, RT_NEP_2017_areal, RT_Nassim_2017_areal, RT_Nfix_2017_areal, RT_Nup_2017_areal,
                       
                       # Nup_2017_MS and Nassim_2017_MS is an NP model
                       RT_ER_2017_MS, RT_GPP_2017_MS, RT_NEP_2017_MS, RT_Nfix_2017_MS, RT_Nup_2017_MS, RT_Nassim_2017_MS) %>% 
                separate(Response, into = c("Year", "Response", "ArealOrMs"), sep = "_", remove = FALSE)

# https://stackoverflow.com/questions/65327289/how-to-represent-a-datatable-in-r-as-a-heatmap
# https://gt.rstudio.com/reference/gtsave.html



## 2015 table ----
as.data.frame(ResponseTable) %>% 
  select(Year, Response, ArealOrMs, MetDate2, Temp_LowN_PerChange:N_Warm_PerChange) %>% 
  filter(Year == "2015") %>% 
  mutate(across(Temp_LowN_PerChange:N_Warm_PerChange, \(x) round(x, digits = 0))) %>% 
  mutate(MetDate2 = ifelse(MetDate2 == "2015-07-15", "Day 1","Day 2")) %>% 
  arrange(Year, Response, ArealOrMs, MetDate2) %>% 
  gt() %>% 
  tab_spanner(
    label = "Temp % change",
    columns = c(Temp_LowN_PerChange, Temp_HighN_PerChange)
  ) %>% 
  cols_label(Temp_LowN_PerChange = "Low nutrient",
             Temp_HighN_PerChange = "High nutrient",
             N_Cold_PerChange = "Cold",
             N_Warm_PerChange = "Warm") %>% 
  tab_spanner(
    label = "Nutrient % change",
    columns = c(N_Cold_PerChange, N_Warm_PerChange)
  ) %>% 
  data_color(columns = 5:8, 
             colors = col_numeric(palette = c("red", "white","blue"),
                                  domain = c(-200,200))) %>% 
  gtsave("04_Tables4MS/ResponseTable2015.html")

## 2016 table ----
as.data.frame(ResponseTable) %>% 
  select(Year, Response, ArealOrMs, MetDate2, Temp_LowN_PerChange:N_Warm_PerChange) %>%
  filter(Year == "2016") %>% 
  mutate(across(Temp_LowN_PerChange:N_Warm_PerChange, \(x) round(x, digits = 0))) %>% 
  mutate(MetDate2 = ifelse(MetDate2 == "2016-07-27", "Day 1","Day 2")) %>% 
  arrange(Year, Response, ArealOrMs, MetDate2) %>% 
  gt() %>% 
  tab_spanner(
    label = "Temp % change",
    columns = c(Temp_LowN_PerChange, Temp_HighN_PerChange)
  ) %>% 
  cols_label(Temp_LowN_PerChange = "Low nutrient",
             Temp_HighN_PerChange = "High nutrient",
             N_Cold_PerChange = "Cold",
             N_Warm_PerChange = "Warm") %>% 
  tab_spanner(
    label = "Nutrient % change",
    columns = c(N_Cold_PerChange, N_Warm_PerChange)
  ) %>% 
  data_color(columns = 5:8, 
             colors = col_numeric(palette = c("red", "white","blue"),
                                  domain = c(-200,200))) %>% 
  gtsave("04_Tables4MS/ResponseTable2016.html")


## 2017 table ----
as.data.frame(ResponseTable) %>% 
  select(Year, Response, ArealOrMs, MetDate2, Temp_LowN_PerChange:N_Warm_PerChange) %>%
  filter(Year == "2017") %>% 
  mutate(across(Temp_LowN_PerChange:N_Warm_PerChange, \(x) round(x, digits = 0))) %>% 
  mutate(MetDate2 = ifelse(MetDate2 == "2017-07-11", "Day 1","Day 2")) %>% 
  arrange(Year, Response, ArealOrMs, MetDate2) %>% 
  gt() %>% 
  tab_spanner(
    label = "Temp % change",
    columns = c(Temp_LowN_PerChange, Temp_HighN_PerChange)
  ) %>% 
  cols_label(Temp_LowN_PerChange = "Low nutrient",
             Temp_HighN_PerChange = "High nutrient",
             N_Cold_PerChange = "Cold",
             N_Warm_PerChange = "Warm") %>% 
  tab_spanner(
    label = "Nutrient % change",
    columns = c(N_Cold_PerChange, N_Warm_PerChange)
  ) %>% 
  data_color(columns = 5:8, 
             colors = col_numeric(palette = c("red", "white","blue"),
                                  domain = c(-200,200))) %>% 
  gtsave("04_Tables4MS/ResponseTable2017.html")


## Table 2 ----
Table2 <- as.data.frame(ResponseTable) %>% 
  select(Year, Response, ArealOrMs, MetDate2, Temp_LowN_PerChange:N_Warm_PerChange) %>% 
  mutate(MetDate2 = case_when(MetDate2 == "2015-07-15" ~ "MD1",
                              MetDate2 == "2015-07-27" ~ "MD2",
                              MetDate2 == "2016-07-27" ~ "MD1",
                              MetDate2 == "2016-08-02" ~ "MD2",
                              MetDate2 == "2017-07-11" ~ "MD1",
                              MetDate2 == "2017-07-22" ~ "MD2")) %>% 
  pivot_wider(id_cols = Year:ArealOrMs, names_from = MetDate2, values_from = Temp_LowN_PerChange:N_Warm_PerChange) %>% 
  mutate(across(Temp_LowN_PerChange_MD1:N_Warm_PerChange_MD2, \(x) round(x, digits = 0))) %>% 
  arrange(Year, Response, ArealOrMs) 

write.csv(Table2, "04_Tables4MS/Table2.csv")

save.image("02b_Script_SavedImages/17_EstimatedResponses_2015.Rdata")
# load("02b_Script_SavedImages/17_EstimatedResponses_2015_Rdat")
