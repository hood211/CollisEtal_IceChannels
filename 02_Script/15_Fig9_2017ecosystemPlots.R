# Generate Fig 9
# Areal and ms responses from 2017
# JMH Jan 2023

# LOAD LIBRARIES ----
library(tidyverse)
library(mgcv)
library(tidymv)
library(RColorBrewer)
library(grid)
library(egg)


# Most likely models ----
## 2017 - mass specific ----
MS_2017_ER_MostLikely <- readRDS("03_Model_RDS/MassSpec_2017_ER_mostlikely.rds")
MS_2017_GPP_MostLikely <- readRDS("03_Model_RDS/MassSpec_2017_GPP_mostlikely.rds")
MS_2017_NEP_MostLikely <-   readRDS("03_Model_RDS/MassSpec_2017_NEP_mostlikely.rds")
MS_2017_Nup_MostLikely <- readRDS("03_Model_RDS/MassSpec_2017_Nup_mostlikely.rds")
MS_2017_Nfix_MostLikely <- readRDS("03_Model_RDS/MassSpec_2017_Nfix_mostlikely.rds")
MS_2017_NAssim_MostLikely <- readRDS("03_Model_RDS/MassSpec_2017_Nass_mostlikely.rds")

## 2017 - areal ----
Areal_2017_ER_MostLikely <- readRDS("03_Model_RDS/Areal_2017_ER_mostlikely.rds")
Areal_2017_GPP_MostLikely <- readRDS("03_Model_RDS/Areal_2017_GPP_mostlikely.rds")
Areal_2017_NEP_MostLikely <- readRDS("03_Model_RDS/Areal_2017_NEP_mostlikely.rds")
Areal_2017_Nup_MostLikely <- readRDS("03_Model_RDS/Areal_2017_Nup_mostlikely.rds")
Areal_2017_Nfix_MostLikely <- readRDS("03_Model_RDS/Areal_2017_Nfix_mostlikely.rds")
Areal_2017_TotAssim_MostLikely <- readRDS("03_Model_RDS/Areal_2017_Nass_mostlikely.rds")
Areal_2017_AFDM_MostLikely <- readRDS("03_Model_RDS/Areal_2017_AFDM_mostlikely.rds")


# standardized ranges ----
  TempRange2017 <- seq(7.9, 21, by = 0.1)
  DINtrt2017 <- c("LowN_T1", "LowN_T2", "HighN_T1", "HighN_T2")
  Ptreats <- c(0, 0.8, 1.6, 3.6, 4.8, 6.5)
  Dates2017 <- c("2017-07-11", "2017-07-22")
  N_uM_date <- c("0.11_2017-07-11", "0.11_2017-07-22", "3.68_2017-07-11", "3.68_2017-07-22")
  P_uM_date <- c("0.36_2017-07-11", "0.36_2017-07-22", "3.94_2017-07-11", "3.94_2017-07-22")
  N_uM <- c("0.11", "3.68")
  P_uM <- c("0.36", "3.94")
  NP <- c("0.31", "0.93", "10.22")

  
# MAKE PREDICTIONS ----
  ## Biomass ----
  
  ## ER -----
    p2017.A.ER <- predict_gam(Areal_2017_ER_MostLikely, values = list(N_uM_date = N_uM_date,
                                                     MeanPre2wksTemp = TempRange2017,
                                                     MetDate2 = Dates2017)) %>% 
                        # separate N_uM_date
                        separate(N_uM_date, into = c("N_uM", "MetDate2b"), sep = "_") %>% 
                        # remove cases where N_uM_date DATE and MetDate2 do not match
                        filter(MetDate2b == MetDate2) %>% 
                        select(-MetDate2b) %>% 
                        mutate(fit.bt = exp(fit)/1000,
                               se.fitL.bt = exp(fit-se.fit)/1000,
                               se.fitU.bt = exp(fit+se.fit)/1000,
                               Basis = "Areal",
                               Res = "ER")

  # This uses P
   p2017.MS.ER <- predict_gam(MS_2017_ER_MostLikely, values = list(P_uM = P_uM,
                                                                    MeanPre2wksTemp = TempRange2017,
                                                                    MetDate2 = Dates2017)) %>%
                        mutate(fit.bt = exp(fit)/1000,
                               se.fitL.bt = exp(fit-se.fit)/1000,
                               se.fitU.bt = exp(fit+se.fit)/1000,
                               Basis = "MS",
                               Res = "ER",
                               N_uM = "0.1-3.7 uM-N",
                               NPratio = "0.3-10.2")%>% 
                       # reorder
                       select("MeanPre2wksTemp", "N_uM", "P_uM", "NPratio", "MetDate2", 
                              "fit", "se.fit", 'fit.bt', "se.fitL.bt", "se.fitU.bt", "Basis", "Res")
     
    
  ## GPP ----
   p2017.A.GPP <- predict_gam(Areal_2017_GPP_MostLikely, values = list(N_uM_date = N_uM_date,
                                                                     MeanPre2wksTemp = TempRange2017,
                                                                     MetDate2 = Dates2017)) %>% 
                         # separate N_uM_date
                         separate(N_uM_date, into = c("N_uM", "MetDate2b"), sep = "_") %>% 
                         # remove cases where N_uM_date DATE and MetDate2 do not match
                         filter(MetDate2b == MetDate2) %>% 
                         select(-MetDate2b) %>% 
                         mutate(fit.bt = exp(fit)/1000,
                                se.fitL.bt = exp(fit-se.fit)/1000,
                                se.fitU.bt = exp(fit+se.fit)/1000,
                                Basis = "Areal",
                                Res = "GPP")
   
  # This uses P-DATE
   p2017.MS.GPP <- predict_gam(MS_2017_GPP_MostLikely, values = list(P_uM_date = P_uM_date,
                                                                     MeanPre2wksTemp = TempRange2017,
                                                                     MetDate2 = Dates2017)) %>%
                       # separate N_uM_date
                       separate(P_uM_date, into = c("P_uM", "MetDate2b"), sep = "_") %>% 
                       # remove cases where P_uM_date DATE and MetDate2 do not match
                       filter(MetDate2b == MetDate2) %>% 
                       select(-MetDate2b) %>% 
                       mutate(fit.bt = exp(fit)/1000,
                              se.fitL.bt = exp(fit-se.fit)/1000,
                              se.fitU.bt = exp(fit+se.fit)/1000,
                              Basis = "MS",
                              Res = "GPP",
                              N_uM = "0.1-3.7 uM-N",
                              NPratio = "0.3-10.2") %>% 
                      # reorder
                      select("MeanPre2wksTemp", "N_uM", "P_uM", "NPratio", "MetDate2", 
                             "fit", "se.fit", 'fit.bt', "se.fitL.bt", "se.fitU.bt", "Basis", "Res")

  
  ## NEP ----
   p2017.A.NEP <- predict_gam(Areal_2017_NEP_MostLikely, values = list(N_uM_date = N_uM_date,
                                                                       MeanPre2wksTemp = TempRange2017,
                                                                       MetDate2 = Dates2017)) %>% 
                    # separate N_uM_date
                    separate(N_uM_date, into = c("N_uM", "MetDate2b"), sep = "_") %>% 
                    # remove cases where N_uM_date DATE and MetDate2 do not match
                    filter(MetDate2b == MetDate2) %>% 
                    select(-MetDate2b) %>% 
                    mutate(fit.bt = exp(fit)/1000,
                           se.fitL.bt = exp(fit-se.fit)/1000,
                           se.fitU.bt = exp(fit+se.fit)/1000,
                           Basis = "Areal",
                           Res = "NEP")
   
  # This uses P
   p2017.MS.NEP <- predict_gam(MS_2017_NEP_MostLikely, values = list(P_uM_date = P_uM_date,
                                                                       MeanPre2wksTemp = TempRange2017,
                                                                       MetDate2 = Dates2017)) %>%
                 # separate N_uM_date
                 separate(P_uM_date, into = c("P_uM", "MetDate2b"), sep = "_") %>% 
                 # remove cases where P_uM_date DATE and MetDate2 do not match
                 filter(MetDate2b == MetDate2) %>% 
                 select(-MetDate2b) %>% 
                 mutate(fit.bt = exp(fit)/1000,
                        se.fitL.bt = exp(fit-se.fit)/1000,
                        se.fitU.bt = exp(fit+se.fit)/1000,
                        Basis = "MS",
                        Res = "NEP",
                        N_uM = "0.1-3.7 uM-N",
                        NPratio = "0.3-10.2") %>% 
                       # reorder
                       select("MeanPre2wksTemp", "N_uM", "P_uM", "NPratio", "MetDate2", 
                              "fit", "se.fit", 'fit.bt', "se.fitL.bt", "se.fitU.bt", "Basis", "Res")

  ## N uptake ----
   p2017.A.Nup <- predict_gam(Areal_2017_Nup_MostLikely, values = list(N_uM = N_uM,
                                                                       MeanPre2wksTemp = TempRange2017,
                                                                       MetDate2 = Dates2017)) %>% 
                     mutate(fit.bt = exp(fit)/1000,
                            se.fitL.bt = exp(fit-se.fit)/1000,
                            se.fitU.bt = exp(fit+se.fit)/1000,
                            Basis = "Areal",
                            Res = "Nup")
   
  
   p2017.MS.Nup <- predict_gam(MS_2017_Nup_MostLikely, values = list(N_uM = N_uM,
                                                                     MeanPre2wksTemp = TempRange2017,
                                                                     MetDate2 = Dates2017)) %>%
                   mutate(fit.bt = exp(fit)/1000,
                          se.fitL.bt = exp(fit-se.fit)/1000,
                          se.fitU.bt = exp(fit+se.fit)/1000,
                          Basis = "MS",
                          Res = "Nup",
                          # N_uM = "0.1-3.7 uM-N",
                          P_uM = "0.4-3.9",
                          NPratio = "0.3-10.2") %>% 
                   # reorder
                   select("MeanPre2wksTemp", "N_uM", "P_uM", "NPratio", "MetDate2", 
                          "fit", "se.fit", 'fit.bt', "se.fitL.bt", "se.fitU.bt", "Basis", "Res")

  ## N fixation ----
   p2017.A.Nfix <- predict_gam(Areal_2017_Nfix_MostLikely, values = list(N_uM = N_uM,
                                                                       MeanPre2wksTemp = TempRange2017,
                                                                       MetDate2 = Dates2017)) %>% 
                   mutate(fit.bt = exp(fit)/1000,
                          se.fitL.bt = exp(fit-se.fit)/1000,
                          se.fitU.bt = exp(fit+se.fit)/1000,
                          Basis = "Areal",
                          Res = "Nfix")
   
  # uses N-date
   p2017.MS.Nfix <- predict_gam(MS_2017_Nfix_MostLikely, values = list(N_uM_date = N_uM_date,
                                                                     MeanPre2wksTemp = TempRange2017,
                                                                     MetDate2 = Dates2017)) %>%
                   # separate N_uM_date
                   separate(N_uM_date, into = c("N_uM", "MetDate2b"), sep = "_") %>% 
                   # remove cases where N_uM_date DATE and MetDate2 do not match
                   filter(MetDate2b == MetDate2) %>% 
                   select(-MetDate2b) %>% 
                   mutate(fit.bt = exp(fit)/1000,
                          se.fitL.bt = exp(fit-se.fit)/1000,
                          se.fitU.bt = exp(fit+se.fit)/1000,
                          Basis = "MS",
                          Res = "Nfix",
                          P_uM = "0.4-3.9",
                          NPratio = "0.3-10.2") %>% 
                 # reorder
                 select("MeanPre2wksTemp", "N_uM", "P_uM", "NPratio", "MetDate2", 
                        "fit", "se.fit", 'fit.bt', "se.fitL.bt", "se.fitU.bt", "Basis", "Res")
   
   
  ## Total N assimilation
   p2017.A.Nass <- predict_gam(Areal_2017_TotAssim_MostLikely, values = list(N_uM = N_uM,
                                                                         MeanPre2wksTemp = TempRange2017,
                                                                         MetDate2 = Dates2017)) %>% 
                   mutate(fit.bt = exp(fit)/1000,
                          se.fitL.bt = exp(fit-se.fit)/1000,
                          se.fitU.bt = exp(fit+se.fit)/1000,
                          Basis = "Areal",
                          Res = "Nassim")
   
   # NP ratio
   p2017.MS.Nass <- predict_gam(MS_2017_NAssim_MostLikely, values = list(N_uM = N_uM,
                                                                       MeanPre2wksTemp = TempRange2017,
                                                                       MetDate2 = Dates2017)) %>%
                 mutate(fit.bt = exp(fit)/1000,
                        se.fitL.bt = exp(fit-se.fit)/1000,
                        se.fitU.bt = exp(fit+se.fit)/1000,
                        Basis = "MS",
                         Res = "Nassim",
                         # N_uM = "0.1-3.7 uM-N",
                         P_uM = "0.4-3.9",
                        NPratio = "0.3-10.2")  %>% 
               # reorder
               select("MeanPre2wksTemp", "N_uM", "P_uM", "NPratio", "MetDate2", 
                      "fit", "se.fit", 'fit.bt', "se.fitL.bt", "se.fitU.bt", "Basis", "Res")
   
# Combine predictions ----
    pComb_Areal <- rbind(p2017.A.ER, 
                         p2017.A.GPP, 
                         p2017.A.NEP, 
                         p2017.A.Nup, 
                         p2017.A.Nfix, 
                         p2017.A.Nass) %>% 
                 mutate(across(c(Res, N_uM, Basis), factor)) %>% 
                 mutate(MetDate2 = ifelse((MetDate2 == "2017-07-11"), "Day1",
                                          ifelse((MetDate2 == "2017-07-22"), "Day2", as.character(MetDate2))),
                        MetDate2 = as.factor(MetDate2)) %>% 
                 mutate(N_uM2 = as.character(N_uM),
                        N_uM2 = fct_recode(N_uM2,
                                           "0.1 µM-N" = "0.11",
                                           "3.7 µM-N" = "3.68"),
                        N_uM2 = fct_relevel(N_uM2, "0.1 µM-N","3.7 µM-N"),
                        Res = fct_relevel(Res, "ER", "GPP", "NEP",
                                          "Nassim", "Nup", "Nfix"))

    pComb_MS <- rbind(p2017.MS.ER,
                      p2017.MS.GPP,
                      p2017.MS.NEP,
                      p2017.MS.Nup,
                      p2017.MS.Nfix,
                      p2017.MS.Nass) %>% 
              mutate(across(c(N_uM, P_uM, NPratio, Basis, Res), factor)) %>% 
              mutate(MetDate2 = ifelse((MetDate2 == "2017-07-11"), "Day1",
                                       ifelse((MetDate2 == "2017-07-22"), "Day2", as.character(MetDate2))),
                     MetDate2 = as.factor(MetDate2),
                     N_uM2 = fct_recode(N_uM, "0 µM-N" = "0.11", "3.6 µM-N" = "3.68"),
                     P_uM2 = fct_recode(P_uM, "0 µM-P" = "0.36", "3.6 µM-P" = "3.94"),
                     # NPratio2 = fct_recode(NPratio, "N:P = 0.3" = "0.31", "N:P = 0.9" = "0.93", "N:P = 10" = "10.22"),
                     Trt = as.factor(ifelse(Res %in% c("GPP", "ER", "NEP"), as.character(P_uM2),
                                            ifelse(Res %in% c("Nup", "Nassim"), as.character(N_uM2),
                                                   ifelse(Res == "Nfix", as.character(N_uM2), "blah")))),
                     Trt = fct_relevel(Trt, "0 µM-N", "3.6 µM-N",
                                            "0 µM-P", "3.6 µM-P"),
                     Res = fct_relevel(Res, "ER", "GPP", "NEP",
                                       "Nassim", "Nup", "Nfix")) 




# Plots ----
    ColorVals.A <- c("black", "#339900")
    linetype.MS <- c("solid", "solid", "dashed", "dashed", "solid", "solid", "solid")
   ColorVals.MS <- c("black", "#339900", 
                      "black", "#339900",
                      "black", "#9900FF", "#FF00CC")

                                                   
    p1 <- ggplot(pComb_Areal %>% 
                   mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                      ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                                                 color = as.factor(N_uM2),
                                                 fill = as.factor(N_uM2))) +
      geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
      geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
      geom_ribbon(color = "transparent", alpha = 0.1) +
      geom_line(size = 1.25) +
      facet_grid(Res ~ MD, scales = "free_y") +
      ylab("") +
      xlab("") +
      scale_color_manual(values = ColorVals.A, name = "Treatment")+
      scale_fill_manual(values = ColorVals.A, name = "Treatment")+
      scale_linetype_manual(values = c("dashed", "solid")) +
      scale_x_continuous(limits = c(5,25), breaks = c(5, 10, 15, 20, 25)) +
      theme(panel.background = element_rect(fill = "white", color = "white"),
            panel.border = element_rect(color = "black", fill = "NA", size = 1),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 60),
            axis.text = element_text(size = 18),
            axis.line = element_line(color = "black", size = 1),
            plot.background = element_rect(fill = "white", color =  "white"),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 24, face = "bold"),
            strip.text.y = element_blank(),
            legend.title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 20),
            legend.key.height = unit(1,"cm")) +
      guides(linetype = "none",
             fill = "none",
             color = "none") 
    
    p2 <- ggplot(pComb_MS %>% 
                   mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                      ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), 
                 aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                               color = Trt,
                               fill = Trt,
                               linetype = Trt)) +
      geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
      geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
      geom_ribbon(color = "transparent", alpha = 0.1) +
      geom_line(size = 1.25) +
      facet_grid(Res ~ MD, scales = "free_y") +
      ylab("") +
      xlab("") +
      scale_color_manual(values = ColorVals.MS, name = "Treatment")+
      scale_fill_manual(values = ColorVals.MS, name = "Treatment")+
      scale_linetype_manual(values = linetype.MS, name = "Treatment")+
      scale_x_continuous(limits = c(5,25), breaks = c(5, 10, 15, 20, 25)) +
      theme(panel.background = element_rect(fill = "white", color = "white"),
            panel.border = element_rect(color = "black", fill = "NA", size = 1),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 60),
            axis.text = element_text(size = 18),
            axis.line = element_line(color = "black", size = 1),
            plot.background = element_rect(fill = "white", color =  "white"),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 24, face = "bold"),
            strip.text.y = element_blank(),
            legend.key.width = unit(1.75, "cm"),
            legend.title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 20),
            legend.key.height = unit(1,"cm")) 
    
    # Print ----
    p1.g <- ggplotGrob(p1)
    p2.g <- ggplotGrob(p2)

    p1.gtf <- gtable_frame(p1.g, width = unit(1, "null"), height = unit(1, "null"))
    p2.gtf <- gtable_frame(p2.g, width = unit(1, "null"), height = unit(1, "null"))

    p12.gtf <- gtable_frame(gtable_cbind(p1.gtf, p2.gtf), width = unit(2,"null"), height = unit(1,"null"))
    
    png("05_Figures4MS/15_FigS13_2017_ArealandMS.png", units = "in", height = 18, width = 18, res = 300)
    grid.newpage()
    grid.draw(p12.gtf)
    
    # Areal labels
    grid.text("Areal ER",
              x = unit(0.02,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("Areal GPP",
              x = unit(0.02,"npc"), y = unit(0.75,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.75,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("Areal NEP",
              x = unit(0.02,"npc"), y = unit(0.59,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.59,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("Areal Total N assim.",
              x = unit(0.02,"npc"), y = unit(0.435,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.435,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("Areal N uptake",
              x = unit(0.02,"npc"), y = unit(0.285,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.285,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(bold(paste("Areal ",N[2]," fixation"))),
              x = unit(0.02,"npc"), y = unit(0.125,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.125,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    

    
    
    
    # Mass-specific labels
    grid.text("MS ER",
              x = unit(0.44,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.46,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("MS GPP",
              x = unit(0.44,"npc"), y = unit(0.75,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.46,"npc"), y = unit(0.75,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("MS NEP",
              x = unit(0.44,"npc"), y = unit(0.59,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.46,"npc"), y = unit(0.59,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("MS Total N assim.",
              x = unit(0.44,"npc"), y = unit(0.435,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.46,"npc"), y = unit(0.435,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text("MS N uptake",
              x = unit(0.44,"npc"), y = unit(0.285,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.46,"npc"), y = unit(0.285,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(bold(paste("MS ",N[2]," fixation"))), 
              x = unit(0.44,"npc"), y = unit(0.125,"npc"), 
              gp=gpar(fontsize = 22, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.46,"npc"), y = unit(0.125,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    

    
    grid.text("Temperature (°C)", 
              x = unit(0.5,"npc"), y = unit(0.015,"npc"), 
              gp=gpar(fontsize = 28, fontface = "bold"))
    dev.off()
    
    
  # save.image
    # save.image("02b_Script_SavedImages/15_Fig4_BM_Met_2015_Rdat")
    # load("02b_Script_SavedImages/15_Fig4_BM_Met_2015_Rdat")
    