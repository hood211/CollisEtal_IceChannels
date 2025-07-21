# Generate Fig 7
# Areal and ms responses from 2016
# JMH Jan 2023

# LOAD LIBRARIES ----
library(tidyverse)
library(mgcv)
library(tidymv)
library(RColorBrewer)
library(grid)
library(egg)

# Most likely models ----
## 2016 - mass specific ----
MS_2016_ER_MostLikely <- readRDS("03_Model_RDS/MS_2016_ER_mostlikely.rds")
MS_2016_GPP_MostLikely <- readRDS("03_Model_RDS/MS_2016_GPP_mostlikely.rds")
MS_2016_NEP_MostLikely <- readRDS("03_Model_RDS/MS_2016_NEP_mostlikely.rds")
# N uptake lost
MS_2016_Nfix_MostLikely <- readRDS("03_Model_RDS/MS_2016_Nfix_mostlikely.rds")

## 2016 - areal ----
Areal_2016_ER_MostLikely <- readRDS("03_Model_RDS/Areal_2016_ER_mostlikely.rds")
Areal_2016_GPP_MostLikely <- readRDS("03_Model_RDS/Areal_2016_GPP_mostlikely.rds")
Areal_2016_NEP_MostLikely <- readRDS("03_Model_RDS/Areal_2016_NEP_mostlikely.rds")
# N uptake lost
Areal_2016_Nfix_MostLikely <- readRDS("03_Model_RDS/Areal_2016_Nfix_mostlikely.rds")



# standardized ranges ----
TempRange2015_16 <- seq(7.9, 25.5, by = 0.1)
DINtreats <- c(0,1.8, 3.6, 7.1, 10.7, 14.3)
Ptreats <- c(0, 0.8, 1.6, 3.6, 4.8, 6.5)
Dates2015 <- c("2015-07-15", "2015-07-27")
Dates2016 <- c("2016-07-27", "2016-08-02")


# MAKE PREDICTIONS ----

## ER -----
p2016.A.ER <- predict_gam(Areal_2016_ER_MostLikely, values = list(NutTrt = Ptreats,
                                                                  MeanPre2wksTemp = TempRange2015_16,
                                                                  MetDate2 = Dates2016)) %>% 
              mutate(fit.bt = exp(fit)/1000,
                     se.fitL.bt = exp(fit-se.fit)/1000,
                     se.fitU.bt = exp(fit+se.fit)/1000,
                     trt = "Ponly",
                     N_uM = 0,
                     Basis = "Areal",
                     Res = "ER")

# no P effect
p2016.MS.ER <- predict_gam(MS_2016_ER_MostLikely, values = list(MeanPre2wksTemp = TempRange2015_16,
                                                               MetDate = Dates2016)) %>% 
  mutate(fit.bt = exp(fit)/1000,
         se.fitL.bt = exp(fit-se.fit)/1000,
         se.fitU.bt = exp(fit+se.fit)/1000,
         trt = "Ponly",
         N_uM = 0,
         P_uM = 0,
         NutTrt = 0,
         Basis = "MS",
         Res = "ER")

## GPP -----
p2016.A.GPP <- predict_gam(Areal_2016_GPP_MostLikely, values = list(NutTrt = Ptreats,
                                                                  MeanPre2wksTemp = TempRange2015_16,
                                                                  MetDate2 = Dates2016)) %>% 
  mutate(fit.bt = exp(fit)/1000,
         se.fitL.bt = exp(fit-se.fit)/1000,
         se.fitU.bt = exp(fit+se.fit)/1000,
         trt = "Ponly",
         N_uM = 0,
         Basis = "Areal",
         Res = "GPP")


p2016.MS.GPP <- predict_gam(MS_2016_GPP_MostLikely, values = list(NutTrt = Ptreats,
                                                                  MeanPre2wksTemp = TempRange2015_16,
                                                                  MetDate = Dates2016)) %>% 
  mutate(fit.bt = exp(fit)/1000,
         se.fitL.bt = exp(fit-se.fit)/1000,
         se.fitU.bt = exp(fit+se.fit)/1000,
         trt = "Ponly",
         N_uM = 0,
         Basis = "MS",
         Res = "GPP")


## NEP -----
p2016.A.NEP <- predict_gam(Areal_2016_NEP_MostLikely, values = list(NutTrt = Ptreats,
                                                                    MeanPre2wksTemp = TempRange2015_16,
                                                                    MetDate2 = Dates2016)) %>% 
  mutate(fit.bt = exp(fit)/1000,
         se.fitL.bt = exp(fit-se.fit)/1000,
         se.fitU.bt = exp(fit+se.fit)/1000,
         trt = "Ponly",
         N_uM = 0,
         Basis = "Areal",
         Res = "NEP")


p2016.MS.NEP <- predict_gam(MS_2016_NEP_MostLikely, values = list(NutTrt = Ptreats,
                                                                  MeanPre2wksTemp = TempRange2015_16,
                                                                  MetDate = Dates2016)) %>% 
  mutate(fit.bt = exp(fit)/1000,
         se.fitL.bt = exp(fit-se.fit)/1000,
         se.fitU.bt = exp(fit+se.fit)/1000,
         trt = "Ponly",
         N_uM = 0,
         Basis = "MS",
         Res = "NEP")

## Nfix -----
# No P effect
p2016.A.Nfix <- predict_gam(Areal_2016_Nfix_MostLikely, values = list(MeanPre2wksTemp = TempRange2015_16,
                                                                    MetDate2 = Dates2016)) %>% 
  mutate(fit.bt = exp(fit)/1000,
         se.fitL.bt = exp(fit-se.fit)/1000,
         se.fitU.bt = exp(fit+se.fit)/1000,
         trt = "Ponly",
         N_uM = 0,
         P_uM = 0,
         NutTrt = 0,
         Basis = "Areal",
         Res = "Nfix")


p2016.MS.Nfix <- predict_gam(MS_2016_Nfix_MostLikely, values = list(MeanPre2wksTemp = TempRange2015_16,
                                                                  MetDate = Dates2016)) %>% 
  mutate(fit.bt = exp(fit)/1000,
         se.fitL.bt = exp(fit-se.fit)/1000,
         se.fitU.bt = exp(fit+se.fit)/1000,
         trt = "Ponly",
         N_uM = 0,
         P_uM = 0,
         NutTrt = 0,
         Basis = "MS",
         Res = "Nfix")

# Combine predictions ----
OrderedColNames <- c("MeanPre2wksTemp", "MetDate2", "NutTrt", "fit", "se.fit", "fit.bt", "se.fitL.bt", "se.fitU.bt", "trt", "N_uM", "Basis", "Res")

  pComb.A <- rbind(p2016.A.ER %>% 
                   select(all_of(OrderedColNames)), 
                 p2016.A.GPP %>% 
                   select(all_of(OrderedColNames)), 
                 p2016.A.NEP %>% 
                   select(all_of(OrderedColNames)), 
                 p2016.A.Nfix %>% 
                   select(all_of(OrderedColNames)))
  
  pComb.MS <- rbind(p2016.MS.ER %>% 
                     select(all_of(OrderedColNames)), 
                   p2016.MS.GPP %>% 
                     select(all_of(OrderedColNames)), 
                   p2016.MS.NEP %>% 
                     select(all_of(OrderedColNames)), 
                   p2016.MS.Nfix %>% 
                     select(all_of(OrderedColNames)))

  pComb <- rbind(pComb.A, pComb.MS) %>% 
            rename(P_uM = NutTrt) %>% 
            mutate(across(c(trt, Res, N_uM, P_uM, Basis), factor)) %>% 
            mutate(MetDate2 = ifelse((MetDate2 == "2016-07-27"), "Day1",
                                     ifelse((MetDate2 == "2016-08-02"), "Day2", as.character(MetDate2))),
                   MetDate2 = as.factor(MetDate2)) %>% 
            mutate(NcompInt = ifelse(P_uM == 0 | P_uM == 3.6,"Comp","NotComp"),
                   NcompInt = as.factor(NcompInt),
                   P_uM2 = fct_recode(P_uM, "0 µM-P" = "0",
                                      "0.8 µM-P" = "0.8",
                                      "1.6 µM-P" = "1.6",
                                      "3.6 µM-P" = "3.6",
                                      "4.8 µM-P" = "4.8",
                                      "6.5 µM-P" = "6.5"),
                   N_uM2 = fct_recode(N_uM, 
                                      "0 µM-" = "0")) %>% 
            mutate(P_uM2 = ifelse(
              (Res == "ER" & MetDate2 == "Day2" & Basis == "Areal") |
                (Res == "ER" & MetDate2 == "Day1" & Basis == "MS") |
                (Res == "ER" & MetDate2 == "Day2" & Basis == "MS") |
                (Res == "GPP" & MetDate2 == "Day2" & Basis == "Areal") |
                (Res == "GPP" & MetDate2 == "Day2" & Basis == "MS") |
                (Res == "NEP" & MetDate2 == "Day2" & Basis == "Areal") |
              (Res == "NEP" & MetDate2 == "Day2" & Basis == "MS") |
                (Res == "Nfix"), "0-6.5 µM-P", as.character(P_uM2)),
              Res = fct_relevel(Res, "ER", "GPP", "NEP", "Nfix"),
              P_uM2 = fct_relevel(P_uM2, "0-6.5 µM-P","0 µM-P",
                                  "0.8 µM-P",
                                  "1.6 µM-P",
                                  "3.6 µM-P",
                                  "4.8 µM-P",
                                  "6.5 µM-P"))

  # Plots ----
  ColorVals <- c("grey60","black", "#C7E9B4", "#7FCDBB", "#339900", "#2C7FB8", "#253494")
  
  
  p1 <- ggplot(pComb %>% 
                 filter(Basis == "Areal") %>% 
                 mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                    ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), 
               aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                                               color = as.factor(P_uM2),
                                               fill = as.factor(P_uM2),
                                               size = NcompInt)) +
    geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
    geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
    geom_ribbon(color = "transparent", alpha = 0.1) +
    geom_line() +
    facet_grid(Res ~ MD, scales = "free_y") +
    ylab("") +
    xlab("") +
    scale_color_manual(values = ColorVals, name = "Treatment")+
    scale_fill_manual(values = ColorVals, name = "Treatment")+
    scale_size_manual(values = c(1.2,0.5)) +
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
           color = "none", 
           size = "none") 
  
  p2 <- ggplot(pComb %>% 
                 filter(Basis == "MS") %>% 
                 mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                    ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                                            color = as.factor(P_uM2),
                                            fill = as.factor(P_uM2),
                                            size = NcompInt)) +
    geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
    geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
    geom_ribbon(color = "transparent", alpha = 0.1) +
    geom_line() +
    facet_grid(Res ~ MD, scales = "free_y") +
    ylab("") +
    xlab("") +
    scale_color_manual(values = ColorVals, name = "Treatment")+
    scale_fill_manual(values = ColorVals, name = "Treatment")+
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_size_manual(values = c(1.2,0.5)) +
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
           size = "none") 
  
  # Print ----
  p1.g <- ggplotGrob(p1)
  p2.g <- ggplotGrob(p2)
  
  p1.gtf <- gtable_frame(p1.g, width = unit(1, "null"), height = unit(1, "null"))
  p2.gtf <- gtable_frame(p2.g, width = unit(1, "null"), height = unit(1, "null"))
  
  p12.gtf <- gtable_frame(gtable_cbind(p1.gtf, p2.gtf), width = unit(2,"null"), height = unit(1,"null"))
  
  png("05_Figures4MS/14_Fig6_2016_ArealandMS.png", units = "in", height = 13, width = 16, res = 300)
  grid.newpage()
  grid.draw(p12.gtf)
  # Areal labels
  grid.text("Areal ER",
            x = unit(0.02,"npc"), y = unit(0.85,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
            x = unit(0.0425,"npc"), y = unit(0.85,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  grid.text("Areal GPP",
            x = unit(0.02,"npc"), y = unit(0.62,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
            x = unit(0.0425,"npc"), y = unit(0.62,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  grid.text("Areal NEP",
            x = unit(0.02,"npc"), y = unit(0.4,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
            x = unit(0.0425,"npc"), y = unit(0.4,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  grid.text(expression(bold(paste("Areal ",N[2]," fixation"))),
            x = unit(0.02,"npc"), y = unit(0.17,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM N ", m^{-2}," ", h^{-1},")")),
            x = unit(0.0425,"npc"), y = unit(0.17,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  
  # MS labels
  grid.text("MS ER",
            x = unit(0.44,"npc"), y = unit(0.85,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
            x = unit(0.465,"npc"), y = unit(0.85,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  grid.text("MS GPP",
            x = unit(0.44,"npc"), y = unit(0.62,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
            x = unit(0.465,"npc"), y = unit(0.62,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  grid.text("MS NEP",
            x = unit(0.44,"npc"), y = unit(0.4,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
            x = unit(0.465,"npc"), y = unit(0.4,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  grid.text(expression(bold(paste("MS ",N[2]," fixation"))),
            x = unit(0.44,"npc"), y = unit(0.17,"npc"), 
            gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
  
  grid.text(expression(paste("(mM N g ", AFDM^{-1}," ", h^{-1},")")),
            x = unit(0.465,"npc"), y = unit(0.17,"npc"), 
            gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
  
  grid.text("Temperature (°C)", 
            x = unit(0.5,"npc"), y = unit(0.015,"npc"), 
            gp=gpar(fontsize = 28, fontface = "bold"))
  
  dev.off()
  
  # save image ----
  # save.image("02b_Script_SavedImages/14_Fig1_BM_Met_2016_Rdat")
  # load("02b_Script_SavedImages/14_Fig1_BM_Met_2016_Rdat")
  