# Generate Fig 3
# Areal biomass  2015-2017
# JMH Jan 2023, Jul 24

# LOAD LIBRARIES ----
library(tidyverse)
library(mgcv)
library(tidymv)
library(RColorBrewer)
library(grid)
library(egg)


# Biomass ----
## Most likely models ----
Areal_2015_AFDM_MostLikely <- readRDS("03_Model_RDS/Areal_2015_AFDM_mostlikely.rds")
Areal_2016_AFDM_MostLikely <- readRDS("03_Model_RDS/Areal_2016_AFDM_mostlikely.rds")
Areal_2017_AFDM_MostLikely <- readRDS("03_Model_RDS/Areal_2017_AFDM_mostlikely.rds")

## standardized ranges ----
TempRange2015_16 <- seq(7.9, 25.5, by = 0.1)
TempRange2017 <- seq(7.9, 21, by = 0.1)
DINtreats <- c(0,1.8, 3.6, 7.1, 10.7, 14.3)
Ptreats <- c(0, 0.8, 1.6, 3.6, 4.8, 6.5)
NP <- c("0.31", "0.93", "10.22")
NPratio_date <- c("0.31_2017-07-11", "0.31_2017-07-22", "0.93_2017-07-11", "0.93_2017-07-22", "10.22_2017-07-11", "10.22_2017-07-22")
Dates2015 <- c("2015-07-15", "2015-07-27")
Dates2016 <- c("2016-07-27", "2016-08-02")
Dates2017 <- c("2017-07-11", "2017-07-22")


## MAKE PREDICTIONS ----
### 2015 -----
p2015.A.AFDM <- predict_gam(Areal_2015_AFDM_MostLikely, values = list(NutTrt = DINtreats,
                                                                  MeanPre2wksTemp = TempRange2015_16,
                                                                  MetDate2 = Dates2015)) %>% 
                              mutate(fit.bt = exp(fit),
                                     se.fitL.bt = exp(fit-se.fit),
                                     se.fitU.bt = exp(fit+se.fit),
                                     trt = "Nonly",
                                     P_uM = 0,
                                     Basis = "Areal") %>% 
                              rename(N_uM = NutTrt) %>% 
                              mutate(across(c(trt, N_uM, P_uM, Basis), factor)) %>% 
                              mutate(MetDate2 = ifelse((MetDate2 == "2015-07-15"), "Day1",
                                                       ifelse((MetDate2 == "2015-07-27"), "Day2", as.character(MetDate2))),
                                     MetDate2 = as.factor(MetDate2)) %>% 
                              mutate(N_uM2 = fct_recode(N_uM,
                                                        "0 µM-N" = "0",
                                                        "1.8 µM-N" = "1.8",
                                                        "3.6 µM-N" = "3.6",
                                                        "7.1 µM-N" = "7.1",
                                                        "10.7 µM-N" = "10.7",
                                                        "14.3 µM-N" = "14.3"))%>% 
                            mutate(NcompInt = ifelse(N_uM == 0 | N_uM == 3.6,"Comp","NotComp"),
                                   NcompInt = as.factor(NcompInt))

### 2016 -----
p2016.A.AFDM <- predict_gam(Areal_2016_AFDM_MostLikely, values = list(MeanPre2wksTemp = TempRange2015_16,
                                                                      MetDate2 = Dates2016)) %>% 
                            mutate(fit.bt = exp(fit),
                                   se.fitL.bt = exp(fit-se.fit),
                                   se.fitU.bt = exp(fit+se.fit),
                                   trt = "Ponly",
                                   N_uM = 0,
                                   P_uM = 0,
                                   Basis = "Areal") %>% 
                            mutate(across(c(trt, N_uM, P_uM, Basis), factor)) %>% 
                            mutate(MetDate2 = ifelse((MetDate2 == "2016-07-27"), "Day1",
                                                     ifelse((MetDate2 == "2016-08-02"), "Day2", as.character(MetDate2))),
                                   MetDate2 = as.factor(MetDate2)) %>% 
                            mutate(P_uM2 = fct_recode(P_uM, "0-6.5 µM-P" = "0"),
                                   N_uM2 = fct_recode(N_uM, 
                                                      "0 µM-N" = "0"))


### 2017 ----
p2017.A.AFDM <- predict_gam(Areal_2017_AFDM_MostLikely, values = list(MeanPre2wksTemp = TempRange2017,
                                                                          MetDate2 = Dates2017)) %>% 
                            mutate(fit.bt = exp(fit),
                                   se.fitL.bt = exp(fit-se.fit),
                                   se.fitU.bt = exp(fit+se.fit),
                                   trt = "N+P",
                                   Basis = "Areal") %>% 
                            mutate(across(c(trt, N_uM, Basis), factor)) %>% 
                            mutate(MetDate2 = ifelse((MetDate2 == "2017-07-11"), "Day1",
                                                     ifelse((MetDate2 == "2017-07-22"), "Day2", as.character(MetDate2))),
                                   MetDate2 = as.factor(MetDate2)) %>% 
                            mutate(N_uM2 = fct_recode(N_uM,
                                                      "0 µM-N" = "0.11",
                                                      "3.6 µM-N" = "3.68"))






                        


## Plots ----
  ColorVals2015 <- c("black", "#FED976", "#339900", "#FD8D3C", "#F03B20", "#BD0026")
  ColorVals2016 <- c("black")
  ColorVals2017 <- c("black", "#339900")
  
  pAFDM_2015 <- ggplot(p2015.A.AFDM %>% 
                         mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                            ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                                               color = as.factor(N_uM2),
                                               fill = as.factor(N_uM2), 
                                               size = NcompInt)) +
    geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
    geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
    geom_ribbon(color = "transparent", alpha = 0.1) +
    geom_line() +
    ylab("") +
    xlab("Temperature (°C)") +
    facet_wrap(vars(MD)) +
    scale_color_manual(values = ColorVals2015, name = "Treatment")+
    scale_fill_manual(values = ColorVals2015, name = "Treatment")+
    scale_size_manual(values = c(1.2,0.5)) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_x_continuous(limits = c(5,25), breaks = c(5,10,15, 20, 25)) +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_rect(color = "black", fill = "NA", size = 1),
          axis.title.x = element_text(size = 32),
          axis.title.y = element_text(size = 60),
          axis.text = element_text(size = 24),
          axis.line = element_line(color = "black", size = 1),
          plot.background = element_rect(fill = "white", color =  "white"),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 32, face = "bold"),
          strip.text.y = element_blank(),
          legend.position = c(0.6,0.66),
          legend.title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 20),
          legend.key.height = unit(1,"cm")) +
    guides(size = "none")

  pAFDM_2016 <-   ggplot(p2016.A.AFDM %>% 
                           mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                              ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                             color = as.factor(P_uM2),
                             fill = as.factor(P_uM2))) +
    geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
    geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
      geom_ribbon(color = "transparent", alpha = 0.25) +
      geom_line(size = 1.25) +
      ylab("") +
      xlab("Temperature (°C)") +
      facet_wrap(vars(MD)) +
      scale_color_manual(values = ColorVals2016, name = "Treatment")+
      scale_fill_manual(values = ColorVals2016, name = "Treatment")+
      scale_linetype_manual(values = c("dashed", "solid")) +
    scale_x_continuous(limits = c(5,25), breaks = c(5,10,15, 20, 25)) +
    scale_y_continuous(limits = c(0,60), breaks = c(0, 20, 40, 60)) +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          panel.border = element_rect(color = "black", fill = "NA", size = 1),
          axis.title.x = element_text(size = 32),
          axis.title.y = element_text(size = 60),
          axis.text = element_text(size = 24),
          axis.line = element_line(color = "black", size = 1),
          plot.background = element_rect(fill = "white", color =  "white"),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 32, face = "bold"),
          strip.text.y = element_blank(),
          legend.position = c(0.6,0.88),
          legend.title = element_text(size = 24, face = "bold"),
          legend.text = element_text(size = 20),
          legend.key.height = unit(1,"cm")) 
  
  pAFDM_2017 <-   ggplot(p2017.A.AFDM %>% 
                           mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                              ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                                           color = N_uM2,
                                           fill = N_uM2)) +
    geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
    geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
      geom_ribbon(color = "transparent", alpha = 0.25) +
      geom_line(size = 1.25) +
      ylab(expression(paste("Biomass (g AFDM ", m^{-2},")"))) +
      xlab("Temperature (°C)") +
      facet_wrap(vars(MD)) +
      scale_color_manual(values = ColorVals2017, name = "Treatment")+
      scale_fill_manual(values = ColorVals2017, name = "Treatment")+
      scale_linetype_manual(values = c("dashed", "solid")) +
      scale_x_continuous(limits = c(5,25), breaks = c(5,10,15, 20, 25)) +
      theme(panel.background = element_rect(fill = "white", color = "white"),
            panel.border = element_rect(color = "black", fill = "NA", size = 1),
            axis.title.x = element_text(size = 32),
            axis.title.y = element_text(size = 22),
            axis.text = element_text(size = 24),
            axis.line = element_line(color = "black", size = 1),
            plot.background = element_rect(fill = "white", color =  "white"),
            strip.background = element_blank(),
            # strip.background = element_rect(fill = "grey", color =  "white"),
            strip.text.x = element_text(size = 32, face = "bold"),
            # strip.text.y = element_text(size = 28, face = "bold"),
            strip.text.y = element_blank(),
            legend.position = c(0.6,0.75),
            legend.title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 20),
            legend.key.height = unit(1,"cm")) 

  # Print ----
  p1.gtf <- gtable_frame(ggplotGrob(pAFDM_2015), width = unit(1, "null"), height = unit(1, "null"))
  p2.gtf <- gtable_frame(ggplotGrob(pAFDM_2016), width = unit(1, "null"), height = unit(1, "null"))
  
  p123.gtf <- gtable_frame(gtable_rbind(p1.gtf, p2.gtf), width = unit(2,"null"), height = unit(1,"null"))
  
  png("05_Figures4MS/16_Fig2_2yrAFDM.png", units = "in", height = 14, width = 16, res = 300)
  grid.newpage()
  grid.draw(p123.gtf)
    grid.text(expression(paste("Biomass (g AFDM ", m^{-2},")")), x = unit(0.04,"npc"), y = unit(0.5,"npc"), gp=gpar(fontsize = 36), rot = 90)
    grid.text("a) N-only", x = unit(0.118, "npc"), y = unit(0.93, "npc"), gp=gpar(fontsize = 36, fontface = "bold"), hjust = 0)
    grid.text("b) P-only", x = unit(0.118, "npc"), y = unit(0.425, "npc"), gp=gpar(fontsize = 36, fontface = "bold"), hjust = 0)
  dev.off()
  
  png("05_Figures4MS/16_FigS10_2017PAFDM.png", units = "in", height = 5, width = 12, res = 300)
  pAFDM_2017
  dev.off()
  
  
# save.image
  # save.image("02b_Script_SavedImages/16_Biomass_PerNfixerPlots.Rdata")
  # load("02b_Script_SavedImages/16_Biomass_PerNfixerPlots_Rdat")
  