# Generate Fig 6
# Areal and ms responses from 2015
# JMH Jan 2023, updated Jun 2023

# LOAD LIBRARIES ----
library(tidyverse)
library(mgcv)
library(tidymv)
library(RColorBrewer)
library(grid)
library(egg)


# Most likely models ----
## 2015 - mass specific ----
MS_2015_ER_MostLikely <- readRDS("03_Model_RDS/MS_2015_ER_mostlikely.rds")
MS_2015_GPP_MostLikely <- readRDS("03_Model_RDS/MS_2015_GPP_mostlikely.rds")
MS_2015_NEP_MostLikely <-   readRDS("03_Model_RDS/MS_2015_NEP_mostlikely.rds")
MS_2015_Nup_MostLikely <- readRDS("03_Model_RDS/MS_2015_Nup_mostlikely.rds")
MS_2015_Nfix_MostLikely <- readRDS("03_Model_RDS/MS_2015_Nfix_mostlikely.rds")
MS_2015_NAssim_MostLikely <- readRDS("03_Model_RDS/MS_2015_NAssim_mostlikely.rds")

## 2015 - areal ----
Areal_2015_ER_MostLikely <- readRDS("03_Model_RDS/Areal_2015_ER_mostlikely.rds")
Areal_2015_GPP_MostLikely <- readRDS("03_Model_RDS/Areal_2015_GPP_mostlikely.rds")
Areal_2015_NEP_MostLikely <- readRDS("03_Model_RDS/Areal_2015_NEP_mostlikely.rds")
Areal_2015_Nup_MostLikely <- readRDS("03_Model_RDS/Areal_2015_Nup_mostlikely.rds")
Areal_2015_Nfix_MostLikely <- readRDS("03_Model_RDS/Areal_2015_Nfix_mostlikely.rds")
Areal_2015_TotAssim_MostLikely <- readRDS("03_Model_RDS/Areal_2015_NTotAssim_mostlikely.rds")
Areal_2015_AFDM_MostLikely <- readRDS("03_Model_RDS/Areal_2015_AFDM_mostlikely.rds")




# standardized ranges ----
  TempRange2015_16 <- seq(7.9, 25.5, by = 0.1)
  DINtreats <- c(0,1.8, 3.6, 7.1, 10.7, 14.3)
  Ptreats <- c(0, 0.8, 1.6, 3.6, 4.8, 6.5)
  Dates2015 <- c("2015-07-15", "2015-07-27")
  Dates2016 <- c("2016-07-27", "2016-08-02")

  
# MAKE PREDICTIONS ----
  ## Biomass ----
  
  ## ER -----
    p2015.A.ER <- predict_gam(Areal_2015_ER_MostLikely, values = list(NutTrt = DINtreats,
                                                     MeanPre2wksTemp = TempRange2015_16,
                                                     MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "Areal")

    p2015.MS.ER0 <- predict_gam(MS_2015_ER_MostLikely, values = list(NutTrt = DINtreats,
                                                                      MeanPre2wksTemp = TempRange2015_16,
                                                                      MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "MS")
    
    TempRange2015_16_For2015MsER <- as.data.frame(cbind(MeanPre2wksTemp = rep(TempRange2015_16, each = 12), 
                                                         MetDate2 = rep(as.character(p2015.MS.ER0$MetDate2), times = length(TempRange2015_16)),
                                                         NutTrt = rep(p2015.MS.ER0$NutTrt, times = length(TempRange2015_16)))) %>% 
      mutate(MeanPre2wksTemp = as.numeric(MeanPre2wksTemp), 
             NutTrt = as.numeric(NutTrt))
    
    p2015.MS.ER <- TempRange2015_16_For2015MsER %>% 
      full_join(p2015.MS.ER0, by = c("MetDate2", "NutTrt"))
    
  ## GPP ----
    p2015.A.GPP <- predict_gam(Areal_2015_GPP_MostLikely, values = list(NutTrt = DINtreats,
                                                                      MeanPre2wksTemp = TempRange2015_16,
                                                                      MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "Areal")
    
    p2015.MS.GPP0 <- predict_gam(MS_2015_GPP_MostLikely, values = list(NutTrt = DINtreats,
                                                                    MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "MS")
    
    TempRange2015_16_For2015MsGpp <- as.data.frame(cbind(MeanPre2wksTemp = rep(TempRange2015_16, each = 12), 
                                                   MetDate2 = rep(as.character(p2015.MS.GPP0$MetDate2), times = length(TempRange2015_16)),
                                                   NutTrt = rep(p2015.MS.GPP0$NutTrt, times = length(TempRange2015_16)))) %>% 
                                      mutate(MeanPre2wksTemp = as.numeric(MeanPre2wksTemp), 
                                             NutTrt = as.numeric(NutTrt))
    
    p2015.MS.GPP <- TempRange2015_16_For2015MsGpp %>% 
                      full_join(p2015.MS.GPP0, by = c("MetDate2", "NutTrt"))
  
  ## NEP ----
    p2015.A.NEP <- predict_gam(Areal_2015_NEP_MostLikely, values = list(NutTrt = DINtreats,
                                                                        MeanPre2wksTemp = TempRange2015_16,
                                                                        MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "Areal")
    
    p2015.MS.NEP0 <- predict_gam(MS_2015_NEP_MostLikely, values = list(NutTrt = DINtreats,
                                                                      MeanPre2wksTemp = TempRange2015_16,
                                                                      MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "MS")
    
    p2015.MS.NEP <- TempRange2015_16_For2015MsGpp %>% 
              full_join(p2015.MS.NEP0, by = c("MetDate2", "NutTrt"))
    
  ## N uptake ----
    p2015.A.Nup <- predict_gam(Areal_2015_Nup_MostLikely, values = list(NutTrt = DINtreats,
                                                                        MeanPre2wksTemp = TempRange2015_16,
                                                                        MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "Areal")
    
    p2015.MS.Nup <- predict_gam(MS_2015_Nup_MostLikely, values = list(NutTrt = DINtreats,
                                                                      MeanPre2wksTemp = TempRange2015_16,
                                                                      MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "MS")
    
  
  ## N fixation ----
    p2015.A.Nfix <- predict_gam(Areal_2015_Nfix_MostLikely, values = list(NutTrt = DINtreats,
                                                                        MeanPre2wksTemp = TempRange2015_16,
                                                                        MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "Areal")
    
    p2015.MS.Nfix <- predict_gam(MS_2015_Nfix_MostLikely, values = list(NutTrt = DINtreats,
                                                                      MeanPre2wksTemp = TempRange2015_16,
                                                                      MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "MS")
  
  ## Total N assimilation
    p2015.A.NAssim <- predict_gam(Areal_2015_TotAssim_MostLikely, values = list(NutTrt = DINtreats,
                                                                          MeanPre2wksTemp = TempRange2015_16,
                                                                          MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "Areal")
    
    p2015.MS.NAssim <- predict_gam(MS_2015_NAssim_MostLikely, values = list(NutTrt = DINtreats,
                                                                        MeanPre2wksTemp = TempRange2015_16,
                                                                        MetDate2 = Dates2015)) %>% 
      mutate(fit.bt = exp(fit)/1000,
             se.fitL.bt = exp(fit-se.fit)/1000,
             se.fitU.bt = exp(fit+se.fit)/1000,
             trt = "Nonly",
             P_uM = 0,
             Basis = "MS")

# Combine predictions ----
    pComb_Areal <- rbind(p2015.A.ER %>% 
                           mutate(Res = "ER"), 
                         p2015.A.GPP %>% 
                           mutate(Res = "GPP"), 
                         p2015.A.NEP %>% 
                           mutate(Res = "NEP"), 
                         p2015.A.Nup %>% 
                           mutate(Res = "Nup"), 
                         p2015.A.Nfix %>% 
                           mutate(Res = "Nfix"), 
                         p2015.A.NAssim %>% 
                           mutate(Res = "NAssim"))
    
    pComb_MS <- rbind(p2015.MS.ER %>% 
                        mutate(Res = "ER"), 
                      p2015.MS.GPP %>% 
                        mutate(Res = "GPP"), 
                      p2015.MS.NEP %>% 
                        mutate(Res = "NEP"), 
                      p2015.MS.Nup %>% 
                        mutate(Res = "Nup"), 
                      p2015.MS.Nfix %>% 
                        mutate(Res = "Nfix"), 
                      p2015.MS.NAssim %>% 
                        mutate(Res = "NAssim"))

    pComb <- rbind(pComb_Areal, pComb_MS) %>% 
                rename(N_uM = NutTrt) %>% 
                mutate(across(c(trt, Res, N_uM, P_uM, Basis), factor)) %>% 
                mutate(MetDate2 = ifelse((MetDate2 == "2015-07-15"), "Day1",
                                        ifelse((MetDate2 == "2015-07-27"), "Day2", as.character(MetDate2))),
                       MetDate2 = as.factor(MetDate2)) %>% 
                mutate(N_uM2 = ifelse(Res == "ER" & MetDate2 == "Day2" & Basis == "Areal",
                                      "0-14.3 µM-N", as.character(N_uM)),
                        N_uM2 = fct_recode(N_uM2,
                                           "0-14.3 µM-N" = "0-14.3 µM-N",
                                           "0 µM-N" = "0",
                                          "1.8 µM-N" = "1.8",
                                          "3.6 µM-N" = "3.6",
                                          "7.1 µM-N" = "7.1",
                                          "10.7 µM-N" = "10.7",
                                          "14.3 µM-N" = "14.3"),
                       N_uM2 = fct_relevel(N_uM2, "0-14.3 µM-N","0 µM-N","1.8 µM-N","3.6 µM-N","7.1 µM-N","10.7 µM-N","14.3 µM-N"),
                       P_uM2 = fct_recode(P_uM, 
                                          "0 µM-P" = "0"),
                       Res = fct_recode(Res, "N up." = "Nup", "N fix." = "Nfix", "N up. + fix." = "NAssim"),
                       Res = fct_relevel(Res, "ER", "GPP", "NEP", "N up. + fix.", "N up.", "N fix."),
                       ) %>% 
                mutate(NcompInt = ifelse(N_uM == 0 | N_uM == 3.6,"Comp","NotComp"),
                       NcompInt = as.factor(NcompInt))
                       
    # putting in dummy row to hack MS legend            
    pComb[dim(pComb)[1]+1,] <- pComb[dim(pComb)[1],]
    pComb[dim(pComb)[1],]$fit <- as.numeric("NA")
    pComb[dim(pComb)[1],]$se.fit <- as.numeric("NA")
    pComb[dim(pComb)[1],]$se.fitL.bt <- as.numeric("NA")
    pComb[dim(pComb)[1],]$se.fitU.bt <- as.numeric("NA")
    pComb[dim(pComb)[1],]$N_uM2 <- "0-14.3 µM-N"
    
# Plots ----
    ColorVals <- c("grey60","black", "#FED976", "#339900", "#FD8D3C", "#F03B20", "#BD0026") ##FED976


                                                   
    p1 <- ggplot(pComb %>% 
                   filter(Basis == "Areal") %>% 
                   mutate(MD = ifelse(MetDate2 == "Day1", "MD1",
                                      ifelse(MetDate2 == "Day2", "MD2", "BLAH"))), aes(ymin = se.fitL.bt, ymax = se.fitU.bt, y = fit.bt, x = MeanPre2wksTemp, 
                                                                                       color = as.factor(N_uM2),
                                                                                       fill = as.factor(N_uM2),
                                                                                       size = NcompInt)) +
      geom_vline(xintercept = 11, linewidth = 0.6, color = "grey") +
      geom_vline(xintercept = 21, linewidth = 0.6, color = "grey") +
      geom_ribbon(color = "transparent", alpha = 0.1) +
      geom_line() + 
      facet_grid(Res ~ MD, scales = "free_y", 
                 switch = "y") +
      ylab("") +
      xlab("") +
      scale_color_manual(values = ColorVals, name = "Treatment")+
      scale_fill_manual(values = ColorVals, name = "Treatment")+
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
            strip.placement = "outside",
            strip.text.x = element_text(size = 20, face = "bold"),
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
                                              color = as.factor(N_uM2),
                                              fill = as.factor(N_uM2),
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
      scale_x_continuous(limits = c(5,25), breaks = c(5, 10, 15, 20, 25)) +
      theme(panel.background = element_rect(fill = "white", color = "white"),
            panel.border = element_rect(color = "black", fill = "NA", size = 1),
            axis.title.x = element_text(size = 22),
            axis.title.y = element_text(size = 60),
            axis.text = element_text(size = 18),
            axis.line = element_line(color = "black", size = 1),
            plot.background = element_rect(fill = "white", color =  "white"),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 20, face = "bold"),
            strip.text.y = element_blank(),
            legend.title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 20),
            legend.key.height = unit(1,"cm")) +
      guides(linetype = "none",
             size = "none") 
    
    # Print ----
    p1.g <- ggplotGrob(p1)
    p2.g <- ggplotGrob(p2)
    
    p1.gtf <- gtable_frame(p1.g, width = unit(1.5, "null"), height = unit(1, "null"))
    p2.gtf <- gtable_frame(p2.g, width = unit(1.5, "null"), height = unit(1, "null"))
    
    p12.gtf <- gtable_frame(gtable_cbind(p1.gtf, p2.gtf), width = unit(3,"null"), height = unit(1,"null"))
    
    # png("05_Figures4MS/13_Fig6_2015_ArealandMS.png", units = "in", height = 18, width = 18, res = 600)
    jpeg("05_Figures4MS/13_Fig5_2015_ArealandMS.jpg", units = "in", height = 18, width = 18, res = 350)
    grid.newpage()
    grid.draw(p12.gtf)
    
    # Areal labels
    grid.text("Areal ER",
              x = unit(0.02,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text("Areal GPP",
              x = unit(0.02,"npc"), y = unit(0.74,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.74,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text("Areal NEP",
              x = unit(0.02,"npc"), y = unit(0.58,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.58,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text("Areal Total N assim.",
              x = unit(0.02,"npc"), y = unit(0.43,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.43,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text("Areal N uptake",
              x = unit(0.02,"npc"), y = unit(0.28,"npc"),  
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.28,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text(expression(bold(paste("Areal ",N[2]," fixation"))),
              x = unit(0.02,"npc"), y = unit(0.12,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N ", m^{-2}," ", h^{-1},")")),
              x = unit(0.04,"npc"), y = unit(0.12,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    

    
    
    # Mass specific
    grid.text("MS ER",
              x = unit(0.45,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.47,"npc"), y = unit(0.9,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    
    grid.text("MS GPP",
              x = unit(0.45,"npc"), y = unit(0.74,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.47,"npc"), y = unit(0.74,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text("MS NEP",
              x = unit(0.45,"npc"), y = unit(0.58,"npc"),
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM C g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.47,"npc"), y = unit(0.58,"npc"),
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text("MS Total N assim.",
              x = unit(0.45,"npc"), y = unit(0.43,"npc"),
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.47,"npc"), y = unit(0.43,"npc"), 
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    
    grid.text("MS N uptake",
              x = unit(0.45,"npc"), y = unit(0.28,"npc"), 
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.47,"npc"), y = unit(0.28,"npc"),
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
   
    grid.text(expression(bold(paste("MS ",N[2]," fixation"))),
              x = unit(0.45,"npc"), y = unit(0.12,"npc"),
              gp=gpar(fontsize = 20, fontface = "bold"), rot = 90)
    
    grid.text(expression(paste("(mM N g ", AFDM^{-1}," ", h^{-1},")")),
              x = unit(0.47,"npc"), y = unit(0.12,"npc"),
              gp=gpar(fontsize = 18, fontface = "bold"), rot = 90)
    

    
    grid.text("Temperature (°C)", 
              x = unit(0.5,"npc"), y = unit(0.015,"npc"), 
              gp=gpar(fontsize = 28, fontface = "bold"))
    dev.off()
    
  # save.image
    # save.image("02b_Script_SavedImages/13_Fig2_BM_Met_2015_Rdat")
    load("02b_Script_SavedImages/13_Fig2_BM_Met_2015_Rdat")
    