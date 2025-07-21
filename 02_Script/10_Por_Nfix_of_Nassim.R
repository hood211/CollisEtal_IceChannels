# Analysis of N-fix/N-assim - 2015 data
# JMH Jan 2023, June 2023

# Libraries ----
# general
library(tidyverse)
# gam models
library(mgcv)
library(mgcViz)
library(MuMIn)
# logit transformation
library(car)
# segmented regression
library(tidygam)

# Data ----
chan <- read.csv("01_Data/IceChan_allyrs_finaldata_imputed_11102022.csv", row.names = 1) %>% 
  mutate(across(c(MetDate, UpDate, NfixDate), as.POSIXct, format = "%m/%d/%Y")) %>% 
  mutate(StartDate = as.POSIXct(StartDate, format = "%Y-%m-%d"),
         # Adding 2015 minimum to Nuptake to allow modeling
         NUp_uM_N_m2_hr = NUp_uM_N_m2_hr,
         TotAssim_uM_N_m2_hr = NUp_uM_N_m2_hr + Nfix_uM_N_m2h,
         PorNfixOfNassim = Nfix_uM_N_m2h/TotAssim_uM_N_m2_hr,
         PorNfixOfNassim = ifelse(PorNfixOfNassim < 0,0,PorNfixOfNassim),
         # remaps some proportations between 0.025 and 0.975
         logit_PorNfixOfNassim = logit(PorNfixOfNassim, percents = FALSE, adjust = 0.001),
         cauch_PorNfixOfNassim = ifelse(tan(pi * (PorNfixOfNassim-0.5)) < -160, -160, tan(pi * (PorNfixOfNassim-0.5))))

chanF <- chan %>% 
  filter(Year == "2015") %>% 
  dplyr::select(channel, MetDate, MeanPre2wksTemp, tempF, N_uM, PorNfixOfNassim, logit_PorNfixOfNassim, cauch_PorNfixOfNassim,
                NUp_uM_N_m2_hr, Nfix_uM_N_m2h, TotAssim_uM_N_m2_hr) %>% 
  mutate(MetDate2 = case_when(MetDate == "2015-07-15" ~ "MD1",
                              MetDate == "2015-07-27" ~ "MD2"),
         MetDate2 = as.factor(MetDate2),
         N_uM_l = log(N_uM))

# GAM fit ----
## model family ----
ModFam <- betar(link = "cauchit",eps= 0.001) #0.0001

## model selection ----
g1_TxN_byDate <- gam(PorNfixOfNassim ~ te(MeanPre2wksTemp, N_uM, bs = "ts", by = MetDate2, k = c(5,5)) + s(MetDate2, bs = "re"),
                     data = chanF,
                     optimizer = "efs",
                     family = ModFam,
                     method = "ML")

g1_TxN_byDate_a <- gam(PorNfixOfNassim ~ te(MeanPre2wksTemp, N_uM, bs = "ts", by = MetDate2, k = c(5,5)),
                       data = chanF,
                       optimizer = "efs",
                       family = ModFam,
                       method = "ML")

g2_TxN <- gam(PorNfixOfNassim ~ te(MeanPre2wksTemp, N_uM, bs = "ts", k = c(5,5)) + s(MetDate2, bs = "re"),
              data = chanF,
              optimizer = "efs",
              family = ModFam,
              method = "ML")

g2_TxN_a <- gam(PorNfixOfNassim ~ te(MeanPre2wksTemp, N_uM, bs = "ts", k = c(5,5)),
                data = chanF,
                optimizer = "efs",
                family = ModFam,
                method = "ML")

g3_TbyDate_p_NbyDate <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, by = MetDate2, k = 5) + s(N_uM, by = MetDate2, k = 6) + s(MetDate2, bs = "re"),
                            data = chanF,
                            optimizer = "efs",
                            family = ModFam,
                            method = "ML")

g3_TbyDate_p_NbyDate_a <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, by = MetDate2, k = 5) + s(N_uM, by = MetDate2, k =6),
                              data = chanF,
                              optimizer = "efs",
                              family = ModFam,
                              method = "ML")

g4_T_p_NbyDate <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, k = 5) + s(N_uM, by = MetDate2, k = 6) + s(MetDate2, bs = "re"),
                      data = chanF,
                      optimizer = "efs",
                      family = ModFam,
                      method = "ML")
g4_T_p_NbyDate_a <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, k = 5) + s(N_uM, by = MetDate2, k = 6),
                        data = chanF,
                        optimizer = "efs",
                        family = ModFam,
                        method = "ML")

g5_TbyDate_p_N <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, by = MetDate2, k = 5) + s(N_uM, k = 6) + s(MetDate2, bs = "re"),
                      data = chanF,
                      optimizer = "efs",
                      family = ModFam,
                      method = "ML")

g5_TbyDate_p_N_a <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, by = MetDate2, k = 5) + s(N_uM, k = 6),
                        data = chanF,
                        optimizer = "efs",
                        family = ModFam,
                        method = "ML")

g6_T_p_N <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, k = 5) + s(N_uM, k = 6) + s(MetDate2, bs = "re"),
                data = chanF,
                optimizer = "efs",
                family = ModFam,
                method = "ML")

g6_T_p_N_a <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, k = 5) + s(N_uM, k = 6),
                  data = chanF,
                  optimizer = "efs",
                  family = ModFam,
                  method = "ML")


g7_T <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, k = 6) + s(MetDate2, bs = "re"),
            data = chanF,
            optimizer = "efs",
            family = ModFam,
            method = "ML")

g7_T_a <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, k = 6),
              data = chanF,
              optimizer = "efs",
              family = ModFam,
              method = "ML")

g8_N <- gam(PorNfixOfNassim ~ s(N_uM, k = 6) + s(MetDate2, bs = "re"),
            data = chanF,
            optimizer = "efs",
            family = ModFam,
            method = "ML")

g8_N_a <- gam(PorNfixOfNassim ~ s(N_uM, k = 6),
              data = chanF,
              optimizer = "efs",
              family = ModFam,
              method = "ML")

g9_interceptOnly <- gam(PorNfixOfNassim ~ 1 + s(MetDate2, bs = "re"),
                        data = chanF,
                        optimizer = "efs",
                        family = ModFam,
                        method = "ML")



betareg_ModelSel <- MuMIn::model.sel(
  g1_TxN_byDate,
  g2_TxN,
  g3_TbyDate_p_NbyDate,
  g4_T_p_NbyDate,
  g5_TbyDate_p_N,
  g6_T_p_N,
  g7_T,
  g8_N,
  g1_TxN_byDate_a,
  g2_TxN_a,
  g3_TbyDate_p_NbyDate_a,
  g4_T_p_NbyDate_a,
  g5_TbyDate_p_N_a,
  g6_T_p_N_a,
  g7_T_a,
  g8_N_a,
  g9_interceptOnly, 
  rank = AICc,
  extra = c(R2 = function(x) summary(x)$r.sq))

betareg_ModelSel

## most likely model
# model g5 has same R2 as model g2, using simpler model
# date random effect isn't doing anything - drop
# https://stats.stackexchange.com/questions/99425/distribution-for-percentage-data

betareg_MostLikely <- gam(PorNfixOfNassim ~ s(MeanPre2wksTemp, by = MetDate2, k = 5) + s(N_uM, k = 6) ,
                          optimizer = "efs",
                          data = chanF,
                          family = ModFam,
                          method = "REML")


summary(betareg_MostLikely)
# this is less than ideal, but the best I can achieve
check(getViz(betareg_MostLikely))  
print(plot(getViz(betareg_MostLikely), allTerms = T), pages = 1)

## Plot ----
# back trans from cauchit
# https://www.rdocumentation.org/packages/VGAM/versions/1.0-3/topics/cauchit
CauchBTfun <- function(x) 0.5 + atan(x)/pi

### N plot ----
chanF_fake <-predict_gam(betareg_MostLikely,
                         values = list(MeanPre2wksTemp = c(9.225614, 24.488832, 12.524347, 15.771314, 19.746454),
                                       N_uM = seq(0.1,14.4, length = 50),
                                       MetDate2 = "MD2"),
                         # backtransform from logit
                         tran_fun = CauchBTfun) %>% 
              mutate(tempF = case_when(MeanPre2wksTemp == 9.225614 ~ "9.2",
                                       MeanPre2wksTemp == 12.524347 ~ "12.5",
                                       MeanPre2wksTemp == 15.771314 ~ "15.7",
                                       MeanPre2wksTemp == 19.746454 ~ "19.7",
                                       MeanPre2wksTemp == 24.488832 ~ "24.5"),
                     tempF = as.factor(tempF),
                     tempF = fct_relevel(tempF, "9.2", "12.5", "15.7", "19.7", "24.5"))

jitter <- position_jitter(width = 0.2, height = 0)
ColorValuesp2 <- c("#0e5ca0", "#81b8d6", "#FED976", "#ef936f", "#bc0019")
p2 <- ggplot() +
  geom_line(data = chanF_fake, aes(y = PorNfixOfNassim, x = N_uM, color = tempF), 
            show.legend = FALSE, linewidth = 1.25) +
  geom_point(data = chanF %>% 
               mutate(tempF2 = case_when(tempF == "A" ~ "9.2",
                                        tempF == "B" ~ "12.5",
                                        tempF == "C" ~ "15.7",
                                        tempF == "D" ~ "19.7",
                                        tempF == "E" ~ "24.5"),
                      tempF2 = as.factor(tempF2),
                      tempF2 = fct_relevel(tempF2, "9.2", "12.5", "15.7", "19.7", "24.5")), 
             aes(y = PorNfixOfNassim, x= N_uM, fill = tempF2, shape = MetDate2),
             position = jitter, size = 5) +
  scale_shape_manual(values = c(21,22)) +
  scale_color_manual(values = ColorValuesp2) +
  scale_linetype_manual(guide = "none") +
  scale_fill_manual(values = ColorValuesp2) +
  guides(fill = guide_legend(override.aes=list(shape=21, size = 6), "Temperature"),
         shape = guide_legend(override.aes=list(size = 6), "Measurement date")) +
  ylab(expression(frac(paste(N[2], " fixation"), "Total N assimilation"))) +
  xlab("µM N") +
  theme_bw()+
  theme(panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_rect(color = "black", fill = "NA", size = 1),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        axis.line = element_line(color = "black", size = 1),
        plot.background = element_rect(fill = "white", color =  "white"),
        strip.background = element_blank(),
        strip.placement = "outside",
        # strip.background = element_rect(fill = "grey", color =  "white"),
        strip.text.x = element_text(size = 20, face = "bold"),
        # strip.text.y = element_text(size = 28, face = "bold"),
        # strip.text.y = element_text(size = 18, face = "bold"),
        strip.text.y = element_blank(),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1,"cm")) 


## Fig 8: data vis plot ----
ColorVals1 <- c("black", "#FED976", "#FEB24C", "#FD8D3C", "#F03B20", "#BD0026")

chanF_p <- chanF %>% 
  mutate(# same as fig 6
    N_uM2 = as.factor(N_uM),
    N_uM2 = fct_recode(N_uM2,
                       "0 µM-N" = "0.11",
                       "1.8 µM-N" = "1.9",
                       "3.6 µM-N" = "3.68",
                       "7.1 µM-N" = "7.25",
                       "10.7 µM-N" = "10.82",
                       "14.3 µM-N" = "14.4"),
    N_uM2 = fct_relevel(N_uM2, "0 µM-N","1.8 µM-N","3.6 µM-N","7.1 µM-N","10.7 µM-N","14.3 µM-N"))



p1 <- ggplot() +
  geom_abline(intercept = 0, slope = 0.3, linetype = "dashed", color = "grey20") +
  geom_abline(intercept = 0, slope = 0.2, linetype = "dashed", color = "grey20") +
  geom_abline(intercept = 0, slope = 0.1, linetype = "dashed", color = "grey20") +
  geom_abline(intercept = 0, slope = 0.01, linetype = "dashed", color = "grey20") +
  geom_point(data = chanF_p, 
             aes(y = Nfix_uM_N_m2h, x = TotAssim_uM_N_m2_hr, 
                 fill = as.factor(N_uM2),
                 shape = MetDate2),
             size = 5)+
  scale_fill_manual(values = ColorVals1, name = "µM N") +
  scale_shape_manual(values = c(21,22), name = "Measurement day") +
  theme_bw()+
  ylab(expression(paste("Areal ", paste(N[2]," fixation (mM N ",m^-2," ",h^-1,")")))) +
  xlab(expression(paste("Areal Total N assim. (mM N ",m^-2," ",h^-1,")"))) +
  xlim(0,4000) +
  ylim(0,450) +
  guides( fill = guide_legend(override.aes=list(shape=21, size = 6))) +
  annotate(geom = "text", x = 1100, y = 450, label = expression(paste(frac(paste(N[2], " fixation"), "Total N assim."), " =")),
           hjust = 1, vjust = 0.5, color = "steelblue4", size = 6) +
  annotate(geom = "text", x = 1400, y = 450, label = "0.3",
           hjust = 1, vjust = 0.5, color = "steelblue4", size = 6) +
  annotate(geom = "text", x = 2125, y = 435, label = "0.2",
           hjust = 1, vjust = 0.5, color = "steelblue4", size = 6) +
  annotate(geom = "text", x = 3250, y = 335, label = "0.1",
           hjust = 1, vjust = 0.5, color = "steelblue4", size = 6) +
  annotate(geom = "text", x = 4000, y = 50, label = "0.01",
           hjust = 1, vjust = 0.5, color = "steelblue4", size = 6) +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        panel.border = element_rect(color = "black", fill = "NA", size = 1),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        axis.line = element_line(color = "black", size = 1),
        plot.background = element_rect(fill = "white", color =  "white"),
        strip.background = element_blank(),
        strip.placement = "outside",
        # strip.background = element_rect(fill = "grey", color =  "white"),
        strip.text.x = element_text(size = 20, face = "bold"),
        # strip.text.y = element_text(size = 28, face = "bold"),
        # strip.text.y = element_text(size = 18, face = "bold"),
        strip.text.y = element_blank(),
        legend.title = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 20),
        legend.key.height = unit(1,"cm")) 

library(cowplot)
png("05_Figures4MS/10_Fig7_porNfix.png", units = "in", height = 14, width = 11, res = 300)
plot_grid(p1,p2, labels = c("a", "b"),
          label_size = 24,
          label_fontface = "bold",
          nrow = 2,
          ncol = 1)
dev.off()


# Segmented regression ----
library(segmented)
# using logit instead of cauchit because logit model has better residuals


seg_lmModel_all.lm <- lm(logit_PorNfixOfNassim ~ N_uM + MeanPre2wksTemp, 
                         data = chanF)


seg_lmModel_all <- segmented(seg_lmModel_all.lm, 
                             seg.Z = ~ N_uM, 
                             psi = list(N_uM = 1.6))

# lm summary
summary(seg_lmModel_all.lm)
plot(seg_lmModel_all.lm)

# seg regression summary
summary(seg_lmModel_all)
confint(seg_lmModel_all)
plot(seg_lmModel_all)

# save/load ----
# save.image("02b_Script_SavedImages/10_Por_Nfix_of_Nassim_Rdat")
# load("02b_Script_SavedImages/10_Por_Nfix_of_Nassim_Rdat")
